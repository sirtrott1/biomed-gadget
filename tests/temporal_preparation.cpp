#include "../src/scientific_temporal_prepare.hpp"
#include <iostream>
#include <rfc/core/builder.hpp>

int main() {
    try {
        using namespace gadget;
        using namespace rfc;
        core::RelationSignature signature;
        const auto type = signature.add_type("declared.multiway");
        core::RelationalComplexBuilder builder(signature);
        for (std::uint64_t k = 10; k < 16; ++k)
            builder.declare_witness_identity(core::IdentityKey{k});
        builder.add_grade_one(core::RelationKey{101}, type,
            {{core::IdentityKey{10}, true, false, Rational{1}, {}, core::OccurrenceKey{1001}},
             {core::IdentityKey{11}, true, true, Rational{2}, {}, core::OccurrenceKey{1002}},
             {core::IdentityKey{12}, true, false, Rational{3}, {}, core::OccurrenceKey{1003}}});
        builder.add_grade_one(core::RelationKey{102}, type,
            {{core::IdentityKey{13}, true, false, Rational{5}, {}, core::OccurrenceKey{1004}},
             {core::IdentityKey{13}, true, true, Rational{5}, {}, core::OccurrenceKey{1005}},
             {core::IdentityKey{14}, false, false, Rational{7}, {}, core::OccurrenceKey{1006}}});
        builder.add_higher(core::Grade{2}, core::RelationKey{201}, type,
            {{1, true, false, Rational{2}, {}, core::OccurrenceKey{2001}}});
        auto source = std::make_shared<const core::RelationalComplex>(builder.freeze());
        auto owner = std::make_shared<const NativeStudy>(source,
            std::vector<field::ExactMetric>{
                field::ExactMetric::diagonal({2, 3, 5, 7, 11, 13}),
                field::ExactMetric::diagonal({17, 19}), field::ExactMetric::diagonal({23})});
        scientific_detail::PreparedTemporalStudy prepared(owner);
        owner.reset(); // The prepared consumer pins the complete admitted owner.
        std::size_t checks = 0;
        for (std::uint16_t k = 0; k < 3; ++k) {
            const auto grade = core::Grade{k};
            const auto n = source->grade(grade).cell_count();
            for (const Rational &lambda : {Rational{1}, Rational{2, 3}, Rational{7, 5}}) {
                auto &plan = prepared.green(grade, lambda);
                if (plan.uses_woodbury() != (k == 0))
                    throw std::runtime_error("source-derived Green factorization selection");
                for (std::size_t c = 0; c <= n; ++c) {
                    Vector x(n);
                    for (std::size_t i = 0; i < n; ++i)
                        x[i] = c == n ? Rational{static_cast<long>(i + 1), 3} : Rational{i == c};
                    auto canonical_input = x;
                    for (auto &value : canonical_input)
                        value.canonicalize();
                    const auto reference = field::exact_green_apply(prepared.study().hodge(grade),
                                                                    lambda, canonical_input);
                    const auto actual = plan.apply(x);
                    if (!reference.exact || actual != reference.value ||
                        plan.matrix().apply(canonical_input) != actual)
                        throw std::runtime_error("complete Green coordinates differ from native action");
                    ++checks;
                }
                if (&prepared.green(grade, lambda) != &plan)
                    throw std::runtime_error("Green preparation was not reused");
            }
        }
        const auto before = prepared.statistics();
        try { (void)prepared.green(core::Grade{0}, Rational{0});
              throw std::runtime_error("zero Green shift accepted"); }
        catch (const std::invalid_argument &) {}
        if (prepared.statistics() != before)
            throw std::runtime_error("failed Green preparation polluted the session cache");
        if (source->grade(core::Grade{1}).incidence_count() != 6 ||
            prepared.sectors(core::Grade{0}).apply(Vector{1, 2, 3, 4, 5, 6})[0].size() != 6)
            throw std::runtime_error("original source occurrences or coordinates lost");
        bool refused = false;
        try { (void)prepared.green(core::Grade{0}, Rational{1}).apply(Vector{1}); }
        catch (const std::invalid_argument &) { refused = true; }
        if (!refused)
            throw std::runtime_error("foreign Green input dimension accepted");
        core::RelationalComplexBuilder empty_builder(signature);
        empty_builder.declare_grade(core::Grade{2});
        auto empty_owner = std::make_shared<const NativeStudy>(empty_builder.freeze(),
            std::vector<field::ExactMetric>(3, field::ExactMetric::identity(0)));
        scientific_detail::PreparedTemporalStudy empty(empty_owner);
        for (std::uint16_t k = 0; k < 3; ++k) {
            auto &plan = empty.green(core::Grade{k}, Rational{1});
            if (!plan.apply({}).empty() || plan.matrix().rows() || plan.matrix().cols())
                throw std::runtime_error("declared empty grade acquired phantom coordinates");
        }
        refused = false;
        try { scientific_detail::PreparedTemporalStudy invalid(std::shared_ptr<const NativeStudy>{}); }
        catch (const std::invalid_argument &) { refused = true; }
        if (!refused)
            throw std::runtime_error("null prepared source owner accepted");
        std::cout << "PASS: " << checks << " full rational Green comparisons, weighted metrics, "
                     "original inactive/cancelling occurrences and owning preparation reuse\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
