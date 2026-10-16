#include <charconv>
#include <chrono>
#include <gadget/io.hpp>
#include <gadget/study.hpp>
#include <iostream>
#include <rfc/core/builder.hpp>

int main(int argc, char **argv) {
    try {
        std::size_t n = 65536;
        if (argc > 1) {
            const std::string s = argv[1];
            const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), n);
            if (error != std::errc{} || end != s.data() + s.size() || n < 2 || n > UINT32_MAX)
                throw std::invalid_argument("coordinate count must fit the native index domain");
        }
        using namespace rfc;
        using gadget::Rational;
        const auto started = std::chrono::steady_clock::now();
        core::RelationSignature signature;
        const auto type = signature.add_type("declared.one_to_many_relation");
        core::RelationalComplexBuilder builder(signature);
        std::vector<core::BoundaryOccurrenceInput> boundary;
        boundary.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            boundary.emplace_back(core::IdentityKey{i + 1}, true, i != 0,
                                  i ? Rational(1, static_cast<unsigned long>(n - 1)) : Rational{1},
                                  std::vector<core::AspectValue>{}, core::OccurrenceKey{i + 1});
        builder.add_grade_one(core::RelationKey{1}, type, boundary);
        auto source = builder.freeze();
        gadget::NativeStudy study(
            std::move(source), {field::ExactMetric::identity(n), field::ExactMetric::identity(1)});
        numeric::RationalVector j(n);
        j[0] = 1;
        const auto prepared = std::chrono::steady_clock::now();
        const auto answer = study.green_action(core::Grade{0}, Rational{1}, j);
        const auto solved = std::chrono::steady_clock::now();
        const auto residual = numeric::sub(
            numeric::add(answer.value, study.hodge(core::Grade{0}).apply(answer.value)), j);
        // Independent rank one inverse: (I + b b*)^-1 j = j - b(b*j)/(1+b*b).
        const Rational tail(1, static_cast<unsigned long>(n - 1));
        const Rational q = 1 + tail;
        const Rational head = 1 - Rational(1) / (1 + q), other = tail / (1 + q);
        if (!answer.exact || !numeric::is_zero(residual))
            throw std::runtime_error("Green residual");
        for (std::size_t i = 0; i < n; ++i)
            if (answer.value[i] != (i ? other : head))
                throw std::runtime_error("independent rank-one identity");
        const auto verified = std::chrono::steady_clock::now();
        auto result = gadget::Json::object();
        result.set("schema", gadget::Json::string("gadget.factored-action.experiment.v1"));
        result.set("coordinates", gadget::Json::integer(n));
        result.set("retained_incidence_occurrences", gadget::Json::integer(n));
        result.set("iterations", gadget::Json::integer(answer.iterations));
        result.set("head", gadget::Json::string(head.get_str()));
        result.set("other", gadget::Json::string(other.get_str()));
        result.set(
            "sum",
            gadget::Json::string(
                numeric::dot(answer.value, numeric::RationalVector(n, Rational{1})).get_str()));
        result.set("complete_field_verified", gadget::Json::boolean(true));
        result.set(
            "output_contract",
            gadget::Json::string("Complete field computed and checked; reported symmetric "
                                 "coordinates and contraction describe this authored fixture"));
        result.set("dense_operator_constructed", gadget::Json::boolean(false));
        result.set("preparation_seconds",
                   gadget::Json::number(std::chrono::duration<double>(prepared - started).count()));
        result.set("green_seconds",
                   gadget::Json::number(std::chrono::duration<double>(solved - prepared).count()));
        result.set("verification_seconds",
                   gadget::Json::number(std::chrono::duration<double>(verified - solved).count()));
        std::cout << result.dump() << '\n';
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
