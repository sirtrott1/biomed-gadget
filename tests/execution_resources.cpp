#include <gadget/execution.hpp>
#include <bit>
#include <iostream>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>

namespace {
using namespace gadget;
using namespace rfc;
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
template <class E, class F> void rejects(F action, const char *message) {
    bool rejected = false;
    try { action(); } catch (const E &) { rejected = true; }
    check(rejected, message);
}
Workspace source(std::size_t ad, std::size_t hd, bool zero = false) {
    Workspace value;
    value.source_digest = "authored original resource source";
    value.genes = {"first", "second"};
    core::RelationSignature signature;
    const auto type = signature.add_type("original binary assertion");
    core::RelationalComplexBuilder builder(signature);
    builder.declare_witness_identity(core::IdentityKey{1}); builder.declare_witness_identity(core::IdentityKey{2});
    const auto count = std::max<std::size_t>({ad, hd, 1});
    for (std::size_t i = 0; i < count; ++i) {
        builder.add_grade_one(core::RelationKey{i + 1}, type,
            {{core::IdentityKey{1}, true, false, Rational{1}, {}, core::OccurrenceKey{2*i + 1}},
             {core::IdentityKey{2}, true, false, Rational{1}, {}, core::OccurrenceKey{2*i + 2}}});
        Assertion assertion;
        assertion.participants = {{0, zero ? Rational{0} : Rational{1}, "original", 2*i + 1, {}},
                                  {1, zero ? Rational{0} : Rational{2}, "original", 2*i + 2, {}}};
        value.assertions.push_back(std::move(assertion));
    }
    value.complex = std::make_shared<const core::RelationalComplex>(builder.freeze());
    for (const auto &[name, width] : std::vector<std::pair<std::string, std::size_t>>{{"AD", ad}, {"HD", hd}}) {
        auto mask = std::make_unique<algebra::SubcomplexMask>(*value.complex);
        for (std::size_t i = 0; i < width; ++i) mask->select({{1}, static_cast<std::uint32_t>(i)});
        mask->close_boundaries(*value.complex);
        value.contexts.emplace(name, std::move(mask));
    }
    return value;
}
struct Resources {
    std::shared_ptr<runtime::HardwareInventory> observed{std::make_shared<runtime::HardwareInventory>()};
    std::shared_ptr<CpuResources> cpu;
    explicit Resources(std::size_t bytes) {
        observed->platform = "authored observed resource control";
        observed->cpu_model = "two declared test CPUs";
        observed->affinity_observed = true;
        observed->allowed_cpus = {{0, 0, 0, {}}, {1, 0, 1, {}}};
        runtime::CpuResourcePolicy policy;
        policy.worker_limit = 2; policy.declared_buffer_byte_limit = bytes;
        cpu = std::make_shared<CpuResources>(*observed, policy, [o = observed] { return *o; });
    }
};
void small() {
    auto original = source(3, 2);
    constexpr std::size_t bytes = 4096;
    Resources resources(bytes);
    std::optional<Execution::GreenWitnesses> retained;
    {
        Execution ad(original, "AD", "rp", resources.cpu), hd(original, "HD", "rp", resources.cpu);
        const auto &gram = ad.relation_gram();
        check(resources.cpu->coordinator().available_declared_buffer_bytes() == bytes - 9*sizeof(double), "cached Gram has no retained charge");
        check(resources.cpu->coordinator().available_workers() == 2, "retained Gram occupies a CPU worker");
        const auto before = resources.cpu->coordinator().available_declared_buffer_bytes();
        check(&gram == &ad.relation_gram() && resources.cpu->coordinator().available_declared_buffer_bytes() == before, "Gram cache duplicated its allocation");
        const auto cross = cross_relation_gram(hd, ad);
        const auto exact_cross = field::boundary_matrix(hd.complex(), {1}).transpose().multiply(field::boundary_matrix(ad.complex(), {1}));
        check(cross.rows() == 2 && cross.cols() == 3, "rectangular Gram was coerced to square");
        for (std::size_t i = 0; i < 2; ++i) for (std::size_t j = 0; j < 3; ++j)
            check(cross(i,j) == exact_cross.get(i,j).get_d(), "cross Gram changed original carrier action");
        check(resources.cpu->coordinator().available_declared_buffer_bytes() == before, "returned cross Gram claims a retained reservation");
        const Vector input{Rational{1}, Rational{"2/7"}};
        Json diagnostic;
        const auto response = ad.green(input, Rational{"1/3"}, 1e-9, diagnostic, &retained);
        const auto expected_system = approximation::DenseMatrix::from_exact(ad.metrics()[1].materialize()).add(gram.scale(Rational{"1/3"}.get_d()));
        for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j)
            check(std::bit_cast<std::uint64_t>(expected_system(i,j)) == std::bit_cast<std::uint64_t>(retained->system(i,j)), "assembly changed binary64 scalar rounding");
        check(retained->solved == binary64_values(expected_system.solve(approximate_values(retained->pulled))), "canonical solve coefficients changed");
        check(retained->residual == numeric::sub(numeric::add(response, numeric::scale(Rational{"1/3"}, ad.hodge().apply(response))), input), "retained residual is a display copy");
        check(retained->relative <= retained->tolerance && retained->absolute == norm(retained->residual), "residual certificate does not reconstruct");
        check(resources.cpu->coordinator().available_workers() == 2 && resources.cpu->coordinator().available_declared_buffer_bytes() == bytes - 18*sizeof(double), "moved retained system lost or duplicated its charge");
        check(diagnostic.at("allocation").at("solve_peak_bytes").text() == std::to_string(18*sizeof(double) + 9*sizeof(double)), "solve buffer scope is inaccurate");
        {
            auto occupied = resources.cpu->coordinator().reserve(2);
            rejects<runtime::CpuResourceBusy>([&] { Json view; ad.green(input, Rational{1}, 1e-9, view); }, "busy Green oversubscribed the caller");
            rejects<runtime::CpuResourceBusy>([&] { (void)cross_relation_gram(hd, ad); }, "busy cross Gram oversubscribed the caller");
        }
        resources.observed->cpu_model = "changed CPU";
        rejects<std::invalid_argument>([&] { (void)ad.relation_gram(); }, "cached execution reused changed observed hardware");
        resources.observed->cpu_model = "two declared test CPUs";
    }
    check(resources.cpu->coordinator().available_declared_buffer_bytes() == bytes - 9*sizeof(double), "retained witness charge ended with Execution");
    retained.reset();
    check(resources.cpu->coordinator().available_declared_buffer_bytes() == bytes && resources.cpu->coordinator().available_workers() == 2, "resource charges survived their buffers");
    Resources limited(9*sizeof(double) - 1);
    {
        Execution ad(original, "AD", "rp", limited.cpu);
        rejects<runtime::ExecutionLimitExceeded>([&] { (void)ad.relation_gram(); }, "insufficient dense memory was admitted");
        check(limited.cpu->coordinator().available_declared_buffer_bytes() == 9*sizeof(double) - 1, "failed allocation leaked a reservation");
    }
    Resources solve_limited(9*sizeof(double) + 18*sizeof(double) + 9*sizeof(double) - 1);
    {
        Execution ad(original, "AD", "rp", solve_limited.cpu); Json view;
        (void)ad.relation_gram();
        rejects<runtime::CpuResourceBusy>([&] { ad.green({1,2}, Rational{1}, 1e-9, view); }, "solve admission omitted its copied dense workspace and live cache");
        check(solve_limited.cpu->coordinator().available_workers() == 2, "failed solve leaked a worker");
    }
    Resources cross_limited(6*sizeof(double) - 1);
    Execution cross_ad(original, "AD", "rp", cross_limited.cpu), cross_hd(original, "HD", "rp", cross_limited.cpu);
    rejects<runtime::ExecutionLimitExceeded>([&] { (void)cross_relation_gram(cross_hd, cross_ad); }, "rectangular buffer admission omitted one axis");
    check(cross_limited.cpu->coordinator().available_workers() == 2, "rectangular refusal leaked a worker");
    auto out_of_range = source(1, 1);
    out_of_range.assertions[0].participants[0].rp = exact_decimal("1e4096");
    Resources range_resources(1024);
    Execution range(out_of_range, "AD", "rp", range_resources.cpu);
    for (int attempt = 0; attempt < 2; ++attempt) {
        rejects<std::invalid_argument>([&] { (void)range.relation_gram(); }, "failed Gram preparation became a valid cache");
        check(range_resources.cpu->coordinator().available_workers() == 2 &&
                  range_resources.cpu->coordinator().available_declared_buffer_bytes() == 1024, "failed Gram preparation retained resources");
    }
    auto empty = source(0, 1);
    Resources empty_resources(1024);
    Execution zero(empty, "AD", "support", empty_resources.cpu), one(empty, "HD", "support", empty_resources.cpu);
    Json view;
    const auto response = zero.green({1,2}, Rational{1}, 1e-9, view, &retained);
    check(response == Vector({1,2}) && retained->system.rows() == 0 && retained->solved.empty() && retained->residual == Vector({0,0}), "empty relation Green did not preserve ambient source");
    check(cross_relation_gram(zero, one).rows() == 0 && cross_relation_gram(one, zero).cols() == 0, "empty rectangular maps lost actual shapes");
    retained.reset();
    auto unrelated = source(3,2); unrelated.source_digest += " independent";
    Execution other(unrelated, "HD", "support", empty_resources.cpu);
    rejects<std::invalid_argument>([&] { (void)cross_relation_gram(one, other); }, "equal dimensions identified independent source carriers");
}
void wide() {
    const std::size_t width = 2049;
    const auto bytes = machine::checked_size_mul(machine::checked_size_mul(width, width), sizeof(double));
    Resources resources(machine::checked_size_add(machine::checked_size_mul(bytes, 4), machine::checked_size_mul(width, 24)));
    auto original = source(width, 1, true);
    {
        Execution ad(original, "AD", "rp", resources.cpu), hd(original, "HD", "rp", resources.cpu);
        const auto &gram = ad.relation_gram();
        check(gram.rows() == width && gram.cols() == width, "legacy ceiling truncated the original relation carrier");
        for (std::size_t i = 0; i < width; ++i) check(gram(i,i) == 0 && ad.source_relations.at(i) == i, "zero columns or original relation identities were omitted");
        const auto cross = cross_relation_gram(ad, hd);
        check(cross.rows() == width && cross.cols() == 1, "wide rectangular Gram refused admitted source");
        Json view; std::optional<Execution::GreenWitnesses> witness;
        const auto response = ad.green({1,2}, Rational{1}, 1e-9, view, &witness);
        check(response == Vector({1,2}) && witness->solved == Vector(width, Rational{0}) && witness->system.rows() == width, "wide zero-supported Green lost its complete solve");
        check(resources.cpu->coordinator().available_workers() == 2, "wide result retained caller capacity");
    }
    check(resources.cpu->coordinator().available_declared_buffer_bytes() == resources.cpu->coordinator().policy().declared_buffer_byte_limit, "wide cached matrices leaked charge");
}
} // namespace
int main() {
    try { small(); wide(); std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n"; }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
