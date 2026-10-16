#include <gadget/execution.hpp>
#include <gadget/snapshot.hpp>
#include <iostream>

namespace {
using namespace gadget;
void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
template<class E, class F> void rejects(F action, const char *message) {
    bool refused = false;
    try { action(); } catch (const E &) { refused = true; }
    check(refused, message);
}
std::shared_ptr<CpuResources> resources(std::size_t bytes) {
    auto hardware = rfc::runtime::detect_hardware();
    rfc::runtime::CpuResourcePolicy policy;
    policy.worker_limit = 1; policy.declared_buffer_byte_limit = bytes;
    return std::make_shared<CpuResources>(hardware, policy,
        [] { return rfc::runtime::detect_hardware(); });
}
}
int main(int argc, char **argv) {
    try {
        if (argc != 3) throw std::invalid_argument("original fixture and destination required");
        fs::create_directories(argv[2]);
        const auto root = private_directory(argv[2], "queries-");
        fs::copy(argv[1], root / "input", fs::copy_options::recursive);
        auto workspace = load_pager(root / "input");
        import_measurements(workspace, root / "input/measurements.json");
        save_snapshot(workspace, root / "snapshot", true);
        const auto source = std::make_shared<const Workspace>(load_snapshot(root / "snapshot"));
        const auto query = Json::parse(
            R"({"operation":"green","context":"AD","realization":"support","grade":0,"field":"AD.G","lambda":"1","mode":"numerical","top_k":3})");
        const auto tiny = resources(1);
        {
            QuerySession limited(source, tiny);
            rejects<rfc::runtime::ExecutionLimitExceeded>([&] { (void)limited.execute(query); },
                "legacy query bypassed chosen memory owner");
            check(tiny->coordinator().available_workers() == 1 &&
                      tiny->coordinator().available_declared_buffer_bytes() == 1,
                  "failed legacy preparation leaked its chosen resources");
        }
        constexpr std::size_t budget = 4096;
        const auto cpu = resources(budget);
        const auto gram_bytes = source->genes.size() * source->genes.size() * sizeof(double);
        {
            QuerySession session(source, cpu);
            const auto first = session.execute(query);
            check(first.at("result").dump() == execute_query(*source, query).at("result").dump(),
                  "chosen resource owner changed original numerical result");
            check(cpu->coordinator().available_workers() == 1 &&
                      cpu->coordinator().available_declared_buffer_bytes() == budget - gram_bytes,
                  "cached original Gram does not retain the chosen byte owner");
            check(session.execute(query).at("result").dump() == first.at("result").dump() &&
                      cpu->coordinator().available_declared_buffer_bytes() == budget - gram_bytes &&
                      session.preparation_count() == 1,
                  "legacy repeat changed result, preparation or memory charge");
            const auto moment = Json::parse(
                R"({"operation":"typed_moment","context":"AD","realization":"support","field":"AD.G","other_field":"AD.E","order":1})");
            check(session.execute(moment).at("result").dump() ==
                      execute_query(*source, moment).at("result").dump(),
                  "single caller exact operation changed its source bound result");
            {
                auto occupied = cpu->coordinator().reserve(1);
                rejects<rfc::runtime::CpuResourceBusy>([&] { (void)session.execute(query); },
                    "warm numerical query bypassed occupied caller");
                rejects<rfc::runtime::CpuResourceBusy>([&] { (void)session.execute(moment); },
                    "exact query bypassed occupied caller");
            }
        }
        check(cpu->coordinator().available_workers() == 1 &&
                  cpu->coordinator().available_declared_buffer_bytes() == budget,
              "destroyed legacy session retains its chosen resource charge");
        std::cout << "PASS: original legacy query results, explicit resource owner, single caller, "
                     "dense refusal, warm cache and final release\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
