#include <gadget/compute.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <iostream>
#include <limits>
#include <thread>

namespace {
using namespace gadget;
namespace native = rfc;
void check(bool ok, const char *message) {
    if (!ok) throw std::logic_error(message);
}
template<class Action> void rejects(Action action) {
    bool refused = false;
    try { action(); } catch (const std::exception &) { refused = true; }
    check(refused, "invalid packed operation admitted");
}
template<class Action> void cancelled(Action action, bool application_exception = false) {
    bool stopped = false;
    try { action(); }
    catch (const native::runtime::ExecutionCancelled &) { stopped = true; }
    catch (const std::runtime_error &error) {
        if (application_exception || std::string_view(error.what()) != "packed execution cancelled")
            throw;
        stopped = true;
    }
    check(stopped, "host stop did not produce the actual cancellation refusal");
}
native::backend::PreparedPackedBinary prepared(std::size_t slots, bool long_word = false) {
    native::composite::AspectSchema schema;
    schema.add("first"); schema.add("second");
    native::composite::AspectPlanes planes(schema, slots);
    for (std::size_t i = 0; i < slots; ++i) {
        planes.set(0, i, i % 3 == 0);
        planes.set(1, i, i % 7 < 2);
    }
    using Op = native::composite::BinaryOpcode;
    native::composite::BinaryProgram program(7);
    program.append({Op::LoadAspect, 0, 0, 0, 0, 0});
    program.append({Op::LoadAspect, 1, 0, 0, 0, 1});
    program.append({Op::And, 2, 0, 1});
    program.append({Op::Or, 3, 0, 1});
    program.append({Op::Xor, 4, 0, 1});
    program.append({Op::Not, 5, 0});
    program.append({Op::Select, 6, 0, 4, 5});
    if (long_word)
        for (std::size_t i = 0; i < 300000; ++i) program.append({Op::Not, 3, 3});
    return native::backend::PreparedPackedBinary(program, planes, {2, 3, 4, 5, 6, 2});
}
Json query(std::size_t workers, std::size_t tile) {
    auto request = Json::object();
    request.set("backend", Json::string("cpu_tiles"));
    request.set("workers", Json::integer(static_cast<std::int64_t>(workers)));
    request.set("tile_words", Json::integer(static_cast<std::int64_t>(tile)));
    return request;
}
std::shared_ptr<CpuResources> resources(const native::runtime::HardwareInventory &hardware,
                                      std::size_t workers, std::size_t bytes,
                                      std::size_t tiles = std::numeric_limits<std::size_t>::max()) {
    native::runtime::CpuResourcePolicy policy;
    policy.worker_limit = workers; policy.declared_buffer_byte_limit = bytes; policy.tile_limit = tiles;
    return std::make_shared<CpuResources>(hardware, policy, [hardware] { return hardware; });
}
void released(CpuResources &resources) {
    const auto &coordinator = resources.coordinator();
    check(coordinator.available_workers() == coordinator.capacity() &&
              coordinator.available_declared_buffer_bytes() == coordinator.policy().declared_buffer_byte_limit,
          "packed completion or refusal retained a cooperating charge");
}
void oracle(const PackedResult &result, std::size_t slots) {
    check(result.outputs.size() == 6, "packed output occurrences lost");
    for (const auto &output : result.outputs) check(output.size() == slots, "packed output domain changed");
    for (std::size_t i = 0; i < slots; ++i) {
        const bool a = i % 3 == 0, b = i % 7 < 2;
        const bool expected[]{a && b, a || b, a != b, !a, a ? (a != b) : !a, a && b};
        for (std::size_t output = 0; output < 6; ++output)
            check(result.outputs[output].get(i) == expected[output], "packed bit differs from independent oracle");
    }
    for (const auto &output : result.outputs)
        if (slots % 64)
            check((output.words().back() >> (slots % 64)) == 0, "packed tail padding exposed an extra coordinate");
    check(result.execution.at("exact_cpu_reference_passed").boolean(), "canonical verification absent");
    const auto report = result.execution.at("cpu_tiles");
    check(report.at("completed_tiles").integer() == report.at("tiles").integer(),
          "packed tiles did not cover the complete source range");
}
} // namespace

int main() {
    try {
        const auto hardware = native::runtime::detect_hardware();
        const auto workers = std::min(std::size_t{8}, hardware.allowed_cpus.size());
        check(workers > 0, "actual CPU affinity unavailable");
        auto host = resources(hardware, workers, 64u << 20);
        std::size_t coordinates = 0;
        for (const std::size_t slots : {0u, 1u, 63u, 64u, 65u, 1025u, 1280193u}) {
            const auto source = prepared(slots);
            const std::array original_planes{source.planes().plane(0).words(), source.planes().plane(1).words()};
            const auto original_program = source.program().instructions();
            PackedOperation operation(source);
            const auto retained_owner = operation.prepared();
            check(retained_owner.same_owner(source), "accessor replaced the actual prepared owner");
            check(retained_owner.program().instructions() == original_program &&
                  retained_owner.outputs() == source.outputs() &&
                  retained_owner.planes().plane(0).words() == original_planes[0] &&
                  retained_owner.planes().plane(1).words() == original_planes[1],
                  "accessor changed ordered program, planes or output occurrences");
            const auto serial = operation.execute(*host, query(1, 1));
            oracle(serial, slots); released(*host);
            const auto parallel = operation.execute(*host, query(workers, 1));
            oracle(parallel, slots); released(*host);
            const auto survivor = [&] {
                PackedOperation temporary(source);
                return temporary.prepared();
            }();
            check(survivor.same_owner(source), "prepared owner expired with operation");
            auto complete = survivor.empty_result(std::numeric_limits<std::size_t>::max());
            survivor.cpu_words(0, survivor.words(), complete);
            for (std::size_t output = 0; output < complete.size(); ++output)
                check(complete[output].words() == serial.outputs[output].words(),
                      "surviving owner changed a complete output occurrence");
            for (std::size_t output = 0; output < serial.outputs.size(); ++output)
                check(serial.outputs[output].words() == parallel.outputs[output].words(),
                      "serial and parallel complete output differ");
            check(std::stoull(parallel.execution.at("declared_output_bytes").text()) ==
                      source.words() * source.outputs().size() * sizeof(std::uint64_t),
                  "packed output occurrence buffers undercounted");
            check(source.planes().plane(0).words() == original_planes[0] &&
                      source.planes().plane(1).words() == original_planes[1] &&
                      source.program().instructions() == original_program,
                  "complete prepared source planes or ordered program were modified");
            if (slots == 1280193) {
                check(parallel.execution.at("cpu_tiles").at("tiles").integer() > 20000,
                      "source did not cross the former tile count");
                const auto broad = operation.execute(*host, query(workers, 20001));
                oracle(broad, slots); released(*host);
                const auto entire = operation.execute(*host, query(workers, std::numeric_limits<std::int64_t>::max()));
                oracle(entire, slots); released(*host);
            }
            coordinates += slots * source.outputs().size();
        }
        const auto source = prepared(4097);
        const auto bytes = source.words() * source.outputs().size() * sizeof(std::uint64_t);
        auto short_budget = resources(hardware, workers, 2 * source.words() * sizeof(std::uint64_t));
        PackedOperation operation(source);
        bool budget = false;
        try { (void)operation.execute(*short_budget, query(1, 1)); }
        catch (const native::runtime::ExecutionLimitExceeded &error) {
            budget = error.resource() == "declared buffer bytes" && error.required() == 2 * bytes;
        }
        check(budget, "full output and verification budget not refused by actual coordinator"); released(*short_budget);
        auto exact_budget = resources(hardware, workers, 2 * bytes);
        oracle(operation.execute(*exact_budget, query(1, 1)), 4097); released(*exact_budget);
        auto tile_budget = resources(hardware, workers, 2 * bytes, 17);
        rejects([&] { (void)operation.execute(*tile_budget, query(1, 1)); }); released(*tile_budget);
        oracle(operation.execute(*tile_budget, query(1, 5)), 4097); released(*tile_budget);
        PackedOperation allocation(source, {0, false, bytes - 1});
        rejects([&] { (void)allocation.execute(*host, query(1, 1)); }); released(*host);
        for (const auto &request : {query(0, 1), query(1, 0), query(host->coordinator().capacity() + 1, 1)})
            rejects([&] { (void)operation.execute(*host, request); });
        released(*host);
        {
            const auto active = host->coordinator().reserve(host->coordinator().capacity());
            rejects([&] { (void)operation.execute(*host, query(1, 1)); });
        }
        released(*host);
        {
            std::array<std::optional<PackedResult>, 2> results;
            std::array<std::exception_ptr, 2> errors;
            std::barrier start(2);
            std::array<std::jthread, 2> calls;
            for (std::size_t i = 0; i < calls.size(); ++i)
                calls[i] = std::jthread([&, i] {
                    start.arrive_and_wait();
                    try { results[i] = operation.execute(*host, query(1, 3)); }
                    catch (...) { errors[i] = std::current_exception(); }
                });
            for (auto &call : calls) call.join();
            for (std::size_t i = 0; i < calls.size(); ++i) {
                if (errors[i]) {
                    try { std::rethrow_exception(errors[i]); }
                    catch (const native::runtime::CpuResourceBusy &) {
                        results[i] = operation.execute(*host, query(1, 3));
                    }
                }
                oracle(results[i].value(), 4097);
            }
            for (std::size_t i = 0; i < source.outputs().size(); ++i)
                check(results[0]->outputs[i].words() == results[1]->outputs[i].words(),
                      "independent concurrent output owners differ");
            released(*host);
        }
        std::stop_source stopped; stopped.request_stop();
        cancelled([&] { (void)operation.execute(*host, query(1, 1), stopped.get_token()); }, true);
        released(*host);
        {
            native::runtime::CpuResourcePolicy policy;
            policy.worker_limit = workers; policy.declared_buffer_byte_limit = 2 * bytes;
            policy.tile_limit = std::numeric_limits<std::size_t>::max();
            std::size_t observations{};
            std::stop_source admission_stop;
            CpuResources admission(hardware, policy, [&] {
                ++observations;
                admission_stop.request_stop();
                return hardware;
            });
            cancelled([&] { (void)operation.execute(admission, query(1, 1), stopped.get_token()); }, true);
            check(observations == 0, "pre-stopped work reached readmission");
            cancelled([&] { (void)operation.execute(admission, query(1, 1), admission_stop.get_token()); }, true);
            check(observations == 1, "stop during readmission was not exercised");
            released(admission);
        }
        {
            auto &parent = host->coordinator();
            const auto available_workers = parent.available_workers();
            const auto available_bytes = parent.available_declared_buffer_bytes();
            {
                auto delegated = native::runtime::CpuResourceCoordinator::delegate(parent.reserve(1, 2 * bytes));
                auto child = host->assigned(delegated);
                cancelled([&] { (void)operation.execute(*child, query(1, 1), stopped.get_token()); }, true);
                released(*child);
                check(parent.available_workers() == available_workers - 1 &&
                          parent.available_declared_buffer_bytes() == available_bytes - 2 * bytes,
                      "cancelled child released its retained parent assignment early");
                oracle(operation.execute(*child, query(1, 3)), 4097);
                released(*child);
                child.reset();
                check(parent.available_workers() == available_workers - 1,
                      "live delegated authority lost its parent charge");
            }
            released(*host);
        }

        const auto slow_source = prepared(65, true);
        PackedOperation slow(slow_source);
        std::stop_source stop;
        std::barrier ready(2);
        std::atomic<bool> observed{};
        auto &coordinator = host->coordinator();
        std::jthread cancel([&](std::stop_token wait_stop) {
            ready.arrive_and_wait();
            while (!wait_stop.stop_requested()) {
                const auto epoch = coordinator.change_epoch();
                if (coordinator.available_workers() < coordinator.capacity()) {
                    observed = true; stop.request_stop(); return;
                }
                coordinator.wait_for_change(epoch, wait_stop);
            }
        });
        ready.arrive_and_wait();
        cancelled([&] { (void)slow.execute(*host, query(1, 1), stop.get_token()); });
        cancel.request_stop(); coordinator.notify_waiters(); cancel.join();
        check(observed && stop.stop_requested(), "active packed CPU reservation/cancellation not exercised"); released(*host);
        oracle(slow.execute(*host, query(1, 1)), 65); released(*host);

        // Authored admission control, separately from the measured physical host.
        auto authored = hardware;
        authored.platform = "authored worker policy fixture";
        authored.allowed_cpus.clear(); authored.cgroups.clear(); authored.caches.clear();
        for (std::size_t i = 0; i < 48; ++i) authored.allowed_cpus.push_back({i, {}, {}, {}});
        auto wide = resources(authored, 48, 64u << 20);
        oracle(operation.execute(*wide, query(48, 1)), 4097); released(*wide);
        std::cout << "PASS: " << coordinates << " complete coordinate occurrences, physical host "
                  << hardware.allowed_cpus.size() << " allowed CPUs, declared 48-worker control, "
                     "full buffers, old tile boundaries, cancellation and retry\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
