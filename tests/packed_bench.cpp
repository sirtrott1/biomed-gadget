#include <algorithm>
#include <gadget/compute.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        using namespace rfc::composite;
        const bool device = argc > 1;
        ComputeOptions options;
        options.allow_software_device = argc > 1 && std::string(argv[1]) == "software";
        // Authored nonclinical throughput fixture, never a biomedical source.
        constexpr std::size_t slots = (std::size_t{1} << 24) + 13;
        AspectSchema schema;
        for (auto i = 0; i < 16; ++i)
            schema.add("fixture-" + std::to_string(i));
        AspectPlanes planes(schema, slots);
        for (std::size_t a = 0; a < 16; ++a)
            for (std::size_t w = 0; w < planes.plane(static_cast<AspectId>(a)).word_count(); ++w)
                planes.plane(static_cast<AspectId>(a))
                    .assign_word(w, std::uint64_t(0x9e3779b97f4a7c15ULL * (w + 1 + a * 17)));
        BinaryProgram program(4);
        program.append({BinaryOpcode::LoadAspect, 0, 0, 0, 0, 0});
        for (std::size_t i = 0; i < 128; ++i) {
            program.append({BinaryOpcode::LoadAspect, 1, 0, 0, 0, static_cast<AspectId>(i % 16)});
            program.append({i % 3 ? BinaryOpcode::Xor : BinaryOpcode::Or, 0, 0, 1});
        }
        PackedOperation operation(rfc::backend::PreparedPackedBinary(program, planes, {0}),
                                  options);
        auto resources = CpuResources::local();
        auto records = Json::array();
        for (auto workers : {1u, 2u, 4u, 8u, 16u}) {
            if (workers > resources->coordinator().capacity())
                continue;
            auto q = Json::object();
            q.set("backend", Json::string("cpu_tiles"));
            q.set("workers", Json::integer(workers));
            q.set("tile_words", Json::integer(1024));
            operation.execute(*resources, q);
            auto samples = Json::array();
            std::vector<std::uint64_t> run_times;
            std::vector<double> complete;
            for (auto i = 0; i < 5; ++i) {
                auto r = operation.execute(*resources, q);
                run_times.push_back(std::stoull(r.execution.at("cpu_run_ns").text()));
                complete.push_back(
                    json_object_get_double(r.execution.at("complete_verified_call_seconds").raw()));
                samples.push(r.execution);
            }
            std::sort(run_times.begin(), run_times.end());
            std::sort(complete.begin(), complete.end());
            auto row = Json::object();
            row.set("workers", Json::integer(workers));
            row.set("median_cpu_run_ns", Json::string(std::to_string(run_times[2])));
            row.set("median_verified_call_seconds", Json::number(complete[2]));
            row.set("samples", samples);
            records.push(row);
        }
        if (device)
            for (const auto &backend : {"vulkan", "cpu_vulkan", "measured"}) {
                auto q = Json::object();
                q.set("backend", Json::string(backend));
                q.set("tile_words", Json::integer(1024));
                if (std::string(backend) == "cpu_vulkan")
                    q.set("workers", Json::integer(4));
                auto r = operation.execute(*resources, q);
                records.push(r.execution);
            }
        auto report = Json::object();
        report.set("schema", Json::string("gadget.packed.benchmark.v1"));
        report.set("scope", Json::string("authored nonclinical packed fixture; kernel timing and "
                                         "complete verified call reported separately"));
        report.set("slots", Json::integer(slots));
        report.set("instructions", Json::integer(program.instruction_count()));
        report.set("hardware", resources->report());
        report.set("samples", records);
        report.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        report.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
        std::cout << report.dump() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
