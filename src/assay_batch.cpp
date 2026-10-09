#include "genomic_detail.hpp"
#include <rfc/runtime/work_graph.hpp>

namespace gadget::assay_detail {
Json query_batch(const AssayWorkspace &workspace, const Json &request,
                 const std::shared_ptr<CpuResources> &resources) {
    members(request, {"schema", "operation", "queries", "concurrent_jobs"});
    (void)resources->readmit();
    const auto concurrency = count(request.at("concurrent_jobs"),
                                   resources->coordinator().capacity(), "concurrent source queries");
    const auto declarations = request.at("queries");
    require(concurrency > 0 && json_object_is_type(declarations.raw(), json_type_array) &&
                declarations.size() > 0,
            "source query batch requires queries and positive admitted concurrency");
    // The native graph retains source-order indices; complete results live
    // under one private owner until all jobs join and publication is assembled.
    const auto stage = private_directory(fs::temp_directory_path(), "gadget-query-batch-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() { std::error_code error; fs::remove_all(path, error); }
    } cleanup{stage};
    std::set<std::string> ids;
    std::vector<rfc::runtime::WorkJob> jobs;
    const auto buffer = resources->coordinator().policy().declared_buffer_byte_limit /
                        std::min(concurrency, declarations.size());
    require(buffer >= sizeof(std::size_t), "native query receipt does not fit resource admission");
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        const auto declaration = declarations.at(i);
        members(declaration, {"id", "query"});
        require(ids.insert(text(declaration.at("id"), 256)).second,
                "duplicate source query job identity");
        const auto query = Json::parse(declaration.at("query").dump());
        require(query.at("operation").text() != "query_batch",
                "nested source query batches are not admitted");
        const auto worker_request =
            query.at("operation").text() == "native_query" ? query.at("query") : query;
        const auto workers = worker_request.has("workers")
                                 ? count(worker_request.at("workers"),
                                         resources->coordinator().capacity(), "query workers")
                                 : 1;
        require(workers > 0 && workers <= resources->coordinator().capacity(),
                "source query worker demand");
        jobs.push_back(rfc::runtime::WorkJob::make<std::size_t>(
            {rfc::machine::checked_size_add(i, 1), 1}, {workers, buffer, sizeof(std::size_t)}, {},
            [&workspace, resources, query, stage, i](const auto &, auto &assigned, std::stop_token stop) {
                if (stop.stop_requested())
                    throw rfc::runtime::ExecutionCancelled{};
                const auto answer = workspace.execute(query, resources->assigned(assigned));
                write_text(stage / std::to_string(i), answer.dump());
                if (stop.stop_requested())
                    throw rfc::runtime::ExecutionCancelled{};
                return i;
            }));
    }
    rfc::runtime::WorkGraphLimits limits;
    limits.jobs = declarations.size();
    limits.edges = 0;
    limits.depth = 1;
    limits.in_flight = concurrency;
    limits.retained_result_bytes = rfc::machine::checked_size_mul(declarations.size(), sizeof(std::size_t));
    const rfc::runtime::WorkGraph graph(std::move(jobs), limits);
    auto work = graph.submit(resources->coordinator());
    work.wait();
    auto answers = Json::array();
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        const auto completion = work.completion(rfc::runtime::WorkKey{i + 1, 1});
        require(*completion.result<std::size_t>() == i, "native result receipt source-order mismatch");
        auto entry = Json::object();
        entry.set("id", declarations.at(i).at("id"));
        entry.set("answer", Json::parse(read_text(stage / std::to_string(i))));
        const auto receipt = completion.receipt();
        entry.set("admitted_workers", cardinality(receipt.demand.workers));
        answers.push(entry);
    }
    auto out = Json::object();
    out.set("results", answers);
    out.set("concurrent_jobs", cardinality(concurrency));
    out.set("dispatcher", Json::string("rfc::runtime::WorkGraph"));
    out.set("publication_order",
            Json::string("declared job order after all owning native jobs join"));
    return out;
}
} // namespace gadget::assay_detail
