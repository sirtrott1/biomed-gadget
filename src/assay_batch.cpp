#include <gadget/build.hpp>
#include "genomic_detail.hpp"
#include <rfc/runtime/work_graph.hpp>

namespace gadget::assay_detail {
NativeRecord query_batch_record(const AssayWorkspace &workspace, const NativeRecord &request,
                               const std::shared_ptr<CpuResources> &resources) {
    using N = NativeRecordValue;
    const auto query_controls = render_record_view(request);
    members(query_controls, {"schema", "operation", "queries", "concurrent_jobs"});
    require(query_controls.at("schema").text() == "gadget.assays.query.v1" &&
                query_controls.at("operation").text() == "query_batch", "source query batch schema");
    (void)resources->readmit();
    const auto concurrency = count(query_controls.at("concurrent_jobs"),
                                   resources->coordinator().capacity(), "concurrent source queries");
    const auto declarations = query_controls.at("queries");
    require(concurrency > 0 && json_object_is_type(declarations.raw(), json_type_array) &&
                declarations.size() > 0,
            "source query batch requires queries and positive admitted concurrency");
    // Complete results retain their private archive owners until all jobs join.
    const auto stage = private_directory(fs::temp_directory_path(), "gadget-query-batch-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() { std::error_code error; fs::remove_all(path, error); }
    } cleanup{stage};
    struct Completion {
        std::size_t index;
        RawByteSource archive;
    };
    std::set<std::string> ids;
    std::map<std::string, std::size_t> summaries;
    std::vector<std::size_t> executions;
    std::vector<rfc::runtime::WorkJob> jobs;
    jobs.reserve(declarations.size());
    executions.reserve(declarations.size());
    const auto buffer = resources->coordinator().policy().declared_buffer_byte_limit /
                        std::min(concurrency, declarations.size());
    require(buffer >= sizeof(Completion), "native query receipt does not fit resource admission");
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        const auto declaration = declarations.at(i);
        members(declaration, {"id", "query"});
        require(ids.insert(text(declaration.at("id"), 256)).second,
                "duplicate source query job identity");
        const auto query = request.at("queries").at(i).at("query");
        const auto query_view = declaration.at("query");
        require(query_view.at("operation").text() != "query_batch",
                "nested source query batches are not admitted");
        require(AssayWorkspace::supports_native_publication(query_view.at("operation").text()),
                "source query child has no native publisher");
        const auto worker_request =
            query_view.at("operation").text() == "native_query" ? query_view.at("query") : query_view;
        const auto workers = worker_request.has("workers")
                                 ? count(worker_request.at("workers"),
                                         resources->coordinator().capacity(), "query workers")
                                 : 1;
        require(workers > 0 && workers <= resources->coordinator().capacity(),
                "source query worker demand");
        // Field summaries read the same immutable workspace. Retain every job
        // occurrence while sharing identical complete requests and publications.
        const bool summary = query_view.at("operation").text() == "native_query" &&
            worker_request.has("operation") && worker_request.at("operation").text() == "field_summary";
        const auto execution = jobs.size();
        if (summary) {
            const auto [previous, inserted] = summaries.emplace(query_view.dump(), execution);
            if (!inserted) {
                executions.push_back(previous->second);
                continue;
            }
        }
        executions.push_back(execution);
        jobs.push_back(rfc::runtime::WorkJob::make<Completion>(
            {rfc::machine::checked_size_add(execution, 1), 1}, {workers, buffer, sizeof(Completion)}, {},
            [&workspace, resources, query, stage, execution](const auto &, auto &assigned, std::stop_token stop) {
                if (stop.stop_requested())
                    throw rfc::runtime::ExecutionCancelled{};
                const auto answer = workspace.execute_record(query, resources->assigned(assigned));
                auto archive = answer.archive_source(1, stage / std::to_string(execution));
                if (stop.stop_requested())
                    throw rfc::runtime::ExecutionCancelled{};
                return Completion{execution, std::move(archive)};
            }));
    }
    rfc::runtime::WorkGraphLimits limits;
    const auto execution_count = jobs.size();
    limits.jobs = execution_count;
    limits.edges = 0;
    limits.depth = 1;
    limits.in_flight = concurrency;
    limits.retained_result_bytes = rfc::machine::checked_size_mul(execution_count, sizeof(Completion));
    const rfc::runtime::WorkGraph graph(std::move(jobs), limits);
    auto work = graph.submit(resources->coordinator());
    work.wait();
    const auto &source = workspace.native_publication();
    std::uint64_t key = 1;
    for (const auto &raw : source.raw_sources()) {
        require(raw.key() < std::numeric_limits<std::uint64_t>::max(),
                "batch raw source identity space exhausted");
        key = std::max(key, raw.key() + 1);
    }
    std::vector<RawByteSource> archives;
    archives.reserve(rfc::machine::checked_size_add(execution_count, 1));
    struct Publication { std::uint64_t key; std::size_t bytes; std::size_t workers; };
    std::vector<Publication> publications;
    publications.reserve(execution_count);
    N::Array answers;
    answers.elements.reserve(declarations.size());
    for (std::size_t i = 0; i < execution_count; ++i) {
        const auto completion = work.completion(rfc::runtime::WorkKey{i + 1, 1});
        const auto retained = completion.result<Completion>();
        require(retained->index == i, "native result receipt source-order mismatch");
        require(key < std::numeric_limits<std::uint64_t>::max(),
                "batch raw source identity space exhausted");
        auto archive = retained->archive.rekey(key);
        const auto size = archive.bytes().size();
        const auto receipt = completion.receipt();
        publications.push_back({key, size, receipt.demand.workers});
        archives.push_back(std::move(archive));
        ++key;
    }
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        const auto &publication = publications.at(executions.at(i));
        answers.elements.emplace_back(N::Object{{
            {"id", N(declarations.at(i).at("id").text())},
            {"child_publication", N(N::ByteSpan{publication.key, 0, publication.bytes, false})},
            {"admitted_workers", N(rfc::numeric::Integer(std::to_string(publication.workers)))}}});
    }
    auto workspace_archive = source.archive_source(key);
    const auto workspace_size = workspace_archive.bytes().size();
    archives.push_back(std::move(workspace_archive));
    N fields(N::Object{{
        {"schema", N(std::string("gadget.assay.batch.fields.v1"))},
        {"results", N(std::move(answers))},
        {"executions", N(rfc::numeric::Integer(std::to_string(execution_count)))},
        {"concurrent_jobs", N(rfc::numeric::Integer(std::to_string(concurrency)))},
        {"dispatcher", N(std::string("rfc::runtime::WorkGraph"))},
        {"publication_order", N(std::string("declared job order after all owning native jobs join"))}}});
    N certificate(N::Object{{
        {"schema", N(std::string("gadget.assay.result.native.v1"))},
        {"query", import_record_view(query_controls)}, {"source_digest", N(workspace.digest())},
        {"result", std::move(fields)},
        {"workspace_publication", N(N::ByteSpan{key, 0, workspace_size, false})},
        {"arithmetic", N(std::string("rational"))},
        {"gadget_sha256", N(std::string(gadget::build::source_sha256()))},
        {"rexgraph_sha256", N(std::string(gadget::build::rexgraph_sha256()))}}});
    auto history = source.publication_history();
    require(history.has_value(), "batch source has no immutable history");
    // Append to the admitted history so temporary declaration cargo is released
    // before the canonical frame is decoded and replayed.
    const auto prepared = NativeRecord::annotate(N(N::Object{}), std::move(*history));
    return retain_operation(prepared, std::move(archives), request, std::move(certificate));
}
Json query_batch_fields_view(const NativeRecord &record) {
    require(record.at("schema").text() == "gadget.assay.batch.fields.v1",
            "native source query batch fields schema");
    auto answers = Json::array();
    const auto entries = record.at("results");
    std::map<std::uint64_t, Json> publications;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto entry = entries.at(i);
        const auto archive = entry.at("child_publication");
        const auto key = archive.byte_source().key();
        auto found = publications.find(key);
        if (found == publications.end()) {
            const auto child = NativeRecord::decode(archive.bytes());
            found = publications.emplace(key, render_record_view(child.at("result"))).first;
        }
        auto answer = Json::object();
        answer.set("id", render_record_view(entry.at("id")));
        answer.set("answer", found->second);
        answer.set("admitted_workers", render_record_view(entry.at("admitted_workers")));
        answers.push(answer);
    }
    auto out = Json::object();
    out.set("results", answers);
    for (const auto *name : {"concurrent_jobs", "dispatcher", "publication_order"})
        out.set(name, render_record_view(record.at(name)));
    if (record.has("executions"))
        out.set("executions", render_record_view(record.at("executions")));
    return out;
}
Json query_batch(const AssayWorkspace &workspace, const Json &request,
                 const std::shared_ptr<CpuResources> &resources) {
    const auto native = query_batch_record(workspace, admit_record_view(request), resources);
    return query_batch_fields_view(native.at("result").at("result"));
}
} // namespace gadget::assay_detail
