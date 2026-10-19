#include "../src/pipeline_result_binding.hpp"
#include "../src/result_answer_witness.hpp"
#include <functional>
#include <iostream>

using namespace gadget;
using namespace gadget::pipeline_detail;
namespace {
std::size_t checks = 0;
void check(bool ok) {
    if (!ok)
        throw std::runtime_error("result binding contract failed");
    ++checks;
}
void refuses(const std::function<void()> &f) {
    try {
        f();
    } catch (const std::exception &) {
        ++checks;
        return;
    }
    throw std::runtime_error("result binding should refuse");
}
Json clone(const Json &j) { return Json::parse(j.dump()); }
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string_view(argv[1]) == "--witness") {
            const auto input=Json::parse(read_text(argv[2]));const fs::path directory=input.at("directory").text();
            const auto reference=result_answer_reference(directory,input.at("result"),input.at("request"),
                input.at("method").text(),input.at("analysis_sha256").text());
            const auto view=result_answer_view(directory,reference);
            if(input.at("method").text()=="pathway_overrepresentation")
                check(view.at("result").at("tests").at(0).at("p_value").text()=="1/3");
            else check(view.at("result").at("features").size()==3);
            auto wrong=Json::parse(reference.dump());
            wrong.at("native_selection").set("analysis_sha256",Json::string("sha256:"+std::string(64,'0')));
            refuses([&]{result_answer_view(directory,wrong);});
            std::cout<<view.dump()<<'\n';return 0;
        }
        if (argc != 2)
            throw std::invalid_argument("test artifact directory required");
        const fs::path root = argv[1];
        fs::create_directories(root);
        const std::string image("native\0archive\xff", 15);
        const auto file = root / "original.rfcxs";
        write_text(file, image);
        const auto hash = "sha256:" + sha256(image);
        check(verified_record_image(file, hash, image.size()) == image);
        refuses([&] { verified_record_image(file, hash, image.size() - 1); });
        refuses(
            [&] { verified_record_image(file, "sha256:" + std::string(64, '0'), image.size()); });
        refuses([&] { verified_record_image(file, hash, 0); });
        auto request = Json::parse(
            R"({"schema":"gadget.result-operation.v1","operation":"pathway_overrepresentation","analysis_binding":"actual_result","analysis_bytes":1048576,"source_record":{"path":"x","sha256":"x","max_bytes":1048576},"query":{}})");
        request.at("source_record").set("path", Json::string(file.string()));
        request.at("source_record").set("sha256", Json::string(hash));
        check(admit_result_operation(request).source_sha256 == hash);
        for (const auto *method : {"pathway_overrepresentation", "feature_evidence_context"}) {
            auto q = clone(request);
            q.set("operation", Json::string(method));
            check(admit_result_operation(q).operation == method);
        }
        for (const auto *member : {"analysis_bytes", "max_bytes"}) {
            auto q = clone(request);
            (std::string(member) == "analysis_bytes" ? q : q.at("source_record"))
                .set(member, Json::integer(0));
            refuses([&] { admit_result_operation(q); });
        }
        {
            auto q = clone(request);
            q.set("operation", Json::string("automatic_efficacy"));
            refuses([&] { admit_result_operation(q); });
        }
        {
            auto q = clone(request);
            q.set("threshold", Json::string("0.05"));
            refuses([&] { admit_result_operation(q); });
        }
        {
            auto q = clone(request);
            q.at("source_record").set("sha256", Json::string("SHA256:" + sha256(image)));
            refuses([&] { admit_result_operation(q); });
        }
        {
            auto q = clone(request);
            q.at("query").set("analysis_sha256", Json::string(hash));
            refuses([&] { admit_result_operation(q); });
        }
        {
            auto q = clone(request);
            q.set("analysis_binding", Json::string("exact"));
            refuses([&] { admit_result_operation(q); });
            q.at("query").set("analysis_sha256", Json::string(hash));
            check(admit_result_operation(q).analysis_binding == "exact");
        }
        auto node = Json::parse(
            R"({"id":"pathway","type":"result","input":"source","result_input":"counts","request":{}})");
        node.set("request", request);
        const std::map<std::string, std::string> declared{
            {"source", "workspace"}, {"counts", "certificate"}, {"derived", "workspace"}};
        admit_result_step(node, declared);
        check(true);
        for (const auto *input : {"source", "derived", "future", "pathway"}) {
            auto n = clone(node);
            n.set("result_input", Json::string(input));
            refuses([&] { admit_result_step(n, declared); });
        }
        for (const auto *member : {"assay", "partition", "normalization_evidence_field"}) {
            auto n = clone(node);
            n.set(member, Json::string("unexpected"));
            refuses([&] { admit_result_step(n, declared); });
        }
        {
            auto n = clone(node);
            n.set("type", Json::string("query"));
            refuses([&] { admit_result_step(n, declared); });
        }
        // The original hash must refuse a source replaced after plan admission.
        write_text(file, "changed");
        refuses([&] { verified_record_image(file, hash, 1024); });
        check(verified_record_image(file, "sha256:" + sha256("changed"), 1024) == "changed");
        std::cout << checks << " pipeline result binding checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
