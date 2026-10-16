#include <gadget/assay.hpp>
#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <iostream>
#include <set>

int main(int argc, char **argv) try {
    if (argc != 4)
        throw std::invalid_argument("usage: gadget_biomedical_collection WORKSPACE REQUEST NEW_RECORD");
    const gadget::fs::path destination(argv[3]);
    if (gadget::fs::exists(destination)) throw std::invalid_argument("collection destination exists");
    const gadget::AssayWorkspace workspace(argv[1]);
    const auto request = gadget::read_record(argv[2]).native();
    if (request.at("schema").text() != "gadget.biomedical.collection.request.v1")
        throw std::invalid_argument("unsupported biomedical collection request");
    const std::set<std::string> required{"schema","id","design","assays","method","sources"};
    const auto names = request.names();
    if (std::set<std::string>(names.begin(),names.end()) != required)
        throw std::invalid_argument("collection request requires its exact declared members");
    std::map<std::string, std::shared_ptr<const gadget::ScientificStudy>> studies;
    const auto selected = request.at("assays");
    if (selected.size() == 0) throw std::invalid_argument("collection requires selected source assays");
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto name = selected.at(i).text();
        if (name.empty() || name.find_first_of("/\\") != std::string::npos || name=="." || name=="..")
            throw std::invalid_argument("assay name cannot escape its collection path");
        if (studies.contains(name)) throw std::invalid_argument("duplicate selected assay");
        const auto study = workspace.study(name);
        auto provenance = gadget::Json::object();
        provenance.set("schema", gadget::Json::string("gadget.biomedical.sources.request.v1"));
        provenance.set("source_digest", gadget::Json::string(study->digest()));
        provenance.set("method", gadget::render_record_view(request.at("method")));
        provenance.set("sources", gadget::render_record_view(request.at("sources")));
        const auto retained = gadget::retain_study_sources(*study,gadget::admit_record_view(provenance));
        auto owner = std::make_shared<const gadget::ScientificStudy>(retained);
        if (owner->digest() != study->digest())
            throw std::runtime_error("original source retention changed the admitted scientific study");
        studies.emplace(name,std::move(owner));
    }
    const auto publication = gadget::AssayWorkspace::publish_studies(request,request.at("id").text(),studies);
    publication.write(destination);
    std::cout << "Retained " << studies.size() << " independent source studies and original byte owners\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
