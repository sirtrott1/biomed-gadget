#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <rfc/runtime/execution_policy.hpp>
#include <algorithm>
#include <iostream>
#include <set>

using namespace gadget;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, "invalid evidence admitted");
}
Json parse(const char *s) { return Json::parse(s); }
void original_sources_contracts(const ScientificStudy &study, const fs::path &root) {
    const std::string bytes("original\0workbook\xff", 18);
    const auto source = root / "original bytes.bin";
    write_text(source, bytes);
    auto request = parse(R"({"schema":"gadget.biomedical.sources.request.v1","source_digest":"","method":"Declared source correspondence; binary provenance remains separate from numerical fields","sources":[{"id":"workbook","path":"","sha256":"","source_url":"offline:authored-original","license":"Authored test input"}]})");
    request.set("source_digest", Json::string(study.digest()));
    request.at("sources").at(0).set("path", Json::string(source.string()));
    request.at("sources").at(0).set("sha256", Json::string("sha256:" + sha256(bytes)));
    const auto original = study.native_publication().encode();
    const auto native_request = admit_record_view(request);
    const auto publication = retain_study_sources(study, native_request);
    const ScientificStudy retained(publication);
    check(retained.digest() == study.digest() && study.native_publication().encode() == original,
          "source provenance changed the original study or digest");
    const auto before = study.field_accession("all.Kd_nM"), after = retained.field_accession("all.Kd_nM");
    check(before.source_namespace() == after.source_namespace() && before.sequence() == after.sequence() &&
              rfc::storage::temporal_endpoint_detail::tensor_equal(
                  *before.tensor_declaration(), *after.tensor_declaration()) &&
              before.carrier()->key == after.carrier()->key,
          "source provenance changed actual numerical field pins");
    publication.write(root / "original-sources.rfcxs");
    const auto reopened = NativeRecord::read(root / "original-sources.rfcxs");
    check(reopened.encode() == publication.encode(), "source provenance archive changed bytes");
    check(NativeRecord::decode(reopened.at("source_provenance").at("request_publication").bytes()).encode() ==
              native_request.encode(), "source provenance lost its original request owner");
    const auto raw = reopened.at("source_provenance").at("sources").at(0).at("bytes").bytes();
    check(std::string(reinterpret_cast<const char *>(raw.data()), raw.size()) == bytes,
          "binary source provenance was rendered or altered");
    rejects([&] { retain_study_sources(retained, admit_record_view(request)); });
    for (const auto name : {"source_digest", "sha256", "id", "path"}) {
        auto bad = Json::parse(request.dump());
        if (std::string(name) == "source_digest") bad.set(name, Json::string("foreign"));
        else if (std::string(name) == "id") bad.at("sources").push(bad.at("sources").at(0));
        else bad.at("sources").at(0).set(name, Json::string("foreign"));
        rejects([&] { retain_study_sources(study, admit_record_view(bad)); });
    }
    auto bad = Json::parse(request.dump());
    bad.set("sources", Json::array());
    rejects([&] { retain_study_sources(study, admit_record_view(bad)); });
    const auto link = root / "source-link";
    fs::create_symlink(source, link);
    bad = Json::parse(request.dump());
    bad.at("sources").at(0).set("path", Json::string(link.string()));
    rejects([&] { retain_study_sources(study, admit_record_view(bad)); });
}
void reported_field_contracts(const EvidenceStore &store, const ScientificStudy &study,
                              const fs::path &root) {
    auto request = parse(R"({"schema":"gadget.biomedical.field.request.v1","source_digest":"","study_digest":"","field":"all.Kd_nM","endpoint":"Kd","units":"nM"})");
    request.set("source_digest", Json::string(store.digest()));
    request.set("study_digest", Json::string(study.digest()));
    const auto binding = store.bind_reported_field(study, admit_record_view(request));
    const auto result = binding.at("result");
    check(result.at("supported").integer() == 3 && result.at("correspondence").size() == 6,
          "reported field lost original complete correspondence");
    const auto rows = result.at("correspondence");
    check(rows.at(0).at("supported").boolean() && !rows.at(1).at("supported").boolean() &&
              !rows.at(2).at("supported").boolean() && rows.at(3).at("supported").boolean() &&
              !rows.at(4).at("supported").boolean() && rows.at(5).at("supported").boolean(),
          "censored, missing, other-unit or zero assertion was substituted");
    const auto selector = result.at("selector");
    const auto pin = study.field_accession("all.Kd_nM");
    check(selector.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
              selector.at("cutoff").integer() == pin.sequence() &&
              selector.at("tensor").integer() == pin.tensor_declaration()->key &&
              selector.at("output_carrier").integer() == *pin.tensor_declaration()->output &&
              selector.at("input_carrier").kind() == NativeRecordKind::Null,
          "reported field selector replaced its native owner/carrier");
    check(binding.source_namespace() == store.native_publication().source_namespace() &&
              binding.source_namespace() != pin.source_namespace(),
          "interpreted and quantitative sources were identified");
    binding.write(root / "field-binding.rfcxs");
    const auto reopened = NativeRecord::read(root / "field-binding.rfcxs");
    const auto field_source = NativeRecord::decode(reopened.at("result").at("field_publication").bytes());
    check(field_source.encode() == study.native_publication().historical_publication(pin.sequence()).encode() &&
              NativeRecord::decode(reopened.at("result").at("study_publication").bytes()).encode() ==
                  study.native_publication().encode(),
          "reported field archive changed its actual study or source");
    const auto actual = field_source.endpoint()->tensors;
    const auto selected = std::find_if(actual.begin(), actual.end(), [&](const auto &tensor) {
        return tensor.key == pin.tensor_declaration()->key;
    });
    check(selected != actual.end() && selected->tensor.state(3) == rfc::numeric::Rational{0} &&
              !selected->tensor.state(1) && !selected->tensor.state(2) && !selected->tensor.state(4),
          "reported field replay replaced supported zero or absent support");
    for (const auto &[name, wrong] : std::vector<std::pair<std::string, std::string>>{
             {"source_digest", "foreign"}, {"study_digest", "foreign"}, {"endpoint", "Ki"},
             {"units", "uM"}, {"field", "all.endpoint"}, {"field", "absent"}}) {
        auto bad = Json::parse(request.dump());
        bad.set(name, Json::string(wrong));
        rejects([&] { store.bind_reported_field(study, admit_record_view(bad)); });
    }
    auto other = Json::parse(request.dump());
    other.set("field", Json::string("all.Ki_nM"));
    other.set("endpoint", Json::string("Ki"));
    check(store.bind_reported_field(study, admit_record_view(other)).at("result").at("supported").integer() == 0,
          "empty endpoint field acquired unsupported values");
    const auto original_map = read_text(root / "store/study/record-map.json");
    auto altered = Json::parse(original_map);
    altered.at(0).set("assertion_key", altered.at(1).at("assertion_key"));
    write_text(root / "store/study/record-map.json", altered.dump());
    auto declaration = Json::parse(read_text(root / "store/study/study.json"));
    const auto documents = declaration.at("documents");
    for (std::size_t i = 0; i < documents.size(); ++i)
        if (documents.at(i).at("id").text() == "native_record_correspondence")
            documents.at(i).set("sha256", Json::string("sha256:" + sha256(altered.dump())));
    write_text(root / "store/study/altered-map-study.json", declaration.dump());
    const ScientificStudy altered_study(root / "store/study/altered-map-study.json");
    auto bad = Json::parse(request.dump());
    bad.set("study_digest", Json::string(altered_study.digest()));
    rejects([&] { store.bind_reported_field(altered_study, admit_record_view(bad)); });
    write_text(root / "store/study/record-map.json", original_map);
}
void proposal_contracts(const EvidenceStore &store, const fs::path &root, const Json &points) {
    auto request = parse(R"({"schema":"gadget.biomedical.context.v1","source_digest":"", "snapshot_sha256":"", "question":"Which original Kd observations support this threshold?", "assertion_keys":[]})");
    request.set("source_digest", Json::string(store.digest()));
    request.set("snapshot_sha256", render_record_view(store.native_publication().at("snapshot").at("snapshot_sha256")));
    for (std::size_t i = 0; i < points.size(); ++i)
        request.at("assertion_keys").push(points.at(i).at("native_assertion_key"));
    const auto native_request = admit_record_view(request);
    const auto context = store.source_context(native_request);
    const auto catalog = context.at("result");
    check(!catalog.at("availability").at("historical_availability_established").boolean(),
          "snapshot acquisition became historical availability");
    check(catalog.at("records").size() == points.size(), "context selection dropped assertions");
    for (std::size_t i = 0; i < points.size(); ++i)
        check(catalog.at("records").at(i).at("record").at("value_token").text() ==
                  points.at(i).at("value_token").text(), "context changed original source token");
    auto proposal = parse(R"({"schema":"gadget.biomedical.proposal.v1","source_digest":"", "context_sha256":"", "question_sha256":"", "action":"query", "query":{"operation":"affinity_evidence","endpoint":"Kd","units":"nM","maximum":"1"}})");
    for (const auto *name : {"source_digest", "context_sha256", "question_sha256"})
        proposal.set(name, render_record_view(catalog.at(name)));
    const auto question = request.at("question").text();
    const auto admit = [&](const Json &input) {
        return store.admit_proposal(native_request, admit_record_view(input), question);
    };
    check(admit(proposal).at("result").at("status").text() == "prospective",
          "prospective query executed or rejected");
    for (const auto *name : {"source_digest", "context_sha256", "question_sha256"}) {
        auto bad = Json::parse(proposal.dump());
        bad.set(name, Json::string("sha256:" + std::string(64, '0')));
        rejects([&] { admit(bad); });
    }
    for (const auto *name : {"units", "endpoint", "operation"}) {
        auto bad = Json::parse(proposal.dump());
        bad.at("query").set(name, Json::string(name == std::string("endpoint") ? "Ki" : "invented"));
        rejects([&] { admit(bad); });
    }
    auto bad = Json::parse(proposal.dump());
    bad.set("available_at", Json::integer(0));
    rejects([&] { admit(bad); });
    rejects([&] { store.admit_proposal(native_request, admit_record_view(proposal), "changed question"); });
    for (const auto *name : {"snapshot_sha256", "source_digest"}) {
        auto invalid_request = Json::parse(request.dump());
        invalid_request.set(name, Json::string("sha256:" + std::string(64, '0')));
        rejects([&] { store.source_context(admit_record_view(invalid_request)); });
    }
    auto invalid_request = Json::parse(request.dump());
    invalid_request.at("assertion_keys").push(points.at(0).at("native_assertion_key"));
    rejects([&] { store.source_context(admit_record_view(invalid_request)); });
    invalid_request = Json::parse(request.dump());
    invalid_request.set("cutoff", Json::integer(0));
    rejects([&] { store.source_context(admit_record_view(invalid_request)); });
    const auto bytes = read_text(root / "snapshot/activity.json");
    auto span = Json::object();
    span.set("document", Json::string("activity"));
    span.set("start", Json::integer(0));
    span.set("end", Json::integer(bytes.size()));
    span.set("text", Json::string(bytes));
    auto evidence = Json::array();
    evidence.push(span);
    auto mention = parse(R"({"id":"Kd-bound","kind":"measurement","assertion_key":"", "polarity":"asserted","interpretation":"Original censored observation, without point conversion.","evidence":[],"endpoint":"Kd","units":"nM","relation":"<","value_token":"2.5"})");
    mention.set("assertion_key", points.at(1).at("native_assertion_key"));
    mention.set("evidence", evidence);
    json_object_object_del(proposal.raw(), "query");
    proposal.set("action", Json::string("annotate"));
    auto mentions = Json::array();
    mentions.push(mention);
    proposal.set("mentions", mentions);
    const auto admitted = admit(proposal);
    check(!admitted.at("result").at("mentions").at(0).at("interpretation_verified").boolean(),
          "annotation interpretation was promoted to evidence");
    for (const auto *name : {"units", "endpoint", "relation", "value_token", "polarity"}) {
        bad = Json::parse(proposal.dump());
        bad.at("mentions").at(0).set(name, Json::string("invented"));
        rejects([&] { admit(bad); });
    }
    for (const auto *name : {"text", "document"}) {
        bad = Json::parse(proposal.dump());
        bad.at("mentions").at(0).at("evidence").at(0).set(name, Json::string("invented"));
        rejects([&] { admit(bad); });
    }
    const auto position = bytes.find("β");
    check(position != std::string::npos, "missing UTF-8 evidence control");
    bad = Json::parse(proposal.dump());
    auto split = bad.at("mentions").at(0).at("evidence").at(0);
    split.set("start", Json::integer(position));
    split.set("end", Json::integer(position + 1));
    split.set("text", Json::string(bytes.substr(position, 1)));
    rejects([&] { admit(bad); });
    mention = parse(R"({"id":"compound","kind":"entity","assertion_key":"", "role":"compound","entity":"chembl:CHEMBL941", "polarity":"uncertain","interpretation":"Reported compound identifier.","evidence":[]})");
    mention.set("assertion_key", points.at(0).at("native_assertion_key"));
    mention.set("evidence", evidence);
    proposal.at("mentions").push(mention);
    check(admit(proposal).at("result").at("mentions").size() == 2, "entity annotation dropped");
    bad = Json::parse(proposal.dump());
    bad.at("mentions").at(1).set("entity", Json::string("pubchem:cid:100"));
    rejects([&] { admit(bad); });
    bad = Json::parse(proposal.dump());
    bad.at("mentions").push(bad.at("mentions").at(1));
    rejects([&] { admit(bad); });
    admitted.write(root / "proposal.rfcxs");
    const auto restored = NativeRecord::read(root / "proposal.rfcxs");
    check(restored.at("result").at("context").at("context_sha256").text() ==
              catalog.at("context_sha256").text(), "proposal replay changed context");
    check(restored.at("result").at("context").at("records").at(1).at("record").at("relation").text() == "<",
          "proposal replay lost original censoring");
    const auto original = store.native_publication().raw_sources();
    const auto retained = restored.raw_sources();
    for (const auto &source : original) {
        const auto found = std::find_if(retained.begin(), retained.end(),
                                       [&](const auto &row) { return row.key() == source.key(); });
        check(found != retained.end() && std::ranges::equal(found->bytes(), source.bytes()),
              "proposal replay changed original provider bytes");
    }
}
void relationship_contracts(const fs::path &root) {
    const auto snapshot = root/"relationship-snapshot";
    fs::copy(root/"snapshot",snapshot,fs::copy_options::recursive);
    fs::remove(snapshot/"snapshot.rfcxs");
    auto manifest = Json::parse(read_text(snapshot/"snapshot.json"));
    auto responses = manifest.at("records");
    const auto add = [&](const char *id,const char *provider,const char *resource,const char *identifier,const Json &body) {
        const auto bytes = body.dump();
        const auto headers = std::string("authored source relationship fixture\n");
        const auto path = std::string(id)+".json", header_path = std::string(id)+".headers";
        write_text(snapshot/path,bytes);write_text(snapshot/header_path,headers);
        auto item = Json::object(), request = Json::object();
        request.set("id",Json::string(id));request.set("provider",Json::string(provider));
        request.set("resource",Json::string(resource));request.set("identifier",Json::string(identifier));
        item.set("request",request);item.set("url",Json::string("offline:authored-source-relationship"));
        item.set("path",Json::string(path));item.set("sha256",Json::string("sha256:"+sha256(bytes)));
        item.set("headers_path",Json::string(header_path));item.set("headers_sha256",Json::string("sha256:"+sha256(headers)));
        responses.push(item);
    };
    add("protein-relationships","uniprot","protein","P00519",parse(R"({
        "primaryAccession":"P00519","uniProtKBCrossReferences":[{"database":"PDB","id":"1TEST"}],
        "comments":[{"commentType":"DISEASE","disease":{"diseaseId":"Authored condition","diseaseAccession":"DI-TEST","diseaseCrossReference":{"database":"MIM","id":"123"}}},
                    {"commentType":"DISEASE","note":{"texts":[{"value":"Authored statement without one condition accession"}]}}]})"));
    add("indication-relationships","chembl","drug_indication","CHEMBL941",parse(R"({"drug_indications":[
        {"molecule_chembl_id":"CHEMBL941","efo_id":"EFO:TEST","efo_term":"Authored condition","mesh_id":"DTEST","mesh_heading":"Authored condition","max_phase_for_ind":"4"},
        {"molecule_chembl_id":"CHEMBL941","efo_id":"EFO:TEST","efo_term":"Authored condition","mesh_id":null,"max_phase_for_ind":"2"},
        {"molecule_chembl_id":"CHEMBL941","efo_id":null,"mesh_id":"DOTHER","mesh_heading":"Other authored condition","max_phase_for_ind":"1"}]})"));
    json_object_object_del(manifest.raw(),"snapshot_sha256");
    manifest.set("snapshot_sha256",Json::string("sha256:"+sha256(manifest.dump())));
    write_text(snapshot/"snapshot.json",manifest.dump());
    import_evidence(snapshot,root/"relationship-store",3);
    rejects([&] { import_evidence(snapshot,root/"relationship-v2",2); });
    const EvidenceStore store(root/"relationship-store");
    check(store.execute(parse(R"({"operation":"coverage"})")).at("result").at("binding_records").integer()==8,
          "source statements became numerical binding observations");
    const auto records = store.native_publication().at("records");
    std::size_t relationships{};
    bool statement=false;
    for (std::size_t i=0;i<records.size();++i) {
        const auto row=records.at(i);
        if (row.at("kind").integer()!=30) continue;
        ++relationships;
        check(row.at("value_token").text().empty() && row.at("endpoint").text().empty(),
              "disease annotation or indication phase became a numerical endpoint");
        statement=statement || row.at("target").text().find("/comments/1")!=std::string::npos;
    }
    check(relationships==5 && statement,"source relationships dropped original duplicate or unassigned condition statements");
    auto query=parse(R"({"operation":"role_join","source_digest":"","left_kind":"binding_observation","left_role":"target_or_right","right_kind":"identifier_cross_reference","right_role":"compound_or_left","activity":"retained","max_matches":100})");
    query.set("source_digest",Json::string(store.digest()));
    check(store.execute(query).at("result").at("selected").integer()==8,"native target role composition lost complete binding occurrences");
    query.set("left_role",Json::string("compound_or_left"));
    query.set("right_kind",Json::string("source_relationship"));
    const auto joined=store.execute_record(admit_record_view(query));
    const auto result=joined.at("result").at("result");
    check(result.at("selected").integer()==18 && result.at("matches").size()==18,
          "native role composition merged distinct source indication rows");
    std::set<std::pair<std::string,std::string>> witnesses;
    for (std::size_t i=0;i<result.at("matches").size();++i) {
        const auto match=result.at("matches").at(i);
        check(match.at("participant").text()=="chembl:CHEMBL941","native join inferred a foreign participant");
        check(match.at("left_position").integer()==0 && match.at("right_position").integer()==0,
              "native join lost source occurrence positions");
        check(witnesses.emplace(match.at("left_occurrence").text(),match.at("right_occurrence").text()).second,
              "native join repeated an occurrence witness");
    }
    joined.write(root/"role-join.rfcxs");
    check(NativeRecord::read(root/"role-join.rfcxs").encode()==joined.encode(),"role composition archive changed its actual identity");
    auto bad=Json::parse(query.dump());bad.set("max_matches",Json::integer(0));
    bool exhausted=false;
    try { (void)store.execute(bad); }
    catch (const rfc::runtime::ExecutionLimitExceeded &error) {
        exhausted=error.resource()=="role join result matches" && error.required()==18 && error.limit()==0;
    }
    check(exhausted,"native role result budget did not report actual complete match demand");
    bad=Json::parse(query.dump());bad.set("source_digest",Json::string("foreign"));
    rejects([&] { store.execute(bad); });
    bad=Json::parse(query.dump());bad.set("activity",Json::string("inferred"));
    rejects([&] { store.execute(bad); });
    bad=Json::parse(query.dump());bad.set("left_role",Json::string("unknown"));
    rejects([&] { store.execute(bad); });
    const auto inspection=store.execute(parse(R"({"operation":"inspect"})")).at("result");
    std::set<std::string> selected;
    for (std::size_t i=0;i<inspection.at("entities").size();++i) {
        auto neighbors=parse(R"({"operation":"neighbors","subject_key":"1"})");
        neighbors.set("subject_key",inspection.at("entities").at(i).at("native_subject_key"));
        const auto rows=store.execute(neighbors).at("result").at("records");
        for (std::size_t j=0;j<rows.size();++j) selected.insert(rows.at(j).at("native_assertion_key").text());
    }
    auto context=parse(R"({"schema":"gadget.biomedical.context.v1","source_digest":"","snapshot_sha256":"","question":"Which retained source role occurrences match?","assertion_keys":[]})");
    context.set("source_digest",Json::string(store.digest()));
    context.set("snapshot_sha256",render_record_view(store.native_publication().at("snapshot").at("snapshot_sha256")));
    for (const auto &key:selected) context.at("assertion_keys").push(Json::string(key));
    const auto native_context=admit_record_view(context);
    const auto catalog=store.source_context(native_context).at("result");
    auto proposal=parse(R"({"schema":"gadget.biomedical.proposal.v1","source_digest":"","context_sha256":"","question_sha256":"","action":"query","query":{}})");
    for (const auto *name:{"source_digest","context_sha256","question_sha256"}) proposal.set(name,render_record_view(catalog.at(name)));
    proposal.set("query",query);
    check(store.admit_proposal(native_context,admit_record_view(proposal),context.at("question").text()).at("result").at("status").text()=="prospective",
          "native role query proposal executed or lost its complete selected context");
    auto partial_keys=Json::array();partial_keys.push(Json::string(*selected.begin()));
    context.set("assertion_keys",partial_keys);
    const auto native_partial_context=admit_record_view(context);
    const auto partial=store.source_context(native_partial_context).at("result");
    proposal.set("context_sha256",render_record_view(partial.at("context_sha256")));
    rejects([&] { store.admit_proposal(native_partial_context,admit_record_view(proposal),context.at("question").text()); });
    fs::create_directory(root/"relationship-file-only");
    fs::copy_file(root/"relationship-store/store.rfcxs",root/"relationship-file-only/store.rfcxs");
    const EvidenceStore cold(root/"relationship-file-only");
    check(cold.native_publication().encode()==store.native_publication().encode() &&
          cold.execute(query).at("result").at("selected").integer()==18,
          "standalone native role composition lost original source or complete witnesses");
}
int main() {
    try {
        const fs::path root = GADGET_BIOMEDICAL_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root / "snapshot");
        auto manifest = Json::object(), responses = Json::array();
        manifest.set("schema", Json::string("gadget.database.snapshot.v1"));
        manifest.set("id", Json::string("authored-binding-fixture"));
        const auto add = [&](const char *id, const char *provider, const char *resource,
                             const Json &body) {
            const auto bytes = body.dump();
            write_text(root / "snapshot" / (std::string(id) + ".json"), bytes);
            write_text(root / "snapshot" / (std::string(id) + ".headers.txt"),
                       "authored offline fixture\n");
            auto r = Json::object(), q = Json::object();
            q.set("id", Json::string(id));
            q.set("provider", Json::string(provider));
            q.set("resource", Json::string(resource));
            q.set("identifier", Json::string("P00519"));
            r.set("request", q);
            r.set("url", Json::string("offline:authored-fixture"));
            r.set("path", Json::string(std::string(id) + ".json"));
            r.set("sha256", Json::string("sha256:" + sha256(bytes)));
            r.set("headers_path", Json::string(std::string(id) + ".headers.txt"));
            r.set("headers_sha256", Json::string("sha256:" + sha256("authored offline fixture\n")));
            responses.push(r);
        };
        auto activities = Json::object(), rows = Json::array();
        for (std::size_t i = 0; i < 6; ++i) {
            auto row = parse(
                R"({"molecule_chembl_id":"CHEMBL941","target_chembl_id":"CHEMBL1862","assay_chembl_id":"CHEMBLTEST","document_chembl_id":"CHEMBLDOC","standard_type":"Kd","standard_relation":"=","standard_units":"nM"})");
            row.set("activity_id", Json::integer(static_cast<std::int64_t>(100 + i)));
            row.set("standard_value", i == 2 ? Json{}
                                             : Json::string(i == 0   ? "1.25"
                                                            : i == 1 ? "2.5"
                                                            : i == 3 ? "0"
                                                            : i == 4 ? "5"
                                                                     : "1.5"));
            if (i == 1)
                row.set("standard_relation", Json::string("<"));
            if (i == 2)
                row.set("standard_relation", Json{});
            if (i == 4)
                row.set("standard_units", Json::string("uM"));
            rows.push(row);
        }
        activities.set("activities", rows);
        activities.set("description", Json::string("β source quote control"));
        add("activity", "chembl", "activity", activities);
        add("binding", "bindingdb", "uniprot",
            parse(
                R"({"getLindsByUniprotsResponse":{"affinities":[{"monomerid":"10","affinity":"<1.000","affinity_type":"Ki","doi":"10.test/example"},{"monomerid":"10","affinity":"1.000","affinity_type":"Ki","doi":"10.test/example"}]}})"));
        add("target", "chembl", "target",
            parse(
                R"({"target_chembl_id":"CHEMBL1862","target_components":[{"accession":"P00519"}]})"));
        add("molecule", "chembl", "molecule",
            parse(
                R"({"molecule_chembl_id":"CHEMBL941","cross_references":[{"xref_src":"PubChem","xref_id":"100"}]})"));
        add("compound", "pubchem", "compound",
            parse(
                R"({"PropertyTable":{"Properties":[{"CID":100,"InChIKey":"AUTHORED-TEST-KEY"}]}})"));
        manifest.set("records", responses);
        manifest.set("snapshot_sha256", Json::string("sha256:" + sha256(manifest.dump())));
        write_text(root / "snapshot/snapshot.json", manifest.dump());
        const auto relative_snapshot =
            capture_database_snapshot(fs::relative(root / "snapshot", fs::current_path()));
        check(relative_snapshot.at("artifacts").size() == responses.size() * 2,
              "relative database capture lost complete response/header bytes");
        const auto retained_response = relative_snapshot.at("artifacts").at(0).at("bytes").bytes();
        const auto original_response = read_text(root / "snapshot/activity.json");
        check(std::string(reinterpret_cast<const char *>(retained_response.data()),
                          retained_response.size()) == original_response,
              "relative database capture changed original provider bytes");
        import_evidence(root / "snapshot", root / "store");
        import_evidence(root / "snapshot", root / "legacy", 1);
        auto legacy_receipt = Json::parse(read_text(root / "legacy/store.json"));
        json_object_object_del(legacy_receipt.raw(), "adapter_version");
        json_object_object_del(legacy_receipt.raw(), "source_digest");
        legacy_receipt.set("source_digest",
                           Json::string("sha256:" + sha256(legacy_receipt.dump())));
        write_text(root / "legacy/store.json", legacy_receipt.dump());
        EvidenceStore legacy(root / "legacy");
        check(legacy.execute(parse(R"({"operation":"coverage"})"))
                      .at("result")
                      .at("binding_records")
                      .integer() == 8,
              "older frozen interpretation changed on load");
        rejects([&] { import_evidence(root / "snapshot", root / "unsupported", 4); });
        EvidenceStore store(root / "store");
        auto coverage = store.execute(parse(R"({"operation":"coverage"})")).at("result");
        const auto native_coverage = store.execute_record(admit_record_view(parse(R"({"operation":"coverage"})")));
        check(native_coverage.at("result").at("result").at("binding_records").kind() ==
                  NativeRecordKind::Integer &&
                  native_coverage.at("result").at("result").at("binding_records").integer() == 8,
              "native coverage did not retain exact provider counts");
        check(coverage.at("binding_records").integer() == 8 &&
                  coverage.at("reported").integer() == 7 && coverage.at("missing").integer() == 1,
              "missing/censoring/duplicates lost");
        check(coverage.at("equality_relation_records").integer() == 5,
              "relation categories coerced");
        const auto selected =
            store
                .execute(parse(
                    R"({"operation":"select","field":"relation","comparison":"equal","value":"<"})"))
                .at("result");
        check(selected.at("selected").integer() == 2, "native exact category selection");
        const auto bounds =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Ki","units":"nM","maximum":"0.5"})"))
                .at("result");
        check(bounds.at("selected").integer() == 0 && bounds.at("contradicted").integer() == 1 &&
                  bounds.at("indeterminate").integer() == 1,
              "censored threshold treated as a point");
        const auto exact =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Ki","units":"nM","maximum":"1"})"))
                .at("result");
        check(exact.at("selected").integer() == 2, "strict bound entailment at cutoff");
        const auto kd =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Kd","units":"nM","maximum":"1"})"))
                .at("result");
        check(kd.at("selected").integer() == 1 && kd.at("contradicted").integer() == 2 &&
                  kd.at("indeterminate").integer() == 2,
              "missing/unit/zero threshold classification");
        auto inspect = store.execute(parse(R"({"operation":"inspect"})"));
        std::string cid, sid, protein;
        const auto entities = inspect.at("result").at("entities");
        for (std::size_t i = 0; i < entities.size(); ++i) {
            const auto j = entities.at(i);
            const auto id = j.at("id").text();
            if (id == "pubchem:cid:100")
                cid = j.at("native_subject_key").text();
            if (id == "pubchem:sid:100")
                sid = j.at("native_subject_key").text();
            if (id == "uniprot:P00519")
                protein = j.at("native_subject_key").text();
        }
        check(!cid.empty() && !sid.empty() && cid != sid, "CID/SID identity collapsed");
        auto neighbors = parse(R"({"operation":"neighbors","subject_key":"1"})");
        neighbors.set("subject_key", Json::string(protein));
        const auto bindings = store.execute(neighbors).at("result").at("records");
        check(bindings.size() == 3, "native role traversal/cross-reference lost");
        auto points =
            store
                .execute(parse(
                    R"({"operation":"select","field":"provider","comparison":"equal","value":"chembl"})"))
                .at("result")
                .at("records");
        auto compare = parse(R"({"operation":"compare","left":"1","right":"2"})");
        proposal_contracts(store, root, points);
        compare.set("left", points.at(0).at("native_assertion_key"));
        compare.set("right", points.at(5).at("native_assertion_key"));
        const auto difference = store.execute(compare).at("result");
        const auto native_difference = store.execute_record(admit_record_view(compare));
        check(native_difference.at("result").at("result").at("reported_difference").rational() ==
                  rfc::numeric::Rational(-1, 4), "native comparison reparsed a numerical view");
        check(native_difference.at("result").at("result").at("reported_ratio").rational() ==
                  rfc::numeric::Rational(5, 6), "native ratio differs from original point tokens");
        native_difference.write(root / "native-comparison.rfcxs");
        auto claims = parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[
            {"path":["reported_difference"],"value_type":"rational","comparison":"equal","value":"-1/4"},
            {"path":["reported_ratio"],"value_type":"rational","comparison":"equal","value":"5/6"}
        ]})");
        claims.set("source_namespace", Json::string(std::to_string(native_difference.source_namespace())));
        claims.set("cutoff", Json::integer(static_cast<std::int64_t>(native_difference.publication_history()->last_sequence())));
        const auto checked_claims = check_certificate_claims_record(native_difference, admit_record_view(claims));
        check(checked_claims.at("result").at("matched").integer() == 2 &&
                  checked_claims.at("result").at("checks").at(0).at("observed").rational() == rfc::numeric::Rational(-1, 4),
              "native evidence claim used a rendered numerical observation");
        check(NativeRecord::decode(checked_claims.at("result").at("source_publication").bytes()).encode() == native_difference.encode(),
              "native evidence claim lost its supplying operation");
        check(NativeRecord::read(root / "native-comparison.rfcxs").at("result").at("result")
                  .at("reported_difference").rational() == rfc::numeric::Rational(-1, 4),
              "native comparison replay lost exact value");
        compare.set("right", points.at(3).at("native_assertion_key"));
        const auto zero_denominator = store.execute_record(admit_record_view(compare));
        check(zero_denominator.at("result").at("result").at("reported_ratio").kind() ==
                  NativeRecordKind::Null &&
                  zero_denominator.at("result").at("result").at("right").at("value_token").text() == "0",
              "undefined ratio lost its supported zero input");
        check(difference.at("reported_difference").at("rational").text() == "-1/4" &&
                  difference.at("reported_ratio").at("rational").text() == "5/6",
              "reported exact arithmetic");
        compare.set("right", points.at(1).at("native_assertion_key"));
        rejects([&] { store.execute(compare); });
        compare.set("right", points.at(4).at("native_assertion_key"));
        rejects([&] { store.execute(compare); });
        auto study = std::make_shared<const ScientificStudy>(root / "store/study/study.json");
        reported_field_contracts(store, *study, root);
        original_sources_contracts(*study, root);
        ScientificSession session(study);
        const auto summary =
            session
                .execute(parse(
                    R"({"operation":"field_summary","field":"all.Kd_nM","support":"reported_only"})"))
                .at("result");
        check(summary.at("reported").integer() == 3 && summary.at("missing").integer() == 3 &&
                  summary.at("sum").at("rational").text() == "11/4",
              "censored threshold or other unit became point/zero lost");
        auto workflow = parse(
            R"({"schema":"gadget.scientific.workflow.v1","id":"evidence-to-math","source_digest":"","evidence_digest":"","limits":{"depth":8,"nodes":16,"model_calls":0,"parallel_workers":2},"definitions":[],"entry":{"id":"root","type":"parallel","workers":2,"steps":[{"id":"evidence","type":"evidence","query":{"operation":"coverage"}},{"id":"math","type":"query","query":{"operation":"support_mask","fields":["all.Kd_nM"],"mode":"all_reported","backend":"cpu_tiles","workers":2,"tile_words":1}}]}})");
        workflow.set("source_digest", Json::string(study->digest()));
        workflow.set("evidence_digest", Json::string(store.digest()));
        const auto receipt =
            execute_workflow(study, workflow, {}, CpuResources::local(), {}, &store);
        check(receipt.at("status").text() == "complete", "native delegated evidence/math workflow");
        auto native_workflow = Json::parse(workflow.dump());
        native_workflow.at("entry").at("steps").at(1).set("query", parse(
            R"({"operation":"linear_model_comparison","responses":["all.Kd_nM"],"terms":[],"controls":[],"selection":[],"support":"joint_reported","metric":"identity","intercept":true})"));
        const auto operation = execute_workflow_record(
            study, admit_record_view(native_workflow), {}, CpuResources::local(), {}, &store);
        check(workflow_record_view(operation).at("status").text() == "complete",
              "owning native workflow did not preserve its actual query results");
        operation.write(root / "complete-operation.rfcxs");
        fs::rename(root / "store", root / "unavailable-store");
        fs::rename(root / "snapshot", root / "unavailable-snapshot");
        const auto restored_operation = NativeRecord::read(root / "complete-operation.rfcxs");
        const auto sources = restored_operation.parent();
        check(sources && sources->at("evidence_digest").text() == store.digest(),
              "workflow lost its actual evidence source binding");
        const auto restored_evidence =
            NativeRecord::decode(sources->at("evidence_publication").bytes());
        check(restored_evidence.encode() == store.native_publication().encode(),
              "workflow changed or omitted its native evidence publication");
        const auto expected_raw = store.native_publication().raw_sources();
        const auto restored_raw = restored_evidence.raw_sources();
        check(expected_raw.size() == restored_raw.size(), "workflow lost original byte owners");
        for (std::size_t i = 0; i < expected_raw.size(); ++i)
            check(expected_raw[i].key() == restored_raw[i].key() &&
                      std::ranges::equal(expected_raw[i].bytes(), restored_raw[i].bytes()),
                  "workflow changed original database response/header bytes");
        fs::rename(root / "unavailable-store", root / "store");
        fs::rename(root / "unavailable-snapshot", root / "snapshot");
        workflow.set("evidence_digest", Json::string(study->digest()));
        rejects([&] { execute_workflow(study, workflow, {}, CpuResources::local(), {}, &store); });
        auto
            media = parse(R"([{"sha256":"sha256:authored-image"}])"),
            observation = parse(
                R"({"schema":"gadget.model.observation.v1","observations":[{"statement":"a white cell","image_index":0}],"limitations":"authored contract fixture"})");
        check(admit_observation(observation, study->digest(), media).at("status").text() ==
                  "prospective",
              "observation promoted or rejected");
        observation.at("observations").at(0).set("image_index", Json::integer(1));
        rejects([&] { admit_observation(observation, study->digest(), media); });
        media.at(0).set("width", Json::integer(100));
        media.at(0).set("height", Json::integer(80));
        observation.at("observations").at(0).set("image_index", Json::integer(0));
        observation.at("observations").at(0).set("region", parse("[0,10,100,80]"));
        const auto regional = admit_observation(observation, study->digest(), media);
        check(regional.at("certificate").at("claims").at(0).at("region").dump() ==
                  parse("[0,10,100,80]").dump() &&
                  regional.at("certificate").at("claims").at(0).at("region_space").text() ==
                      "original_image_pixel_edges",
              "image region lost original coordinates or became a measurement");
        for (const auto *invalid : {"[0,0,101,80]", "[0,0,100,81]", "[10,0,10,20]",
                                    "[20,0,10,20]", "[-1,0,10,20]", "[0,0,10.5,20]"}) {
            observation.at("observations").at(0).set("region", parse(invalid));
            rejects([&] { admit_observation(observation, study->digest(), media); });
        }
        observation.at("observations").at(0).set("region", parse("[0,0,10,20]"));
        observation.at("observations").at(0).set("image_index", parse("null"));
        rejects([&] { admit_observation(observation, study->digest(), media); });
        const auto before = store.execute(parse(R"({"operation":"coverage"})")).dump();
        coverage.set("missing", Json::integer(0));
        check(before == store.execute(parse(R"({"operation":"coverage"})")).dump(),
              "public JSON mutated source");
        relationship_contracts(root);
        auto corrupt = read_text(root / "store/evidence.rfcxi");
        corrupt.back() ^= 1;
        write_text(root / "store/evidence.rfcxi", corrupt);
        rejects([&] { EvidenceStore damaged(root / "store"); });
        write_text(root / "snapshot/activity.json", "{}");
        rejects([&] { import_evidence(root / "snapshot", root / "corrupt"); });
        check(!fs::exists(root / "corrupt"), "failed import published output");
        std::cout << "PASS: retained binding evidence, native queries, censoring, missing/zero, "
                     "CID/SID, exact comparison, source binding, delegated workflow and "
                     "prospective observations\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
