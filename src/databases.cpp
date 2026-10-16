#include "scientific_detail.hpp"
#include <array>
#include <chrono>
#include <ctime>
#include <gadget/biomedical.hpp>
#include <rfc/machine/checked_size.hpp>
#include <thread>
#if GADGET_HAS_PUBLIC_DATABASES
#include <curl/curl.h>
#endif

namespace gadget {
namespace {
using namespace scientific_detail;
struct Provider {
    const char *id, *docs, *license;
};
constexpr Provider providers[] = {
    {"chembl",
     "https://chembl.gitbook.io/chembl-interface-documentation/web-services/"
     "chembl-data-web-services",
     "https://chembl.github.io/faq-is-there-license-agreement-i-need/"},
    {"uniprot", "https://www.uniprot.org/help/api_queries", "https://www.uniprot.org/help/license"},
    {"pubchem", "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
     "https://pubchem.ncbi.nlm.nih.gov/docs/usage-and-privacy"},
    {"pdb", "https://www.rcsb.org/docs/programmatic-access/web-apis-overview",
     "https://www.wwpdb.org/about/usage-policies"},
    {"bindingdb", "https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp",
     "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp"}};
std::string token(const Json &j, std::size_t maximum = source_byte_limit()) {
    const auto value = text(j, maximum);
    require(!value.empty() && std::all_of(value.begin(), value.end(),
                                          [](unsigned char c) {
                                              return (c >= 'A' && c <= 'Z') ||
                                                     (c >= 'a' && c <= 'z') ||
                                                     (c >= '0' && c <= '9') || c == '_' || c == '-';
                                          }),
            "database identifier must be an explicit ASCII token");
    return value;
}
[[maybe_unused]] std::string url(const Json &item) {
    members(item, {"id", "provider", "resource", "identifier"},
            {"limit", "offset", "endpoint", "cutoff_nm"});
    (void)token(item.at("id"));
    const auto provider = item.at("provider").text(), resource = item.at("resource").text();
    const auto id = token(item.at("identifier"));
    if (provider == "chembl") {
        require(resource == "activity" || resource == "assay" || resource == "target" ||
                    resource == "molecule" || resource == "status" || resource == "drug_indication",
                "unsupported ChEMBL resource");
        require(!item.has("cutoff_nm"), "cutoff applies only to BindingDB");
        if (resource != "activity" && resource != "drug_indication") {
            require(!item.has("limit") && !item.has("offset") && !item.has("endpoint"),
                    "record fetch has no page filters");
            return "https://www.ebi.ac.uk/chembl/api/data/" + resource +
                   (resource == "status" ? ".json" : "/" + id + ".json");
        }
        require(id.starts_with("CHEMBL"), "paged ChEMBL requests require a molecule identifier");
        require(resource=="activity" || !item.has("endpoint"), "drug indication pages have no binding endpoint filter");
        const auto limit =
            item.has("limit") ? count(item.at("limit"), 1000, "API page limit") : 100;
        // ChEMBL imposes the 1000-row page size. It does not impose an
        // application-wide observation count or a million-record offset.
        const auto offset = item.has("offset")
                                ? count(item.at("offset"), std::numeric_limits<std::size_t>::max(),
                                        "API page offset")
                                : 0;
        require(limit > 0, "empty API page limit");
        auto result =
            "https://www.ebi.ac.uk/chembl/api/data/" + resource + ".json?molecule_chembl_id=" + id +
            "&limit=" + std::to_string(limit) + "&offset=" + std::to_string(offset);
        if (item.has("endpoint")) {
            const auto e = token(item.at("endpoint"));
            require(e == "Kd" || e == "Ki" || e == "IC50" || e == "EC50",
                    "unsupported measurement endpoint");
            result += "&standard_type=" + e;
        }
        return result;
    }
    require(!item.has("limit") && !item.has("offset") && !item.has("endpoint"),
            "unsupported record filters");
    if (provider == "bindingdb") {
        require(resource == "uniprot", "BindingDB profile requires explicit UniProt target");
        require(item.has("cutoff_nm"),
                "BindingDB finite profile requires declared affinity cutoff in nM");
        const auto cutoff = count(item.at("cutoff_nm"), std::numeric_limits<std::size_t>::max(),
                                  "BindingDB cutoff");
        require(cutoff > 0, "BindingDB cutoff must be positive");
        return "https://www.bindingdb.org/rest/getLigandsByUniprots?uniprot=" + id +
               "&cutoff=" + std::to_string(cutoff) + "&response=application/json";
    }
    require(!item.has("cutoff_nm"), "cutoff applies only to BindingDB");
    if (provider == "uniprot" && resource == "protein")
        return "https://rest.uniprot.org/uniprotkb/" + id + ".json";
    if (provider == "pubchem" && resource == "compound") {
        require(
            std::all_of(id.begin(), id.end(), [](unsigned char c) { return c >= '0' && c <= '9'; }),
            "PubChem requires a CID");
        return "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + id +
               "/property/InChIKey,MolecularFormula/JSON";
    }
    if (provider == "pdb" && resource == "entry")
        return "https://data.rcsb.org/rest/v1/core/entry/" + id;
    if (provider == "pdb" && resource == "polymer") {
        const auto split = id.find('_');
        require(split != std::string::npos && split && split + 1 < id.size(),
                "PDB polymer requires entry_entity");
        return "https://data.rcsb.org/rest/v1/core/polymer_entity/" + id.substr(0, split) + "/" +
               id.substr(split + 1);
    }
    throw std::invalid_argument("unsupported database provider/resource");
}
#if GADGET_HAS_PUBLIC_DATABASES
struct Response {
    std::string body, headers;
    std::size_t body_limit{}, header_limit{};
    bool overflow{};
};
size_t body(char *data, size_t size, size_t count, void *user) noexcept {
    auto &out = *static_cast<Response *>(user);
    try {
        const auto n = rfc::machine::checked_size_mul(size, count);
        if (n > out.body_limit - out.body.size()) {
            out.overflow = true;
            return 0;
        }
        out.body.append(data, n);
        return n;
    } catch (...) {
        out.overflow = true;
        return 0;
    }
}
size_t header(char *data, size_t size, size_t count, void *user) noexcept {
    auto &out = *static_cast<Response *>(user);
    try {
        const auto n = rfc::machine::checked_size_mul(size, count);
        if (n > out.header_limit - out.headers.size()) {
            out.overflow = true;
            return 0;
        }
        out.headers.append(data, n);
        return n;
    } catch (...) {
        out.overflow = true;
        return 0;
    }
}
Response retrieve(const std::string &address, const Json &transport, std::size_t remaining) {
    static const auto initialized = curl_global_init(CURL_GLOBAL_DEFAULT);
    require(initialized == CURLE_OK, "CURL initialization failed");
    std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> handle(curl_easy_init(), curl_easy_cleanup);
    require(static_cast<bool>(handle), "CURL allocation failed");
    Response out;
    out.body_limit = std::min(
        count(transport.at("max_body_bytes"), source_byte_limit(), "response bytes"), remaining);
    out.header_limit = count(transport.at("max_header_bytes"), source_byte_limit(), "header bytes");
    curl_easy_setopt(handle.get(), CURLOPT_URL, address.c_str());
    curl_easy_setopt(handle.get(), CURLOPT_PROTOCOLS_STR, "https");
    curl_easy_setopt(handle.get(), CURLOPT_FOLLOWLOCATION, 0L);
    curl_easy_setopt(handle.get(), CURLOPT_CONNECTTIMEOUT,
                     static_cast<long>(transport.at("connect_timeout_seconds").integer()));
    curl_easy_setopt(handle.get(), CURLOPT_TIMEOUT,
                     static_cast<long>(transport.at("request_timeout_seconds").integer()));
    curl_easy_setopt(handle.get(), CURLOPT_NOSIGNAL, 1L);
    curl_easy_setopt(handle.get(), CURLOPT_USERAGENT, "GADGET/0.1 public-evidence-retrieval");
    curl_easy_setopt(handle.get(), CURLOPT_WRITEFUNCTION, body);
    curl_easy_setopt(handle.get(), CURLOPT_WRITEDATA, &out);
    curl_easy_setopt(handle.get(), CURLOPT_HEADERFUNCTION, header);
    curl_easy_setopt(handle.get(), CURLOPT_HEADERDATA, &out);
    const auto status = curl_easy_perform(handle.get());
    if (status != CURLE_OK)
        throw std::runtime_error(std::string("database transport failed: ") +
                                 curl_easy_strerror(status));
    long http = 0;
    curl_easy_getinfo(handle.get(), CURLINFO_RESPONSE_CODE, &http);
    require(http == 200 && !out.body.empty() && !out.overflow,
            "database did not return a complete HTTP 200 response within declared resources");
    (void)Json::parse(out.body);
    return out;
}
#endif
} // namespace
NativeRecord capture_database_snapshot(const fs::path &directory) {
    const auto path=directory/"snapshot.rfcxs";
    if(fs::is_regular_file(path))return NativeRecord::read(path);
    const auto original=read_text(directory/"snapshot.json");const auto manifest=Json::parse(original);
    using N=NativeRecordValue;std::vector<RawByteSource> sources;
    sources.emplace_back(1,std::span<const std::byte>{reinterpret_cast<const std::byte *>(original.data()),original.size()});
    N::Array artifacts;const auto records=manifest.at("records");std::set<std::string> paths;
    const auto add=[&](const Json &record) {
        const auto name=record.at("path").text();if(!paths.insert(name).second)return;
        const auto bytes=scientific_detail::local_artifact(directory,record);const auto key=sources.size()+1;
        sources.emplace_back(key,std::span<const std::byte>{reinterpret_cast<const std::byte *>(bytes.data()),bytes.size()});
        artifacts.elements.emplace_back(N::Object{{{"path",N(name)},{"bytes",N(N::ByteSpan{key,0,bytes.size(),false})}}});
    };
    for(std::size_t i=0;i<records.size();++i) {
        const auto row=records.at(i);add(row);auto headers=Json::object();
        headers.set("path",row.at("headers_path"));headers.set("sha256",row.at("headers_sha256"));add(headers);
    }
    return NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.database.snapshot.native.v1"))},
        {"snapshot",import_record_view(manifest)},{"declaration",N(N::ByteSpan{1,0,original.size(),false})},
        {"artifacts",N(std::move(artifacts))}}}),new_native_record_namespace(),std::move(sources));
}
Json plan_database_fetch(const Json &request) {
    const auto spec = Json::parse(request.dump());
    members(spec, {"schema", "id", "requests"}, {"transport"});
    require(spec.at("schema").text() == "gadget.database.fetch.v1", "unsupported database request");
    (void)token(spec.at("id"));
    const auto requests = spec.at("requests");
    source_array(requests, "database requests");
    require(requests.size(), "database requests require explicit records");
    auto transport = Json::object();
    const auto supplied = spec.has("transport") ? spec.at("transport") : Json::object();
    members(supplied, {},
            {"max_body_bytes", "max_header_bytes", "max_total_body_bytes",
             "connect_timeout_seconds", "request_timeout_seconds"});
    for (const auto *name : {"max_body_bytes", "max_header_bytes", "max_total_body_bytes"}) {
        const auto n = supplied.has(name) ? count(supplied.at(name), source_byte_limit(), name)
                                          : source_byte_limit();
        require(n > 0, "transport byte resources must be positive");
        transport.set(name, cardinality(n));
    }
    for (const auto &[name, fallback] : std::array<std::pair<const char *, std::size_t>, 2>{
             {{"connect_timeout_seconds", 10}, {"request_timeout_seconds", 30}}}) {
        const auto n = supplied.has(name)
                           ? count(supplied.at(name), std::numeric_limits<long>::max(), name)
                           : fallback;
        require(n > 0, "transport timeout must be positive");
        transport.set(name, cardinality(n));
    }
    auto result = Json::object(), addresses = Json::array();
    std::set<std::string> ids;
    for (std::size_t i = 0; i < requests.size(); ++i) {
        require(ids.insert(requests.at(i).at("id").text()).second,
                "duplicate database request identity");
        auto item = Json::object();
        item.set("request", requests.at(i));
        item.set("url", Json::string(url(requests.at(i))));
        addresses.push(item);
    }
    result.set("schema", Json::string("gadget.database.request.v1"));
    result.set("request", spec);
    result.set("requests", addresses);
    result.set("transport", transport);
    return seal(result, "request_sha256");
}
Json database_capabilities() {
    auto result = Json::object(), list = Json::array();
    for (const auto &p : providers) {
        auto item = Json::object();
        item.set("provider", Json::string(p.id));
        item.set("documentation", Json::string(p.docs));
        item.set("data_terms", Json::string(p.license));
        item.set("retrieval_compiled", Json::boolean(GADGET_HAS_PUBLIC_DATABASES));
        list.push(item);
    }
    result.set("providers", list);
    result.set("drugbank", Json::string("licensed access adapter is not implemented"));
    result.set("scope",
               Json::string(
                   "explicit HTTPS requests with declared transport resources; original JSON and "
                   "response headers; no "
                   "automatic pagination, whole-database completeness or clinical recommendation"));
    return result;
}
Json fetch_databases(const Json &request, const fs::path &destination) {
#if GADGET_HAS_PUBLIC_DATABASES
    const auto admitted = plan_database_fetch(request);
    const auto spec = admitted.at("request");
    const auto id = token(spec.at("id"));
    const auto requests = spec.at("requests");
    const auto transport = admitted.at("transport");
    const auto total_limit = count(transport.at("max_total_body_bytes"), source_byte_limit(),
                                   "aggregate response bytes");
    require(!fs::exists(destination), "database snapshot destination already exists");
    fs::create_directories(destination.parent_path());
    require(fs::create_directory(destination), "database snapshot creation conflict");
    try {
        auto manifest = Json::object(), records = Json::array();
        manifest.set("schema", Json::string("gadget.database.snapshot.v1"));
        manifest.set("request", spec);
        std::size_t total = 0;
        for (std::size_t i = 0; i < requests.size(); ++i) {
            // Below PubChem's five requests/s policy, also gentle to other providers.
            if (i)
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
            const auto address = admitted.at("requests").at(i).at("url").text();
            const auto received = retrieve(address, transport, total_limit - total);
            total = rfc::machine::checked_size_add(total, received.body.size());
            const auto rid = requests.at(i).at("id").text();
            write_text(destination / (rid + ".json"), received.body);
            write_text(destination / (rid + ".headers.txt"), received.headers);
            auto record = Json::object();
            record.set("request", requests.at(i));
            record.set("url", Json::string(address));
            record.set("path", Json::string(rid + ".json"));
            record.set("sha256", Json::string("sha256:" + sha256(received.body)));
            record.set("headers_path", Json::string(rid + ".headers.txt"));
            record.set("headers_sha256", Json::string("sha256:" + sha256(received.headers)));
            record.set("received_unix_seconds", Json::string(std::to_string(std::time(nullptr))));
            record.set("scope", Json::string("one explicit API response; page/cutoff selection "
                                             "retained; arrival time is not experimental time"));
            for (const auto &p : providers)
                if (requests.at(i).at("provider").text() == p.id) {
                    record.set("documentation", Json::string(p.docs));
                    record.set("data_terms", Json::string(p.license));
                }
            records.push(record);
        }
        manifest.set("records", records);
        manifest.set("id", Json::string(id));
        manifest = seal(manifest, "snapshot_sha256");
        write_text(destination / "snapshot.json", manifest.dump() + "\n");
        capture_database_snapshot(destination).write(destination/"snapshot.rfcxs");
        return manifest;
    } catch (...) {
        fs::remove_all(destination);
        throw;
    }
#else
    (void)request;
    (void)destination;
    throw std::invalid_argument(
        "public database retrieval is disabled; enable GADGET_ENABLE_PUBLIC_DATABASES or import a "
        "frozen snapshot offline");
#endif
}
} // namespace gadget
