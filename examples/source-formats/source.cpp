#include "bindingdb.hpp"
#include "gtex.hpp"
#include "pubchem_evidence.hpp"
#include <gadget/record_view.hpp>
#include <gadget/result_operations.hpp>
#include <gadget/pipeline.hpp>
#include <iostream>
#include <fstream>
#include <limits>

using namespace gadget;
using namespace gadget::examples;
namespace {
Json number(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("input coordinate exceeds the external control representation");
    return Json::integer(static_cast<std::int64_t>(value));
}
std::size_t limit(const Json &request, const char *name) {
    const auto n = request.at("limits").at(name).integer();
    if (n <= 0 || static_cast<std::uint64_t>(n) > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("source reader limit must be positive and addressable");
    return static_cast<std::size_t>(n);
}
Json text(std::string_view value) {return Json::string(std::string(value));}
void emit(const Json &value) {
    std::cout << value.dump_compact() << '\n';
    if (!std::cout) throw std::runtime_error("source input consumer closed");
}
Json strings(std::span<const std::string_view> values) {
    auto out = Json::array(); for (const auto value : values) out.push(text(value)); return out;
}
RawByteSource source(const fs::path &path, std::size_t maximum, std::uint64_t key = 1) {
    if (fs::file_size(path) > maximum) throw std::length_error("original source byte admission");
    auto bytes = RawByteSource::read(key, path);
    if (bytes.bytes().size() > maximum) throw std::length_error("original source changed size");
    return bytes;
}
std::string_view view(const RawByteSource &source) {
    const auto bytes = source.bytes();
    return {reinterpret_cast<const char *>(bytes.data()), bytes.size()};
}
void gtex(const RawByteSource &bytes, const Json &q) {
    const GtexLimits limits{limit(q,"source_bytes"),limit(q,"rows"),limit(q,"columns"),
        limit(q,"record_bytes"),limit(q,"field_bytes"),limit(q,"cells")};
    GtexDeclaration declaration;
    const auto layout=q.at("layout").text(), feature=q.at("feature").text(), measurement=q.at("measurement").text();
    if (layout == "gct12") declaration.layout=GtexLayout::Gct12;
    else if (layout == "portal_tissue13") declaration.layout=GtexLayout::PortalTissue13;
    else throw std::invalid_argument("unknown declared GCT layout");
    if (feature == "gene") declaration.feature=GtexFeature::Gene;
    else if (feature == "transcript") declaration.feature=GtexFeature::Transcript;
    else throw std::invalid_argument("unknown declared feature kind");
    if (measurement == "integer_counts") declaration.measurement=GtexMeasurement::IntegerCounts;
    else if (measurement == "estimated_counts") declaration.measurement=GtexMeasurement::EstimatedCounts;
    else if (measurement == "tpm") declaration.measurement=GtexMeasurement::Tpm;
    else throw std::invalid_argument("unknown declared measurement kind");
    const auto samples=source(q.at("sample_attributes").text(),limits.source_bytes,2);
    const auto donors=source(q.at("subject_phenotypes").text(),limits.source_bytes,3);
    const auto attributes=read_gtex_metadata(view(samples),GtexMetadataKind::SampleAttributes,limits);
    const auto phenotypes=read_gtex_metadata(view(donors),GtexMetadataKind::SubjectPhenotypes,limits);
    bool header{};
    const auto summary=visit_gtex_gct(view(bytes),declaration,limits,[&](const GtexRowView &row) {
        if (!header) {
            const auto joined=join_gtex_samples(row.header.samples,attributes,phenotypes);
            auto h=Json::object(); h.set("event",text("header")); h.set("rows",number(row.header.features));
            h.set("columns",number(row.header.samples.size()));
            h.set("sample_columns",strings(attributes.columns));h.set("donor_columns",strings(phenotypes.columns));emit(h);
            for (std::size_t i=0;i<joined.samples.size();++i) {
                const auto &sample=joined.samples[i];
                auto item=Json::object(); item.set("event",text("sample"));item.set("index",number(i));
                item.set("id",text(sample.sample));item.set("donor",text(sample.donor));item.set("tissue",text(sample.tissue));
                item.set("repeated_donor",Json::boolean(sample.repeated_donor));
                item.set("sample_fields",strings(sample.sample_attributes->fields));
                item.set("donor_fields",strings(sample.subject_phenotypes->fields));
                item.set("sample_offset",number(sample.sample_attributes->source.byte_offset));
                item.set("sample_length",number(sample.sample_attributes->source.byte_length));
                item.set("donor_offset",number(sample.subject_phenotypes->source.byte_offset));
                item.set("donor_length",number(sample.subject_phenotypes->source.byte_length));emit(item);
            }
            header=true;
        }
        auto item=Json::object(),values=Json::array();item.set("event",text("row"));item.set("index",number(row.ordinal));
        item.set("id",text(row.feature));item.set("description",text(row.description));
        item.set("offset",number(row.source.byte_offset));item.set("length",number(row.source.byte_length));
        for (const auto &value:row.values) values.push(value.observed?text(value.token):Json{});
        item.set("values",values);emit(item);
    });
    if (!header || !summary.rows) throw std::invalid_argument("empty measurement carrier");
}
void bindingdb(const fs::path &path, const Json &q) {
    BindingDBBulkLimits limits{limit(q,"source_bytes"),limit(q,"rows"),limit(q,"columns"),limit(q,"field_bytes")};
    if(q.has("allow_unterminated_final_line"))limits.allow_unterminated_final_line=q.at("allow_unterminated_final_line").boolean();
    const auto length=fs::file_size(path);
    if(length>limits.source_bytes)throw std::length_error("original source byte admission");
    limits.expected_source_bytes=static_cast<std::size_t>(length);
    std::ifstream input(path,std::ios::binary);
    if(!input)throw std::runtime_error("original source cannot be opened");
    const auto summary=visit_bindingdb_bulk(input,limits,[](const auto &){});
    input.clear();input.seekg(0);
    if(!input)throw std::runtime_error("original source cannot be rewound");
    if(!summary.rows)throw std::invalid_argument("empty original evidence carrier");
    const auto header=[&](const BindingDBBulkHeader &source_header) {
            auto h=Json::object(), endpoints=Json::array();h.set("event",text("header"));h.set("rows",number(summary.rows));
            h.set("source_columns",strings(source_header.columns));
            for(std::size_t column=0;column<source_header.columns.size();++column) {
                const auto name=source_header.columns[column];
                if(name!="Ki (nM)"&&name!="Kd (nM)"&&name!="IC50 (nM)"&&name!="EC50 (nM)")continue;
                auto endpoint=Json::object();endpoint.set("column",number(column));endpoint.set("endpoint",text(name.substr(0,name.size()-5)));
                endpoints.push(endpoint);
            }
            h.set("columns",number(endpoints.size()));
            h.set("endpoints",endpoints);
            if(source_header.chain_start) h.set("chain_start",number(*source_header.chain_start));
            emit(h);
    };
    const auto converted=visit_bindingdb_bulk(input,limits,[&](const BindingDBBulkRowView &row) {
        auto item=Json::object(),values=Json::array(),comparisons=Json::array();
        item.set("event",text("row"));item.set("index",number(row.ordinal));item.set("offset",number(row.byte_offset));item.set("length",number(row.byte_length));
        item.set("fields",strings(row.fields));
        for (const auto &affinity:row.affinities) {
            values.push(affinity.relation?text(affinity.number_token):Json{});
            comparisons.push(affinity.relation?Json::integer(static_cast<std::int64_t>(*affinity.relation)+1):Json{});
        }
        item.set("values",values);item.set("comparisons",comparisons);emit(item);
    },{},header);
    if(converted.rows!=summary.rows||fs::file_size(path)!=length)
        throw std::runtime_error("original source changed during conversion");
}
}
int main(int argc,char **argv) {
    try {
        if (argc<4 || argc>6) throw std::invalid_argument("usage: gadget_source FORMAT SOURCE REQUEST [OUTPUT [CONCISE]]");
        const std::string format=argv[1];
        if(format=="pubchem-query" || format=="pubchem-result" || format=="pubchem-pipeline") {
            if((format=="pubchem-query" && argc!=5)||(format!="pubchem-query" && argc!=6))
                throw std::invalid_argument("native evidence mode requires a request and new output");
            const EvidenceStore evidence(NativeRecord::read(argv[2]),pubchem_evidence_adapter());
            const auto request=read_record(argv[3]).native();
            const auto output=argv[argc-1];
            if(fs::exists(output))throw std::invalid_argument("native evidence output already exists");
            if(format=="pubchem-query") {
                const auto result=evidence.execute_record(request);result.write(output);
                std::cout<<render_record_view(result.at("result")).dump()<<'\n';
            } else if(format=="pubchem-result") {
                const auto analysis=read_analysis_record(argv[4],request);
                const auto result=execute_result_record(analysis,request,&evidence);result.write(output);
                std::cout<<render_record_view(result.at("result")).dump()<<'\n';
            } else {
                const auto workspace=std::make_shared<const AssayWorkspace>(argv[4]);
                const auto control=render_record_view(request);
                const auto plan=control.at("schema").text()=="gadget.pipeline.default.v1"
                    ? request.derive(import_record_view(default_pipeline(*workspace,control))) : request;
                const auto result=execute_pipeline_record(workspace,plan,output,{},{},&evidence);
                std::cout<<render_record_view(result.at("result")).dump()<<'\n';
            }
            return 0;
        }
        const auto q=Json::parse(read_text(argv[3]));
        if(format=="bindingdb"&&argc==4) {bindingdb(argv[2],q);return 0;}
        const auto bytes=source(argv[2],limit(q,"source_bytes"));
        if (format=="gtex" && argc==4) gtex(bytes,q);
        else if (format=="pubchem" && argc>=5) {
            using N=NativeRecordValue;std::vector<RawByteSource> original{bytes};
            N::Object root{{{"schema",N(std::string("gadget.pubchem.sources.native.v1"))},
                {"assay",N(N::ByteSpan{1,0,bytes.bytes().size(),false})}}};
            if (argc==6) {
                auto concise=source(argv[5],limit(q,"source_bytes"),2);
                root.members.emplace_back("concise",N(N::ByteSpan{2,0,concise.bytes().size(),false}));original.push_back(std::move(concise));
            }
            if (fs::exists(argv[4])) throw std::invalid_argument("native evidence output already exists");
            const auto input=NativeRecord::create(N(std::move(root)),new_native_record_namespace(),std::move(original));
            const auto result=admit_pubchem_evidence(input,admit_record_view(q));result.write(argv[4]);
        } else throw std::invalid_argument("unsupported source format or arguments");
        return 0;
    } catch (const std::exception &error) {std::cerr<<error.what()<<'\n';return 1;}
}
