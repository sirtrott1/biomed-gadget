#include "../examples/source-formats/gtex.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

using namespace gadget::examples;
namespace {
std::size_t checks{};
void check(bool ok, const char *message) {
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void refuses(F action) {
    bool rejected{};
    try {
        action();
    } catch (const std::invalid_argument &) {
        rejected = true;
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    check(rejected, "invalid GTEx source accepted");
}
const std::string gene = "#1.2\n2\t2\nName\tDescription\tGTEX-AA-1\tGTEX-BB-2\n"
                         "ENSG0001.1\tgene one\t0\tNA\nENSG0002.4\tgene two\t12\t3\n";
const std::string sample_bytes = "SAMPID\tSMTSD\tSMRDTTL\tAGE\tunknown\n"
                                 "GTEX-AA-1\tBrain\t100\tNA\t\n"
                                 "GTEX-AA-3\tLiver\tNA\t\tNA\n"
                                 "GTEX-BB-2\tBrain\t200\t\t\n";
const std::string subject_bytes = "SUBJID\tAGE\tSEX\nGTEX-AA\t60-69\t1\nGTEX-BB\t40-49\t2\n";
std::string changed(std::string source, const std::string &from, const std::string &to) {
    const auto at = source.find(from);
    if (at == std::string::npos)
        throw std::runtime_error("test replacement absent");
    source.replace(at, from.size(), to);
    return source;
}
GtexSummary inspect(const std::string &source, GtexDeclaration declaration = {},
                    GtexLimits limits = {}) {
    return visit_gtex_gct(source, declaration, limits, [](const GtexRowView &) {});
}
std::string value_source(const std::string &value) {
    return "#1.2\n1\t1\nName\tDescription\tGTEX-AA-1\nENSG0001.1\tgene\t" + value + "\n";
}
void contracts() {
    std::size_t seen{};
    const auto summary = visit_gtex_gct(gene, {}, {}, [&](const GtexRowView &row) {
        check(row.ordinal == seen++, "row order");
        check(row.header.samples == std::vector<std::string_view>{"GTEX-AA-1", "GTEX-BB-2"},
              "sample identities");
        check(row.source.line == row.ordinal + 4, "source line");
        check(row.source.original_bytes ==
                  std::string_view(gene).substr(row.source.byte_offset, row.source.byte_length),
              "source span");
        check(row.source.original_bytes.back() == '\n', "source newline retained");
        if (!row.ordinal) {
            check(row.feature == "ENSG0001.1" && row.description == "gene one",
                  "feature description");
            check(row.values[0].observed && row.values[0].token == "0", "supported zero");
            check(!row.values[1].observed && row.values[1].token == "NA", "missing value");
            check(row.header.version_record.original_bytes == "#1.2\n", "version source record");
        }
    });
    check(summary.rows == 2 && summary.samples == 2 && summary.observed == 3 &&
              summary.missing == 1 && summary.observed_zeros == 1,
          "support census");
    visit_gtex_gct(changed(gene, "ENSG0001.1", "ENSG0001.1_PAR_Y"), {}, {},
                   [&](const GtexRowView &row) {
                       if (!row.ordinal)
                           check(row.feature == "ENSG0001.1_PAR_Y",
                                 "Y pseudoautosomal source identity collapsed");
                   });
    refuses([&] { inspect(changed(gene, "ENSG0001.1", "ENSG0001.1_PAR_X")); });
    auto crlf = gene;
    for (std::size_t at = 0; (at = crlf.find('\n', at)) != std::string::npos; at += 2)
        crlf.insert(at, "\r");
    visit_gtex_gct(crlf, {}, {}, [&](const GtexRowView &row) {
        check(row.source.original_bytes.ends_with("\r\n"), "CRLF source retained");
    });
    check(inspect(gene.substr(0, gene.size() - 1)).rows == 2, "unterminated final row");
    for (const std::string token :
         {"0", "000", "0.000", "0e999999999999999999999", "+12", "12.0", "120e-1", ".0", "1.",
          "1.0e+0", "1e20", "10000000000000000000000000000000000000000",
          "0e-999999999999999999999"})
        check(inspect(value_source(token)).observed == 1, "exact integer token refused");
    for (const std::string token :
         {"",     " ",     "-1",  "-0",   "NaN",    "Infinity",
          "nan",  "1,000", "1.2", "1e-1", "1.01e1", "1e-99999999999999999999",
          ".",    "+",     "1e",  "1e+",  "1e--1",  "1.2.3",
          "1\r2", "1 2",   " 1",  "1 ",   "N/A",    "na"})
        refuses([&] { inspect(value_source(token)); });
    for (const std::string token : {"1.2", "1e-10", ".5", "0.0001"})
        check(inspect(value_source(token),
                      {GtexLayout::Gct12, GtexFeature::Gene, GtexMeasurement::Tpm})
                      .observed == 1,
              "fractional TPM refused");
    const std::string transcript =
        "#1.2\n1\t1\ntranscript_id\tgene_id\tGTEX-AA-1\nENST0001.3\tENSG0001.1\t1.5\n";
    check(inspect(transcript,
                  {GtexLayout::Gct12, GtexFeature::Transcript, GtexMeasurement::EstimatedCounts})
                  .observed == 1,
          "estimated transcript count");
    check(inspect(transcript, {GtexLayout::Gct12, GtexFeature::Transcript, GtexMeasurement::Tpm})
                  .observed == 1,
          "transcript TPM");
    refuses([&] {
        inspect(transcript,
                {GtexLayout::Gct12, GtexFeature::Transcript, GtexMeasurement::IntegerCounts});
    });
    refuses([&] {
        inspect(gene, {GtexLayout::Gct12, GtexFeature::Gene, GtexMeasurement::EstimatedCounts});
    });
    refuses([&] {
        inspect(changed(transcript, "ENSG0001.1", "gene symbol"),
                {GtexLayout::Gct12, GtexFeature::Transcript, GtexMeasurement::EstimatedCounts});
    });
    const std::string portal = "#1.3\n2\t4\t0\t0\nid\tName\tDescription\tGTEX-AA-1\tGTEX-BB-2\n"
                               "0\tENSG0001.1\tgene one\t0\tNA\n1\tENSG0002.4\tgene two\t12\t3\n";
    const GtexDeclaration portal_declaration{GtexLayout::PortalTissue13, GtexFeature::Gene,
                                             GtexMeasurement::IntegerCounts};
    check(inspect(portal, portal_declaration).samples == 2,
          "portal annotation columns became samples");
    refuses([&] { inspect(portal); });
    for (const auto &pair :
         std::vector<std::pair<std::string, std::string>>{{"4\t0\t0", "4\t1\t0"},
                                                          {"4\t0\t0", "2\t0\t0"},
                                                          {"id\tName", "id\tgene"},
                                                          {"1\tENSG0002", "2\tENSG0002"}})
        refuses([&] { inspect(changed(portal, pair.first, pair.second), portal_declaration); });
    refuses([&] {
        inspect(portal,
                {GtexLayout::PortalTissue13, GtexFeature::Transcript, GtexMeasurement::Tpm});
    });
    for (const auto &pair : std::vector<std::pair<std::string, std::string>>{
             {"#1.2", "#1.3"},
             {"2\t2\n", "3\t2\n"},
             {"2\t2\n", "1\t2\n"},
             {"2\t2\n", "0\t2\n"},
             {"2\t2\n", "2\t0\n"},
             {"2\t2\n", "2.0\t2\n"},
             {"2\t2\n", "184467440737095516160\t2\n"},
             {"Name\tDescription", "id\tDescription"},
             {"GTEX-BB-2", "GTEX-AA-1"},
             {"ENSG0002.4", "ENSG0001.1"},
             {"ENSG0002.4", "symbol"},
             {"ENSG0002.4", "ENSG0002"},
             {"ENSG0002.4", "ENSG0002.x"},
             {"gene two\t12\t3", "gene two\t12"},
             {"gene two\t12\t3", "gene two\t12\t3\t4"}})
        refuses([&] { inspect(changed(gene, pair.first, pair.second)); });
    refuses([&] { inspect(gene + "\n"); });
    refuses([&] { inspect(gene + std::string(1, '\0')); });
    refuses([&] { inspect(""); });
    refuses([&] { visit_gtex_gct(gene, {}, {}, {}); });
    for (const auto member :
         {&GtexLimits::source_bytes, &GtexLimits::rows, &GtexLimits::columns,
          &GtexLimits::record_bytes, &GtexLimits::field_bytes, &GtexLimits::cells}) {
        GtexLimits limits;
        limits.*member = 0;
        refuses([&] { inspect(gene, {}, limits); });
    }
    GtexLimits limits;
    limits.source_bytes = gene.size() - 1;
    refuses([&] { inspect(gene, {}, limits); });
    limits = {};
    limits.rows = 1;
    refuses([&] { inspect(gene, {}, limits); });
    limits = {};
    limits.columns = 1;
    refuses([&] { inspect(gene, {}, limits); });
    limits = {};
    limits.cells = 3;
    refuses([&] { inspect(gene, {}, limits); });
    limits = {};
    limits.record_bytes = 5;
    refuses([&] { inspect(gene, {}, limits); });
    limits = {};
    limits.field_bytes = 5;
    refuses([&] { inspect(gene, {}, limits); });
    std::stop_source cancellation;
    cancellation.request_stop();
    refuses([&] {
        visit_gtex_gct(gene, {}, {}, [](const GtexRowView &) {}, cancellation.get_token());
    });
    std::stop_source mid_visit;
    std::size_t provisional{};
    refuses([&] {
        visit_gtex_gct(
            gene, {}, {},
            [&](const GtexRowView &) {
                ++provisional;
                mid_visit.request_stop();
            },
            mid_visit.get_token());
    });
    check(provisional == 1, "cancelled provisional callback count");
    std::size_t late_callbacks{};
    refuses([&] {
        visit_gtex_gct(changed(gene, "gene two\t12", "gene two\t-1"), {}, {},
                       [&](const GtexRowView &) { ++late_callbacks; });
    });
    check(late_callbacks == 1, "late invalid row callbacks must be provisional");
    const auto samples = read_gtex_metadata(sample_bytes, GtexMetadataKind::SampleAttributes, {});
    const auto subjects =
        read_gtex_metadata(subject_bytes, GtexMetadataKind::SubjectPhenotypes, {});
    check(samples.rows.at("GTEX-AA-1").fields[3] == "NA" &&
              samples.rows.at("GTEX-AA-1").fields[4].empty(),
          "NA and empty metadata preserved");
    check(subjects.rows.at("GTEX-AA").fields[1] == "60-69", "age range coerced");
    check(samples.rows.at("GTEX-AA-1").source.original_bytes == "GTEX-AA-1\tBrain\t100\tNA\t\n",
          "metadata source bytes");
    const std::vector<std::string_view> selected{"GTEX-AA-1", "GTEX-AA-3", "GTEX-BB-2"};
    const auto joined = join_gtex_samples(selected, samples, subjects);
    check(joined.samples.size() == 3 && joined.distinct_donors == 2 && joined.repeated_donors == 1,
          "repeated donor structure lost");
    check(joined.samples[0].repeated_donor && joined.samples[1].repeated_donor &&
              !joined.samples[2].repeated_donor,
          "repeated donor flags");
    check(joined.samples[0].sample_attributes->fields[3] == "NA" &&
              joined.samples[0].subject_phenotypes->fields[1] == "60-69",
          "same-name metadata overwritten");
    check(joined.samples[0].tissue == "Brain" && joined.samples[1].tissue == "Liver",
          "tissue detail");
    refuses([&] { joined.require_independent_donors(); });
    const std::vector<std::string_view> independent{"GTEX-AA-1", "GTEX-BB-2"};
    join_gtex_samples(independent, samples, subjects).require_independent_donors();
    ++checks;
    for (const std::vector<std::string_view> &invalid :
         {std::vector<std::string_view>{"GTEX-AA-1", "GTEX-AA-1"},
          {"GTEX-CC-1"},
          {"sample"},
          {"GTEX--1"},
          {"GTEX-AA-"}})
        refuses([&] { join_gtex_samples(invalid, samples, subjects); });
    refuses([&] { join_gtex_samples(independent, subjects, samples); });
    refuses([&] { join_gtex_samples(independent, samples, subjects, cancellation.get_token()); });
    for (const auto &pair : std::vector<std::pair<std::string, std::string>>{
             {"SAMPID", "ID"},
             {"SMTSD", "AGE"},
             {"GTEX-BB-2", "GTEX-AA-1"},
             {"GTEX-BB-2", ""},
             {"GTEX-BB-2\tBrain\t200\t\t", "GTEX-BB-2\tBrain\t200"}})
        refuses([&] {
            read_gtex_metadata(changed(sample_bytes, pair.first, pair.second),
                               GtexMetadataKind::SampleAttributes, {});
        });
    refuses([&] { read_gtex_metadata(sample_bytes, GtexMetadataKind::SubjectPhenotypes, {}); });
    refuses([&] {
        read_gtex_metadata(sample_bytes, GtexMetadataKind::SampleAttributes, {},
                           cancellation.get_token());
    });
    limits = {};
    limits.rows = 2;
    refuses([&] { read_gtex_metadata(sample_bytes, GtexMetadataKind::SampleAttributes, limits); });
    limits = {};
    limits.cells = 14;
    refuses([&] { read_gtex_metadata(sample_bytes, GtexMetadataKind::SampleAttributes, limits); });
    const auto no_tissue_bytes = changed(sample_bytes, "GTEX-AA-1\tBrain", "GTEX-AA-1\tNA");
    const auto no_tissue =
        read_gtex_metadata(no_tissue_bytes, GtexMetadataKind::SampleAttributes, {});
    refuses([&] { join_gtex_samples(independent, no_tissue, subjects); });
    const auto no_subject_bytes = changed(subject_bytes, "GTEX-AA", "GTEX-CC");
    const auto no_subject =
        read_gtex_metadata(no_subject_bytes, GtexMetadataKind::SubjectPhenotypes, {});
    refuses([&] { join_gtex_samples(independent, samples, no_subject); });
    std::cout << "PASS: " << checks << " GTEx source-format checks\n";
}
std::string read(const char *path) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    if (!stream || stream.tellg() < 0 ||
        stream.tellg() > static_cast<std::streamoff>(GtexLimits{}.source_bytes))
        throw std::invalid_argument("GTEx test input unavailable or exceeds byte limit");
    stream.seekg(0);
    return {std::istreambuf_iterator<char>(stream), {}};
}
void join_dump(int argc, char **argv) {
    if (argc != 5)
        throw std::invalid_argument("--join SAMPLE_IDS SAMPLE_TSV SUBJECT_TSV");
    const auto ids_source = read(argv[2]), sample_source = read(argv[3]),
               subject_source = read(argv[4]);
    std::vector<std::string_view> ids;
    for (std::size_t start = 0; start < ids_source.size();) {
        auto end = ids_source.find('\n', start);
        if (end == std::string::npos)
            end = ids_source.size();
        ids.push_back(std::string_view(ids_source).substr(start, end - start));
        start = end == ids_source.size() ? end : end + 1;
    }
    const auto samples = read_gtex_metadata(sample_source, GtexMetadataKind::SampleAttributes, {});
    const auto subjects =
        read_gtex_metadata(subject_source, GtexMetadataKind::SubjectPhenotypes, {});
    const auto joined = join_gtex_samples(ids, samples, subjects);
    bool independent = true;
    try {
        joined.require_independent_donors();
    } catch (const std::invalid_argument &) {
        independent = false;
    }
    std::cout << "J\t" << joined.samples.size() << '\t' << joined.distinct_donors << '\t'
              << joined.repeated_donors << '\t' << independent << '\n';
    for (const auto &sample : joined.samples)
        std::cout << "S\t" << sample.sample << '\t' << sample.donor << '\t' << sample.tissue << '\t'
                  << sample.repeated_donor << '\n';
}
void dump(int argc, char **argv) {
    if (argc != 8)
        throw std::invalid_argument("--dump SOURCE {gct12|portal13} {gene|transcript} "
                                    "{integer_counts|estimated_counts|tpm} SAMPLE_TSV SUBJECT_TSV");
    const std::string_view layout = argv[3], feature = argv[4], measurement = argv[5];
    if ((layout != "gct12" && layout != "portal13") ||
        (feature != "gene" && feature != "transcript") ||
        (measurement != "integer_counts" && measurement != "estimated_counts" &&
         measurement != "tpm"))
        throw std::invalid_argument("GTEx dump declaration");
    const GtexDeclaration declaration{
        layout == "gct12" ? GtexLayout::Gct12 : GtexLayout::PortalTissue13,
        feature == "gene" ? GtexFeature::Gene : GtexFeature::Transcript,
        measurement == "integer_counts"     ? GtexMeasurement::IntegerCounts
        : measurement == "estimated_counts" ? GtexMeasurement::EstimatedCounts
                                            : GtexMeasurement::Tpm};
    const auto source = read(argv[2]), sample_source = read(argv[6]),
               subject_source = read(argv[7]);
    const auto samples = read_gtex_metadata(sample_source, GtexMetadataKind::SampleAttributes, {});
    const auto subjects =
        read_gtex_metadata(subject_source, GtexMetadataKind::SubjectPhenotypes, {});
    const auto summary = visit_gtex_gct(source, declaration, {}, [&](const GtexRowView &row) {
        if (!row.ordinal) {
            std::cout << "H\t" << row.header.features << '\t' << row.header.declared_columns;
            for (const auto id : row.header.samples)
                std::cout << '\t' << id;
            std::cout << '\n' << "AH";
            for (const auto name : samples.columns)
                std::cout << '\t' << name;
            std::cout << '\n' << "PH";
            for (const auto name : subjects.columns)
                std::cout << '\t' << name;
            std::cout << '\n';
            const auto joined = join_gtex_samples(row.header.samples, samples, subjects);
            for (const auto &sample : joined.samples) {
                std::cout << "S\t" << sample.sample << '\t' << sample.donor << '\t' << sample.tissue
                          << '\t' << sample.repeated_donor << '\t'
                          << sample.sample_attributes->source.byte_offset << '\t'
                          << sample.sample_attributes->source.byte_length << '\t'
                          << sample.subject_phenotypes->source.byte_offset << '\t'
                          << sample.subject_phenotypes->source.byte_length << '\n';
                std::cout << 'A';
                for (const auto field : sample.sample_attributes->fields)
                    std::cout << '\t' << field;
                std::cout << '\n' << 'P';
                for (const auto field : sample.subject_phenotypes->fields)
                    std::cout << '\t' << field;
                std::cout << '\n';
            }
        }
        std::cout << "R\t" << row.ordinal << '\t' << row.source.line << '\t'
                  << row.source.byte_offset << '\t' << row.source.byte_length << '\t' << row.feature
                  << '\t' << row.description;
        for (const auto value : row.values)
            std::cout << '\t' << (value.observed ? '1' : '0') << ':' << value.token;
        std::cout << '\n';
    });
    std::cout << "Z\t" << summary.rows << '\t' << summary.samples << '\t' << summary.observed
              << '\t' << summary.missing << '\t' << summary.observed_zeros << '\n';
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc > 1 && std::string_view(argv[1]) == "--dump")
            dump(argc, argv);
        else if (argc > 1 && std::string_view(argv[1]) == "--join")
            join_dump(argc, argv);
        else if (argc == 1)
            contracts();
        else
            throw std::invalid_argument("unknown GTEx test command");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
