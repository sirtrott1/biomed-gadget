#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
#include <cstdlib>
#include <cstring>
#include <gadget/record_view.hpp>
#include <limits>
#include <rfc/core/builder.hpp>
#if GADGET_HAS_HTSLIB
#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#endif

namespace gadget::assay_detail {
#if GADGET_HAS_HTSLIB
namespace {
struct HtsCloser {
    void operator()(htsFile *p) const {
        if (p)
            hts_close(p);
    }
};
struct HeaderCloser {
    void operator()(sam_hdr_t *p) const { sam_hdr_destroy(p); }
};
struct BamCloser {
    void operator()(bam1_t *p) const { bam_destroy1(p); }
};
struct VcfHeaderCloser {
    void operator()(bcf_hdr_t *p) const { bcf_hdr_destroy(p); }
};
struct VcfCloser {
    void operator()(bcf1_t *p) const { bcf_destroy(p); }
};
struct KString {
    kstring_t value{0, 0, nullptr};
    ~KString() { std::free(value.s); }
};
// HTSlib takes ownership of this separate, already hash-checked byte image.
std::unique_ptr<htsFile, HtsCloser> open_image(const ResourceImage &r) {
    require(!r.bytes.empty(), "empty HTS source");
    auto *buffer = static_cast<char *>(std::malloc(r.bytes.size()));
    if (!buffer)
        throw std::bad_alloc();
    std::memcpy(buffer, r.bytes.data(), r.bytes.size());
    auto *stream = hopen("mem:retained", "r:", buffer, r.bytes.size());
    require(stream != nullptr,
            "HTSlib retained memory stream"); // hopen owns buffer even on failure
    auto *file = hts_hopen(stream, "retained", "r");
    if (!file) {
        hclose_abruptly(stream);
        throw std::invalid_argument("HTSlib could not decode retained source");
    }
    std::unique_ptr<htsFile, HtsCloser> owner(file);
    const auto format = r.declaration.at("format").text();
    const auto actual = hts_get_format(file)->format;
    require((format == "sam" && actual == sam) || (format == "bam" && actual == bam) ||
                (format == "cram" && actual == cram) || (format == "vcf" && actual == vcf) ||
                (format == "bcf" && actual == bcf),
            "declared HTS format differs from decoded bytes");
    if (actual == bam || actual == cram || actual == bcf)
        require(hts_check_EOF(file) == 1, "binary HTS source lacks a complete end marker");
    return owner;
}
std::map<std::string, std::uint64_t> contig_lengths(const Json &reference) {
    (void)reference_record(reference);
    std::map<std::string, std::uint64_t> lengths;
    for (std::size_t i = 0; i < reference.at("contigs").size(); ++i) {
        const auto c = reference.at("contigs").at(i);
        lengths.emplace(c.at("id").text(), static_cast<std::uint64_t>(c.at("length").integer()));
    }
    return lengths;
}
// CRAM is decoded only against a caller-supplied, complete reference image.
// A private directory holds HTSlib's required FASTA index; source paths/UR tags
// and ambient reference caches are never selected as input references.
struct ReferenceImage {
    fs::path root;
    explicit ReferenceImage(const ResourceImage &r) {
        require(!r.reference_bytes.empty() &&
                    "sha256:" + sha256(r.reference_bytes) ==
                        r.declaration.at("reference").at("sequence_sha256").text(),
                "CRAM requires the hash-bound reference FASTA");
        const auto expected = contig_lengths(r.declaration.at("reference"));
        auto s = Json::object();
        s.set("format", Json::string("fasta"));
        const auto inspected = inspect_resource(r.reference_bytes, s);
        std::map<std::string, std::uint64_t> found;
        for (std::size_t i = 0; i < inspected.at("records").size(); ++i) {
            const auto c = inspected.at("records").at(i);
            found.emplace(c.at("id").text(), static_cast<std::uint64_t>(c.at("length").integer()));
        }
        require(found == expected, "CRAM reference sequence dictionaries differ");
        root = private_directory(fs::temp_directory_path(), "gadget-reference-");
        try {
            write_text(root / "reference.fa", r.reference_bytes);
            require(fai_build((root / "reference.fa").c_str()) == 0, "reference FASTA index");
        } catch (...) {
            fs::remove_all(root);
            throw;
        }
    }
    ~ReferenceImage() {
        std::error_code ec;
        fs::remove_all(root, ec);
    }
};
Json alignment_json(const Alignment &r, std::size_t i) {
    auto out = Json::object(), blocks = Json::array(), junctions = Json::array();
    out.set("source_record_index", cardinality(i));
    out.set("query", Json::string(r.query));
    out.set("read_group", Json::string(r.group));
    out.set("flags", cardinality(r.flags));
    out.set("mapping_quality", cardinality(r.mapq));
    out.set("contig", Json::string(r.contig));
    out.set("strand", Json::string(r.strand));
    out.set("start", cardinality(r.start));
    out.set("end", cardinality(r.end));
    out.set("query_length", cardinality(r.query_length));
    for (auto [a, b] : r.blocks) {
        auto p = Json::array();
        p.push(cardinality(a));
        p.push(cardinality(b));
        blocks.push(p);
    }
    for (auto [a, b] : r.junctions) {
        auto p = Json::array();
        p.push(cardinality(a));
        p.push(cardinality(b));
        junctions.push(p);
    }
    out.set("aligned_blocks", blocks);
    out.set("splice_junctions", junctions);
    out.set("decoded_sam", Json::string(r.decoded));
    return out;
}
struct VariantCall {
    std::vector<std::optional<std::int32_t>> alleles;
    std::vector<bool> phased;
    bool complete{};
};
struct VariantRecord {
    std::size_t source_record_index{};
    std::uint64_t start{}, reference_end{};
    std::string contig, id, decoded_vcf;
    std::vector<std::string> alleles;
    std::vector<VariantCall> calls;
    bool pass{};
};
struct DecodedVariants {
    std::vector<std::string> samples;
    std::vector<VariantRecord> records;
};
DecodedVariants decode_variants(const ResourceImage &r) {
    const auto lengths = contig_lengths(r.declaration.at("reference"));
    auto file = open_image(r);
    std::unique_ptr<bcf_hdr_t, VcfHeaderCloser> header(bcf_hdr_read(file.get()));
    require(header != nullptr, "VCF/BCF header");
    for (int i = 0; i < header->nhrec; ++i) {
        const auto *record = header->hrec[i];
        if (record->type != BCF_HL_CTG)
            continue;
        std::string id, length;
        for (int k = 0; k < record->nkeys; ++k) {
            if (std::string(record->keys[k]) == "ID")
                id = record->vals[k];
            if (std::string(record->keys[k]) == "length")
                length = record->vals[k];
        }
        require(lengths.contains(id) &&
                    (length.empty() || unsigned_token(length) == lengths.at(id)),
                "variant contig header differs from declared reference");
    }
    const auto ns = bcf_hdr_nsamples(header.get());
    require(ns >= 0, "variant sample count");
    DecodedVariants decoded;
    auto &samples = decoded.samples;
    auto &records = decoded.records;
    std::set<std::string> ids;
    for (int i = 0; i < ns; ++i) {
        const std::string id(header->samples[i]);
        require(ids.insert(id).second, "duplicate variant sample");
        samples.push_back(id);
    }
    std::unique_ptr<bcf1_t, VcfCloser> record(bcf_init());
    require(record != nullptr, "variant record allocation");
    int status;
    std::size_t genotypes = 0;
    while ((status = bcf_read(file.get(), header.get(), record.get())) == 0) {
        require(record->errcode == 0 && bcf_unpack(record.get(), BCF_UN_ALL) == 0,
                "invalid VCF/BCF record");
        require(record->rid >= 0 && record->pos >= 0 && record->n_allele > 0,
                "variant coordinates/alleles");
        const std::string contig(bcf_hdr_id2name(header.get(), record->rid));
        const std::string ref(record->d.allele[0]);
        require(lengths.contains(contig) &&
                    static_cast<std::uint64_t>(record->pos) < lengths.at(contig) &&
                    ref.size() <= lengths.at(contig) - static_cast<std::uint64_t>(record->pos),
                "variant REF span outside reference");
        VariantRecord item;
        auto &alleles = item.alleles;
        auto &calls = item.calls;
        item.source_record_index = records.size();
        item.contig = contig;
        item.start = record->pos;
        item.reference_end = record->pos + ref.size();
        item.id = record->d.id ? record->d.id : ".";
        for (int i = 0; i < record->n_allele; ++i)
            alleles.emplace_back(record->d.allele[i]);
        KString formatted;
        require(vcf_format1(header.get(), record.get(), &formatted.value) == 0, "VCF rendering");
        item.decoded_vcf.assign(formatted.value.s, formatted.value.l);
        item.pass = split(item.decoded_vcf, '\t').at(6) == "PASS";
        int32_t *gt = nullptr;
        int capacity = 0;
        const int n = bcf_get_genotypes(header.get(), record.get(), &gt, &capacity);
        std::unique_ptr<int32_t, decltype(&std::free)> owner(gt, &std::free);
        require(n >= 0 || n == -1 || n == -3, "invalid genotype GT encoding/type");
        require(n <= 0 || (ns > 0 && n % ns == 0), "genotype storage shape");
        const int stride = n > 0 ? n / ns : 0;
        for (int s = 0; s < ns; ++s) {
            VariantCall call;
            auto &copies = call.alleles;
            auto &phased = call.phased;
            bool complete = n > 0;
            for (int k = 0; k < stride; ++k) {
                const auto g = gt[s * stride + k];
                if (g == bcf_int32_vector_end)
                    break;
                if (bcf_gt_is_missing(g)) {
                    copies.push_back(std::nullopt);
                    complete = false;
                } else {
                    const auto index = bcf_gt_allele(g);
                    require(index >= 0 && index < record->n_allele,
                            "GT allele outside REF/ALT dictionary");
                    copies.push_back(index);
                }
                if (k > 0)
                    phased.push_back(bcf_gt_is_phased(g));
                genotypes = rfc::machine::checked_size_add(genotypes, 1);
            }
            if (copies.size() == 0)
                complete = false;
            call.complete = complete;
            calls.push_back(std::move(call));
        }
        records.push_back(std::move(item));
    }
    require(status == -1, "truncated or invalid variant stream");
    return decoded;
}
Json variants_view(const DecodedVariants &decoded) {
    auto out = Json::object(), samples = Json::array(), records = Json::array();
    for (const auto &name : decoded.samples)
        samples.push(Json::string(name));
    for (const auto &r : decoded.records) {
        auto item = Json::object(), alleles = Json::array(), calls = Json::array();
        item.set("source_record_index", cardinality(r.source_record_index));
        item.set("contig", Json::string(r.contig));
        item.set("start", cardinality(r.start));
        item.set("reference_end", cardinality(r.reference_end));
        item.set("id", Json::string(r.id));
        item.set("decoded_vcf", Json::string(r.decoded_vcf));
        for (const auto &allele : r.alleles)
            alleles.push(Json::string(allele));
        item.set("alleles", alleles);
        for (std::size_t s = 0; s < r.calls.size(); ++s) {
            const auto &c = r.calls[s];
            auto call = Json::object(), copies = Json::array(), phased = Json::array();
            for (const auto &allele : c.alleles)
                copies.push(allele ? cardinality(*allele) : Json{});
            for (auto phase : c.phased)
                phased.push(Json::boolean(phase));
            call.set("sample", Json::string(decoded.samples[s]));
            call.set("allele_indices", copies);
            call.set("phased_separators", phased);
            call.set("complete", Json::boolean(c.complete));
            calls.push(call);
        }
        item.set("genotypes", calls);
        records.push(item);
    }
    out.set("samples", samples);
    out.set("records", records);
    out.set("record_count", cardinality(records.size()));
    return out;
}
Json read_variants(const ResourceImage &r) { return variants_view(decode_variants(r)); }

} // namespace

Alignments decode_alignments(const ResourceImage &r, bool render_records) {
    const auto lengths = contig_lengths(r.declaration.at("reference"));
    auto file = open_image(r);
    std::unique_ptr<sam_hdr_t, HeaderCloser> header(sam_hdr_read(file.get()));
    require(header != nullptr && sam_hdr_nref(header.get()) > 0, "alignment reference header");
    for (int i = 0; i < sam_hdr_nref(header.get()); ++i) {
        const std::string name(sam_hdr_tid2name(header.get(), i));
        require(lengths.contains(name) && sam_hdr_tid2len(header.get(), i) >= 0 &&
                    static_cast<std::uint64_t>(sam_hdr_tid2len(header.get(), i)) ==
                        lengths.at(name),
                "alignment @SQ differs from declared reference");
    }
    std::unique_ptr<ReferenceImage> reference;
    if (r.declaration.at("format").text() == "cram") {
        reference = std::make_unique<ReferenceImage>(r);
        require(hts_set_fai_filename(file.get(), (reference->root / "reference.fa").c_str()) == 0,
                "CRAM reference binding");
    }
    Alignments out;
    out.header = sam_hdr_str(header.get());
    std::unique_ptr<bam1_t, BamCloser> record(bam_init1());
    require(record != nullptr, "alignment record allocation");
    int status;
    std::size_t pieces = 0;
    while ((status = sam_read1(file.get(), header.get(), record.get())) >= 0) {
        const auto &c = record->core;
        require(c.tid >= -1 && c.tid < sam_hdr_nref(header.get()) && c.mtid >= -1 &&
                    c.mtid < sam_hdr_nref(header.get()) && c.l_qseq >= 0,
                "alignment reference indices/query length");
        require((c.flag & ~0xfffu) == 0 && !((c.flag & BAM_FREAD1) && (c.flag & BAM_FREAD2)) &&
                    (!(c.flag & (BAM_FREAD1 | BAM_FREAD2)) || (c.flag & BAM_FPAIRED)),
                "invalid alignment flags/mate identity");
        Alignment a;
        a.flags = c.flag;
        a.mapq = c.qual;
        a.query = bam_get_qname(record.get());
        require(!a.query.empty(), "alignment query name");
        a.query_length = c.l_qseq;
        a.strand = (c.flag & BAM_FREVERSE) ? "-" : "+";
        const auto *group = bam_aux_get(record.get(), "RG");
        if (group) {
            const auto *name = bam_aux2Z(group);
            require(name != nullptr, "RG tag must be text");
            a.group = name;
            KString rg;
            require(sam_hdr_find_line_id(header.get(), "RG", "ID", name, &rg.value) == 0,
                    "alignment RG tag lacks header dictionary");
        }
        if (render_records) {
            KString decoded;
            require(sam_format1(header.get(), record.get(), &decoded.value) >= 0, "SAM rendering");
            a.decoded.assign(decoded.value.s, decoded.value.l);
        }
        if (!(c.flag & BAM_FUNMAP)) {
            require(c.tid >= 0 && c.pos >= 0 && c.n_cigar > 0,
                    "mapped record lacks coordinate/CIGAR");
            a.contig = sam_hdr_tid2name(header.get(), c.tid);
            a.start = c.pos;
            std::uint64_t position = a.start, consumed_query = 0;
            for (std::uint32_t i = 0; i < c.n_cigar; ++i) {
                const auto cig = bam_get_cigar(record.get())[i];
                const auto op = bam_cigar_op(cig), length = bam_cigar_oplen(cig);
                require(length > 0 && op <= BAM_CDIFF, "invalid CIGAR operation");
                if (bam_cigar_type(op) & 1)
                    consumed_query += length;
                if (bam_cigar_type(op) & 2) {
                    require(position <= lengths.at(a.contig) &&
                                length <= lengths.at(a.contig) - position,
                            "alignment CIGAR outside declared contig");
                    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
                        a.blocks.emplace_back(position, position + length);
                    if (op == BAM_CREF_SKIP)
                        a.junctions.emplace_back(position, position + length);
                    position += length;
                }
                pieces = rfc::machine::checked_size_add(pieces, 1);
            }
            require(c.l_qseq == 0 || consumed_query == static_cast<std::uint64_t>(c.l_qseq),
                    "CIGAR query length differs from SEQ");
            a.end = position;
        } else {
            a.contig = "*";
            a.strand = ".";
        }
        out.records.push_back(std::move(a));
    }
    require(status == -1, "truncated or invalid alignment stream");
    // Closing the decoder precedes deleting its private reference index.
    require(hts_close(file.release()) == 0, "alignment decoder completion");
    return out;
}
std::string hts_codec_version() { return hts_version(); }
Json inspect_hts(const ResourceImage &r, bool summary) {
    const auto format = r.declaration.at("format").text();
    Json out;
    if (format == "vcf" || format == "bcf")
        out = read_variants(r);
    else {
        const auto decoded = decode_alignments(r, !summary);
        out = Json::object();
        out.set("header", Json::string(decoded.header));
        if (!summary) {
            auto records = Json::array();
            for (std::size_t i = 0; i < decoded.records.size(); ++i)
                records.push(alignment_json(decoded.records[i], i));
            out.set("records", records);
        }
        out.set("record_count", cardinality(decoded.records.size()));
    }
    out.set("format", r.declaration.at("format"));
    out.set("reference", r.declaration.at("reference"));
    out.set("source_sha256", Json::string("sha256:" + sha256(r.bytes)));
    out.set("htslib_version", Json::string(hts_version()));
    out.set("record_identity", Json::string("Source file hash and zero-based decoded record index; "
                                            "binary offsets are not text spans"));
    return out;
}
Json genotype_analysis(const ResourceImage &r, const Json &request,
                       std::optional<AssayFieldsSource> *retained) {
    members(request, {"schema", "operation", "resource", "filter_policy"});
    const auto filter = request.at("filter_policy").text();
    require(filter == "all" || filter == "pass", "variant filter policy");
    const auto decoded = decode_variants(r);
    using N = NativeRecordValue;
    using namespace rfc::storage;
    const auto integer = [](std::size_t x) { return N(rfc::numeric::Integer(std::to_string(x))); };
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("upstream_variant_record",
                                         rfc::core::BoundaryPolicy::DeclaredExact, 1, 0);
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(rfc::core::Grade{1});
    for (std::size_t s = 0; s < decoded.samples.size(); ++s)
        builder.declare_witness_identity(s);
    for (std::size_t v = 0; v < decoded.records.size(); ++v) {
        std::vector<rfc::core::BoundaryOccurrenceInput> occurrences;
        for (std::size_t s = 0; s < decoded.samples.size(); ++s)
            occurrences.emplace_back(s, true, false, Rational(1),
                                     std::vector<rfc::core::AspectValue>{});
        builder.add_grade_one(rfc::core::RelationKey{v}, type, occurrences);
    }
    TemporalEndpoint endpoint(builder.freeze());
    const auto ns = decoded.samples.size(), nr = decoded.records.size();
    endpoint.carriers.push_back({1,
                                 TemporalCarrierKind::SourceGrade,
                                 rfc::core::Grade{0},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(ns)});
    endpoint.carriers.push_back({2,
                                 TemporalCarrierKind::SourceGrade,
                                 rfc::core::Grade{1},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(nr)});
    Vector complete, ploidy, starts, ends, allele_values, phase_values;
    std::vector<std::size_t> allele_support, phase_support;
    for (const auto &record : decoded.records)
        for (const auto &call : record.calls) {
            complete.emplace_back(call.complete ? 1 : 0);
            ploidy.emplace_back(call.alleles.size());
            starts.emplace_back(allele_values.size());
            for (std::size_t k = 0; k < call.alleles.size(); ++k) {
                const auto coordinate = allele_values.size();
                allele_values.emplace_back(call.alleles[k].value_or(0));
                phase_values.emplace_back(k > 0 && call.phased[k - 1] ? 1 : 0);
                if (call.alleles[k])
                    allele_support.push_back(coordinate);
                if (k > 0)
                    phase_support.push_back(coordinate);
            }
            ends.emplace_back(allele_values.size());
        }
    endpoint.carriers.push_back({3,
                                 TemporalCarrierKind::Explicit,
                                 {},
                                 {},
                                 1,
                                 rfc::numeric::IndexDomain::linear(complete.size())});
    endpoint.carriers.push_back({4,
                                 TemporalCarrierKind::Explicit,
                                 {},
                                 {},
                                 2,
                                 rfc::numeric::IndexDomain::linear(allele_values.size())});
    SourceGradeFields store(endpoint);
    const auto column = [&](std::string name, std::uint64_t carrier, const Vector &values,
                            rfc::numeric::StructuralSupport support) {
        return store.on_carrier(
            std::move(name), carrier,
            SourceGradeFields::tensor(
                values.size(), [&](auto i) { return values[i]; }, std::move(support)),
            TemporalTensorRole::Field);
    };
    N::Object calls;
    for (const auto &[name, values] :
         std::vector<std::pair<std::string, const Vector *>>{{"complete", &complete},
                                                             {"ploidy", &ploidy},
                                                             {"slot_begin", &starts},
                                                             {"slot_end", &ends}})
        calls.members.emplace_back(name,
                                   column("call." + name, 3, *values,
                                          rfc::numeric::StructuralSupport::full(values->size())));
    calls.members.emplace_back("alleles", column("call.allele_index", 4, allele_values,
                                                 rfc::numeric::StructuralSupport::packed(
                                                     allele_values.size(), allele_support)));
    calls.members.emplace_back("phase", column("call.phased_separator", 4, phase_values,
                                               rfc::numeric::StructuralSupport::packed(
                                                   phase_values.size(), phase_support)));
    N::Array samples, records, declarations, allele_fields;
    for (const auto &name : decoded.samples)
        samples.elements.emplace_back(name);
    for (std::size_t v = 0; v < nr; ++v) {
        const auto &record = decoded.records[v];
        N::Array alleles;
        for (const auto &allele : record.alleles)
            alleles.elements.emplace_back(allele);
        records.elements.emplace_back(N::Object{{{"contig", N(record.contig)},
                                                 {"id", N(record.id)},
                                                 {"alleles", N(std::move(alleles))},
                                                 {"decoded_vcf", N(record.decoded_vcf)},
                                                 {"pass", N(record.pass)}}});
        if (filter == "pass" && !record.pass)
            continue;
        for (std::size_t a = 1; a < record.alleles.size(); ++a) {
            Vector counts(ns), observed_counts;
            std::vector<std::size_t> support;
            for (std::size_t s = 0; s < ns; ++s) {
                const auto &call = record.calls[s];
                if (!call.complete)
                    continue;
                Vector indicator;
                for (const auto &index : call.alleles)
                    indicator.emplace_back(*index == static_cast<std::int32_t>(a) ? 1 : 0);
                counts[s] = rfc::numeric::dot(indicator, Vector(indicator.size(), 1));
                support.push_back(s);
                observed_counts.push_back(counts[s]);
            }
            const auto id = "record:" + std::to_string(v) + ":allele:" + std::to_string(a);
            const auto tensor_key = store.column(
                id, 0, ns, [&](auto s) { return counts[s]; },
                rfc::numeric::StructuralSupport::packed(ns, support));
            declarations.elements.emplace_back(
                N::Object{{{"id", N(id)},
                           {"key", tensor_key},
                           {"quantity", N(std::string("called_allele_copy_number"))},
                           {"units", N(std::string("allele_copies"))},
                           {"method", N(std::string("Retained upstream allele calls"))}}});
            allele_fields.elements.emplace_back(N::Object{
                {{"id", N(id)},
                 {"key", tensor_key},
                 {"source_record_index", integer(v)},
                 {"allele_index", integer(a)},
                 {"observed_genotypes", integer(support.size())},
                 {"observed_copy_total",
                  N(rfc::numeric::dot(observed_counts, Vector(observed_counts.size(), 1)))}}});
        }
    }
    const auto source_start = store.column(
        "record.start", 1, nr, [&](auto v) { return Rational(decoded.records[v].start); });
    const auto source_end = store.column("record.reference_end", 1, nr, [&](auto v) {
        return Rational(decoded.records[v].reference_end);
    });
    N declaration(N::Object{
        {{"schema", N(std::string("gadget.genotype.fields.v1"))},
         {"fields", N(std::move(declarations))},
         {"samples", N(samples)},
         {"cell_ids", N(N::Array{{N(std::move(samples)), N{}}})},
         {"records", N(std::move(records))},
         {"calls", N(std::move(calls))},
         {"record_start", source_start},
         {"record_end", source_end},
         {"allele_fields", N(std::move(allele_fields))},
         {"request", import_record_view(request)},
         {"reference", import_record_view(r.declaration.at("reference"))},
         {"scope", N(std::string("Retained upstream genotype calls. Partial or missing GT is "
                                 "absent from the copy-number field; zero requires a complete "
                                 "call. No pathogenicity or genotype likelihood is inferred."))}}});
    if (retained) {
        retained->emplace(AssayFieldsSource{std::move(endpoint), std::move(declaration), {}});
        return Json::object();
    }
    const auto publication = NativeRecord::create_from_endpoint(
        std::move(declaration), new_native_record_namespace(), {}, std::move(endpoint));
    return genotype_fields_view(publication);
}

#else
std::string hts_codec_version() {
    throw std::invalid_argument("HTS input requires GADGET_ENABLE_HTSLIB");
}
Alignments decode_alignments(const ResourceImage &, bool) {
    throw std::invalid_argument("HTS input requires GADGET_ENABLE_HTSLIB");
}
Json inspect_hts(const ResourceImage &, bool) {
    throw std::invalid_argument("HTS input requires GADGET_ENABLE_HTSLIB");
}
Json genotype_analysis(const ResourceImage &, const Json &, std::optional<AssayFieldsSource> *) {
    throw std::invalid_argument("HTS input requires GADGET_ENABLE_HTSLIB");
}
#endif
Json genotype_fields_view(const NativeRecord &record) {
    const auto value = [&](const NativeRecord &key, std::size_t i) {
        return retained_tensor(record, key).tensor.state(i);
    };
    const auto coordinate = [&](const NativeRecord &key, std::size_t i) {
        const auto state = value(key, i);
        require(state.has_value() && *state >= 0 && state->get_den() == 1,
                "genotype coordinate column");
        return static_cast<std::size_t>(state->get_num().get_ui());
    };
    const auto samples = record.at("samples"), source_records = record.at("records"),
               calls = record.at("calls");
    auto out = Json::object(), records = Json::array(), fields = Json::array();
    out.set("samples", render_record_view(samples));
    for (std::size_t v = 0; v < source_records.size(); ++v) {
        auto item = Json::object();
        for (const auto *name : {"contig", "id", "alleles", "decoded_vcf"})
            item.set(name, render_record_view(source_records.at(v).at(name)));
        item.set("source_record_index", cardinality(v));
        item.set("start", cardinality(coordinate(record.at("record_start"), v)));
        item.set("reference_end", cardinality(coordinate(record.at("record_end"), v)));
        auto genotypes = Json::array();
        for (std::size_t s = 0; s < samples.size(); ++s) {
            const auto call_index = rfc::machine::checked_size_add(
                rfc::machine::checked_size_mul(v, samples.size()), s);
            const auto begin = coordinate(calls.at("slot_begin"), call_index),
                       end = coordinate(calls.at("slot_end"), call_index);
            auto call = Json::object(), alleles = Json::array(), phase = Json::array();
            for (std::size_t k = begin; k < end; ++k) {
                const auto allele = value(calls.at("alleles"), k);
                alleles.push(allele ? cardinality(allele->get_num().get_ui()) : Json{});
                if (k > begin) {
                    const auto separator = value(calls.at("phase"), k);
                    require(separator.has_value(), "genotype phase support");
                    phase.push(Json::boolean(*separator != 0));
                }
            }
            call.set("sample", Json::string(samples.at(s).text()));
            call.set("allele_indices", alleles);
            call.set("phased_separators", phase);
            call.set("complete",
                     Json::boolean(value(calls.at("complete"), call_index).value() != 0));
            genotypes.push(call);
        }
        item.set("genotypes", genotypes);
        records.push(item);
    }
    for (std::size_t f = 0; f < record.at("allele_fields").size(); ++f) {
        const auto field = record.at("allele_fields").at(f);
        auto item = Json::object(), values = Json::array();
        for (const auto *name : {"id", "source_record_index", "allele_index", "observed_genotypes",
                                 "observed_copy_total"})
            item.set(name, render_record_view(field.at(name)));
        const auto v = field.at("source_record_index").integer().get_ui();
        const auto &tensor = retained_tensor(record, field.at("key")).tensor;
        for (std::size_t s = 0; s < samples.size(); ++s) {
            const auto c = rfc::machine::checked_size_add(
                rfc::machine::checked_size_mul(v, samples.size()), s);
            const auto state = tensor.state(s);
            auto sample = Json::object();
            sample.set("sample", Json::string(samples.at(s).text()));
            sample.set("complete", Json::boolean(state.has_value()));
            sample.set("ploidy", cardinality(coordinate(calls.at("ploidy"), c)));
            if (state.has_value())
                sample.set("copy_number", Json::string(state->get_str()));
            values.push(sample);
        }
        item.set("values", values);
        fields.push(item);
    }
    out.set("records", records);
    out.set("record_count", cardinality(source_records.size()));
    out.set("allele_fields", fields);
    out.set("quantity", Json::string("called_allele_copy_number"));
    out.set("units", Json::string("allele_copies"));
    out.set("reference", render_record_view(record.at("reference")));
    out.set("scope", Json::string(record.at("scope").text()));
    return out;
}
} // namespace gadget::assay_detail
