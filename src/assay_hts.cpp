#include "genomic_detail.hpp"
#include <cstdlib>
#include <cstring>
#include <limits>
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
Json read_variants(const ResourceImage &r) {
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
    auto samples = Json::array(), records = Json::array();
    std::set<std::string> ids;
    for (int i = 0; i < ns; ++i) {
        const std::string id(header->samples[i]);
        require(ids.insert(id).second, "duplicate variant sample");
        samples.push(Json::string(id));
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
        auto item = Json::object(), alleles = Json::array(), calls = Json::array();
        item.set("source_record_index", cardinality(records.size()));
        item.set("contig", Json::string(contig));
        item.set("start", cardinality(record->pos));
        item.set("reference_end", cardinality(record->pos + ref.size()));
        item.set("id", Json::string(record->d.id ? record->d.id : "."));
        for (int i = 0; i < record->n_allele; ++i)
            alleles.push(Json::string(record->d.allele[i]));
        item.set("alleles", alleles);
        KString formatted;
        require(vcf_format1(header.get(), record.get(), &formatted.value) == 0, "VCF rendering");
        item.set("decoded_vcf", Json::string(std::string(formatted.value.s, formatted.value.l)));
        int32_t *gt = nullptr;
        int capacity = 0;
        const int n = bcf_get_genotypes(header.get(), record.get(), &gt, &capacity);
        std::unique_ptr<int32_t, decltype(&std::free)> owner(gt, &std::free);
        require(n >= 0 || n == -1 || n == -3, "invalid genotype GT encoding/type");
        require(n <= 0 || (ns > 0 && n % ns == 0), "genotype storage shape");
        const int stride = n > 0 ? n / ns : 0;
        for (int s = 0; s < ns; ++s) {
            auto call = Json::object(), copies = Json::array(), phased = Json::array();
            call.set("sample", samples.at(s));
            bool complete = n > 0;
            for (int k = 0; k < stride; ++k) {
                const auto g = gt[s * stride + k];
                if (g == bcf_int32_vector_end)
                    break;
                if (bcf_gt_is_missing(g)) {
                    copies.push(Json{});
                    complete = false;
                } else {
                    const auto index = bcf_gt_allele(g);
                    require(index >= 0 && index < record->n_allele,
                            "GT allele outside REF/ALT dictionary");
                    copies.push(cardinality(index));
                }
                if (k > 0)
                    phased.push(Json::boolean(bcf_gt_is_phased(g)));
                genotypes = rfc::machine::checked_size_add(genotypes, 1);
            }
            if (copies.size() == 0)
                complete = false;
            call.set("allele_indices", copies);
            call.set("phased_separators", phased);
            call.set("complete", Json::boolean(complete));
            calls.push(call);
        }
        item.set("genotypes", calls);
        records.push(item);
    }
    require(status == -1, "truncated or invalid variant stream");
    auto out = Json::object();
    out.set("samples", samples);
    out.set("records", records);
    out.set("record_count", cardinality(records.size()));
    return out;
}
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
Json genotype_analysis(const ResourceImage &r, const Json &request) {
    members(request, {"schema", "operation", "resource", "filter_policy"});
    const auto filter = request.at("filter_policy").text();
    require(filter == "all" || filter == "pass", "variant filter policy");
    auto decoded = read_variants(r), fields = Json::array();
    std::size_t entries = 0;
    for (std::size_t v = 0; v < decoded.at("records").size(); ++v) {
        const auto record = decoded.at("records").at(v);
        const auto cells = split(record.at("decoded_vcf").text(), '\t');
        if (filter == "pass" && cells.at(6) != "PASS")
            continue;
        for (std::size_t allele = 1; allele < record.at("alleles").size(); ++allele) {
            auto item = Json::object(), values = Json::array();
            item.set("id", Json::string("record:" + std::to_string(v) +
                                        ":allele:" + std::to_string(allele)));
            item.set("source_record_index", cardinality(v));
            item.set("allele_index", cardinality(allele));
            Vector supported_counts;
            std::size_t observed = 0;
            for (std::size_t s = 0; s < decoded.at("samples").size(); ++s) {
                const auto call = record.at("genotypes").at(s);
                auto value = Json::object();
                value.set("sample", call.at("sample"));
                value.set("complete", call.at("complete"));
                value.set("ploidy", cardinality(call.at("allele_indices").size()));
                if (call.at("complete").boolean()) {
                    Vector indicator;
                    for (std::size_t k = 0; k < call.at("allele_indices").size(); ++k)
                        indicator.emplace_back(call.at("allele_indices").at(k).integer() ==
                                                       static_cast<std::int64_t>(allele)
                                                   ? 1
                                                   : 0);
                    const auto count = rfc::numeric::dot(indicator, Vector(indicator.size(), 1));
                    value.set("copy_number", Json::string(count.get_str()));
                    supported_counts.push_back(count);
                    ++observed;
                }
                values.push(value);
                entries = rfc::machine::checked_size_add(entries, 1);
            }
            item.set("values", values);
            item.set("observed_genotypes", cardinality(observed));
            item.set(
                "observed_copy_total",
                Json::string(rfc::numeric::dot(supported_counts, Vector(supported_counts.size(), 1))
                                 .get_str()));
            fields.push(item);
        }
    }
    decoded.set("allele_fields", fields);
    decoded.set("quantity", Json::string("called_allele_copy_number"));
    decoded.set("units", Json::string("allele_copies"));
    decoded.set("reference", r.declaration.at("reference"));
    decoded.set("scope",
                Json::string("Retained upstream genotype calls. Partial or missing GT is absent "
                             "from the copy-number field; zero requires a complete call. No "
                             "pathogenicity or genotype likelihood is inferred."));
    return decoded;
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
Json genotype_analysis(const ResourceImage &, const Json &) {
    throw std::invalid_argument("HTS input requires GADGET_ENABLE_HTSLIB");
}
#endif
} // namespace gadget::assay_detail
