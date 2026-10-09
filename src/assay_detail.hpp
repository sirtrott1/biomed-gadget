#pragma once
#include "scientific_detail.hpp"
#include <gadget/assay.hpp>
#include <rfc/numeric/sparse_tensor.hpp>

namespace gadget::assay_detail {
using namespace scientific_detail;
struct Axis {
    std::string bytes;
    std::vector<std::string> ids;
    Json records = Json::array();
    std::vector<std::pair<std::size_t, std::size_t>> spans;
};
struct Matrix {
    Axis features, observations;
    std::vector<rfc::numeric::SparseRationalTensor> columns;
    std::string bytes;
};
inline bool is_count_quantity(const std::string &quantity) {
    return quantity == "read_count" || quantity == "umi_count" || quantity == "fragment_count";
}
struct Interval {
    std::string contig, id, strand;
    std::string feature_type;
    std::map<std::string, std::vector<std::string>> attributes;
    std::uint64_t start{}, end{};
    Rational weight{1};
    std::size_t source_start{}, source_end{};
};
Axis read_axis(const std::string &bytes);
Matrix read_matrix(const fs::path &root, const Json &spec);
Matrix read_h5ad(const fs::path &root, const Json &spec);
std::vector<Interval> read_intervals(const std::string &bytes, const Json &resource);
struct SequenceSpan {
    std::size_t start{}, end{}, id_start{}, id_length{}, length{};
};
std::vector<SequenceSpan> read_sequence_spans(std::string_view bytes, const std::string &format);
Json inspect_resource(const std::string &bytes, const Json &resource,
                      const std::string &reference_bytes = {});
Json interval_analysis(const std::string &bytes, const Json &resource, const Json &request);
Json salmon_classes(const std::string &bytes);
Json motif_analysis(const std::string &bytes, const Json &resource, const Json &request);
Json reference_record(const Json &value);
Json span(const std::string &document, const std::string &bytes, std::size_t start,
          std::size_t end);
Json file_record(const std::string &path, const std::string &bytes);
std::vector<std::string> split(const std::string &s, char delimiter);
std::vector<std::pair<std::size_t, std::size_t>> lines(const std::string &bytes);
std::uint64_t unsigned_token(const std::string &s);
unsigned iupac_mask(char c);

} // namespace gadget::assay_detail
