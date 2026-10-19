#pragma once
#include "../scientific_result_fields.hpp"
#include <gadget/io.hpp>
#include <stop_token>

namespace gadget::statistics::pathway_detail {
struct Source {
    std::string id, release, identifier_namespace, sha256;
    std::string mapping_id, mapping_release, mapping_origin, feature_namespace, mapping_sha256;
    std::string_view gmt, mapping;
};
struct Parameters {
    std::string analysis_sha256, gene_sets_sha256, feature_namespace;
    Json features_path, feature_id_path;
    std::vector<std::size_t> background, selected;
    std::vector<std::string> set_ids;
    std::string ambiguous, unmapped, collisions, correction;
    std::size_t min_size{}, max_size{}, source_bytes{}, mapping_rows{}, membership_cells{}, sets{},
        features{}, field_bytes{}, kernel_work{}, archive_bytes{}, declared_buffer_bytes{};
};
Parameters admit_request(const Json &request);
// Native admission supplies exact identity strings selected through the actual
// NativeRecord, plus immutable GMT/mapping byte views. No numerical receipt or
// rendered molecular result is an input to this application/kernel boundary.
scientific_detail::ResultValue overrepresentation(const std::vector<std::string> &feature_ids, const Source &source,
                        const Json &request, scientific_detail::ResultFields &fields, std::stop_token stop = {});
} // namespace gadget::statistics::pathway_detail
