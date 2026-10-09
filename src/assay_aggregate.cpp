#include "genomic_detail.hpp"
#include <rfc/database/native_tensor_window.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/numeric/prepared_sparse_map.hpp>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::assay_detail {
namespace {
Json aggregate_impl(const Matrix &matrix, const Json &declaration, const Json &request,
                    const rfc::runtime::CpuResourceCoordinator::Reservation &reservation,
                    const std::vector<NativeAssayPartition> &publications = {}) {
    const bool native = !publications.empty();
    members(request,
            {"schema", "operation", "assay", "feature_groups", "observation_groups", "map_method",
             "map_evidence", "support_policy"},
            {"count_policy", "workers"});
    require(request.at("support_policy").text() == "require_all_nonzero_inputs",
            "aggregation support policy");
    (void)text(request.at("map_method"), source_byte_limit());
    (void)text(request.at("map_evidence"), source_byte_limit());
    const auto workers =
        request.has("workers")
            ? count(request.at("workers"), reservation.workers(), "aggregation workers")
            : 1;
    require(workers > 0, "aggregation workers must be positive");
    require(workers <= reservation.workers(), "aggregation exceeds the native CPU reservation");
    const auto map = [](const Axis &axis, const Json &groups) {
        require(json_object_is_type(groups.raw(), json_type_array) && groups.size() > 0,
                "aggregation requires declared groups");
        std::map<std::string, std::size_t> source;
        for (std::size_t i = 0; i < axis.ids.size(); ++i)
            source.emplace(axis.ids[i], i);
        rfc::numeric::ExactMatrix m(groups.size(), axis.ids.size());
        std::set<std::string> targets;
        for (std::size_t i = 0; i < groups.size(); ++i) {
            const auto g = groups.at(i);
            members(g, {"id", "members"});
            require(targets.insert(text(g.at("id"), source_byte_limit())).second,
                    "duplicate aggregate group");
            require(g.at("members").size() > 0, "empty aggregate group");
            std::set<std::string> seen;
            for (std::size_t j = 0; j < g.at("members").size(); ++j) {
                const auto member = g.at("members").at(j);
                members(member, {"id", "coefficient"});
                const auto id = text(member.at("id"), source_byte_limit());
                require(source.contains(id) && seen.insert(id).second,
                        "unknown/duplicate aggregate member");
                m.set(i, source.at(id),
                      exact_rational(text(member.at("coefficient"), source_byte_limit())));
            }
        }
        return m;
    };
    const auto feature_groups = request.at("feature_groups"),
               observation_groups = request.at("observation_groups");
    const auto feature_map = map(matrix.features, feature_groups),
               observation_map = map(matrix.observations, observation_groups);
    const bool count_sum = request.has("count_policy");
    if (count_sum) {
        require(request.at("count_policy").text() == "disjoint_count_sum" &&
                    is_count_quantity(declaration.at("quantity").text()),
                "count sums require a read-count or UMI-count source");
        // Admission of measurement meaning, not a second aggregation kernel.
        // The declared native maps must preserve feature identities and count
        // disjoint observation units once. General linear maps keep their own
        // derived quantity below.
        std::vector<std::size_t> feature_inputs(feature_groups.size()),
            observation_inputs(observation_groups.size());
        std::set<std::size_t> used_observations;
        for (const auto &[index, coefficient] : feature_map.entries()) {
            require(coefficient == 1 && feature_groups.at(index.first).at("id").text() ==
                                            matrix.features.ids.at(index.second),
                    "count sums must preserve feature identity with unit coefficients");
            ++feature_inputs.at(index.first);
        }
        require(std::all_of(feature_inputs.begin(), feature_inputs.end(),
                            [](auto n) { return n == 1; }),
                "count sums require one source feature per output feature");
        for (const auto &[index, coefficient] : observation_map.entries()) {
            require(coefficient == 1 && used_observations.insert(index.second).second,
                    "count sums require disjoint observation groups with unit coefficients");
            ++observation_inputs.at(index.first);
        }
        require(std::all_of(observation_inputs.begin(), observation_inputs.end(),
                            [](auto n) { return n > 0; }),
                "empty count-sum observation group");
        // read_matrix() already admits the count integer domain. Reuse
        // that source contract rather than rescanning the count values here.
    }
    // Restrict numerical inputs to the native maps' nonzero dependencies.
    // These are execution coordinates, not a new source or inferred quotient.
    struct SelectedMap {
        rfc::numeric::ExactMatrix action;
        std::vector<std::size_t> inputs;
        std::vector<std::vector<std::size_t>> dependencies;
        explicit SelectedMap(const rfc::numeric::ExactMatrix &full) {
            std::set<std::size_t> used;
            for (const auto &[index, coefficient] : full.entries()) {
                (void)coefficient;
                used.insert(index.second);
            }
            inputs.assign(used.begin(), used.end());
            rfc::numeric::ExactMatrix inclusion(full.cols(), inputs.size());
            for (std::size_t i = 0; i < inputs.size(); ++i)
                inclusion.set(inputs[i], i, 1);
            action = full.multiply(inclusion);
            dependencies.resize(full.rows());
            for (const auto &[index, coefficient] : action.entries()) {
                (void)coefficient;
                dependencies[index.first].push_back(index.second);
            }
        }
    };
    const SelectedMap features_selected(feature_map), observations_selected(observation_map);
    (void)rfc::machine::checked_size_mul(feature_groups.size(), observation_groups.size());
    const auto feature_first_size =
        rfc::machine::checked_size_mul(feature_groups.size(), observations_selected.inputs.size());
    const auto observation_first_size =
        rfc::machine::checked_size_mul(observation_groups.size(), features_selected.inputs.size());
    const auto whole_sum = [](const rfc::numeric::ExactMatrix &m) {
        if (m.rows() != 1 || m.nonzeros() != m.cols() || !m.cols())
            return false;
        return std::all_of(m.entries().begin(), m.entries().end(),
                           [](const auto &e) { return e.second == 1; });
    };
    const bool native_features = native && declaration.at("carrier").text() == "features";
    const bool sparse_reduction =
        native && whole_sum(native_features ? feature_map : observation_map);
    const bool features_first =
        native ? native_features : feature_first_size <= observation_first_size;
    const bool fused_sparse_maps = native && !sparse_reduction;
    const auto &first = features_first ? features_selected : observations_selected;
    const auto &second = features_first ? observations_selected : features_selected;
    std::vector<rfc::storage::PinnedNativeAccession> pins;
    std::vector<std::shared_ptr<const rfc::numeric::SparseRationalTensor>> field_tensors;
    std::vector<std::size_t> packet_offsets(publications.size());
    std::vector<std::pair<std::size_t, std::size_t>> source_coordinates;
    auto bindings = Json::array();
    std::size_t required_packets = 0;
    const bool count_quantity = is_count_quantity(declaration.at("quantity").text());
    if (native) {
        const auto &axis = native_features ? matrix.features : matrix.observations;
        const auto &field_axis = native_features ? matrix.observations : matrix.features;
        const auto &field_inputs =
            native_features ? observations_selected.inputs : features_selected.inputs;
        const auto &carrier_inputs =
            native_features ? features_selected.inputs : observations_selected.inputs;
        std::size_t covered = 0;
        std::set<std::uint64_t> namespaces;
        for (const auto &packet : publications) {
            require(packet.owner && packet.begin == covered && packet.end > packet.begin &&
                        packet.end <= axis.ids.size(),
                    "native assay packet range/coverage");
            require(namespaces.insert(packet.owner->native_history()->source_namespace()).second,
                    "native assay packet namespace repeated");
            covered = packet.end;
        }
        require(covered == axis.ids.size(), "native assay packets must cover the original axis");
        std::vector<unsigned char> required(publications.size());
        for (std::size_t i = 0, p = 0; i < carrier_inputs.size(); ++i) {
            const auto coordinate = carrier_inputs[i];
            while (coordinate >= publications.at(p).end)
                ++p;
            source_coordinates.emplace_back(p, coordinate - publications[p].begin);
            required[p] = !field_inputs.empty();
        }
        const auto dictionary_hash = "sha256:" + sha256(axis.bytes);
        for (std::size_t p = 0; p < publications.size(); ++p) {
            const auto &packet = publications[p];
            const auto endpoint = packet.owner->endpoint();
            const auto meaning = packet.owner->declaration();
            const auto extent = packet.end - packet.begin;
            std::map<std::uint64_t, std::string> labels;
            for (std::size_t i = 0; i < meaning.at("cells").size(); ++i) {
                const auto cell = meaning.at("cells").at(i);
                if (cell.at("grade").integer() == 0)
                    labels.emplace(key(cell.at("key")), cell.at("id").text());
            }
            const auto &arena = endpoint->source.grade(rfc::core::Grade{0});
            require(arena.cell_count() == extent, "native reduction carrier size");
            for (std::size_t i = 0; i < extent; ++i)
                require(labels.at(arena.identities.at(i)) == axis.ids[packet.begin + i],
                        "native reduction source identity/order");
            auto binding = Json::object();
            binding.set("partition", cardinality(p));
            binding.set("begin", cardinality(packet.begin));
            binding.set("end", cardinality(packet.end));
            binding.set("axis_dictionary_sha256", Json::string(dictionary_hash));
            binding.set("source_digest", Json::string(packet.owner->digest()));
            binding.set(
                "source_namespace",
                Json::string(std::to_string(packet.owner->native_history()->source_namespace())));
            binding.set("cutoff",
                        Json::string(std::to_string(packet.owner->native_history()->cutoff())));
            binding.set("required_native_fields",
                        cardinality(required[p] ? field_inputs.size() : 0));
            bindings.push(binding);
            packet_offsets[p] = field_tensors.size();
            if (!required[p])
                continue;
            ++required_packets;
            std::map<std::string, Json> fields;
            for (std::size_t i = 0; i < meaning.at("fields").size(); ++i) {
                const auto field = meaning.at("fields").at(i);
                fields.emplace(field.at("id").text(), field);
            }
            std::map<std::uint64_t, const rfc::storage::TemporalTensor *> tensors;
            std::map<std::uint64_t, const rfc::storage::TemporalCarrier *> carriers;
            for (const auto &t : endpoint->tensors)
                tensors.emplace(t.key, &t);
            for (const auto &c : endpoint->carriers)
                carriers.emplace(c.key, &c);
            for (const auto input : field_inputs) {
                const auto id = "all." + field_axis.ids.at(input);
                require(fields.contains(id), "native reduction field identity");
                const auto field = fields.at(id);
                require(field.at("value_type").text() == "quantitative",
                        "native aggregation requires quantitative fields");
                for (const auto *name : {"quantity", "units", "channel", "method"})
                    require(field.at(name).text() == declaration.at(name).text(),
                            "native reduction measurement meaning");
                const auto pin = packet.owner->field_accession(id);
                const auto tensor_owner = pin.tensor();
                const auto &tensor = *tensor_owner;
                require(tensor.domain().extents() == std::vector<std::size_t>{extent},
                        "native reduction field domain");
                const auto tensor_key = key(field.at("key"));
                require(tensors.contains(tensor_key), "native reduction tensor key");
                const auto &t = *tensors.at(tensor_key);
                require(t.output && t.interpretation == field_axis.ids.at(input),
                        "native reduction field accession");
                require(carriers.contains(*t.output), "native reduction carrier key");
                const auto &carrier = *carriers.at(*t.output);
                require(carrier.kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                            carrier.grade.value == 0,
                        "native reduction source grade");
                pins.push_back(pin);
                field_tensors.push_back(tensor_owner);
            }
        }
    }
    rfc::numeric::ExactMatrix packet_sum(1, publications.size());
    if (sparse_reduction && publications.size() > 1)
        for (std::size_t p = 0; p < publications.size(); ++p)
            packet_sum.set(0, p, 1);
    // A retained source field supplies all selected coordinates. Validate the
    // original count domain once per required field, also when the action reads
    // across fields rather than along the stored carrier. Arithmetic below stays
    // in the canonical matrix/window actions.
    rfc::runtime::CpuTileExecutor domain_executor({workers, workers, field_tensors.size()});
    const auto domain_report = domain_executor.run_tiles(
        count_quantity ? field_tensors.size() : 0,
        [&](std::size_t s, std::size_t, const rfc::runtime::CpuTileStop &stop) {
            for (const auto &entry : field_tensors.at(s)->values().entries()) {
                stop.check();
                require(entry.value >= 0 && entry.value.get_den() == 1,
                        "native count aggregation integer domain");
            }
        });
    std::vector<std::size_t> sparse_visits(second.inputs.size());
    auto out = Json::object(), features = Json::array(), observations = Json::array(),
         rows = Json::array();
    for (std::size_t i = 0; i < feature_groups.size(); ++i)
        features.push(feature_groups.at(i).at("id"));
    for (std::size_t i = 0; i < observation_groups.size(); ++i)
        observations.push(observation_groups.at(i).at("id"));
    std::vector<Json> output_slices(first.action.rows());
    rfc::runtime::CpuTileRunReport first_report, second_report;
    if (fused_sparse_maps) {
        // The declared carrier map acts on actual sparse field owners. The
        // cross-field map accumulates its outputs without retaining a dense
        // feature-by-observation intermediate or a copied source complex.
        const auto &carrier_map = native_features ? feature_map : observation_map;
        std::vector<rfc::numeric::ExactMatrix> local_maps;
        for (const auto &packet : publications)
            local_maps.emplace_back(carrier_map.rows(), packet.end - packet.begin);
        std::vector<std::size_t> ends;
        for (const auto &packet : publications)
            ends.push_back(packet.end);
        for (const auto &[index, coefficient] : carrier_map.entries()) {
            const auto packet = static_cast<std::size_t>(
                std::upper_bound(ends.begin(), ends.end(), index.second) - ends.begin());
            local_maps.at(packet).set(index.first, index.second - publications.at(packet).begin,
                                      coefficient);
        }
        std::vector<rfc::numeric::PreparedSparseExactMap> carrier_actions;
        for (auto &map : local_maps)
            carrier_actions.emplace_back(std::move(map));
        const rfc::numeric::PreparedSparseExactMap field_action(second.action);
        std::vector<rfc::numeric::PreparedSparseExactMap::FieldAccumulator> scratch;
        for (std::size_t slot = 0; slot < std::min(workers, second.inputs.size()); ++slot)
            scratch.push_back(field_action.field_accumulator(carrier_map.rows()));
        rfc::runtime::CpuTileExecutor field_executor({workers, workers, second.inputs.size()});
        first_report = field_executor.run_tiles(
            second.inputs.size(),
            [&](std::size_t field, std::size_t slot, const rfc::runtime::CpuTileStop &stop) {
                const rfc::numeric::SparseMapCheck check = [&] { stop.check(); };
                if (publications.size() == 1 && carrier_actions[0].matrix().nonzeros()) {
                    const auto &tensor = *field_tensors.at(packet_offsets[0] + field);
                    sparse_visits[field] = tensor.stored_nonzeros();
                    scratch.at(slot).add_column(field, carrier_actions[0].apply(tensor, check),
                                                check);
                    return;
                }
                rfc::numeric::SparseExactFieldSum total(carrier_map.rows());
                for (std::size_t packet = 0; packet < publications.size(); ++packet) {
                    if (!carrier_actions[packet].matrix().nonzeros())
                        continue;
                    const auto &tensor = *field_tensors.at(packet_offsets[packet] + field);
                    sparse_visits[field] = rfc::machine::checked_size_add(sparse_visits[field],
                                                                          tensor.stored_nonzeros());
                    total.add(carrier_actions[packet].apply(tensor, check), Rational{1}, check);
                }
                scratch.at(slot).add_column(field, total.finish(check), check);
            });
        auto total = field_action.field_accumulator(carrier_map.rows());
        for (const auto &part : scratch)
            total.merge(part);
        std::vector<rfc::numeric::SparseRationalTensor> fields;
        fields.reserve(second.action.rows());
        for (std::size_t column = 0; column < second.action.rows(); ++column)
            fields.push_back(total.finish_column(column));
        rfc::runtime::CpuTileExecutor output_executor({workers, workers, first.action.rows()});
        second_report = output_executor.run_tiles(
            first.action.rows(),
            [&](std::size_t row, std::size_t, const rfc::runtime::CpuTileStop &stop) {
                auto slice = Json::array();
                for (const auto &field : fields) {
                    stop.check();
                    const auto value = field.state(row);
                    slice.push(Json::string(value ? value->get_str() : "NA"));
                }
                output_slices[row] = std::move(slice);
            });
    } else {
        std::vector<Vector> transformed(second.inputs.size());
        std::vector<std::vector<unsigned char>> observed(second.inputs.size());
        rfc::runtime::CpuTileExecutor first_executor({workers, workers, second.inputs.size()});
        first_report = first_executor.run_tiles(
            second.inputs.size(),
            [&](std::size_t s, std::size_t, const rfc::runtime::CpuTileStop &stop) {
                stop.check();
                if (sparse_reduction) {
                    Vector sums(publications.size(), 0);
                    bool complete = true;
                    for (std::size_t p = 0; p < publications.size(); ++p) {
                        const auto &pin = pins.at(packet_offsets[p] + s);
                        const rfc::runtime::WorkAvailability availability{pin.source_namespace(),
                                                                          pin.sequence()};
                        const auto lease = source_execution_lease("assay-native-unit-reduction");
                        const auto window = rfc::database::NativeTensorWindow::prepare(
                            pin, {0}, {publications[p].end - publications[p].begin}, availability,
                            lease);
                        const auto summary = window.summarize(availability, lease);
                        const bool supported = summary.supported_count() ==
                                               publications[p].end - publications[p].begin;
                        sums[p] = supported && summary.sum() ? *summary.sum() : Rational{0};
                        complete = complete && supported;
                        sparse_visits[s] = rfc::machine::checked_size_add(
                            sparse_visits[s], summary.visited_entries());
                        stop.check();
                    }
                    transformed[s] =
                        publications.size() == 1 ? std::move(sums) : packet_sum.apply(sums);
                    observed[s] = std::vector<unsigned char>{static_cast<unsigned char>(complete)};
                    stop.check();
                    return;
                }
                Vector values(first.inputs.size(), 0);
                std::vector<unsigned char> supported(first.inputs.size(), false);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    const auto observation = features_first ? second.inputs[s] : first.inputs[i];
                    const auto feature = features_first ? first.inputs[i] : second.inputs[s];
                    std::optional<Rational> v;
                    if (native) {
                        const auto [packet, coordinate] = source_coordinates.at(
                            native_features ? (features_first ? i : s) : (features_first ? s : i));
                        const auto field =
                            native_features ? (features_first ? s : i) : (features_first ? i : s);
                        v = field_tensors.at(packet_offsets[packet] + field)->state(coordinate);
                    } else {
                        v = matrix.columns.at(observation).state(feature);
                    }
                    if (v) {
                        values[i] = *v;
                        supported[i] = true;
                    }
                }
                transformed[s] = first.action.apply(values);
                auto &complete = observed[s];
                complete.assign(first.action.rows(), true);
                for (std::size_t f = 0; f < complete.size(); ++f)
                    for (auto input : first.dependencies[f])
                        if (!supported[input])
                            complete[f] = false;
                stop.check();
            });
        rfc::runtime::CpuTileExecutor second_executor({workers, workers, first.action.rows()});
        second_report = second_executor.run_tiles(
            first.action.rows(),
            [&](std::size_t f, std::size_t, const rfc::runtime::CpuTileStop &stop) {
                stop.check();
                Vector values(second.inputs.size(), 0);
                for (std::size_t s = 0; s < values.size(); ++s)
                    if (observed[s][f])
                        values[s] = transformed[s][f];
                const auto output = second.action.apply(values);
                std::vector<unsigned char> complete(second.action.rows(), true);
                for (std::size_t s = 0; s < complete.size(); ++s)
                    for (auto input : second.dependencies[s])
                        if (!observed[input][f])
                            complete[s] = false;
                auto slice = Json::array();
                for (std::size_t s = 0; s < output.size(); ++s)
                    slice.push(Json::string(complete[s] ? output[s].get_str() : "NA"));
                output_slices[f] = std::move(slice);
                stop.check();
            });
    }
    for (std::size_t f = 0; f < feature_groups.size(); ++f) {
        if (features_first)
            rows.push(output_slices[f]);
        else {
            auto row = Json::array();
            for (std::size_t s = 0; s < observation_groups.size(); ++s)
                row.push(output_slices[s].at(f));
            rows.push(row);
        }
    }
    const auto &feature_report = fused_sparse_maps || features_first ? first_report : second_report;
    const auto &observation_report =
        fused_sparse_maps || !features_first ? first_report : second_report;
    out.set("feature_ids", features);
    out.set("observation_ids", observations);
    out.set("rows", rows);
    out.set("quantity", Json::string(count_sum ? declaration.at("quantity").text()
                                               : "declared_linear_aggregate_of_" +
                                                     declaration.at("quantity").text()));
    if (count_sum) {
        out.set("count_policy", request.at("count_policy"));
        out.set("count_contract_verified", Json::boolean(true));
    }
    out.set("units", declaration.at("units"));
    out.set("feature_map_nonzeros", cardinality(feature_map.nonzeros()));
    out.set("observation_map_nonzeros", cardinality(observation_map.nonzeros()));
    auto execution = Json::object();
    execution.set("owner", Json::string("rfc::runtime::CpuTileExecutor"));
    execution.set("map_owner", Json::string("rfc::numeric::ExactMatrix"));
    execution.set("reservation_scope",
                  Json::string("source loading, map preparation and both joined stages"));
    execution.set("requested_workers", cardinality(workers));
    execution.set("first_map", Json::string(features_first ? "features" : "observations"));
    execution.set("feature_stage_worker_slots", cardinality(feature_report.worker_slots));
    execution.set("observation_stage_worker_slots", cardinality(observation_report.worker_slots));
    execution.set("feature_stage_completed_tiles", cardinality(feature_report.completed_tiles));
    execution.set("observation_stage_completed_tiles",
                  cardinality(observation_report.completed_tiles));
    execution.set("required_source_features", cardinality(features_selected.inputs.size()));
    execution.set("required_source_observations", cardinality(observations_selected.inputs.size()));
    execution.set("intermediate_coordinates",
                  cardinality(fused_sparse_maps ? 0
                                                : (features_first ? feature_first_size
                                                                  : observation_first_size)));
    if (native) {
        execution.set("source_owner", Json::string("rfc::storage::PinnedNativeAccession"));
        execution.set("source_owner_reused", Json::boolean(true));
        execution.set("raw_matrix_reparsed", Json::boolean(false));
        execution.set("required_native_fields", cardinality(field_tensors.size()));
        execution.set("count_domain_completed_fields", cardinality(domain_report.completed_tiles));
        execution.set("count_domain_worker_slots", cardinality(domain_report.worker_slots));
        execution.set("reservation_scope",
                      Json::string("source binding, count admission and joined native actions"));
        execution.set("source_partition_count", cardinality(publications.size()));
        execution.set("required_native_partitions", cardinality(required_packets));
        execution.set("source_bindings", bindings);
        if (publications.size() == 1)
            execution.set("source_digest", Json::string(publications[0].owner->digest()));
    }
    if (fused_sparse_maps) {
        execution.set("sparse_map_owner", Json::string("rfc::numeric::PreparedSparseExactMap"));
        execution.set("output_worker_slots", cardinality(second_report.worker_slots));
        execution.set("output_completed_tiles", cardinality(second_report.completed_tiles));
        std::size_t visited = 0;
        for (const auto count : sparse_visits)
            visited = rfc::machine::checked_size_add(visited, count);
        execution.set("bound_sparse_entries", cardinality(visited));
    }
    if (sparse_reduction) {
        execution.set("reduction_owner", Json::string("rfc::database::NativeTensorWindow"));
        if (publications.size() > 1)
            execution.set("reduction_composition_owner", Json::string("rfc::numeric::ExactMatrix"));
        std::size_t visited = 0;
        for (auto n : sparse_visits)
            visited = rfc::machine::checked_size_add(visited, n);
        execution.set("visited_sparse_entries", cardinality(visited));
    }
    out.set("execution", execution);
    out.set("scope",
            Json::string("Exact declared linear feature and observation maps; unsupported nonzero "
                         "input dependencies propagate absence. No learned grouping, implicit "
                         "averaging, imputation or transfer of topology is performed."));
    return out;
}
} // namespace
Json aggregate_analysis(const Matrix &matrix, const Json &declaration, const Json &request,
                        const rfc::runtime::CpuResourceCoordinator::Reservation &reservation) {
    return aggregate_impl(matrix, declaration, request, reservation);
}
Json native_aggregate_analysis(
    const Axis &features, const Axis &observations,
    const std::vector<NativeAssayPartition> &partitions, const Json &declaration,
    const Json &request, const rfc::runtime::CpuResourceCoordinator::Reservation &reservation) {
    require(!partitions.empty(), "native aggregation requires admitted source owners");
    Matrix axes;
    axes.features = features;
    axes.observations = observations;
    return aggregate_impl(axes, declaration, request, reservation, partitions);
}
} // namespace gadget::assay_detail
