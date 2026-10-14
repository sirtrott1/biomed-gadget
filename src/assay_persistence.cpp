#include "genomic_detail.hpp"
#include "scientific_result_fields.hpp"
#include <rfc/algebra/materialize.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/geometry/chain_map.hpp>
#include <rfc/topology/integral_homology.hpp>

namespace gadget::assay_detail {
NativeRecord relational_persistence_record(const ScientificStudy &study,
                                           const NativeRecord &native_request) {
    using N = NativeRecordValue;
    using V = scientific_detail::ResultValue;
    using Fields = scientific_detail::ResultFields;
    const auto request = render_record_view(native_request);
    Fields publication(study.native_publication());
    using namespace rfc;
    using numeric::ExactMatrix;
    members(request,
            {"schema", "operation", "assay", "partition", "levels", "births", "filtration_meaning",
             "evidence"},
            {"lattice_snapshots"});
    require(!request.has("lattice_snapshots") ||
                json_object_is_type(request.at("lattice_snapshots").raw(), json_type_boolean),
            "lattice snapshots must be boolean");
    const bool lattice_snapshots =
        request.has("lattice_snapshots") && request.at("lattice_snapshots").boolean();
    (void)text(request.at("filtration_meaning"), 16384);
    (void)text(request.at("evidence"), 16384);
    const auto endpoint = study.endpoint();
    const auto &source = endpoint->source;
    const auto level_values = request.at("levels"), entries = request.at("births");
    require(level_values.size() > 0, "persistence requires levels");
    Vector levels;
    for (std::size_t i = 0; i < level_values.size(); ++i) {
        levels.push_back(exact_rational(text(level_values.at(i))));
        require(i == 0 || levels[i] > levels[i - 1], "persistence levels must increase strictly");
    }
    std::size_t total = 0;
    std::vector<std::vector<std::optional<std::size_t>>> births;
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        const auto n = source.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count();
        total = machine::checked_size_add(total, n);
        births.emplace_back(n);
    }
    require(entries.size() == total, "persistence requires all source cells");
    auto coordinates = V::array();
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto e = entries.at(i);
        members(e, {"grade", "index", "level"});
        const auto k = count(e.at("grade"), source.grade_count() - 1, "persistence grade");
        require(!births[k].empty(), "persistence empty grade coordinate");
        const auto cell = count(e.at("index"), births[k].size() - 1, "persistence cell index");
        const auto level = count(e.at("level"), levels.size() - 1, "persistence birth index");
        require(!births[k][cell], "duplicate persistence cell birth");
        births[k][cell] = level;
        const auto grade = core::Grade{static_cast<std::uint16_t>(k)};
        const auto &arena = source.grade(grade);
        auto c = V::object();
        c.set("grade", V::count(k));
        c.set("index", V::count(cell));
        c.set("level", V::count(level));
        c.set("source_identity",
              V::string(std::to_string(k == 0 ? arena.identities[cell]
                                                 : arena.relation_identities[cell].value)));
        coordinates.push(c);
    }
    std::vector<algebra::MaterializedSubcomplex> stages;
    for (std::size_t t = 0; t < levels.size(); ++t) {
        algebra::SubcomplexMask mask(source);
        for (std::size_t k = 0; k < births.size(); ++k)
            for (std::size_t i = 0; i < births[k].size(); ++i) {
                require(births[k][i].has_value(), "missing persistence birth");
                if (*births[k][i] <= t)
                    mask.select({core::Grade{static_cast<std::uint16_t>(k)},
                                 static_cast<std::uint32_t>(i)});
            }
        // Includes inactive retained incidences too: the declaration, not just
        // its realized nonzero boundary, must be closed at every level.
        require(mask.boundary_closed(source),
                "persistence birth order violates retained boundary closure");
        stages.push_back(algebra::materialize_subcomplex(source, mask));
    }
    const auto boundary = [](const core::RelationalComplex &c, std::size_t k) {
        if (k == 0)
            return ExactMatrix(0, c.grade(core::Grade{0}).cell_count());
        if (k == c.grade_count())
            return ExactMatrix(c.grade(core::Grade{static_cast<std::uint16_t>(k - 1)}).cell_count(),
                               0);
        return field::materialize_boundary_matrix(c, core::Grade{static_cast<std::uint16_t>(k)});
    };
    for (const auto &stage : stages)
        for (std::size_t k = 2; k < source.grade_count(); ++k)
            require(
                boundary(stage.complex, k - 1).multiply(boundary(stage.complex, k)).nonzeros() == 0,
                "persistence requires a chain complex at every level");
    const auto matrix = [](Fields &owner, const ExactMatrix &value,
                           std::uint64_t input, std::uint64_t output) {
        auto reference = owner.matrix(value, input, output);
        reference.set("display", V::string("matrix_triplets"));
        return reference;
    };
    const auto integer_matrix = [&](Fields &owner, const numeric::IntegerMatrix &value,
                                   std::uint64_t input, std::uint64_t output) {
        ExactMatrix exact(value.rows(), value.cols());
        for (std::size_t i = 0; i < value.rows(); ++i)
            for (std::size_t j = 0; j < value.cols(); ++j)
                if (value.get(i, j) != 0) exact.set(i, j, Rational(value.get(i, j)));
        return matrix(owner, exact, input, output);
    };
    const auto integer_values = [](const auto &values) {
        auto out = V::array();
        for (const auto &value : values) out.push(V(N(value)));
        return out;
    };
    struct StageBinding {
        NativeRecord owner;
        std::vector<std::uint64_t> grades, cycles, boundaries;
        std::vector<std::uint64_t> map_axes, cycle_axes;
    };
    std::vector<StageBinding> owners;
    std::vector<std::vector<ExactMatrix>> cycles, boundaries;
    auto stage_records = V::array();
    for (std::size_t t = 0; t < stages.size(); ++t) {
        cycles.emplace_back();
        boundaries.emplace_back();
        storage::TemporalEndpoint stage_endpoint(stages[t].complex);
        const auto base = NativeRecord::create_from_endpoint(N(N::Object{{
            {"schema", N(std::string("gadget.persistence.stage.source.v1"))},
            {"original_source_namespace", N(std::to_string(study.native_publication().source_namespace()))},
            {"original_source_cutoff", N(numeric::Integer(std::to_string(study.native_publication().publication_history()->last_sequence())))},
            {"original_bytes_context", N(std::string("complete parent persistence publication"))}
        }}), new_native_record_namespace(), {}, std::move(stage_endpoint));
        Fields fields(base);
        std::vector<std::uint64_t> grade_carriers, cycle_carriers, boundary_carriers;
        for (std::size_t k = 0; k < source.grade_count(); ++k)
            grade_carriers.push_back(fields.source_grade(k));
        auto stage_view = V::object(), grades = V::array();
        stage_view.set("level", V::string(level_values.at(t).text()));
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const auto lower = boundary(stages[t].complex, k),
                       upper = boundary(stages[t].complex, k + 1);
            cycles[t].push_back(lower.kernel_frame());
            boundaries[t].push_back(upper.image_frame());
            auto g = V::object();
            std::vector<std::size_t> ambient;
            for (auto i : stages[t].to_ambient[k])
                ambient.push_back(i);
            g.set("grade", V::count(k));
            cycle_carriers.push_back(fields.axis(cycles[t][k].cols()));
            boundary_carriers.push_back(fields.axis(boundaries[t][k].cols()));
            const auto lower_carrier = k == 0 ? fields.axis(0) : grade_carriers[k - 1];
            const auto upper_carrier = k + 1 == source.grade_count() ? fields.axis(0) : grade_carriers[k + 1];
            g.set("ambient_coordinates", fields.indices(ambient, grade_carriers[k]));
            g.set("boundary", matrix(fields, lower, grade_carriers[k], lower_carrier));
            g.set("cycle_frame", matrix(fields, cycles[t][k], cycle_carriers[k], grade_carriers[k]));
            g.set("boundary_frame", matrix(fields, boundaries[t][k], boundary_carriers[k], grade_carriers[k]));
            g.set("betti", V::count(cycles[t][k].cols() - boundaries[t][k].cols()));
            auto binding = V::object(), from_ambient = V::array();
            binding.set("grade", V::count(k));
            binding.set("source_carrier", V::count(grade_carriers[k]));
            binding.set("cycle_carrier", V::count(cycle_carriers[k]));
            binding.set("boundary_carrier", V::count(boundary_carriers[k]));
            binding.set("upper_boundary", matrix(fields, upper, upper_carrier, grade_carriers[k]));
            for (const auto i : stages[t].from_ambient[k])
                from_ambient.push(i == algebra::RelationalMorphism::absent ? V{} : V::count(i));
            binding.set("from_ambient", std::move(from_ambient));
            if (lattice_snapshots) {
                const auto integral = topology::integral_homology(
                    stages[t].complex, core::Grade{static_cast<std::uint16_t>(k)});
                auto lattice = V::object(), torsion = V::array(), scales = V::array();
                for (const auto &x : integral.torsion_invariants)
                    torsion.push(V::string(x.get_str()));
                for (const auto &x : stages[t]
                                         .complex.grade(core::Grade{static_cast<std::uint16_t>(k)})
                                         .cell_scales)
                    scales.push(V::string(x.get_str()));
                lattice.set("free_rank", V::count(integral.free_rank));
                lattice.set("torsion_invariants", torsion);
                lattice.set("retained_cell_scales", scales);
                auto witness = V::object(), smith = V::object();
                const auto integral_cycles = fields.axis(integral.cycle_rank);
                const auto smith_rows = fields.axis(integral.boundary_smith.diagonal.rows());
                const auto smith_columns = fields.axis(integral.boundary_smith.diagonal.cols());
                witness.set("coefficient_domain", V::string("Z; integer coefficients retained exactly in rational tensors"));
                witness.set("chain_rank", V::count(integral.chain_rank));
                witness.set("cycle_rank", V::count(integral.cycle_rank));
                witness.set("boundary_rank", V::count(integral.boundary_rank));
                witness.set("free_rank", V::count(integral.free_rank));
                witness.set("torsion_invariants", integer_values(integral.torsion_invariants));
                witness.set("cycle_basis", integer_matrix(fields, integral.cycle_basis,
                            integral_cycles, grade_carriers[k]));
                witness.set("boundaries_in_cycle_coordinates", integer_matrix(fields,
                            integral.boundaries_in_cycle_coordinates, upper_carrier, integral_cycles));
                const auto integral_lower = topology::integral_boundary_matrix(stages[t].complex,
                            core::Grade{static_cast<std::uint16_t>(k)});
                const auto integral_upper = k + 1 == source.grade_count() ? numeric::IntegerMatrix(lower.cols(), 0)
                    : topology::integral_boundary_matrix(stages[t].complex,
                            core::Grade{static_cast<std::uint16_t>(k + 1)});
                witness.set("lower_boundary", integer_matrix(fields, integral_lower, grade_carriers[k], lower_carrier));
                witness.set("upper_boundary", integer_matrix(fields, integral_upper, upper_carrier, grade_carriers[k]));
                smith.set("diagonal", integer_matrix(fields, integral.boundary_smith.diagonal, smith_columns, smith_rows));
                smith.set("left", integer_matrix(fields, integral.boundary_smith.left, integral_cycles, smith_rows));
                smith.set("right", integer_matrix(fields, integral.boundary_smith.right, smith_columns, upper_carrier));
                smith.set("rank", V::count(integral.boundary_smith.rank));
                smith.set("invariant_factors", integer_values(integral.boundary_smith.invariant_factors));
                witness.set("boundary_smith", std::move(smith));
                binding.set("integral_homology", std::move(witness));
                g.set("integral_homology", lattice);
            }
            fields.bind(std::move(binding));
            grades.push(std::move(g));
        }
        stage_view.set("grades", std::move(grades));
        auto owner = fields.finish(std::move(stage_view));
        auto binding = V::object(), axes = V::array();
        binding.set("stage", V::count(t));
        binding.set("source_namespace", V::string(std::to_string(owner.source_namespace())));
        binding.set("cutoff", V::count(owner.publication_history()->last_sequence()));
        binding.set("publication", publication.publication(owner));
        std::vector<std::uint64_t> map_axes, cycle_axes;
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            map_axes.push_back(publication.axis(stages[t].to_ambient[k].size()));
            cycle_axes.push_back(publication.axis(cycles[t][k].cols()));
            auto axis = V::object();
            axis.set("grade", V::count(k));
            axis.set("source_carrier", V::count(grade_carriers[k]));
            axis.set("cycle_carrier", V::count(cycle_carriers[k]));
            axis.set("map_axis", V::count(map_axes[k]));
            axis.set("cycle_axis", V::count(cycle_axes[k]));
            ExactMatrix lift(source.grade({static_cast<std::uint16_t>(k)}).cell_count(), stages[t].to_ambient[k].size());
            for (std::size_t i = 0; i < stages[t].to_ambient[k].size(); ++i) lift.set(stages[t].to_ambient[k][i], i, 1);
            axis.set("ambient_lift", matrix(publication, lift, map_axes[k], publication.source_grade(k)));
            axes.push(std::move(axis));
        }
        binding.set("axes", std::move(axes));
        publication.bind(std::move(binding));
        stage_records.push(publication.publication(owner));
        owners.push_back({std::move(owner), std::move(grade_carriers), std::move(cycle_carriers),
                          std::move(boundary_carriers), std::move(map_axes), std::move(cycle_axes)});
    }
    const auto l = stages.size();
    std::vector<std::vector<std::vector<numeric::Integer>>> ranks(
        source.grade_count(),
        std::vector<std::vector<numeric::Integer>>(l, std::vector<numeric::Integer>(l)));
    auto maps = V::array();
    for (std::size_t i = 0; i < l; ++i)
        for (std::size_t j = i; j < l; ++j) {
            std::vector<geometry::LinearGradeMap> grade_maps;
            std::vector<ExactMatrix> matrices;
            for (std::size_t k = 0; k < source.grade_count(); ++k) {
                geometry::LinearGradeMap grade_map(stages[i].to_ambient[k].size());
                matrices.emplace_back(stages[j].to_ambient[k].size(),
                                      stages[i].to_ambient[k].size());
                for (std::size_t c = 0; c < grade_map.size(); ++c) {
                    const auto target = stages[j].from_ambient[k].at(stages[i].to_ambient[k][c]);
                    require(target != algebra::RelationalMorphism::absent,
                            "persistence inclusion identity");
                    grade_map[c].push_back({target, Rational{1}});
                    matrices[k].set(target, c, Rational{1});
                }
                grade_maps.push_back(std::move(grade_map));
            }
            require(geometry::RelationalChainMap(stages[i].complex, stages[j].complex,
                                                 std::move(grade_maps))
                        .is_chain_map(),
                    "native persistence chain-map certificate");
            for (std::size_t k = 0; k < source.grade_count(); ++k) {
                const auto mapped = matrices[k].multiply(cycles[i][k]);
                const auto joint = ExactMatrix::horizontal_concat(boundaries[j][k], mapped);
                ranks[k][i][j] = numeric::Integer(std::to_string(joint.rank() - boundaries[j][k].cols()));
                auto m = V::object();
                m.set("grade", V::count(k));
                m.set("from", V::count(i));
                m.set("to", V::count(j));
                m.set("chain_map_verified", V::boolean(true));
                m.set("inclusion", matrix(publication, matrices[k], owners[i].map_axes[k], owners[j].map_axes[k]));
                m.set("mapped_cycles", matrix(publication, mapped, owners[i].cycle_axes[k], owners[j].map_axes[k]));
                m.set("homology_image_rank", V(N(ranks[k][i][j])));
                auto binding = V::object();
                binding.set("grade", V::count(k));
                binding.set("from_stage", V::count(i));
                binding.set("to_stage", V::count(j));
                binding.set("joint_frame", matrix(publication, joint, publication.axis(joint.cols()), owners[j].map_axes[k]));
                binding.set("homology_image_rank", V(N(ranks[k][i][j])));
                publication.bind(std::move(binding));
                maps.push(m);
            }
        }
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        const auto size = machine::checked_size_mul(l, l);
        std::vector<std::size_t> support;
        for (std::size_t i = 0; i < l; ++i)
            for (std::size_t j = i; j < l; ++j)
                support.push_back(machine::checked_size_add(machine::checked_size_mul(i, l), j));
        auto binding = V::object();
        binding.set("grade", V::count(k));
        binding.set("level_count", V::count(l));
        binding.set("rank_invariant", Fields::staged(SourceGradeFields::tensor(size,
            [&](auto i) { return Rational(ranks[k][i / l][i % l]); },
            numeric::StructuralSupport::sparse(size, std::move(support))), {},
            publication.axis(size), "vector"));
        publication.bind(std::move(binding));
    }
    auto barcodes = V::array();
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        struct Interval {
            std::size_t birth, death;
            numeric::Integer multiplicity;
        };
        std::vector<Interval> bars;
        const auto rank = [&](std::optional<std::size_t> i, std::size_t j) {
            return !i || j == l ? numeric::Integer(0) : ranks[k][*i][j];
        };
        for (std::size_t b = 0; b < l; ++b)
            for (std::size_t d = b + 1; d <= l; ++d) {
                const auto previous = b == 0 ? std::nullopt : std::optional<std::size_t>(b - 1);
                const numeric::Integer mu = rank(b, d - 1) - rank(previous, d - 1) -
                                            rank(b, d) + rank(previous, d);
                require(mu >= 0, "negative persistence interval multiplicity");
                if (mu == 0)
                    continue;
                bars.push_back({b, d, mu});
                auto bar = V::object();
                bar.set("grade", V::count(k));
                bar.set("birth_index", V::count(b));
                bar.set("birth", V::string(level_values.at(b).text()));
                if (d != l) {
                    bar.set("death_index", V::count(d));
                    bar.set("death", V::string(level_values.at(d).text()));
                }
                bar.set("survives_final_level", V::boolean(d == l));
                bar.set("multiplicity", V(N(mu)));
                barcodes.push(bar);
            }
        for (std::size_t i = 0; i < l; ++i)
            for (std::size_t j = i; j < l; ++j) {
                numeric::Integer reconstructed = 0;
                for (const auto &bar : bars)
                    if (bar.birth <= i && j < bar.death)
                        reconstructed += bar.multiplicity;
                require(reconstructed == ranks[k][i][j], "barcode rank-invariant reconstruction");
            }
    }
    auto out = V::object();
    out.set("coefficient_field", V::string("Q"));
    out.set("metric", V::string("not required for chain homology"));
    out.set("source_digest", V::string(study.digest()));
    out.set("coordinates", coordinates);
    out.set("stages", stage_records);
    out.set("inclusions", maps);
    out.set("barcode", barcodes);
    out.set("rank_reconstruction_verified", V::boolean(true));
    out.set("endpoint_semantics", V::string("Birth included, death excluded; surviving bars are "
                                               "censored at the last supplied level"));
    out.set(
        "coefficient_semantics",
        V::string("Rational chain homology; torsion requires a separately declared lattice"));
    if (lattice_snapshots)
        out.set("lattice_snapshot_semantics",
                V::string("Free rank and torsion of the retained lattice at every level; the "
                             "barcode and inclusion ranks remain over Q"));
    auto retained = publication.finish(std::move(out));
    std::uint64_t key = 1;
    for (const auto &bytes : retained.raw_sources())
        key = machine::checked_size_add(std::max(key - 1, bytes.key()), 1);
    auto original = study.native_publication().archive_source(key);
    const auto size = original.bytes().size();
    return retained.derive(N(N::Object{{
        {"schema", N(std::string("gadget.persistence.fields.v1"))},
        {"fields", retained.unbound_value()},
        {"study_publication", N(N::ByteSpan{key, 0, size, false})}
    }}), {std::move(original)});
}
Json persistence_fields_view(const NativeRecord &record) {
    return scientific_detail::result_fields_view(record.at("fields"));
}
Json relational_persistence(const ScientificStudy &study, const Json &request) {
    return persistence_fields_view(relational_persistence_record(study, admit_record_view(request)));
}
} // namespace gadget::assay_detail
