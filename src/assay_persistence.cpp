#include "genomic_detail.hpp"
#include <rfc/algebra/materialize.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/geometry/chain_map.hpp>
#include <rfc/topology/integral_homology.hpp>

namespace gadget::assay_detail {
Json relational_persistence(const ScientificStudy &study, const Json &request) {
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
    auto coordinates = Json::array();
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
        auto c = Json::parse(e.dump());
        c.set("source_identity",
              Json::string(std::to_string(k == 0 ? arena.identities[cell]
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
    const auto matrix_json = [](const ExactMatrix &matrix) {
        auto out = Json::object(), values = Json::array();
        out.set("rows", cardinality(matrix.rows()));
        out.set("columns", cardinality(matrix.cols()));
        for (const auto &[ij, value] : matrix.entries()) {
            auto e = Json::array();
            e.push(cardinality(ij.first));
            e.push(cardinality(ij.second));
            e.push(Json::string(value.get_str()));
            values.push(e);
        }
        out.set("entries", values);
        return out;
    };
    std::vector<std::vector<ExactMatrix>> cycles, boundaries;
    auto stage_records = Json::array();
    for (std::size_t t = 0; t < stages.size(); ++t) {
        cycles.emplace_back();
        boundaries.emplace_back();
        auto s = Json::object(), grades = Json::array();
        s.set("level", level_values.at(t));
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const auto lower = boundary(stages[t].complex, k),
                       upper = boundary(stages[t].complex, k + 1);
            cycles[t].push_back(lower.kernel_frame());
            boundaries[t].push_back(upper.image_frame());
            auto g = Json::object(), ambient = Json::array();
            for (auto i : stages[t].to_ambient[k])
                ambient.push(cardinality(i));
            g.set("grade", cardinality(k));
            g.set("ambient_coordinates", ambient);
            g.set("boundary", matrix_json(lower));
            g.set("cycle_frame", matrix_json(cycles[t][k]));
            g.set("boundary_frame", matrix_json(boundaries[t][k]));
            g.set("betti", cardinality(cycles[t][k].cols() - boundaries[t][k].cols()));
            if (lattice_snapshots) {
                const auto integral = topology::integral_homology(
                    stages[t].complex, core::Grade{static_cast<std::uint16_t>(k)});
                auto lattice = Json::object(), torsion = Json::array(), scales = Json::array();
                for (const auto &x : integral.torsion_invariants)
                    torsion.push(Json::string(x.get_str()));
                for (const auto &x : stages[t]
                                         .complex.grade(core::Grade{static_cast<std::uint16_t>(k)})
                                         .cell_scales)
                    scales.push(Json::string(x.get_str()));
                lattice.set("free_rank", cardinality(integral.free_rank));
                lattice.set("torsion_invariants", torsion);
                lattice.set("retained_cell_scales", scales);
                g.set("integral_homology", lattice);
            }
            grades.push(g);
        }
        s.set("grades", grades);
        stage_records.push(s);
    }
    const auto l = stages.size();
    std::vector<std::vector<std::vector<std::int64_t>>> ranks(
        source.grade_count(),
        std::vector<std::vector<std::int64_t>>(l, std::vector<std::int64_t>(l)));
    auto maps = Json::array();
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
                ranks[k][i][j] = static_cast<std::int64_t>(joint.rank() - boundaries[j][k].cols());
                auto m = Json::object();
                m.set("grade", cardinality(k));
                m.set("from", cardinality(i));
                m.set("to", cardinality(j));
                m.set("chain_map_verified", Json::boolean(true));
                m.set("inclusion", matrix_json(matrices[k]));
                m.set("mapped_cycles", matrix_json(mapped));
                m.set("homology_image_rank", Json::integer(ranks[k][i][j]));
                maps.push(m);
            }
        }
    auto barcodes = Json::array();
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        struct Interval {
            std::size_t birth, death;
            std::int64_t multiplicity;
        };
        std::vector<Interval> bars;
        const auto rank = [&](std::int64_t i, std::size_t j) {
            return i < 0 || j == l ? std::int64_t{0} : ranks[k][i][j];
        };
        for (std::size_t b = 0; b < l; ++b)
            for (std::size_t d = b + 1; d <= l; ++d) {
                const auto mu = rank(b, d - 1) - rank(static_cast<std::int64_t>(b) - 1, d - 1) -
                                rank(b, d) + rank(static_cast<std::int64_t>(b) - 1, d);
                require(mu >= 0, "negative persistence interval multiplicity");
                if (mu == 0)
                    continue;
                bars.push_back({b, d, mu});
                auto bar = Json::object();
                bar.set("grade", cardinality(k));
                bar.set("birth_index", cardinality(b));
                bar.set("birth", level_values.at(b));
                if (d != l) {
                    bar.set("death_index", cardinality(d));
                    bar.set("death", level_values.at(d));
                }
                bar.set("survives_final_level", Json::boolean(d == l));
                bar.set("multiplicity", Json::integer(mu));
                barcodes.push(bar);
            }
        for (std::size_t i = 0; i < l; ++i)
            for (std::size_t j = i; j < l; ++j) {
                std::int64_t reconstructed = 0;
                for (const auto &bar : bars)
                    if (bar.birth <= i && j < bar.death)
                        reconstructed += bar.multiplicity;
                require(reconstructed == ranks[k][i][j], "barcode rank-invariant reconstruction");
            }
    }
    auto out = Json::object();
    out.set("coefficient_field", Json::string("Q"));
    out.set("metric", Json::string("not required for chain homology"));
    out.set("source_digest", Json::string(study.digest()));
    out.set("coordinates", coordinates);
    out.set("stages", stage_records);
    out.set("inclusions", maps);
    out.set("barcode", barcodes);
    out.set("rank_reconstruction_verified", Json::boolean(true));
    out.set("endpoint_semantics", Json::string("Birth included, death excluded; surviving bars are "
                                               "censored at the last supplied level"));
    out.set(
        "coefficient_semantics",
        Json::string("Rational chain homology; torsion requires a separately declared lattice"));
    if (lattice_snapshots)
        out.set("lattice_snapshot_semantics",
                Json::string("Free rank and torsion of the retained lattice at every level; the "
                             "barcode and inclusion ranks remain over Q"));
    return out;
}
} // namespace gadget::assay_detail
