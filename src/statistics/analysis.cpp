#include <gadget/build.hpp>
#include "../scientific_result_fields.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <gadget/statistics/execution.hpp>
#include <gadget/snapshot.hpp>
#include <set>
#include <iostream>
#include <numeric>
#include <rfc/statistics/field_moments.hpp>
#include <sstream>

namespace gadget::statistics {
namespace {
std::ostringstream stream() {
    std::ostringstream out;
    out.precision(17);
    return out;
}
std::string channels(const Workspace &w, const std::string &context, std::size_t gene) {
    std::string result;
    for (char c : std::string("GET"))
        if (w.fields.at(context + "." + c)[gene] != 0) {
            if (!result.empty())
                result += '+';
            result += c;
        }
    return result;
}
Vector sparse_values(const rfc::numeric::SparseRationalVector &value) {
    Vector out;
    out.reserve(value.nonzeros());
    for (const auto &e : value.entries())
        out.push_back(e.value);
    return out;
}
Rational sparse_inner(const rfc::numeric::SparseRationalVector &a,
                      const rfc::numeric::SparseRationalVector &b) {
    // Extract common coordinates; the arithmetic is canonical numeric::dot.
    Vector x, y;
    std::size_t i = 0, j = 0;
    while (i < a.nonzeros() && j < b.nonzeros()) {
        const auto &p = a.entries()[i];
        const auto &q = b.entries()[j];
        if (p.index < q.index)
            ++i;
        else if (q.index < p.index)
            ++j;
        else {
            x.push_back(p.value);
            y.push_back(q.value);
            ++i;
            ++j;
        }
    }
    return rfc::numeric::dot(x, y);
}
Rational restricted_pairing(const rfc::approximation::DenseMatrix &gram,
                            const rfc::numeric::SparseRationalVector &left,
                            const rfc::numeric::SparseRationalVector &right) {
    // An explicit coordinate selection, not a second contraction kernel.
    // Keep native sorted coface indices and evaluate with the canonical matrix.
    rfc::approximation::DenseMatrix selected(left.nonzeros(), right.nonzeros());
    rfc::approximation::Vector x, y;
    x.reserve(left.nonzeros());
    y.reserve(right.nonzeros());
    for (const auto &entry : left.entries())
        x.push_back(entry.value.get_d());
    for (const auto &entry : right.entries())
        y.push_back(entry.value.get_d());
    for (std::size_t i = 0; i < left.nonzeros(); ++i)
        for (std::size_t j = 0; j < right.nonzeros(); ++j)
            selected(i, j) = gram(left.entries()[i].index, right.entries()[j].index);
    return rfc::numeric::dot(binary64_values(x), binary64_values(selected.apply(y)));
}
using Fields = scientific_detail::ResultFields;
using V = scientific_detail::ResultValue;
using N = NativeRecordValue;
using Matrix = rfc::numeric::ExactMatrix;
V approximate(double value) {
    if (!std::isfinite(value)) throw std::runtime_error("analysis binary64 range exceeded");
    return V(N(N::ExactNumber{Rational(value), true}));
}
Matrix exact_coefficients(const rfc::approximation::DenseMatrix &value) {
    Matrix out(value.rows(), value.cols());
    for (std::size_t i = 0; i < value.rows(); ++i)
        for (std::size_t j = 0; j < value.cols(); ++j) {
            if (!std::isfinite(value(i, j))) throw std::runtime_error("analysis matrix range exceeded");
            if (value(i, j) != 0) out.set(i, j, Rational(value(i, j)));
        }
    return out;
}
Vector complete(const rfc::numeric::SparseRationalVector &value) {
    Vector out(value.size(), Rational{0});
    for (const auto &entry : value.entries()) out.at(entry.index) = entry.value;
    return out;
}
struct AnalysisFields {
    NativeRecord base;
    Fields fields;
    std::uint64_t genes, relations, channels;
    V view{V::object()}, moments{V::array()}, solves{V::array()}, exact_checks{V::array()};
    Matrix pulls, coefficients;
    std::optional<V> system;
    explicit AnalysisFields(const Execution &op, const NativeRecord &original)
        : base(NativeRecord::create_from_endpoint(N(N::Object{{
              {"schema", N(std::string("gadget.analysis.context.source.v1"))},
              {"original_source_namespace", N(std::to_string(original.source_namespace()))},
              {"original_source_cutoff", N(rfc::numeric::Integer(std::to_string(original.publication_history()->last_sequence())))},
              {"original_bytes_context", N(std::string("complete parent analysis publication"))}}}),
              new_native_record_namespace(), {}, rfc::storage::TemporalEndpoint(op.complex()))),
          fields(base), genes(fields.source_grade(0)), relations(fields.source_grade(1)),
          channels(fields.axis(3)), pulls(op.boundary().domain_size(), op.boundary().codomain_size()),
          coefficients(pulls.rows(), pulls.cols()) {
        view.set("context", V::string(op.context));
        view.set("realization", V::string(op.realization));
        view.set("source_digest", V::string(op.source_digest()));
        view.set("source_base_namespace", V::string(std::to_string(base.source_namespace())));
        view.set("source_base_cutoff", V::count(base.publication_history()->last_sequence()));
        view.set("gene_carrier", V::count(genes));
        view.set("relation_carrier", V::count(relations));
        view.set("channel_carrier", V::count(channels));
        auto labels = V::array();
        for (char c : std::string("GET")) labels.push(V::string(std::string(1, c)));
        view.set("channels", std::move(labels));
        std::vector<std::size_t> ambient;
        for (auto id : op.source_relations) ambient.push_back(id);
        view.set("source_relations", fields.indices(ambient, relations));
        ambient.clear();
        for (std::size_t i = 0; i < pulls.cols(); ++i) ambient.push_back(i);
        view.set("source_genes", fields.indices(ambient, genes));
        view.set("boundary", fields.matrix(rfc::field::boundary_matrix(op.complex(), {1}), relations, genes));
        view.set("gene_metric", fields.matrix(op.metrics()[0].materialize(), genes, genes));
        view.set("relation_metric", fields.matrix(op.metrics()[1].materialize(), relations, relations));
        view.set("arithmetic", V::string("rational actions with explicitly lifted binary64 Gram and solve coefficients"));
        view.set("metric", V::string("M0=identity; M1=exact boundary quadrance; zero-column convention=1"));
    }
    V green(const Execution::GreenWitnesses &witness, const Vector &source, const Vector &response) {
        auto result = V::object();
        result.set("source", fields.vector(source, genes));
        result.set("pulled", fields.vector(witness.pulled, relations));
        result.set("solved", fields.vector(witness.solved, relations));
        result.set("response", fields.vector(response, genes));
        result.set("residual", fields.vector(witness.residual, genes));
        if (!system) system = fields.matrix(exact_coefficients(witness.system), relations, relations);
        result.set("system", *system);
        auto diagnostic = V::object();
        diagnostic.set("backend", V::string("rfc.approximation.DenseMatrix.solve / relation-space Green"));
        diagnostic.set("exact", V::boolean(false));
        diagnostic.set("relative_residual", approximate(witness.relative));
        diagnostic.set("absolute_residual", approximate(witness.absolute));
        diagnostic.set("tolerance", approximate(witness.tolerance));
        diagnostic.set("relation_dimension", V::count(pulls.rows()));
        diagnostic.set("metric", V::string("M0=identity; M1=exact boundary quadrance; zero-column convention=1"));
        diagnostic.set("lambda", V(N(N::ExactNumber{Rational{1}, true})));
        result.set("diagnostic", diagnostic);
        solves.push(std::move(result));
        return diagnostic;
    }
    void gene_column(std::size_t gene, const rfc::numeric::SparseRationalVector &pulled,
                     const rfc::numeric::SparseRationalVector &coefficient) {
        for (const auto &entry : pulled.entries()) pulls.set(entry.index, gene, entry.value);
        for (const auto &entry : coefficient.entries()) coefficients.set(entry.index, gene, entry.value);
    }
    NativeRecord finish(Execution &op) {
        view.set("relation_gram", fields.matrix(exact_coefficients(op.relation_gram()), relations, relations));
        view.set("moments", std::move(moments));
        view.set("green_solves", std::move(solves));
        view.set("gene_pulls", fields.matrix(pulls, genes, relations));
        view.set("gene_adjoint", fields.matrix(coefficients, genes, relations));
        view.set("exact_checks", std::move(exact_checks));
        return fields.finish(std::move(view));
    }
};
void matrix_csv(const fs::path &path, const std::array<std::array<Rational, 3>, 3> &matrix,
                bool normalized) {
    auto out = stream();
    out << "channel,G,E,T\n";
    for (std::size_t i = 0; i < 3; ++i) {
        out << "GET"[i];
        for (std::size_t j = 0; j < 3; ++j) {
            double x = matrix[i][j].get_d();
            if (normalized) {
                double denominator = std::sqrt(matrix[i][i].get_d() * matrix[j][j].get_d());
                x = denominator > 0 ? x / denominator : 0;
            }
            out << ',' << x;
        }
        out << '\n';
    }
    write_text(path, out.str());
}
} // namespace
void analyze(const Workspace &w, const fs::path &requested) {
    if (fs::exists(requested))
        throw std::invalid_argument("output destination already exists: " + requested.string());
    const auto destination =
        fs::absolute(requested).parent_path() / (requested.filename().string() + ".partial");
    fs::create_directories(destination.parent_path());
    if (!fs::create_directory(destination))
        throw std::invalid_argument("analysis staging directory already exists: " +
                                    destination.string());
    const auto started = std::chrono::steady_clock::now();
    fs::create_directories(destination / "tables");
    const auto original = capture_workspace(w);
    Fields native(original);
    const auto genes = native.source_grade(0);
    auto retained_results = V::object(), retained_contexts = V::array(), retained_cross = V::array(), retained_sweeps = V::array();
    auto native_certificates = V::array();
    auto certificates = Json::array();
    write_text(destination / "manifests/source.json", w.manifest().dump() + "\n");
    auto metrics = Json::object();
    auto native_performance = V::object();
    const auto timing = [&](std::string name, double value) {
        metrics.set(name, Json::number(value));
        native_performance.set(std::move(name), approximate(value));
    };
    for (const auto &realization : {"support", "rp"}) {
        std::cerr << "Preparing canonical " << realization << " contexts\n";
        Execution ad(w, "AD", realization), hd(w, "HD", realization);
        AnalysisFields ad_fields(ad, original), hd_fields(hd, original);
        for (auto *op : {&ad, &hd}) {
            auto &retained = op == &ad ? ad_fields : hd_fields;
            auto prefix = op->context == "AD" ? "alzheimers" : "huntingtons";
            std::string base = std::string(prefix) + "_" + realization;
            std::array<Vector, 3> fields, pulled, powers;
            for (std::size_t i = 0; i < 3; ++i) {
                fields[i] = w.fields.at(op->context + "." + "GET"[i]);
                pulled[i] = op->boundary().transpose_apply(fields[i]);
                powers[i] = op->metrics()[1].solve(pulled[i]);
            }
            auto original_fields = V::array(), original_pulls = V::array();
            for (std::size_t i = 0; i < 3; ++i) {
                original_fields.push(retained.fields.vector(fields[i], retained.genes));
                original_pulls.push(retained.fields.vector(pulled[i], retained.relations));
            }
            retained.view.set("source_fields", std::move(original_fields));
            retained.view.set("pulled_fields", std::move(original_pulls));
            for (std::size_t power = 0; power <= 4; ++power) {
                std::array<std::array<Rational, 3>, 3> matrix;
                for (std::size_t i = 0; i < 3; ++i)
                    for (std::size_t j = 0; j < 3; ++j)
                        matrix[i][j] = power == 0 ? op->hodge().metric().inner(fields[i], fields[j])
                                                  : rfc::numeric::dot(pulled[i], powers[j]);
                matrix_csv(destination / "tables" /
                               (base + "_typed_moment_p" + std::to_string(power) + ".csv"),
                           matrix, false);
                matrix_csv(destination / "tables" /
                               (base + "_typed_moment_cosine_p" + std::to_string(power) + ".csv"),
                           matrix, true);
                auto cert = Json::object();
                cert.set("operation", Json::string("typed_moment"));
                cert.set("context", Json::string(op->context));
                cert.set("realization", Json::string(realization));
                cert.set("order", Json::integer(power));
                cert.set("exact", Json::boolean(power <= 1));
                cert.set(
                    "backend",
                    Json::string(
                        power <= 1
                            ? "rfc exact boundary/metric/dot"
                            : "rfc binary64 relation Gram action with exact source and metric"));
                auto rows = Json::array();
                for (const auto &row : matrix) {
                    auto cells = Json::array();
                    for (const auto &x : row)
                        cells.push(exact_number(x));
                    rows.push(cells);
                }
                cert.set("values", rows);
                certificates.push(cert);
                auto native_cert = V::object();
                native_cert.set("operation", V::string("typed_moment"));
                native_cert.set("context", V::string(op->context));
                native_cert.set("realization", V::string(realization));
                native_cert.set("order", V::count(power));
                native_cert.set("exact", V::boolean(power <= 1));
                native_cert.set("backend", V::string(cert.at("backend").text()));
                auto native_rows = V::array();
                Matrix full(3, 3), cosine(3, 3);
                for (std::size_t i = 0; i < 3; ++i) {
                    auto row = V::array();
                    for (std::size_t j = 0; j < 3; ++j) {
                        row.push(native.scalar(matrix[i][j]));
                        full.set(i, j, matrix[i][j]);
                        const auto denominator = std::sqrt(matrix[i][i].get_d() * matrix[j][j].get_d());
                        cosine.set(i, j, Rational(denominator > 0 ? matrix[i][j].get_d() / denominator : 0));
                    }
                    native_rows.push(std::move(row));
                }
                native_cert.set("values", std::move(native_rows));
                native_certificates.push(std::move(native_cert));
                auto moment = V::object(), coefficients = V::array();
                moment.set("order", V::count(power));
                moment.set("exact", V::boolean(power <= 1));
                moment.set("matrix", retained.fields.matrix(full, retained.channels, retained.channels));
                moment.set("cosine", retained.fields.matrix(cosine, retained.channels, retained.channels));
                for (const auto &value : powers) coefficients.push(retained.fields.vector(value, retained.relations));
                moment.set("relation_coefficients", std::move(coefficients));
                retained.moments.push(std::move(moment));
                if (power > 0 && power < 4)
                    for (auto &value : powers)
                        value = op->metrics()[1].solve(
                            binary64_values(op->relation_gram().apply(approximate_values(value))));
            }
            auto diagnostics = Json::array();
            std::array<Vector, 3> responses;
            for (std::size_t i = 0; i < 3; ++i) {
                Json diagnostic;
                std::optional<Execution::GreenWitnesses> witness;
                responses[i] = op->green(fields[i], Rational{1}, 1e-9, diagnostic, &witness);
                retained.green(*witness, fields[i], responses[i]);
                diagnostic.set("channel", Json::string(std::string(1, "GET"[i])));
                diagnostics.push(diagnostic);
            }
            std::array<std::array<Rational, 3>, 3> gram;
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    gram[i][j] = rfc::numeric::dot(responses[i], responses[j]);
            Matrix response_gram(3, 3);
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j) response_gram.set(i, j, gram[i][j]);
            retained.view.set("response_gram", retained.fields.matrix(response_gram, retained.channels, retained.channels));
            matrix_csv(destination / "tables" / (base + "_green_response_gram.csv"), gram, false);
            write_text(destination / "manifests" / (base + "_green.json"),
                       diagnostics.dump() + "\n");
        }
        if (std::string(realization) == "support")
            for (const auto &source_name : {"AD_GET", "HD_GET", "UNION_GET"}) {
                Vector source(w.genes.size(), Rational{0});
                for (char channel : std::string("GET")) {
                    auto field =
                        std::string(source_name) == "UNION_GET"
                            ? union_field(w, std::string(1, channel))
                            : w.fields.at(std::string(source_name).substr(0, 2) + "." + channel);
                    for (std::size_t i = 0; i < source.size(); ++i)
                        if (field[i] != 0)
                            source[i] = 1;
                }
                Json da, dh, dr;
                std::optional<Execution::GreenWitnesses> wa, wh, wr;
                auto ga = ad.green(source, Rational{1}, 1e-9, da, &wa),
                     gh = hd.green(source, Rational{1}, 1e-9, dh, &wh);
                auto native_da = ad_fields.green(*wa, source, ga), native_dh = hd_fields.green(*wh, source, gh);
                wa.reset(); wh.reset();
                auto delta = rfc::numeric::sub(gh, ga);
                auto action = rfc::numeric::sub(hd.hodge().apply(ga), ad.hodge().apply(ga));
                const auto response = hd.green(action, Rational{1}, 1e-9, dr, &wr);
                auto native_dr = hd_fields.green(*wr, action, response);
                wr.reset();
                auto predicted = rfc::numeric::scale(Rational{-1}, response);
                auto residual = rfc::numeric::sub(delta, predicted);
                auto cert = Json::object();
                cert.set("operation", Json::string("cross_context_green_identity"));
                cert.set("source", Json::string(source_name));
                cert.set("relative_residual",
                         Json::number(norm(residual) / std::max(1e-30, norm(delta))));
                cert.set(
                    "formula",
                    Json::string("GH-GA = -lambda GH (LH-LA) GA on common ambient gene carrier"));
                cert.set("solves", Json::array());
                cert.at("solves").push(da);
                cert.at("solves").push(dh);
                cert.at("solves").push(dr);
                certificates.push(cert);
                auto native_cert = V::object();
                native_cert.set("operation", V::string("cross_context_green_identity"));
                native_cert.set("source", V::string(source_name));
                native_cert.set("relative_residual", approximate(norm(residual) / std::max(1e-30, norm(delta))));
                native_cert.set("formula", V::string(cert.at("formula").text()));
                auto native_solves = V::array();
                native_solves.push(std::move(native_da)); native_solves.push(std::move(native_dh)); native_solves.push(std::move(native_dr));
                native_cert.set("solves", std::move(native_solves));
                native_certificates.push(std::move(native_cert));
                auto cross = V::object();
                cross.set("source_name", V::string(source_name));
                cross.set("source", native.vector(source, genes));
                cross.set("response_AD", native.vector(ga, genes));
                cross.set("response_HD", native.vector(gh, genes));
                cross.set("delta", native.vector(delta, genes));
                cross.set("operator_difference_action", native.vector(action, genes));
                cross.set("predicted", native.vector(predicted, genes));
                cross.set("residual", native.vector(residual, genes));
                retained_cross.push(std::move(cross));
                std::vector<std::size_t> order(w.genes.size());
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(), [&](auto a, auto b) {
                    auto av = std::abs(delta[a].get_d()), bv = std::abs(delta[b].get_d());
                    return av == bv ? w.genes[a] < w.genes[b] : av > bv;
                });
                auto out = stream();
                out << "rank,gene,AD_GET,HD_GET,source,response_AD,response_HD,delta_response\n";
                std::size_t rank = 0;
                for (auto i : order)
                    out << ++rank << ',' << csv_cell(w.genes[i]) << ',' << channels(w, "AD", i)
                        << ',' << channels(w, "HD", i) << ',' << source[i].get_d() << ','
                        << ga[i].get_d() << ',' << gh[i].get_d() << ',' << delta[i].get_d() << '\n';
                write_text(destination / "tables" /
                               (std::string("cross_context_green_field_") + source_name + ".csv"),
                           out.str());
            }
        std::cerr << "Evaluating declared numerical " << realization
                  << " gene operator differentials\n";
        auto cross_gram = cross_relation_gram(hd, ad);
        const auto &aa = ad.relation_gram();
        const auto &hh = hd.relation_gram();
        auto out = stream();
        out << "gene,AD_GET,HD_GET,in_AD_subcomplex,in_HD_subcomplex,mu1_AD,mu2_AD,mu1_HD,mu2_HD,"
               "cross_AD_HD,delta_L_squared,delta_L_norm\n";
        std::vector<std::pair<Rational, std::uint32_t>> ranking;
        ranking.reserve(w.genes.size());
        std::vector<double> first_ad(w.genes.size()), first_hd(w.genes.size());
        std::array<Vector, 8> sweep;
        for (auto &column : sweep) column.resize(w.genes.size());
        auto phase = std::chrono::steady_clock::now();
        for (std::uint32_t i = 0; i < w.genes.size(); ++i) {
            auto basis = rfc::numeric::SparseRationalVector::basis(w.genes.size(), i);
            auto ad_pulled = ad.boundary().transpose_apply(basis);
            auto hd_pulled = hd.boundary().transpose_apply(basis);
            auto ad_coeff = ad.metrics()[1].solve(ad_pulled.view());
            auto hd_coeff = hd.metrics()[1].solve(hd_pulled.view());
            auto qa = restricted_pairing(aa, ad_coeff, ad_coeff);
            auto qh = restricted_pairing(hh, hd_coeff, hd_coeff);
            auto cross = restricted_pairing(cross_gram, hd_coeff, ad_coeff);
            ad_fields.gene_column(i, ad_pulled, ad_coeff);
            hd_fields.gene_column(i, hd_pulled, hd_coeff);
            Rational squared = qa + qh - 2 * cross;
            sweep[6][i] = squared;
            squared.canonicalize();
            if (squared < 0) {
                if (std::abs(squared.get_d()) > 1e-10 * std::max(1.0, Rational(qa + qh).get_d()))
                    throw std::runtime_error(
                        "negative numerical squared differential beyond roundoff budget");
                squared = 0;
            }
            ranking.emplace_back(squared, i);
            auto ma = sparse_inner(ad_pulled, ad_coeff), mh = sparse_inner(hd_pulled, hd_coeff);
            sweep[0][i] = ma; sweep[1][i] = qa; sweep[2][i] = mh; sweep[3][i] = qh;
            sweep[4][i] = cross; sweep[5][i] = squared; sweep[7][i] = Rational(std::sqrt(squared.get_d()));
            first_ad[i] = ma.get_d();
            first_hd[i] = mh.get_d();
            out << csv_cell(w.genes[i]) << ',' << channels(w, "AD", i) << ','
                << channels(w, "HD", i) << ','
                << w.contexts.at("AD")->contains({rfc::core::Grade{0}, i}) << ','
                << w.contexts.at("HD")->contains({rfc::core::Grade{0}, i}) << ',' << ma.get_d()
                << ',' << qa.get_d() << ',' << mh.get_d() << ',' << qh.get_d() << ','
                << cross.get_d() << ',' << squared.get_d() << ',' << std::sqrt(squared.get_d())
                << '\n';
            if ((i + 1) % 2000 == 0)
                std::cerr << realization << " differential " << i + 1 << '/' << w.genes.size()
                          << "\n";
        }
        write_text(destination / "tables" /
                       (std::string("gene_role_differential_") + realization + ".csv"),
                   out.str());
        timing(
            std::string(realization) + "_sweep_seconds",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count());
        const auto exact_checks_started = std::chrono::steady_clock::now();
        std::sort(ranking.begin(), ranking.end(), [&](const auto &a, const auto &b) {
            return a.first == b.first ? w.genes[a.second] < w.genes[b.second] : a.first > b.first;
        });
        auto top = Json::array();
        auto cofaces = stream();
        cofaces << "gene_rank,gene,context,channel_rank,external_pag_id,retained_id,relation_name,"
                   "arity,Bdagger_coefficient\n";
        for (std::size_t rank = 0; rank < std::min<std::size_t>(30, ranking.size()); ++rank) {
            auto [squared, gene] = ranking[rank];
            auto basis = rfc::numeric::SparseRationalVector::basis(w.genes.size(), gene);
            auto exact_ad = ad.hodge().apply(basis.view()),
                 exact_hd = hd.hodge().apply(basis.view());
            auto exact_a = sparse_values(exact_ad), exact_h = sparse_values(exact_hd);
            Rational verified = rfc::numeric::dot(exact_a, exact_a) +
                                rfc::numeric::dot(exact_h, exact_h) -
                                2 * sparse_inner(exact_ad, exact_hd);
            verified.canonicalize();
            double error =
                std::abs(Rational(verified - squared).get_d()) / std::max(1.0, verified.get_d());
            if (error > 1e-10)
                throw std::runtime_error("exact top-gene differential check failed");
            auto row = Json::object();
            row.set("rank", Json::integer(rank + 1));
            row.set("gene", Json::string(w.genes[gene]));
            row.set("delta_L_squared", exact_number(verified));
            row.set("delta_L_norm", Json::number(std::sqrt(verified.get_d())));
            row.set("exact", Json::boolean(true));
            row.set("numerical_relative_error", Json::number(error));
            row.set("provenance", w.provenance(gene));
            top.push(row);
            for (auto *op : {&ad, &hd}) {
                auto &retained = op == &ad ? ad_fields : hd_fields;
                auto exact = V::object();
                exact.set("gene", V::count(gene));
                exact.set("rank", V::count(rank + 1));
                exact.set("hodge_action", retained.fields.vector(complete(op == &ad ? exact_ad : exact_hd), retained.genes));
                exact.set("verified_squared_difference", retained.fields.scalar(verified));
                exact.set("relative_error", approximate(error));
                retained.exact_checks.push(std::move(exact));
            }
            for (auto *op : {&ad, &hd}) {
                auto coefficient =
                    op->boundary().adjoint(basis.view(), op->metrics()[0], op->metrics()[1]);
                auto entries = coefficient.entries();
                std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) {
                    Rational av = abs(a.value), bv = abs(b.value);
                    return av == bv ? a.index < b.index : av > bv;
                });
                for (std::size_t k = 0; k < std::min<std::size_t>(12, entries.size()); ++k) {
                    const auto &entry = entries[k];
                    auto id = op->source_relations.at(entry.index);
                    const auto &assertion = w.assertions[id];
                    std::string name;
                    if (assertion.observations.contains(op->context))
                        name = assertion.observations.at(op->context).at("Name").text();
                    cofaces << rank + 1 << ',' << csv_cell(w.genes[gene]) << ',' << op->context
                            << ',' << k + 1 << ',' << csv_cell(assertion.external) << ','
                            << csv_cell(assertion.id) << ',' << csv_cell(name) << ','
                            << assertion.participants.size() << ',' << entry.value.get_d() << '\n';
                }
            }
        }
        write_text(
            destination / "tables" /
                (std::string("top_context_differential_") + realization + "_coface_channels.csv"),
            cofaces.str());
        auto rank_csv = stream();
        rank_csv << "rank,gene,AD_GET,HD_GET,mu1_AD,mu1_HD,delta_L_norm\n";
        for (std::size_t rank = 0; rank < std::min<std::size_t>(150, ranking.size()); ++rank) {
            auto [squared, gene] = ranking[rank];
            rank_csv << rank + 1 << ',' << csv_cell(w.genes[gene]) << ',' << channels(w, "AD", gene)
                     << ',' << channels(w, "HD", gene) << ',' << first_ad[gene] << ','
                     << first_hd[gene] << ',' << std::sqrt(squared.get_d()) << '\n';
        }
        write_text(destination / "tables" /
                       (std::string("top_context_differential_") + realization + ".csv"),
                   rank_csv.str());
        write_text(destination / "manifests" / (std::string("top_") + realization + ".json"),
                   top.dump() + "\n");
        timing(std::string(realization) + "_exact_checks_and_explanations_seconds",
                    std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                                               exact_checks_started)
                                     .count());
        auto retained_sweep = V::object();
        retained_sweep.set("realization", V::string(realization));
        const std::array<std::string, 8> names{"mu1_AD", "mu2_AD", "mu1_HD", "mu2_HD", "cross_AD_HD", "delta_L_squared", "unclamped_delta_L_squared", "delta_L_norm"};
        for (std::size_t c = 0; c < names.size(); ++c) retained_sweep.set(names[c], native.vector(sweep[c], genes));
        const auto ad_owner = ad_fields.finish(ad), hd_owner = hd_fields.finish(hd);
        retained_contexts.push(native.publication(ad_owner)); retained_contexts.push(native.publication(hd_owner));
        auto binding = V::object();
        binding.set("realization", V::string(realization));
        binding.set("left_context", native.publication(hd_owner)); binding.set("right_context", native.publication(ad_owner));
        const auto left = native.axis(cross_gram.rows()), right = native.axis(cross_gram.cols());
        binding.set("left_relation_axis", V::count(left)); binding.set("right_relation_axis", V::count(right));
        binding.set("left_namespace", V::string(std::to_string(hd_owner.source_namespace())));
        binding.set("right_namespace", V::string(std::to_string(ad_owner.source_namespace())));
        binding.set("left_cutoff", V::count(hd_owner.publication_history()->last_sequence()));
        binding.set("right_cutoff", V::count(ad_owner.publication_history()->last_sequence()));
        binding.set("left_context_relation_carrier", V::count(hd_fields.relations));
        binding.set("right_context_relation_carrier", V::count(ad_fields.relations));
        std::vector<std::size_t> left_coordinates, right_coordinates;
        for (std::size_t i = 0; i < cross_gram.rows(); ++i) left_coordinates.push_back(i);
        for (std::size_t i = 0; i < cross_gram.cols(); ++i) right_coordinates.push_back(i);
        binding.set("left_coordinates", native.indices(left_coordinates, left));
        binding.set("right_coordinates", native.indices(right_coordinates, right));
        binding.set("cross_relation_gram", native.matrix(exact_coefficients(cross_gram), right, left));
        native.bind(std::move(binding));
        retained_sweeps.push(std::move(retained_sweep));
        timing(
            std::string(realization) + "_differential_seconds",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count());
    }
    auto summary = Json::object();
    summary.set("schema", Json::string("gadget.analysis.v1"));
    summary.set("source_digest", Json::string(w.source_digest));
    if (w.snapshot_origin.has("study_digest"))
        summary.set("snapshot", w.snapshot_origin);
    if (w.store_origin.has("root"))
        summary.set("store", w.store_origin);
    summary.set("rexgraph_headers_sha256", Json::string(gadget::build::rexgraph_sha256()));
    summary.set("gadget_sources_sha256", Json::string(gadget::build::source_sha256()));
    summary.set(
        "numerical_contract",
        Json::string("p0/p1 moments exact; higher moments and full gene sweep use declared "
                     "binary64 relation Gram actions; top 30 per realization independently checked "
                     "by canonical exact sparse Hodge with relative error budget 1e-10; every "
                     "Green solve checks its canonical source residual <=1e-9"));
    summary.set("audit", w.audit);
    summary.set("certificates", certificates);
    timing(
        "analysis_seconds",
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
    summary.set("performance", metrics);
    summary.set("interpretation",
                Json::string("Operator differences and GET-channel coupling describe this selected "
                             "source. Disease-selected PAGER inputs are not independent validation "
                             "or evidence of causality. No statistical significance is inferred."));
    write_text(destination / "manifests/analysis.json", summary.dump() + "\n");
    auto identity = stream();
    identity << "external_pag_id,retained_id,in_AD,in_HD,arity,content_hash\n";
    for (const auto &assertion : w.assertions)
        identity << csv_cell(assertion.external) << ',' << csv_cell(assertion.id) << ','
                 << assertion.contexts.contains("AD") << ',' << assertion.contexts.contains("HD")
                 << ',' << assertion.participants.size() << ','
                 << assertion.id.substr(assertion.id.find('@') + 1) << '\n';
    write_text(destination / "tables/primary_relation_identity_audit.csv", identity.str());
    auto completion = Json::object();
    completion.set("schema", Json::string("gadget.artifacts.v1"));
    completion.set("source_digest", Json::string(w.source_digest));
    completion.set("rexgraph_headers_sha256", Json::string(gadget::build::rexgraph_sha256()));
    completion.set("gadget_sources_sha256", Json::string(gadget::build::source_sha256()));
    auto artifacts = Json::array();
    std::vector<fs::path> paths;
    for (const auto &entry : fs::recursive_directory_iterator(destination))
        if (entry.is_regular_file())
            paths.push_back(entry.path());
    std::sort(paths.begin(), paths.end());
    for (const auto &path : paths) {
        auto item = Json::object();
        item.set("path", Json::string(path.lexically_relative(destination).generic_string()));
        item.set("sha256", Json::string("sha256:" + sha256_file(path)));
        item.set("bytes", Json::integer(fs::file_size(path)));
        artifacts.push(item);
    }
    auto native_summary = V::object();
    native_summary.set("schema", V::string("gadget.analysis.v1"));
    native_summary.set("source_digest", V::string(w.source_digest));
    if (w.snapshot_origin.has("study_digest")) native_summary.set("snapshot", V(render_record_view(w.snapshot_origin.native())));
    if (w.store_origin.has("root")) native_summary.set("store", V(render_record_view(w.store_origin.native())));
    native_summary.set("rexgraph_headers_sha256", V::string(gadget::build::rexgraph_sha256()));
    native_summary.set("gadget_sources_sha256", V::string(gadget::build::source_sha256()));
    native_summary.set("numerical_contract", V::string(summary.at("numerical_contract").text()));
    native_summary.set("audit", V(w.audit.native().unbound_value()));
    native_summary.set("certificates", std::move(native_certificates));
    native_summary.set("performance", std::move(native_performance));
    native_summary.set("interpretation", V::string(summary.at("interpretation").text()));
    retained_results.set("contexts", std::move(retained_contexts));
    retained_results.set("cross_context", std::move(retained_cross));
    retained_results.set("gene_sweeps", std::move(retained_sweeps));
    retained_results.set("roundoff_relative_budget", approximate(1e-10));
    retained_results.set("exact_check_selection", V::count(30));
    native.bind(std::move(retained_results));
    std::uint64_t next = native.next_source_key();
    std::vector<RawByteSource> output_bytes;
    N::Array outputs;
    for(const auto &path:paths) {
        const auto key = next;
        next = rfc::machine::checked_size_add(next, 1);
        auto bytes=rfc::storage::read_file_exact(path,source_byte_limit());const auto count=bytes.size();
        output_bytes.push_back(RawByteSource::own(key,std::move(bytes)));
        outputs.elements.emplace_back(N::Object{{{"path",N(path.lexically_relative(destination).generic_string())},
            {"bytes",N(N::ByteSpan{key,0,count,false})}}});
    }
    native.retain_sources(std::move(output_bytes));
    const auto publication = native.finish_wrapped(std::move(native_summary), N::Object{{
        {"schema", N(std::string("gadget.analysis.native.v1"))},
        {"outputs", N(std::move(outputs))}}}, "summary");
    publication.write(destination / "analysis.rfcxs");
    auto binary_artifact=Json::object();
    binary_artifact.set("path",Json::string("analysis.rfcxs"));binary_artifact.set("bytes",Json::integer(fs::file_size(destination/"analysis.rfcxs")));
    binary_artifact.set("sha256",Json::string("sha256:"+sha256_file(destination/"analysis.rfcxs")));artifacts.push(binary_artifact);
    completion.set("artifacts", artifacts);
    write_text(destination / "manifests/artifacts.json", completion.dump() + "\n");
    if (fs::exists(requested))
        throw std::runtime_error("output destination appeared before commit");
    fs::rename(destination, requested);
}
} // namespace gadget::statistics
