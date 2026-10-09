#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <gadget/execution.hpp>
#include <gadget/snapshot.hpp>
#include <set>
#include <iostream>
#include <numeric>
#include <rfc/statistics/field_moments.hpp>
#include <sstream>

namespace gadget {
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
    auto certificates = Json::array();
    write_text(destination / "manifests/source.json", w.manifest().dump() + "\n");
    auto metrics = Json::object();
    for (const auto &realization : {"support", "rp"}) {
        std::cerr << "Preparing canonical " << realization << " contexts\n";
        Execution ad(w, "AD", realization), hd(w, "HD", realization);
        for (auto *op : {&ad, &hd}) {
            auto prefix = op->context == "AD" ? "alzheimers" : "huntingtons";
            std::string base = std::string(prefix) + "_" + realization;
            std::array<Vector, 3> fields, pulled, powers;
            for (std::size_t i = 0; i < 3; ++i) {
                fields[i] = w.fields.at(op->context + "." + "GET"[i]);
                pulled[i] = op->boundary().transpose_apply(fields[i]);
                powers[i] = op->metrics()[1].solve(pulled[i]);
            }
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
                if (power > 0 && power < 4)
                    for (auto &value : powers)
                        value = op->metrics()[1].solve(
                            binary64_values(op->relation_gram().apply(approximate_values(value))));
            }
            auto diagnostics = Json::array();
            std::array<Vector, 3> responses;
            for (std::size_t i = 0; i < 3; ++i) {
                Json diagnostic;
                responses[i] = op->green(fields[i], Rational{1}, 1e-9, diagnostic);
                diagnostic.set("channel", Json::string(std::string(1, "GET"[i])));
                diagnostics.push(diagnostic);
            }
            std::array<std::array<Rational, 3>, 3> gram;
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    gram[i][j] = rfc::numeric::dot(responses[i], responses[j]);
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
                auto ga = ad.green(source, Rational{1}, 1e-9, da),
                     gh = hd.green(source, Rational{1}, 1e-9, dh);
                auto delta = rfc::numeric::sub(gh, ga);
                auto action = rfc::numeric::sub(hd.hodge().apply(ga), ad.hodge().apply(ga));
                auto predicted =
                    rfc::numeric::scale(Rational{-1}, hd.green(action, Rational{1}, 1e-9, dr));
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
            Rational squared = qa + qh - 2 * cross;
            squared.canonicalize();
            if (squared < 0) {
                if (std::abs(squared.get_d()) > 1e-10 * std::max(1.0, Rational(qa + qh).get_d()))
                    throw std::runtime_error(
                        "negative numerical squared differential beyond roundoff budget");
                squared = 0;
            }
            ranking.emplace_back(squared, i);
            auto ma = sparse_inner(ad_pulled, ad_coeff), mh = sparse_inner(hd_pulled, hd_coeff);
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
        metrics.set(
            std::string(realization) + "_sweep_seconds",
            Json::number(
                std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()));
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
        metrics.set(std::string(realization) + "_exact_checks_and_explanations_seconds",
                    Json::number(std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                                               exact_checks_started)
                                     .count()));
        metrics.set(
            std::string(realization) + "_differential_seconds",
            Json::number(
                std::chrono::duration<double>(std::chrono::steady_clock::now() - phase).count()));
    }
    auto summary = Json::object();
    summary.set("schema", Json::string("gadget.analysis.v1"));
    summary.set("source_digest", Json::string(w.source_digest));
    if (w.snapshot_origin.has("study_digest"))
        summary.set("snapshot", w.snapshot_origin);
    if (w.store_origin.has("root"))
        summary.set("store", w.store_origin);
    summary.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    summary.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
    summary.set(
        "numerical_contract",
        Json::string("p0/p1 moments exact; higher moments and full gene sweep use declared "
                     "binary64 relation Gram actions; top 30 per realization independently checked "
                     "by canonical exact sparse Hodge with relative error budget 1e-10; every "
                     "Green solve checks its canonical source residual <=1e-9"));
    summary.set("audit", w.audit);
    summary.set("certificates", certificates);
    metrics.set(
        "analysis_seconds",
        Json::number(
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()));
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
    completion.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    completion.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
    auto artifacts = Json::array();
    std::vector<fs::path> paths;
    for (const auto &entry : fs::recursive_directory_iterator(destination))
        if (entry.is_regular_file())
            paths.push_back(entry.path());
    std::sort(paths.begin(), paths.end());
    for (const auto &path : paths) {
        auto bytes = read_text(path);
        auto item = Json::object();
        item.set("path", Json::string(path.lexically_relative(destination).generic_string()));
        item.set("sha256", Json::string("sha256:" + sha256(bytes)));
        item.set("bytes", Json::integer(bytes.size()));
        artifacts.push(item);
    }
    // Preserve the complete admitted complex and source bytes with the exact
    // certificates and original bytes of every generated compatibility export.
    const auto publication=capture_workspace(w);
    using N=NativeRecordValue;
    std::set<std::uint64_t> used;
    for(const auto &source:publication.raw_sources())used.insert(source.key());
    std::uint64_t next=1;std::vector<RawByteSource> output_bytes;N::Array outputs;
    for(const auto &path:paths) {
        while(used.contains(next)) {
            if(next==std::numeric_limits<std::uint64_t>::max())throw std::overflow_error("analysis byte identity space exhausted");
            ++next;
        }
        const auto key=next;used.insert(key);
        auto bytes=rfc::storage::read_file_exact(path,source_byte_limit());const auto count=bytes.size();
        output_bytes.push_back(RawByteSource::own(key,std::move(bytes)));
        outputs.elements.emplace_back(N::Object{{{"path",N(path.lexically_relative(destination).generic_string())},
            {"bytes",N(N::ByteSpan{key,0,count,false})}}});
    }
    publication.derive(N(N::Object{{{"schema",N(std::string("gadget.analysis.native.v1"))},
        {"summary",import_record_view(summary)},{"outputs",N(std::move(outputs))}}}),std::move(output_bytes)).write(destination/"analysis.rfcxs");
    const auto binary=read_text(destination/"analysis.rfcxs");auto binary_artifact=Json::object();
    binary_artifact.set("path",Json::string("analysis.rfcxs"));binary_artifact.set("bytes",Json::integer(binary.size()));
    binary_artifact.set("sha256",Json::string("sha256:"+sha256(binary)));artifacts.push(binary_artifact);
    completion.set("artifacts", artifacts);
    write_text(destination / "manifests/artifacts.json", completion.dump() + "\n");
    if (fs::exists(requested))
        throw std::runtime_error("output destination appeared before commit");
    fs::rename(destination, requested);
}
} // namespace gadget
