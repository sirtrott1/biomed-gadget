#include "fixture.hpp"
#include <iostream>

namespace {
gadget::Json count_json(std::size_t n) {
    return gadget::Json::integer(static_cast<std::int64_t>(n));
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::invalid_argument("Usage: gadget_temporal_fixture NEW_DIRECTORY");
        using namespace gadget;
        using namespace rfc;
        using namespace rfc::storage;
        const fs::path root = argv[1];
        example::scientific_study(root);
        ScientificStudy admitted(root / "study.json");
        auto declaration = admitted.declaration();
        TemporalEndpoint endpoint = *admitted.endpoint();
        // Three admitted subcomplexes retain original relation identities; only
        // the middle endpoint includes the declared grade two relation.
        for (auto &block : endpoint.local_blocks)
            if (block.block_key == 1)
                std::erase_if(block.cells, [](const auto &cell) { return cell.grade.value == 2; });
        auto terminal = endpoint.local_blocks.front();
        terminal.block_key = 3;
        endpoint.local_blocks.push_back(terminal);
        auto contexts = declaration.at("contexts");
        auto final = Json::parse(contexts.at(0).dump());
        final.set("id", Json::string("final"));
        final.set("block", Json::string("3"));
        final.at("time").set("position", Json::string("5"));
        contexts.push(final);
        for (std::size_t i = 0; i < contexts.size(); ++i)
            contexts.at(i).at("metrics").at(1).set("kind", Json::string("identity"));
        // Replace diagonal declarations completely; identity has no values.
        for (std::size_t i = 0; i < contexts.size(); ++i) {
            auto metrics = contexts.at(i).at("metrics");
            auto m = Json::object();
            m.set("kind", Json::string("identity"));
            json_object_array_put_idx(metrics.raw(), 1, json_object_get(m.raw()));
        }
        const std::string bytes =
            "context\tR101\tR102\tR103\tweights\n"
            "baseline\t1\t0\t2\t1,3,2\nfollowup\t3\t1\t6\t2,1,4\nfinal\t13\t0\t-7\t3,1,2\n"
            "J0 columns R101,R103; rows R101,R102,R103: [[1,0],[0,1],[1,1]].\n"
            "J1 columns R101,R102,R103; rows R101,R103: [[1,0,1],[0,1,-1]].\n"
            "Grade zero uses retained identity; grade two uses retained identity.\n"
            "Authored algebra fixture, not clinical observations.\n";
        write_text(root / "temporal-measurements.txt", bytes);
        auto e = Json::object();
        e.set("document", Json::string("temporal"));
        e.set("start", Json::integer(0));
        e.set("end", count_json(bytes.size()));
        e.set("text", Json::string(bytes));
        auto evidence = Json::array();
        evidence.push(e);
        auto doc = Json::object();
        doc.set("id", Json::string("temporal"));
        doc.set("path", Json::string("temporal-measurements.txt"));
        doc.set("sha256", Json::string("sha256:" + sha256(bytes)));
        doc.set("format", Json::string("utf8"));
        declaration.at("documents").push(doc);
        std::size_t key = 5000;
        for (const auto &[name, values, weights] :
             {std::tuple{"baseline", Vector{1, 0, 2}, Vector{1, 3, 2}},
              std::tuple{"followup", Vector{3, 1, 6}, Vector{2, 1, 4}},
              std::tuple{"final", Vector{13, 0, -7}, Vector{3, 1, 2}}}) {
            for (const auto &[suffix, x] :
                 {std::pair{"x", values}, std::pair{"weights", weights}}) {
                TemporalTensor tensor;
                tensor.key = ++key;
                tensor.role = TemporalTensorRole::Field;
                tensor.output = 2;
                tensor.interpretation = std::string(name) + "." + suffix;
                std::vector<numeric::SparseRationalEntry> nonzero;
                for (std::size_t i = 0; i < x.size(); ++i)
                    if (x[i] != 0)
                        nonzero.push_back({i, x[i]});
                tensor.tensor = {numeric::IndexDomain::linear(x.size()),
                                 numeric::StructuralSupport::full(x.size()),
                                 numeric::SparseRationalVector(x.size(), std::move(nonzero))};
                endpoint.tensors.push_back(tensor);
                auto f = Json::object();
                f.set("id", Json::string(tensor.interpretation));
                f.set("key", Json::string(std::to_string(tensor.key)));
                f.set("context", Json::string(name));
                f.set("channel", Json::string("observation"));
                f.set("quantity",
                      Json::string(std::string(suffix) == "x" ? "authored_relation_amplitude"
                                                              : "authored_multiplier"));
                f.set("units", Json::string("dimensionless"));
                f.set("method", Json::string("Authored exact temporal amplitudes"));
                f.set("value_type", Json::string("quantitative"));
                f.set("evidence", evidence);
                declaration.at("fields").push(f);
            }
        }
        auto archive = encode_native_series(NativeTemporalSeries(endpoint, 9876));
        const std::string native(reinterpret_cast<const char *>(archive.data()), archive.size());
        write_text(root / "study.rfcxs", native);
        declaration.at("native").set("sha256", Json::string("sha256:" + sha256(native)));
        write_text(root / "study.json", declaration.dump() + "\n");
        auto query = Json::parse(
            R"({"operation":"temporal_field_history","endpoints":[{"context":"baseline","field":"baseline.x"},{"context":"followup","field":"followup.x"},{"context":"final","field":"final.x"}],"correspondences":[],"word":[{"kind":"diagonal_field","fields":["baseline.weights","followup.weights","final.weights"]},{"kind":"hodge"},{"kind":"boundary"}],"lambda":"1"})");
        for (std::size_t t = 0; t < 2; ++t) {
            auto maps = Json::array();
            for (std::size_t grade = 0; grade < 3; ++grade) {
                auto map = Json::object();
                map.set("grade", count_json(grade));
                map.set("kind", Json::string(grade == 1 ? "explicit" : "retained_identity"));
                map.set(
                    "meaning",
                    Json::string("Authored correspondence on the indicated retained source grade"));
                map.set("evidence", evidence);
                if (grade == 1) {
                    auto entries = Json::array();
                    const auto rows =
                        t ? std::vector<unsigned>{101, 103} : std::vector<unsigned>{101, 102, 103};
                    const auto cols =
                        t ? std::vector<unsigned>{101, 102, 103} : std::vector<unsigned>{101, 103};
                    const std::vector<std::vector<int>> values =
                        t ? std::vector<std::vector<int>>{{1, 0, 1}, {0, 1, -1}}
                          : std::vector<std::vector<int>>{{1, 0}, {0, 1}, {1, 1}};
                    for (std::size_t i = 0; i < rows.size(); ++i)
                        for (std::size_t j = 0; j < cols.size(); ++j)
                            if (values[i][j]) {
                                auto entry = Json::object();
                                entry.set("row",
                                          Json::string("assertion:R" + std::to_string(rows[i])));
                                entry.set("column",
                                          Json::string("assertion:R" + std::to_string(cols[j])));
                                entry.set("value", Json::string(std::to_string(values[i][j])));
                                entries.push(entry);
                            }
                    map.set("entries", entries);
                }
                maps.push(map);
            }
            query.at("correspondences").push(maps);
        }
        write_text(root / "history-query.json", query.dump() + "\n");
        std::cout << (root / "study.json").string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
