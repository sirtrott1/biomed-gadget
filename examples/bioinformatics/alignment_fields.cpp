#include <algorithm>
#include <fstream>
#include <gadget/io.hpp>
#include <gadget/native_record.hpp>
#include <iostream>
#include <map>

int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("usage: verify_alignment RECORD OUTPUT_PREFIX");
        const auto publication = gadget::NativeRecord::read(argv[1]);
        const auto record = publication.at("result").at("result");
        using Tensor = rfc::numeric::SparseRationalTensor;
        const auto fields = [&](const gadget::NativeRecord &columns) {
            std::map<std::string, const Tensor *> output;
            for (const auto &name : columns.names()) {
                const auto key = std::stoull(columns.at(name).integer().get_str());
                const auto &tensors = publication.endpoint()->tensors;
                const auto it = std::find_if(tensors.begin(), tensors.end(),
                                             [&](const auto &t) { return t.key == key; });
                if (it == tensors.end())
                    throw std::invalid_argument("missing source field");
                output.emplace(name, &it->tensor);
            }
            return output;
        };
        if (record.at("schema").text() == "gadget.feature.count.fields.v1") {
            const auto quantitative = fields(record.at("fields"));
            std::ofstream output(std::string(argv[2]) + "-bounds.tsv");
            const auto candidates = record.at("candidates");
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                const auto lo = quantitative.at("lower")->state(i),
                           hi = quantitative.at("upper")->state(i);
                output << candidates.at(i).text() << '\t' << (lo ? lo->get_str() : "NA") << '\t'
                       << (hi ? hi->get_str() : "NA") << '\n';
            }
            output.close();
            output.open(std::string(argv[2]) + "-classes.tsv");
            const auto &grade = publication.endpoint()->source.grade({1});
            for (std::size_t i = 0; i < grade.cell_count(); ++i) {
                output << i << '\t' << quantitative.at("class_count")->state(i)->get_str();
                const auto [begin, end] = grade.incidence_range(i);
                for (auto j = begin; j < end; ++j)
                    output << '\t' << candidates.at(grade.lower_index[j]).text();
                output << '\n';
            }
            output.close();
            output.open(std::string(argv[2]) + "-segments.tsv");
            const auto segments = record.at("segments");
            std::map<std::string, const Tensor *> columns;
            for (const auto &name : segments.names())
                columns.emplace(name, &segments.at(name).tensor());
            const auto number = [&](const auto &name, auto i) {
                return std::stoull(columns.at(name)->state(i)->get_str());
            };
            const auto names = record.at("decoded_query_names").bytes();
            const auto groups = record.at("read_groups");
            for (std::size_t i = 0; i < std::stoull(record.at("segment_count").integer().get_str());
                 ++i) {
                const auto start = number("query_name_offset", i),
                           length = number("query_name_length", i);
                if (start > names.size() || length > names.size() - start)
                    throw std::invalid_argument("segment query span");
                output << i << '\t'
                       << std::string_view(reinterpret_cast<const char *>(names.data() + start),
                                           length)
                       << '\t' << groups.at(number("read_group", i)).text() << '\t'
                       << number("mate", i) << '\t' << number("class_index", i);
                const auto begin = number("source_record_offset", i),
                           count = number("source_record_count", i);
                for (auto j = begin; j < begin + count; ++j)
                    output << '\t' << number("source_records", j);
                output << '\n';
            }
            output.close();
            output.open(std::string(argv[2]) + "-overlaps.tsv");
            const auto overlaps = record.at("overlaps");
            const auto rows = std::stoull(record.at("overlap_count").integer().get_str());
            for (std::size_t i = 0; i < rows; ++i) {
                output << overlaps.at("source_record_index").tensor().state(i)->get_str() << '\t';
                const auto feature =
                    std::stoull(overlaps.at("feature").tensor().state(i)->get_str());
                output << candidates.at(feature).text();
                for (const auto *name : {"overlap_bases", "aligned_bases", "eligible"})
                    output << '\t' << overlaps.at(name).tensor().state(i)->get_str();
                if (overlaps.has("junction_compatible"))
                    output << '\t'
                           << overlaps.at("junction_compatible").tensor().state(i)->get_str();
                output << '\n';
            }
            for (const auto *name : {"original_alignments", "original_annotation"}) {
                const auto raw = record.at(name).bytes();
                std::cout << name << '\t'
                          << gadget::sha256(std::string_view(
                                 reinterpret_cast<const char *>(raw.data()), raw.size()))
                          << '\t' << raw.size() << '\n';
            }
            return 0;
        }
        if (record.at("schema").text() != "gadget.alignment.fields.v1")
            throw std::invalid_argument("alignment fields schema");
        const auto records = fields(record.at("record_fields"));
        const auto relations = fields(record.at("relation_fields"));
        const auto label = [&](const auto &columns, std::string_view field,
                               std::string_view dictionary, std::size_t i) {
            const auto value = columns.at(std::string(field))->state(i).value();
            return record.at(dictionary).at(std::stoull(value.get_str())).text();
        };
        const auto value = [](const auto &columns, std::string_view field, std::size_t i) {
            const auto v = columns.at(std::string(field))->state(i);
            return v ? v->get_str() : "NA";
        };
        const auto names = record.at("decoded_query_names").bytes();
        std::ofstream output(std::string(argv[2]) + "-records.tsv");
        for (std::size_t i = 0; i < publication.endpoint()->source.grade({0}).cell_count(); ++i) {
            const auto begin = std::stoull(value(records, "query_name_offset", i));
            const auto length = std::stoull(value(records, "query_name_length", i));
            if (begin > names.size() || length > names.size() - begin)
                throw std::invalid_argument("query name span");
            output << i << '\t'
                   << std::string_view(reinterpret_cast<const char *>(names.data() + begin), length)
                   << '\t' << label(records, "read_group", "read_groups", i) << '\t'
                   << value(records, "flags", i) << '\t' << value(records, "mapping_quality", i)
                   << '\t' << label(records, "contig", "contigs", i) << '\t'
                   << label(records, "strand", "strands", i) << '\t' << value(records, "start", i)
                   << '\t' << value(records, "end", i) << '\t' << value(records, "query_length", i)
                   << '\t' << value(records, "selected", i) << '\n';
        }
        output.close();
        output.open(std::string(argv[2]) + "-relations.tsv");
        const auto &grade = publication.endpoint()->source.grade({1});
        for (std::size_t i = 0; i < grade.cell_count(); ++i) {
            output << i << '\t' << value(relations, "kind", i) << '\t'
                   << label(relations, "contig", "contigs", i) << '\t'
                   << label(relations, "strand", "strands", i) << '\t'
                   << value(relations, "start", i) << '\t' << value(relations, "end", i) << '\t'
                   << value(relations, "selected_record_count", i);
            const auto [begin, end] = grade.incidence_range(i);
            for (auto slot = begin; slot < end; ++slot)
                output << '\t' << grade.lower_index[slot];
            output << '\n';
        }
        const auto raw = record.at("original").bytes();
        const auto raw_hash =
            "sha256:" + gadget::sha256(std::string_view(reinterpret_cast<const char *>(raw.data()),
                                                        raw.size()));
        if (raw_hash != record.at("original_sha256").text())
            throw std::invalid_argument("original hash");
        std::cout << raw_hash << '\t' << raw.size() << '\t' << record.at("record_count").integer()
                  << '\t' << record.at("block_relation_count").integer() << '\t'
                  << record.at("junction_relation_count").integer() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
