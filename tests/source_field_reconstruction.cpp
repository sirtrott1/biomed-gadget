#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <algorithm>
#include <iostream>
#include <set>

int main(int argc, char **argv) try {
    using namespace gadget;
    std::size_t checks{};
    const auto check = [&](bool value, const char *message) {
        ++checks;
        if (!value) throw std::runtime_error(message);
    };
    check(argc == 3, "usage: source_field_reconstruction ARCHIVE EXPECTED");
    const auto original = NativeRecord::read(argv[1]);
    const auto restored = NativeRecord::decode(original.encode());
    check(restored.encode() == original.encode(), "complete native archive roundtrip differs");
    const ScientificStudy study(restored);
    const auto expected = read_record(argv[2]).native();
    const auto catalog = study.catalog();
    for (std::size_t r = 0; r < expected.at("patient_ids").size(); ++r)
        check(study.cell_identifier({0}, r) == expected.at("patient_ids").at(r).text(),
              "original patient identity or source row order differs");
    std::set<std::string> selected;
    for (std::size_t f = 0; f < expected.at("fields").size(); ++f) {
        const auto field = expected.at("fields").at(f);
        const auto name = field.at("id").text();
        check(selected.insert(name).second, "repeated expected field");
        const auto pin = study.field_accession(name);
        const auto &tensor = *pin.tensor_declaration();
        check(std::to_string(pin.source_namespace()) == expected.at("namespace").text(),
              "native source namespace differs");
        check(tensor.tensor.domain().logical_size() == field.at("values").size(),
              "native field omitted original source rows");
        for (std::size_t r = 0; r < field.at("values").size(); ++r) {
            const auto value = field.at("values").at(r);
            const auto present = value.kind() != NativeRecordKind::Null;
            check(tensor.tensor.support().contains(r) == present, "native missing support differs");
            if (present)
                check(tensor.tensor.state(r) == rfc::numeric::Rational(value.text()),
                      "complete native coefficient differs from original source reading");
        }
        bool found{};
        for (std::size_t i = 0; i < catalog.at("fields").size(); ++i) {
            const auto metadata = catalog.at("fields").at(i);
            if (metadata.at("id").text() != name) continue;
            found = true;
            check(metadata.at("value_type").text() == field.at("value_type").text() &&
                      metadata.at("units").text() == field.at("units").text(),
                  "source field type or units differ");
            if (field.has("codebook")) {
                constexpr std::string_view prefix = "Lexical category codebook (no ordinal or linear meaning): ";
                const auto method = metadata.at("method").text();
                check(method.starts_with(prefix), "categorical field omitted its lexical codebook");
                const auto codebook = Json::parse(method.substr(prefix.size()));
                check(static_cast<std::size_t>(json_object_object_length(codebook.raw())) ==
                          field.at("codebook").names().size(),
                      "category codebook changed label count");
                for (const auto &label : field.at("codebook").names())
                    check(codebook.at(label).text() == field.at("codebook").at(label).text(),
                          "category codebook changed an original spelling or code");
            }
        }
        check(found, "native catalog omitted original source field");
    }
    const auto raw = restored.raw_sources();
    for (std::size_t i = 0; i < expected.at("originals").size(); ++i) {
        const auto path = expected.at("originals").at(i).text();
        const auto bytes = RawByteSource::read(1, path);
        check(std::any_of(raw.begin(), raw.end(), [&](const auto &owner) {
            return std::ranges::equal(owner.bytes(), bytes.bytes());
        }), "native publication omitted or changed complete original bytes");
    }
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks
              << ",\"original_fields\":" << selected.size() << "}\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
