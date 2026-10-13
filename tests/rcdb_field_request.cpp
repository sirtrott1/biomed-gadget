#include <gadget/native_record.hpp>
#include <algorithm>
#include <iostream>

namespace {
using N = gadget::NativeRecordValue;
N number(std::uint64_t value) {
    return N(rfc::numeric::Integer(std::to_string(value)));
}
N carrier(const std::optional<std::uint64_t> &value) {
    return value ? number(*value) : N{};
}
N coordinates(const std::vector<std::size_t> &values) {
    N::Array result;
    for (const auto value : values)
        result.elements.push_back(number(value));
    return N(std::move(result));
}
}

int main(int argc, char **argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("expected native source, tensor key and fresh request directory");
        const auto source = gadget::NativeRecord::read(argv[1]);
        const auto key = std::stoull(argv[2]);
        const auto endpoint = source.endpoint();
        const auto field = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                                       [&](const auto &value) { return value.key == key; });
        if (field == endpoint->tensors.end() || std::filesystem::exists(argv[3]))
            throw std::invalid_argument("actual tensor unavailable or request directory exists");
        std::filesystem::create_directories(argv[3]);
        N::Object request{{{"schema", N(std::string("gadget.rcdb.field.query.v1"))},
            {"source_namespace", N(std::to_string(source.source_namespace()))},
            {"cutoff", number(endpoint->sequence)}, {"tensor", number(key)},
            {"input_carrier", carrier(field->input)}, {"output_carrier", carrier(field->output)},
            {"interpretation", N(field->interpretation)}}};
        const auto root = std::filesystem::path(argv[3]);
        gadget::NativeRecord::create(N(request), gadget::new_native_record_namespace()).write(root/"full.rfcxs");
        auto lower = std::vector<std::size_t>(field->tensor.domain().rank(), 0);
        const auto upper = field->tensor.domain().extents();
        if (!upper.empty() && upper.front())
            lower.front() = upper.front()-1;
        request.members.emplace_back("window", N(N::Object{{
            {"lower", coordinates(lower)}, {"upper", coordinates(upper)}}}));
        request.members.emplace_back("execution", N(N::Object{{
            {"name", N(std::string("complete-source-final-row"))},
            {"records", number(10000000)}, {"rational_slots", number(10000000)},
            {"sparse_entries", number(10000000)}, {"scalar_work", number(128000000)},
            {"integer_bits", number(1048576)}}}));
        gadget::NativeRecord::create(N(std::move(request)), gadget::new_native_record_namespace()).write(root/"window.rfcxs");
        std::cout << "{\"status\":\"PASS\",\"tensor\":" << key
                  << ",\"source_cutoff\":" << endpoint->sequence
                  << ",\"raw_owners\":" << source.raw_sources().size()
                  << ",\"complete_support\":" << field->tensor.support().count()
                  << ",\"complete_nonzeros\":" << field->tensor.stored_nonzeros() << "}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
