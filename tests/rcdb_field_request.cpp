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
void write_requests(const gadget::NativeRecord &source,
                    const rfc::storage::TemporalTensor &field,
                    const std::filesystem::path &root, bool refusal_controls = false) {
    std::filesystem::create_directories(root);
    N::Object request{{{"schema", N(std::string("gadget.rcdb.field.query.v1"))},
        {"source_namespace", N(std::to_string(source.source_namespace()))},
        {"cutoff", number(source.endpoint()->sequence)}, {"tensor", number(field.key)},
        {"input_carrier", carrier(field.input)}, {"output_carrier", carrier(field.output)},
        {"interpretation", N(field.interpretation)}}};
    const auto write = [&](const char *name, const N::Object &value) {
        gadget::NativeRecord::create(N(value), gadget::new_native_record_namespace()).write(root/name);
    };
    write("full.rfcxs", request);
    if (refusal_controls) {
        auto future = request;
        rfc::numeric::Integer future_cutoff = rfc::numeric::Integer(std::to_string(source.endpoint()->sequence)) + 1;
        future.members.at(2).second = N(std::move(future_cutoff));
        write("future.rfcxs", future);
        auto foreign = request;
        foreign.members.at(1).second = N(source.source_namespace() == 1 ? std::string("2") : std::string("1"));
        write("foreign.rfcxs", foreign);
        auto wrong_carrier = request;
        wrong_carrier.members.at(4).second = field.input ? N{} : number(1);
        write("wrong-carrier.rfcxs", wrong_carrier);
    }
    auto lower = std::vector<std::size_t>(field.tensor.domain().rank(), 0);
    const auto upper = field.tensor.domain().extents();
    if (!upper.empty() && upper.front()) lower.front() = upper.front()-1;
    request.members.emplace_back("window", N(N::Object{{
        {"lower", coordinates(lower)}, {"upper", coordinates(upper)}}}));
    request.members.emplace_back("execution", N(N::Object{{
        {"name", N(std::string("complete-source-final-row"))},
        {"records", number(10000000)}, {"rational_slots", number(10000000)},
        {"sparse_entries", number(10000000)}, {"scalar_work", number(128000000)},
        {"integer_bits", number(1048576)}}}));
    write("window.rfcxs", request);
}
}

int main(int argc, char **argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("expected native source, tensor key and fresh request directory");
        if (std::string_view(argv[1]) == "--all-fields") {
            if (std::filesystem::exists(argv[3]))
                throw std::invalid_argument("request directory exists");
            const auto source = gadget::NativeRecord::read(argv[2]);
            const std::filesystem::path root(argv[3]);
            std::filesystem::create_directories(root);
            N::Array fields;
            for (const auto &field : source.endpoint()->tensors) {
                write_requests(source, field, root/std::to_string(field.key), true);
                fields.elements.push_back(number(field.key));
            }
            gadget::NativeRecord::create(N(N::Object{{
                {"schema", N(std::string("gadget.rcdb.field.test_manifest.v1"))},
                {"source_namespace", N(std::to_string(source.source_namespace()))},
                {"cutoff", number(source.endpoint()->sequence)}, {"fields", N(std::move(fields))}}}),
                gadget::new_native_record_namespace()).write(root/"manifest.rfcxs");
            std::cout << "{\"status\":\"PASS\",\"fields\":" << source.endpoint()->tensors.size() << "}\n";
            return 0;
        }
        const auto source = gadget::NativeRecord::read(argv[1]);
        const auto key = std::stoull(argv[2]);
        const auto endpoint = source.endpoint();
        const auto field = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                                       [&](const auto &value) { return value.key == key; });
        if (field == endpoint->tensors.end() || std::filesystem::exists(argv[3]))
            throw std::invalid_argument("actual tensor unavailable or request directory exists");
        write_requests(source, *field, argv[3]);
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
