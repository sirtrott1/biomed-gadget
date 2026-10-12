#include <gadget/native_record.hpp>
#include <gadget/record_view.hpp>
#include <rfc/database/native_tensor_window.hpp>
#include <charconv>
#include <iostream>
#include <limits>
#include <set>

namespace {
using N = gadget::NativeRecordValue;
using R = gadget::NativeRecord;
void require(bool value, const char *message) {
    if (!value) throw std::invalid_argument(message);
}
std::uint64_t integer(const R &value) {
    require(value.kind() == gadget::NativeRecordKind::Integer,
            "native selector requires exact integer selectors");
    const auto input = value.integer();
    require(input >= 0 && input <= rfc::numeric::Integer(std::to_string(std::numeric_limits<std::uint64_t>::max())),
            "native selector integer is outside its host representation");
    return std::stoull(input.get_str());
}
std::uint64_t source_namespace(const R &value) {
    const auto text = value.text();
    std::uint64_t result{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), result);
    require(!text.empty() && error == std::errc{} && end == text.data() + text.size() && result,
            "native selector requires an unsigned source namespace");
    return result;
}
N number(std::uint64_t value) { return N(rfc::numeric::Integer(std::to_string(value))); }
N nullable(const std::optional<std::uint64_t> &value) {
    return value ? number(*value) : N{};
}
N array(const std::vector<std::size_t> &values) {
    N::Array output;
    for (const auto value : values) output.elements.push_back(number(value));
    return N(std::move(output));
}
std::vector<std::size_t> coordinates(const R &values) {
    require(values.kind() == gadget::NativeRecordKind::Array,
            "native window coordinates require an array");
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < values.size(); ++i)
        result.push_back(rfc::machine::checked_size(integer(values.at(i))));
    return result;
}
void members(const R &record, const std::set<std::string> &required,
             const std::set<std::string> &optional = {}) {
    require(record.kind() == gadget::NativeRecordKind::Object, "native selector requires an object");
    for (const auto &name : required) require(record.has(name), "native selector member missing");
    for (const auto &name : record.names())
        require(required.contains(name) || optional.contains(name), "native selector member unknown");
}
void carrier(const R &value, const std::optional<std::uint64_t> &expected) {
    require(expected ? value.kind() == gadget::NativeRecordKind::Integer && integer(value) == *expected
                     : value.kind() == gadget::NativeRecordKind::Null,
            "native tensor carrier selection mismatch");
}
} // namespace

int main(int argc, char **argv) {
    try {
        std::map<std::string, std::filesystem::path> options;
        for (int i = 1; i < argc; ++i) {
            const std::string name(argv[i]);
            require((name == "--input" || name == "--request" || name == "--out") &&
                        i + 1 < argc && !options.contains(name),
                    "usage: gadget_native_field_query --input RFCXS --request FILE --out NEW_RFCXS");
            options.emplace(name, argv[++i]);
        }
        require(options.size() == 3, "native selector requires input, request and output");
        require(!std::filesystem::exists(options.at("--out")), "native selector destination exists");
        const auto publication = R::read(options.at("--input"));
        const auto request = gadget::read_record(options.at("--request")).native();
        members(request, {"schema", "source_namespace", "cutoff", "tensor", "input_carrier",
                          "output_carrier", "interpretation"}, {"window", "execution"});
        require(request.at("schema").text() == "gadget.rcdb.field.query.v1",
                "native field query schema mismatch");
        const auto ns = source_namespace(request.at("source_namespace"));
        require(ns == publication.source_namespace(), "native field source namespace mismatch");
        const auto cutoff = integer(request.at("cutoff"));
        const auto selected = publication.historical_publication(cutoff);
        const auto endpoint = selected.endpoint();
        const auto key = integer(request.at("tensor"));
        const auto found = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                                       [&](const auto &tensor) { return tensor.key == key; });
        require(found != endpoint->tensors.end(), "native selected tensor is unavailable");
        carrier(request.at("input_carrier"), found->input);
        carrier(request.at("output_carrier"), found->output);
        require(request.at("interpretation").text() == found->interpretation,
                "native tensor interpretation selection mismatch");
        rfc::storage::NativeAccessionRequest accession;
        accession.kind = rfc::storage::NativeAccessionKind::Tensor;
        accession.subject = key;
        const auto pin = selected.accession_directory()->pin(accession);
        require(pin.source_namespace() == ns && pin.sequence() == cutoff && pin.endpoint() == endpoint,
                "native tensor accession owner mismatch");
        auto lower = std::vector<std::size_t>(found->tensor.domain().rank(), 0);
        auto upper = found->tensor.domain().extents();
        auto supported = found->tensor.support().count(), nonzeros = found->tensor.stored_nonzeros();
        require(request.has("window") == request.has("execution"),
                "native window requires explicit execution limits");
        if (request.has("window")) {
            const auto bounds = request.at("window"), execution = request.at("execution");
            members(bounds, {"lower", "upper"});
            members(execution, {"name", "records", "rational_slots", "sparse_entries", "scalar_work", "integer_bits"});
            lower = coordinates(bounds.at("lower")); upper = coordinates(bounds.at("upper"));
            rfc::runtime::ExecutionHostPolicy policy;
            policy.name = execution.at("name").text();
            policy.limits = {rfc::machine::checked_size(integer(execution.at("records"))),
                             rfc::machine::checked_size(integer(execution.at("rational_slots"))),
                             rfc::machine::checked_size(integer(execution.at("sparse_entries"))),
                             rfc::machine::checked_size(integer(execution.at("scalar_work"))),
                             rfc::machine::checked_size(integer(execution.at("integer_bits")))};
            const auto window = rfc::database::NativeTensorWindow::prepare(
                pin, lower, upper, {ns, cutoff}, rfc::runtime::ExecutionLease(policy));
            require(window.source_publication() == endpoint, "native window changed source owner");
            supported = window.supported_count(); nonzeros = window.nonzero_count();
        }
        // The receipt's source is administrative. The unchanged selected prefix
        // retains the mathematical fields, domains, carriers and original bytes.
        auto original = selected.archive_source(1);
        const auto bytes = original.bytes().size();
        N receipt(N::Object{{{"schema", N(std::string("gadget.rcdb.field.selection.v1"))},
            {"source_namespace", N(std::to_string(ns))}, {"cutoff", number(cutoff)},
            {"cutoff_kind", N(std::string("publication_sequence"))},
            {"tensor", number(key)}, {"input_carrier", nullable(found->input)},
            {"output_carrier", nullable(found->output)},
            {"interpretation", N(found->interpretation)},
            {"source_domain", array(found->tensor.domain().extents())},
            {"lower", array(lower)}, {"upper", array(upper)},
            {"supported_coordinates", number(supported)}, {"nonzero_entries", number(nonzeros)},
            {"selected_publication", N(N::ByteSpan{1, 0, bytes, false})}}});
        rfc::core::RelationalComplexBuilder builder(rfc::core::RelationSignature{});
        builder.declare_witness_identity(1);
        gadget::retain_operation(rfc::storage::TemporalEndpoint(builder.freeze()),
            gadget::new_native_record_namespace(), {std::move(original)}, request,
            std::move(receipt)).write(options.at("--out"));
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
