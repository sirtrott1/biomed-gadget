#include <gadget/data_field.hpp>
#include <gadget/io.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <algorithm>
#include <iostream>

namespace {
using namespace rfc;
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void refuses(F operation) {
    bool refused = false;
    try { operation(); } catch (const std::exception &) { refused = true; }
    check(refused, "invalid native data import accepted");
}
gadget::RawByteSource source(bool complete) {
    core::RelationSignature signature;
    const auto type = signature.add_type("original-row");
    core::RelationalComplexBuilder builder(signature);
    for (std::uint64_t row = 0; row < 3; ++row)
        builder.add_grade_one(core::RelationKey{row + 1}, type,
            {{row + 1, true, true, numeric::Rational{1}, {}, core::OccurrenceKey{2 * row + 1}},
             {4, true, false, numeric::Rational{1}, {}, core::OccurrenceKey{2 * row + 2}}});
    storage::TemporalEndpoint endpoint(builder.freeze());
    const auto domain = numeric::IndexDomain::linear(6);
    endpoint.carriers.push_back({71, storage::TemporalCarrierKind::Explicit, core::Grade{0}, {}, 71, domain});
    storage::TemporalTensor field;
    field.key = 72; field.role = storage::TemporalTensorRole::Field; field.output = 71;
    field.tensor = {domain, complete ? numeric::StructuralSupport::full(6) :
        numeric::StructuralSupport::sparse(6, {0, 2, 3, 4, 5}),
        numeric::SparseRationalVector(6, {{2, numeric::Rational{7, 3}}, {3, numeric::Rational{-5}},
                                         {5, numeric::Rational{9}}})};
    endpoint.tensors.push_back(std::move(field));
    const auto bytes = storage::encode_native_series(storage::NativeTemporalSeries(std::move(endpoint), 801, gadget::native_data_limits()));
    return gadget::RawByteSource(1, bytes);
}
gadget::Json coordinates(const gadget::RawByteSource &bytes, bool complete) {
    auto result = gadget::Json::parse(R"({"schema":"gadget.native-data.field.v1","shape":[3,2],"rows":[5,2,0],"columns":[4,1],"row_labels":["a","a","b"],"column_labels":["x","y"],"observed_count":5})");
    auto binding = gadget::Json::object();
    binding.set("namespace", gadget::Json::string("801"));
    binding.set("sha256", gadget::Json::string(gadget::sha256(std::string_view(
        reinterpret_cast<const char *>(bytes.bytes().data()), bytes.bytes().size()))));
    result.set("native_source", binding);
    if (complete) result.set("observed_count", gadget::Json::integer(6));
    return result;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 3 && (std::string_view(argv[1]) == "--read" || std::string_view(argv[1]) == "--apply" ||
                          std::string_view(argv[1]) == "--summary")) {
            const auto record = gadget::NativeRecord::read(argv[2]);
            const auto endpoint = record.endpoint();
            const auto field = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                                            [](const auto &value) { return value.key == 72; });
            check(field != endpoint->tensors.end(), "selected native field missing");
            auto result = gadget::Json::object(), values = gadget::Json::array();
            if (std::string_view(argv[1]) == "--summary") {
                result.set("coordinates", gadget::Json::integer(field->tensor.domain().logical_size()));
                result.set("observed", gadget::Json::integer(field->tensor.support().count()));
                result.set("nonzeros", gadget::Json::integer(field->tensor.values().entries().size()));
                const auto first = field->tensor.state(0), last = field->tensor.state(field->tensor.domain().logical_size() - 1);
                result.set("first", first ? gadget::Json::string(first->get_str()) : gadget::Json{});
                result.set("last", last ? gadget::Json::string(last->get_str()) : gadget::Json{});
                std::cout << result.dump() << '\n';
                return 0;
            }
            if (std::string_view(argv[1]) == "--apply") {
                check(field->tensor.support().count() == field->tensor.domain().logical_size(),
                      "exact matrix consumer requires observed support");
                for (const auto &value : storage::temporal_tensor_matrix(*field).apply(numeric::RationalVector{2, 3}))
                    values.push(gadget::Json::string(value.get_str()));
            } else {
                for (std::size_t i = 0; i < field->tensor.domain().logical_size(); ++i) {
                    const auto value = field->tensor.state(i);
                    values.push(value ? gadget::Json::string(value->get_str()) : gadget::Json{});
                }
            }
            result.set("values", values);
            result.set("coordinates", gadget::render_record_view(record));
            std::cout << result.dump() << '\n';
            return 0;
        }
        for (bool complete : {false, true}) {
            const auto bytes = source(complete);
            const auto control = coordinates(bytes, complete);
            const auto record = gadget::import_data_field(bytes, gadget::admit_record_view(control), 802);
            const auto cold = gadget::NativeRecord::decode(record.encode());
            const auto &field = *std::find_if(cold.endpoint()->tensors.begin(), cold.endpoint()->tensors.end(),
                                            [](const auto &value) { return value.key == 72; });
            check(field.tensor.domain().extents() == std::vector<std::size_t>{3, 2}, "native field shape changed");
            check(field.tensor.state(0) == numeric::Rational{0} &&
                      (complete ? field.tensor.state(1) == numeric::Rational{0} : !field.tensor.state(1)) &&
                      field.tensor.state(2) == numeric::Rational{7, 3}, "zero, absence or exact coefficient changed");
            if (complete) {
                const auto matrix = storage::temporal_tensor_matrix(field);
                const auto product = matrix.apply(numeric::RationalVector{2, 3});
                check(product == numeric::RationalVector{0, numeric::Rational{-31, 3}, 27},
                      "canonical exact matrix consumer differs from original coefficients");
            }
            auto bad = gadget::Json::parse(control.dump());
            bad.set("observed_count", gadget::Json::integer(1));
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(bad), 803); });
            bad = gadget::Json::parse(control.dump());
            bad.at("native_source").set("namespace", gadget::Json::string("900"));
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(bad), 803); });
            bad = gadget::Json::parse(control.dump());
            bad.at("native_source").set("sha256", gadget::Json::string(std::string(64, '0')));
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(bad), 803); });
            bad = gadget::Json::parse(control.dump());
            bad.set("shape", gadget::Json::parse("[2,3]"));
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(bad), 803); });
            bad = gadget::Json::parse(control.dump());
            bad.set("rows", gadget::Json::parse("[0,0,0]"));
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(bad), 803); });
            refuses([&] { (void)gadget::import_data_field(bytes, gadget::admit_record_view(control), 0); });
        }
        std::cout << "Native data field source, support, coordinates and exact consumer PASS\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
