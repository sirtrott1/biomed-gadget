#include <gadget/native_record.hpp>
#include <gadget/record_view.hpp>
#include <rfc/database/native_tensor_window.hpp>
#include <rfc/storage/rfcx.hpp>
#include <iostream>

namespace {
void check(bool value, const char *message) { if (!value) throw std::runtime_error(message); }
}
int main(int argc, char **argv) {
    try {
        check(argc == 4, "source, selection and request required");
        const auto source = gadget::NativeRecord::read(argv[1]);
        const auto selection = gadget::NativeRecord::read(argv[2]);
        const auto request = gadget::read_record(argv[3]).native();
        const auto receipt = selection.at("result");
        check(receipt.at("schema").text() == "gadget.rcdb.field.selection.v1", "selection schema");
        const auto cutoff = request.at("cutoff").integer().get_ui();
        const auto original = source.historical_publication(cutoff);
        const auto archive = receipt.at("selected_publication").bytes();
        const auto encoded = original.encode();
        check(archive.size() == encoded.size() && std::equal(archive.begin(), archive.end(), encoded.begin()),
              "complete selected original prefix bytes changed");
        const auto reopened = gadget::NativeRecord::decode(archive);
        check(reopened.source_namespace() == original.source_namespace() &&
              reopened.endpoint()->sequence == cutoff, "selected owner or cutoff changed");
        check(reopened.endpoint()->carriers == original.endpoint()->carriers &&
              reopened.endpoint()->tensors.size() == original.endpoint()->tensors.size() &&
              rfc::storage::encode_rfcx(reopened.endpoint()->source).bytes ==
                  rfc::storage::encode_rfcx(original.endpoint()->source).bytes,
              "original source, occurrence order or carriers changed");
        for (std::size_t i = 0; i < original.endpoint()->tensors.size(); ++i)
            check(rfc::storage::temporal_endpoint_detail::tensor_equal(
                reopened.endpoint()->tensors[i], original.endpoint()->tensors[i]),
                "complete original field/support changed");
        const auto tensor_key = receipt.at("tensor").integer().get_ui();
        rfc::storage::NativeAccessionRequest address;
        address.kind = rfc::storage::NativeAccessionKind::Tensor; address.subject = tensor_key;
        const auto pin = reopened.accession_directory()->pin(address);
        const auto tensor = pin.tensor();
        if (tensor_key == 888) {
            check(tensor->state(0) == rfc::numeric::Rational{0} && !tensor->state(1) &&
                      tensor->state(2) == rfc::numeric::Rational{-5, 13},
                  "supported zero, absence or negative exact coefficient changed");
        } else if (tensor_key == 891) {
            check(tensor->state(0) == rfc::numeric::Rational{-7, 5}, "source grade field value changed");
        }
        check(receipt.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
              receipt.at("cutoff").integer() == cutoff, "receipt points to another owner");
        const auto lower = receipt.at("lower"), upper = receipt.at("upper");
        std::size_t supported{}, nonzeros{}, selected_count{1};
        std::vector<std::size_t> lo, hi;
        for (std::size_t i = 0; i < lower.size(); ++i) {
            lo.push_back(lower.at(i).integer().get_ui()); hi.push_back(upper.at(i).integer().get_ui());
            selected_count *= hi.back() - lo.back();
        }
        for (std::size_t index = 0; index < tensor->domain().logical_size(); ++index) {
            bool included = true;
            for (std::size_t axis = 0; axis < lo.size(); ++axis) {
                const auto coordinate = index / tensor->domain().strides()[axis] % tensor->domain().extents()[axis];
                included = included && coordinate >= lo[axis] && coordinate < hi[axis];
            }
            if (included) {
                const auto value = tensor->state(index);
                supported += value.has_value(); nonzeros += value && *value != 0;
            }
        }
        check(receipt.at("supported_coordinates").integer() == supported &&
              receipt.at("nonzero_entries").integer() == nonzeros,
              "window receipt differs from complete independent coordinate enumeration");
        check(selection.endpoint()->source.grade({0}).cell_count() == 1,
              "administrative receipt mislabeled mathematical source");
        std::cout << "PASS: retained complete original source, exact field/carriers, support and "
                     "declared window (" << selected_count << " coordinates)\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
