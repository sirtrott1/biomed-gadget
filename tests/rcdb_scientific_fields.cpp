#include <gadget/native_record.hpp>
#include <rfc/storage/temporal_endpoint.hpp>
#include <algorithm>
#include <iostream>

namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void same_bytes(const gadget::NativeRecord &left, const gadget::NativeRecord &right,
                const char *message) {
    const auto a = left.bytes(), b = right.bytes();
    require(std::ranges::equal(a, b), message);
}
}

int main(int argc, char **argv) {
    try {
        require(argc == 3, "expected original and replayed native operations");
        const auto original = gadget::NativeRecord::read(argv[1]);
        const auto replayed = gadget::NativeRecord::read(argv[2]);
        const auto a = original.endpoint(), b = replayed.endpoint();
        require(a->tensors.size() == b->tensors.size(), "replay changed complete native tensor family");
        std::size_t supported = 0, nonzeros = 0;
        for (std::size_t i = 0; i < a->tensors.size(); ++i) {
            require(rfc::storage::temporal_endpoint_detail::tensor_equal(a->tensors[i], b->tensors[i]),
                    "replay changed native domain, support, coefficients, role or carriers");
            supported += a->tensors[i].tensor.support().count();
            nonzeros += a->tensors[i].tensor.stored_nonzeros();
        }
        same_bytes(original.at("request_publication"), replayed.at("request_publication"),
                   "replay changed the original native request owner");
        same_bytes(original.at("result").at("source_study_publication"),
                   replayed.at("result").at("source_study_publication"),
                   "replay changed or extended the complete original study archive");
        require(original.raw_sources().size() == replayed.raw_sources().size(),
                "replay added original archive owners");
        std::cout << "{\"status\":\"PASS\",\"complete_tensors\":" << a->tensors.size()
                  << ",\"supported_coordinates\":" << supported << ",\"nonzero_entries\":" << nonzeros
                  << ",\"original_request_and_study_bytes_equal\":true}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
