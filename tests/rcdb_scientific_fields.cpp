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
struct Counts {
    std::size_t tensors{}, supported{}, nonzeros{}, plans{}, steps{};
};
void compare(const gadget::NativeRecord &original, const gadget::NativeRecord &replayed, Counts &counts) {
    const auto a = original.endpoint(), b = replayed.endpoint();
    require(a->tensors.size() == b->tensors.size(), "replay changed complete native tensor family");
    for (std::size_t i = 0; i < a->tensors.size(); ++i) {
        require(rfc::storage::temporal_endpoint_detail::tensor_equal(a->tensors[i], b->tensors[i]),
                "replay changed native domain, support, coefficients, role or carriers");
        counts.tensors = rfc::machine::checked_size_add(counts.tensors, 1);
        counts.supported = rfc::machine::checked_size_add(counts.supported, a->tensors[i].tensor.support().count());
        counts.nonzeros = rfc::machine::checked_size_add(counts.nonzeros, a->tensors[i].tensor.stored_nonzeros());
    }
    same_bytes(original.at("request_publication"), replayed.at("request_publication"),
               "replay changed the original native request owner");
    require(original.raw_sources().size() == replayed.raw_sources().size(),
            "replay added original archive owners");
    const auto result = original.at("result"), actual = replayed.at("result");
    require(result.at("schema").text() == actual.at("schema").text(), "replay changed result schema");
    if (result.at("schema").text() == "gadget.scientific.plan.fields.v1") {
        const auto steps = result.at("steps"), other = actual.at("steps");
        require(steps.size() > 0 && steps.size() == other.size(), "replay changed complete plan step count");
        require(result.at("source_digest").text() == actual.at("source_digest").text(),
                "replay changed plan source identity");
        counts.plans = rfc::machine::checked_size_add(counts.plans, 1);
        for (std::size_t i = 0; i < steps.size(); ++i) {
            require(steps.at(i).at("id").text() == other.at(i).at("id").text(),
                    "replay changed original plan step order or identity");
            counts.steps = rfc::machine::checked_size_add(counts.steps, 1);
            compare(gadget::NativeRecord::decode(steps.at(i).at("publication").bytes()),
                    gadget::NativeRecord::decode(other.at(i).at("publication").bytes()), counts);
        }
    } else {
        require(result.at("schema").text() == "gadget.scientific.result.fields.v1",
                "unsupported scientific replay result");
        same_bytes(result.at("source_study_publication"), actual.at("source_study_publication"),
                   "replay changed or extended the complete original study archive");
    }
}
}

int main(int argc, char **argv) {
    try {
        require(argc == 3, "expected original and replayed native operations");
        const auto original = gadget::NativeRecord::read(argv[1]);
        const auto replayed = gadget::NativeRecord::read(argv[2]);
        Counts counts;
        compare(original, replayed, counts);
        std::cout << "{\"status\":\"PASS\",\"complete_tensors\":" << counts.tensors
                  << ",\"supported_coordinates\":" << counts.supported << ",\"nonzero_entries\":" << counts.nonzeros
                  << ",\"plans\":" << counts.plans << ",\"ordered_steps\":" << counts.steps
                  << ",\"original_request_and_study_bytes_equal\":true}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
