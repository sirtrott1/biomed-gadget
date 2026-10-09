#pragma once
#include <gadget/native_record.hpp>
#include <rfc/machine/checked_size.hpp>

namespace gadget::assay_detail {
struct SourceGradeFields {
    rfc::storage::TemporalEndpoint &endpoint;
    template <class Value>
    NativeRecordValue column(std::string name, std::uint16_t grade, std::size_t size, Value value) {
        std::vector<rfc::numeric::SparseRationalEntry> entries;
        for (std::size_t i = 0; i < size; ++i) {
            rfc::numeric::Rational x = value(i);
            if (x != 0)
                entries.push_back({i, std::move(x)});
        }
        rfc::storage::TemporalTensor tensor;
        tensor.key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
        tensor.role = rfc::storage::TemporalTensorRole::Field;
        tensor.output = grade + 1;
        tensor.interpretation = std::move(name);
        tensor.tensor = {rfc::numeric::IndexDomain::linear(size),
                         rfc::numeric::StructuralSupport::full(size),
                         rfc::numeric::SparseRationalVector(size, std::move(entries))};
        const auto key = tensor.key;
        endpoint.tensors.push_back(std::move(tensor));
        return NativeRecordValue(rfc::numeric::Integer(std::to_string(key)));
    }
    NativeRecordValue field(std::string name, std::uint16_t grade, const Vector &values) {
        return column(std::move(name), grade, values.size(), [&](auto i) { return values[i]; });
    }
};
} // namespace gadget::assay_detail
