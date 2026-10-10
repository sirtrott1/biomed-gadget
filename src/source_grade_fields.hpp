#pragma once
#include <gadget/native_record.hpp>
#include <rfc/machine/checked_size.hpp>

namespace gadget::assay_detail {
struct SourceGradeFields {
    rfc::storage::TemporalEndpoint &endpoint;
    template <class Value>
    NativeRecordValue column(std::string name, std::uint16_t grade, std::size_t size, Value value) {
        return column(std::move(name), grade, size, std::move(value),
                      rfc::numeric::StructuralSupport::full(size));
    }
    template <class Value>
    NativeRecordValue column(std::string name, std::uint16_t grade, std::size_t size, Value value,
                             rfc::numeric::StructuralSupport support) {
        auto tensor_value = tensor(size, std::move(value), std::move(support));
        rfc::storage::TemporalTensor tensor;
        tensor.key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
        tensor.role = rfc::storage::TemporalTensorRole::Field;
        tensor.output = grade + 1;
        tensor.interpretation = std::move(name);
        tensor.tensor = std::move(tensor_value);
        const auto key = tensor.key;
        endpoint.tensors.push_back(std::move(tensor));
        return NativeRecordValue(rfc::numeric::Integer(std::to_string(key)));
    }
    template <class Value>
    static rfc::numeric::SparseRationalTensor tensor(std::size_t size, Value value,
                                                     rfc::numeric::StructuralSupport support) {
        if (support.size() != size)
            throw std::invalid_argument("source column support shape");
        std::vector<rfc::numeric::SparseRationalEntry> entries;
        support.for_each_active([&](std::size_t i) {
            rfc::numeric::Rational x = value(i);
            if (x != 0)
                entries.push_back({i, std::move(x)});
        });
        return {rfc::numeric::IndexDomain::linear(size), std::move(support),
                rfc::numeric::SparseRationalVector(size, std::move(entries))};
    }
    NativeRecordValue field(std::string name, std::uint16_t grade, const Vector &values) {
        return column(std::move(name), grade, values.size(), [&](auto i) { return values[i]; });
    }
};
} // namespace gadget::assay_detail
