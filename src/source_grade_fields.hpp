#pragma once
#include <algorithm>
#include <cstdint>
#include <gadget/native_record.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail {
inline const rfc::storage::TemporalTensor &retained_tensor(const NativeRecord &record,
                                                           const NativeRecord &key) {
    const auto integer = key.integer();
    if (integer <= 0 || !mpz_fits_ulong_p(integer.get_mpz_t()))
        throw std::invalid_argument("invalid retained tensor key");
    const auto id = integer.get_ui();
    for (const auto &tensor : record.endpoint()->tensors)
        if (tensor.key == id)
            return tensor;
    throw std::invalid_argument("retained tensor key is absent");
}
struct SourceGradeFields {
    rfc::storage::TemporalEndpoint &endpoint;
    std::uint64_t next_key{1};
    explicit SourceGradeFields(rfc::storage::TemporalEndpoint &source) : endpoint(source) {
        for (const auto &tensor : endpoint.tensors)
            next_key = std::max<std::uint64_t>(next_key, rfc::machine::checked_size_add(tensor.key, 1));
    }
    std::uint64_t explicit_carrier(std::size_t size) {
        std::uint64_t key = 1;
        for (const auto &carrier : endpoint.carriers)
            key = std::max<std::uint64_t>(key, rfc::machine::checked_size_add(carrier.key, 1));
        endpoint.carriers.push_back({key, rfc::storage::TemporalCarrierKind::Explicit,
                                     {}, {}, key, rfc::numeric::IndexDomain::linear(size)});
        return key;
    }
    NativeRecordValue vector(std::string name, std::uint64_t carrier, const rfc::numeric::RationalVector &values) {
        return on_carrier(std::move(name), carrier,
                          tensor(values.size(), [&](auto i) { return values[i]; },
                                 rfc::numeric::StructuralSupport::full(values.size())));
    }
    NativeRecordValue insert(rfc::storage::TemporalTensor tensor) {
        const auto key = next_key;
        next_key = rfc::machine::checked_size_add(next_key, 1);
        tensor.key = key;
        endpoint.tensors.push_back(std::move(tensor));
        return NativeRecordValue(rfc::numeric::Integer(std::to_string(key)));
    }
    NativeRecordValue
    on_carrier(std::string name, std::uint64_t carrier, rfc::numeric::SparseRationalTensor value,
               rfc::storage::TemporalTensorRole role = rfc::storage::TemporalTensorRole::Field) {
        rfc::storage::TemporalTensor field;
        field.role = role;
        field.output = carrier;
        field.interpretation = std::move(name);
        field.tensor = std::move(value);
        return insert(std::move(field));
    }
    NativeRecordValue
    map(std::string name, std::uint64_t input, std::uint64_t output,
        const rfc::numeric::ExactMatrix &value,
        rfc::storage::TemporalTensorRole role = rfc::storage::TemporalTensorRole::Observation,
        rfc::storage::TemporalTensorOwner owner = {}) {
        rfc::storage::TemporalTensor field;
        field.role = role;
        field.owner = owner;
        field.input = input;
        field.output = output;
        field.interpretation = std::move(name);
        field.tensor = rfc::storage::temporal_matrix_tensor(value);
        return insert(std::move(field));
    }
    template <class Value>
    NativeRecordValue column(std::string name, std::uint16_t grade, std::size_t size, Value value) {
        return column(std::move(name), grade, size, std::move(value),
                      rfc::numeric::StructuralSupport::full(size));
    }
    template <class Value>
    NativeRecordValue column(std::string name, std::uint16_t grade, std::size_t size, Value value,
                             rfc::numeric::StructuralSupport support) {
        auto tensor_value = tensor(size, std::move(value), std::move(support));
        return on_carrier(std::move(name), grade + 1, std::move(tensor_value));
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
    NativeRecordValue field(std::string name, std::uint16_t grade, const rfc::numeric::RationalVector &values) {
        return column(std::move(name), grade, values.size(), [&](auto i) { return values[i]; });
    }
    static NativeRecordValue matrix(const rfc::numeric::ExactMatrix &value) {
        return NativeRecordValue(rfc::storage::temporal_matrix_tensor(value));
    }
};
} // namespace gadget::assay_detail
