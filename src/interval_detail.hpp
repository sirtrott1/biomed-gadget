#pragma once
#include "assay_detail.hpp"
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail::interval_detail {
using Frame = std::pair<std::string, std::string>;
struct Record {
    Interval interval;
    std::size_t sample{}, index{};
};
struct Atom {
    Frame frame;
    std::uint64_t start{}, end{};
    std::vector<std::size_t> records;
};
struct Overlay {
    std::vector<Atom> atoms;
    std::size_t memberships{};
    explicit Overlay(const std::vector<Record> &records, const std::string &strand_policy);
    rfc::numeric::ExactMatrix coverage(std::size_t records) const;
    Vector lengths() const;
};
// Connected components of a declared threshold on the complete physical union.
Json threshold(const std::vector<Atom> &atoms, const Vector &lengths, const Vector &values,
               const Rational &threshold, bool report_atoms,
               NativeRecordValue *retained = nullptr);
} // namespace gadget::assay_detail::interval_detail
