#include "interval_detail.hpp"
#include <iostream>
#include <stdexcept>
using namespace gadget::assay_detail;
int main() {
    // Every long record covers the entire finely divided union. The memberships
    // exceed the former ceiling; no source or output membership is sampled.
    std::vector<interval_detail::Record> records;
    constexpr std::size_t long_count = 101, short_count = 20001;
    for (std::size_t i = 0; i < long_count; ++i) {
        Interval r;
        r.contig = "chr";
        r.strand = "+";
        r.end = short_count;
        records.push_back({r, 0, i});
    }
    for (std::size_t i = 0; i < short_count; ++i) {
        Interval r;
        r.contig = "chr";
        r.strand = "+";
        r.start = i;
        r.end = i + 1;
        records.push_back({r, 0, long_count + i});
    }
    const interval_detail::Overlay geometry(records, "separate");
    if (geometry.atoms.size() != short_count ||
        geometry.memberships != short_count * (long_count + 1))
        throw std::runtime_error("complete interval overlay shape");
    const auto coverage = geometry.coverage(records.size());
    const auto result = coverage.apply(gadget::Vector(records.size(), 1));
    for (const auto &q : result)
        if (q != long_count + 1)
            throw std::runtime_error("native exact coverage");
    if (rfc::numeric::dot(geometry.lengths(), result) != short_count * (long_count + 1))
        throw std::runtime_error("exact complete integral");
    std::cout << "20001 atoms and " << geometry.memberships
              << " complete native memberships PASS\n";
}
