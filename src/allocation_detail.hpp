#pragma once
#include "genomic_detail.hpp"
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail::allocation_detail {
using rfc::numeric::ExactMatrix;
using rfc::numeric::Integer;
Json vector_json(const Vector &values);
Json matrix_json(const ExactMatrix &matrix);
bool satisfies(const Rational &value, const std::string &sense, const Rational &rhs);

struct Problem {
    Json declaration, names = Json::array();
    std::string domain, direction;
    Vector lower, upper, rhs;
    std::vector<std::string> senses;
    ExactMatrix constraints, objective;
    Integer declared_box{1};
    explicit Problem(const Json &request);
};
struct Block {
    std::vector<std::size_t> parameters, variables, inequalities;
    ExactMatrix embedding, normals, objective;
    Vector offset, rhs, lower, upper;
    std::vector<std::string> senses;
};
struct Factors {
    ExactMatrix equalities, kernel, inequalities, grouping, indicator, objective;
    Vector equality_rhs, offset, rhs;
    std::vector<std::string> senses;
    std::vector<std::size_t> free_variables;
    std::vector<Block> blocks;
    std::size_t equality_rank{};
    Json equality_origins = Json::array(), inequality_origins = Json::array();
    Json carrier, contradiction = Json::object();
    bool consistent{true};
    explicit Factors(const Problem &problem);
    Json certificate(const Problem &problem) const;
};
} // namespace gadget::assay_detail::allocation_detail
