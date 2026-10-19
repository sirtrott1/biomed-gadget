#pragma once
#include <gadget/io.hpp>
#include <rfc/database/native_interpretation.hpp>
#include <stop_token>

namespace gadget::feature_context_detail {
struct Limits {
    std::size_t source_bytes, features, mapping_rows, assertions, participants, matches, work, output_bytes;
};
struct Request {
    std::string analysis_sha256, mapping_sha256, evidence_digest, question;
    Json feature_path, feature_id_path;
    std::vector<std::uint64_t> assertions;
    Limits limits;
};
class Work {
    std::size_t limit_, used_{};
    std::stop_token stop_;
  public:
    explicit Work(std::size_t limit,std::stop_token stop={}):limit_(limit),stop_(stop) {}
    void take(std::size_t n=1);
    void check() const;
    std::size_t used() const {return used_;}
};
using Roles = std::map<std::string, std::uint64_t>;
struct Participant {
    std::string identity,role;
    std::uint64_t occurrence{};
    std::size_t position{};
    bool active{};
};
struct Assertion {std::uint64_t key{},kind{};std::vector<Participant> participants;};
struct MappingRow {
    std::size_t source_record{},begin{},length{};
    std::string feature,identity,role;
};
struct Match {std::size_t mapping_row{},assertion{},participant{};};
struct Feature {std::string id,status;std::size_t candidates{};std::vector<std::size_t> rows;std::vector<Match> matches;};
struct Join {std::vector<Feature> features;std::vector<std::string> row_status;};
Request admit_request(const Json &request);
Roles roles(const rfc::database::InterpretationSchema &schema);
std::vector<MappingRow> admit_mapping(const Json &declaration,const std::string &source,
    const Roles &roles,const Limits &limits,Work &work);
std::vector<Assertion> project_assertions(const rfc::database::InterpretedSnapshot &source,
    const Request &request,const Roles &roles,Work &work);
Join join_features(const std::vector<std::string> &features,const std::vector<MappingRow> &mapping,
    const std::vector<Assertion> &assertions,const Limits &limits,Work &work);
} // namespace gadget::feature_context_detail
