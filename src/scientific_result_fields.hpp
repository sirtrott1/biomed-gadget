#pragma once
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>

namespace gadget::scientific_detail {
class ResultValue {
    NativeRecordValue value_;
  public:
    ResultValue() = default;
    explicit ResultValue(NativeRecordValue value) : value_(std::move(value)) {}
    ResultValue(const Json &control) : value_(import_record_view(control)) {}
    static ResultValue object() {return ResultValue(NativeRecordValue(NativeRecordValue::Object{}));}
    static ResultValue array() {return ResultValue(NativeRecordValue(NativeRecordValue::Array{}));}
    static ResultValue string(std::string value) {return ResultValue(NativeRecordValue(std::move(value)));}
    static ResultValue boolean(bool value) {return ResultValue(NativeRecordValue(value));}
    static ResultValue count(std::size_t value) {
        return ResultValue(NativeRecordValue(rfc::numeric::Integer(std::to_string(value))));
    }
    void set(std::string name, ResultValue value) {
        auto &items = std::get<NativeRecordValue::Object>(value_.data).members;
        for (auto &[key, item] : items) if (key == name) {item = std::move(value.value_); return;}
        items.emplace_back(std::move(name), std::move(value.value_));
    }
    void push(ResultValue value) {
        std::get<NativeRecordValue::Array>(value_.data).elements.push_back(std::move(value.value_));
    }
    ResultValue at(std::string_view name) const {
        for (const auto &[key, value] : std::get<NativeRecordValue::Object>(value_.data).members)
            if (key == name) return ResultValue(value);
        throw std::invalid_argument("scientific result member absent");
    }
    ResultValue at(std::size_t index) const {
        return ResultValue(std::get<NativeRecordValue::Array>(value_.data).elements.at(index));
    }
    std::string text() const {return std::get<std::string>(value_.data);}
    std::size_t count_value() const {
        const auto &value = std::get<rfc::numeric::Integer>(value_.data);
        if (value < 0 || !mpz_fits_ulong_p(value.get_mpz_t()) ||
            value > rfc::numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())))
            throw std::invalid_argument("native result count exceeds machine addressability");
        return static_cast<std::size_t>(value.get_ui());
    }
    std::size_t size() const {return std::get<NativeRecordValue::Array>(value_.data).elements.size();}
    NativeRecordValue release() && {return std::move(value_);}
};

class ResultFields {
    static rfc::storage::TemporalEndpoint preparation_endpoint(const NativeRecord &source) {
        const auto original = source.endpoint();
        rfc::storage::TemporalEndpoint endpoint(original->source, original->sequence);
        endpoint.carriers = original->carriers;
        return endpoint;
    }
    NativeRecord source_;
    rfc::storage::TemporalEndpoint endpoint_;
    std::size_t original_fields_, original_carriers_;
    assay_detail::SourceGradeFields fields_;
    std::uint64_t scalar_;
    std::vector<RawByteSource> publications_;
    struct ArchivedPublication { NativeRecord owner; ResultValue reference; };
    std::map<std::pair<const void *, std::uint64_t>, ArchivedPublication> archived_;
    ResultValue bindings_{ResultValue::array()};
    ResultValue reference(NativeRecordValue key, std::string display,
                          std::optional<NativeRecordValue> coordinates = {}) {
        auto value = ResultValue::object();
        value.set("schema", ResultValue::string("gadget.scientific.tensor.view.v1"));
        value.set("tensor", ResultValue(std::move(key)));
        value.set("display", ResultValue::string(std::move(display)));
        if (coordinates) value.set("coordinates", ResultValue(std::move(*coordinates)));
        return value;
    }
  public:
    explicit ResultFields(NativeRecord source)
        : source_(std::move(source)), endpoint_(preparation_endpoint(source_)),
          original_fields_(0), original_carriers_(endpoint_.carriers.size()),
          fields_(endpoint_), scalar_(fields_.explicit_carrier(1)) {
        // Existing declarations and numerical fields stay with their immutable
        // source owner. Preparation only allocates the new result fields.
        for (const auto &tensor : source_.endpoint()->tensors)
            fields_.next_key = std::max(fields_.next_key, rfc::machine::checked_size_add(tensor.key, 1));
    }
    std::uint64_t axis(std::size_t size) {return fields_.explicit_carrier(size);}
    std::uint64_t stalk(rfc::storage::StableCellKey cell, std::size_t size) {
        if (!rfc::storage::StableSourceIndex(endpoint_.source).cell(cell))
            throw std::invalid_argument("scientific stalk source cell absent");
        const auto key = fields_.explicit_carrier(size);
        auto &carrier = endpoint_.carriers.back();
        carrier.kind = rfc::storage::TemporalCarrierKind::CellStalk;
        carrier.grade = cell.grade;
        carrier.cell = cell;
        return key;
    }
    std::uint64_t source_grade(std::size_t grade, std::optional<std::uint64_t> preferred = {});
    std::uint64_t next_source_key() const {
        std::uint64_t key = 1;
        for (const auto &bytes : source_.raw_sources())
            key = std::max(key, rfc::machine::checked_size_add(bytes.key(), 1));
        for (const auto &bytes : publications_)
            key = std::max(key, rfc::machine::checked_size_add(bytes.key(), 1));
        return key;
    }
    ResultValue publication(const NativeRecord &owner) {
        if (!owner.is_root() || !owner.publication_history())
            throw std::invalid_argument("scientific archive requires its actual root publication");
        const auto identity = std::pair<const void *, std::uint64_t>{owner.endpoint().get(),
            owner.publication_history()->last_sequence()};
        const auto previous = archived_.find(identity);
        if (previous != archived_.end()) return previous->second.reference;
        const auto key = next_source_key();
        auto bytes = owner.archive_source(key);
        auto out = ResultValue::object();
        out.set("schema", ResultValue::string("gadget.scientific.publication.view.v1"));
        out.set("publication", ResultValue(NativeRecordValue(
            NativeRecordValue::ByteSpan{key, 0, bytes.bytes().size(), false})));
        out.set("path", ResultValue::array());
        publications_.push_back(std::move(bytes));
        archived_.emplace(identity, ArchivedPublication{owner, out});
        return out;
    }
    static ResultValue member(ResultValue publication, std::initializer_list<std::string> names) {
        auto path = publication.at("path");
        for (const auto &name : names) path.push(ResultValue::string(name));
        publication.set("path", std::move(path));
        return publication;
    }
    static ResultValue member(ResultValue publication, std::size_t index) {
        auto path = publication.at("path");
        path.push(ResultValue::count(index));
        publication.set("path", std::move(path));
        return publication;
    }
    static ResultValue overlay(ResultValue base, ResultValue additions) {
        auto value = ResultValue::object();
        value.set("schema", ResultValue::string("gadget.scientific.view.overlay.v1"));
        value.set("base", std::move(base));
        value.set("additions", std::move(additions));
        return value;
    }
    ResultValue stored(std::uint64_t key, std::string display) {
        const auto original = source_.endpoint();
        const auto matches = [&](const auto &tensor) {return tensor.key == key;};
        if (std::none_of(original->tensors.begin(), original->tensors.end(), matches) &&
            std::none_of(endpoint_.tensors.begin(), endpoint_.tensors.end(), matches))
            throw std::invalid_argument("scientific source tensor absent");
        return reference(NativeRecordValue(rfc::numeric::Integer(std::to_string(key))),
                         std::move(display));
    }
    ResultValue indices(const std::vector<std::size_t> &values, std::uint64_t carrier) {
        Vector coefficients;
        for (auto value : values) coefficients.emplace_back(std::to_string(value));
        return reference(fields_.vector("scientific_coordinates", carrier, coefficients), "indices");
    }
    ResultValue scalar(const Rational &value) {
        return reference(fields_.vector("scientific_scalar", scalar_, Vector{value}), "rational");
    }
    // Workers stage actual native tensors. Publication consumes them in result
    // order after all workers join, preserving deterministic field identities.
    static ResultValue staged(rfc::numeric::SparseRationalTensor tensor,
                              std::optional<std::uint64_t> input, std::uint64_t output,
                              std::string display,
                              std::optional<NativeRecordValue> coordinates = {},
                              std::optional<rfc::storage::TemporalTensorRole> role = {}) {
        auto value = ResultValue::object();
        value.set("schema", ResultValue::string("gadget.scientific.tensor.stage.v1"));
        value.set("tensor", ResultValue(NativeRecordValue(std::move(tensor))));
        if (input) value.set("input", ResultValue::count(*input));
        value.set("output", ResultValue::count(output));
        value.set("display", ResultValue::string(std::move(display)));
        if (coordinates) value.set("coordinates", ResultValue(std::move(*coordinates)));
        if (role) value.set("role", ResultValue::count(static_cast<std::size_t>(*role)));
        return value;
    }
    ResultValue stage_scalar(const Rational &value) const {
        return stage_vector(Vector{value}, scalar_, "rational");
    }
    static ResultValue stage_vector(const Vector &value, std::uint64_t carrier,
                                    std::string display = "vector") {
        return staged(assay_detail::SourceGradeFields::tensor(value.size(),
            [&](auto i) {return value[i];}, rfc::numeric::StructuralSupport::full(value.size())),
            {}, carrier, std::move(display));
    }
    static ResultValue stage_cell_vector(const Vector &value, std::uint64_t carrier,
                                         rfc::storage::StableCellKey cell,
                                         std::string display = "vector") {
        auto result = stage_vector(value, carrier, std::move(display));
        result.set("cell_grade", ResultValue::count(cell.grade.value));
        result.set("cell_key", ResultValue::count(cell.key));
        return result;
    }
    ResultValue cell_vector(const Vector &value, std::uint64_t carrier,
                            rfc::storage::StableCellKey cell, std::string display = "vector") {
        return ResultValue(consume_staged(std::move(stage_cell_vector(value, carrier, cell,
            std::move(display))).release()));
    }
    static ResultValue stage_matrix(const rfc::numeric::ExactMatrix &value,
                                    std::uint64_t input, std::uint64_t output,
                                    std::optional<rfc::storage::TemporalTensorRole> role = {}) {
        return staged(rfc::storage::temporal_matrix_tensor(value), input, output, "sparse_matrix", {}, role);
    }
    static ResultValue stage_selected(const Vector &value, std::uint64_t carrier,
                                     std::size_t size, const std::vector<std::size_t> &rows) {
        if (value.size() != rows.size()) throw std::invalid_argument("staged selected field shape");
        std::vector<rfc::numeric::SparseRationalEntry> entries;
        Vector coordinates;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            coordinates.emplace_back(std::to_string(rows[i]));
            if (value[i] != 0) entries.push_back({rows[i], value[i]});
        }
        rfc::numeric::SparseRationalTensor tensor(rfc::numeric::IndexDomain::linear(size),
            rfc::numeric::StructuralSupport::sparse(size, rows),
            rfc::numeric::SparseRationalVector(size, std::move(entries)));
        return staged(std::move(tensor), {}, carrier, "vector", NativeRecordValue(
            assay_detail::SourceGradeFields::tensor(coordinates.size(),
                [&](auto i) {return coordinates[i];},
                rfc::numeric::StructuralSupport::full(coordinates.size()))));
    }
    ResultValue vector(const Vector &value, std::uint64_t carrier) {
        return reference(fields_.vector("scientific_vector", carrier, value), "vector");
    }
    ResultValue matrix(const rfc::numeric::ExactMatrix &value,
                       std::uint64_t input, std::uint64_t output) {
        return reference(fields_.map("scientific_map", input, output, value), "sparse_matrix");
    }
    NativeRecordValue observation(rfc::numeric::SparseRationalTensor value,
                                   std::uint64_t input, std::uint64_t output) {
        rfc::storage::TemporalTensor tensor;
        tensor.role = rfc::storage::TemporalTensorRole::Observation;
        tensor.input = input;
        tensor.output = output;
        tensor.interpretation = "declared_finite_statistics";
        tensor.tensor = std::move(value);
        return fields_.insert(std::move(tensor));
    }
    ResultValue element(NativeRecordValue key, std::size_t index) {
        auto value = reference(std::move(key), "rational");
        value.set("index", ResultValue::count(index));
        return value;
    }
    NativeRecordValue coordinates(const std::vector<std::size_t> &rows, std::uint64_t carrier) {
        Vector values;
        for (auto row : rows) values.emplace_back(std::to_string(row));
        return fields_.vector("source_coordinates", carrier, values);
    }
    ResultValue selected(const Vector &value, std::uint64_t carrier,
                          const std::vector<std::size_t> &rows, const NativeRecordValue &coordinates) {
        if (rows.size() != value.size()) throw std::invalid_argument("selected scientific field shape");
        const auto original = std::find_if(endpoint_.carriers.begin(), endpoint_.carriers.end(),
            [&](const auto &item) {return item.key == carrier;});
        if (original == endpoint_.carriers.end() ||
            original->kind != rfc::storage::TemporalCarrierKind::SourceGrade)
            throw std::invalid_argument("selected scientific field requires its original source carrier");
        const auto size = original->domain.logical_size();
        std::vector<rfc::numeric::SparseRationalEntry> entries;
        for (std::size_t i = 0; i < rows.size(); ++i)
            if (value[i] != 0) entries.push_back({rows[i], value[i]});
        rfc::numeric::SparseRationalTensor tensor(rfc::numeric::IndexDomain::linear(size),
            rfc::numeric::StructuralSupport::sparse(size, rows),
            rfc::numeric::SparseRationalVector(size, std::move(entries)));
        return reference(fields_.on_carrier("selected_scientific_field", carrier, std::move(tensor)),
                         "vector", coordinates);
    }
    void bind(ResultValue binding) {bindings_.push(std::move(binding));}
    void retain_sources(std::vector<RawByteSource> sources) {
        const auto original = source_.raw_sources();
        for (auto &bytes : sources) {
            const auto found = std::find_if(original.begin(), original.end(),
                [&](const auto &item) { return item.key() == bytes.key(); });
            if (found != original.end()) {
                if (!found->shares_owner(bytes)) throw std::invalid_argument("scientific source byte identity collision");
                continue;
            }
            const auto retained = std::find_if(publications_.begin(), publications_.end(),
                [&](const auto &item) { return item.key() == bytes.key(); });
            if (retained != publications_.end()) {
                if (!retained->shares_owner(bytes)) throw std::invalid_argument("scientific retained byte identity collision");
            } else {
                publications_.push_back(std::move(bytes));
            }
        }
    }
    NativeRecord finish(ResultValue view) {
        return finish_wrapped(std::move(view), {}, {});
    }
    NativeRecord finish_wrapped(ResultValue view, NativeRecordValue::Object wrapper,
                                std::string member) {
        view = ResultValue(consume_staged(std::move(view).release()));
        bindings_ = ResultValue(consume_staged(std::move(bindings_).release()));
        auto root = ResultValue::object();
        root.set("schema", ResultValue::string("gadget.scientific.fields.v1"));
        root.set("view", std::move(view));
        root.set("bindings", std::move(bindings_));
        auto original_bytes = ResultValue::array();
        for (const auto &bytes : source_.raw_sources())
            original_bytes.push(ResultValue(NativeRecordValue(
                NativeRecordValue::ByteSpan{bytes.key(), 0, bytes.bytes().size(), false})));
        root.set("source_bytes", std::move(original_bytes));
        std::vector<rfc::storage::TemporalTensor> fields;
        for (auto i = original_fields_; i < endpoint_.tensors.size(); ++i)
            fields.push_back(std::move(endpoint_.tensors[i]));
        std::vector<rfc::storage::TemporalCarrier> carriers;
        for (auto i = original_carriers_; i < endpoint_.carriers.size(); ++i)
            carriers.push_back(std::move(endpoint_.carriers[i]));
        auto value = std::move(root).release();
        if (!member.empty()) {
            wrapper.members.emplace_back(std::move(member), std::move(value));
            value = NativeRecordValue(std::move(wrapper));
        } else if (!wrapper.members.empty()) {
            throw std::invalid_argument("scientific field wrapper requires its selected member");
        }
        // Complete immutable archives now own the child publications. Their
        // decoded preparation caches end before the final publication allocates.
        archived_.clear();
        return NativeRecord::publish_fields(source_, new_native_record_namespace(), std::move(fields),
            std::move(value), std::move(publications_), std::move(carriers));
    }
  private:
    NativeRecordValue consume_staged(NativeRecordValue value);
};

// Each derived row axis has an explicit lift to its original source grade.
// Independent grades and coefficient axes are never identified by dimension.
struct ResultFieldContext {
    ResultFields &publication;
    struct Grade {
        std::uint64_t source, rows;
        std::size_t size;
        std::vector<std::size_t> coordinates;
        NativeRecordValue coordinate_field;
    };
    std::vector<Grade> grades;
    ResultFieldContext(ResultFields &fields, const NativeRecord &source,
                       const std::vector<std::vector<std::uint32_t>> &coordinates,
                       std::optional<std::pair<std::size_t, std::uint64_t>> field_carrier = {});
    ResultValue field(const Vector &value, std::size_t grade) const;
};
Json result_tensor_view(const NativeRecord &record);
Json result_fields_view(const NativeRecord &record);
NativeRecord retain_scientific_plan(const NativeRecord &source, const NativeRecord &plan,
    const std::string &source_digest, const std::vector<std::pair<std::string, NativeRecord>> &steps);
} // namespace gadget::scientific_detail
