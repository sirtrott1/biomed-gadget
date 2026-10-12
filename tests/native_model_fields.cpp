#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <gadget/resources.hpp>
#include <gadget/workflow.hpp>
#include <rfc/storage/rfcx.hpp>
#include <iostream>

using namespace gadget;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
const rfc::storage::TemporalTensor *field(const NativeRecord &owner, NativeRecord descriptor) {
    const auto key = (descriptor.kind() == NativeRecordKind::Object ? descriptor.at("tensor") : descriptor).integer();
    for (const auto &tensor : owner.endpoint()->tensors) if (key == tensor.key) return &tensor;
    throw std::runtime_error("retained field absent");
}
Rational state(const rfc::numeric::SparseRationalTensor &tensor, std::size_t i) {
    check(tensor.support().contains(i), "required coefficient absent");
    return tensor.state(i).value();
}
NativeRecordValue &member(NativeRecordValue &value, const std::string &name) {
    for (auto &[key, item] : std::get<NativeRecordValue::Object>(value.data).members)
        if (key == name) return item;
    throw std::runtime_error("fixture member absent");
}
void verify_model(const NativeRecord &publication, const NativeRecord &original) {
    const auto fields = publication.at("result").at("fields");
    check(fields.at("schema").text() == "gadget.scientific.fields.v1", "model has no native owner");
    check(rfc::storage::encode_rfcx(fields.endpoint()->source).bytes ==
          rfc::storage::encode_rfcx(original.endpoint()->source).bytes, "original complex changed");
    const auto view = fields.at("view");
    for (std::size_t s = 0; s < view.at("strata").size(); ++s) {
        const auto stratum = view.at("strata").at(s);
        if (stratum.at("status").text() != "completed") continue;
        const auto rows = stratum.at("coordinates");
        const auto &metric = *field(fields, stratum.at("metric_diagonal"));
        check(metric.tensor.support().count() == rows.size(), "selected support replaced by ambient support");
        for (std::size_t m = 0; m < stratum.at("models").size(); ++m) {
            const auto model = stratum.at("models").at(m);
            const auto &design = *field(fields, model.at("design"));
            const auto &frame = *field(fields, model.at("independent_frame"));
            const auto &kernel = *field(fields, model.at("coefficient_kernel"));
            const auto columns = design.tensor.domain().extents()[1];
            const auto free = kernel.tensor.domain().extents()[1];
            check(design.output == frame.output && design.input != frame.input &&
                  design.input == kernel.output && kernel.input != kernel.output,
                  "parameter and frame axes merged by dimension");
            check(columns == model.at("columns").size() && free == model.at("coefficient_nullity").integer(),
                  "complete coefficient family missing");
            for (std::size_t i = 0; i < rows.size(); ++i)
                for (std::size_t k = 0; k < free; ++k) {
                    Rational sum{0};
                    for (std::size_t j = 0; j < columns; ++j)
                        sum += state(design.tensor, i * columns + j) * state(kernel.tensor, j * free + k);
                    check(sum == 0, "coefficient family changes represented field");
                }
            for (std::size_t r = 0; r < model.at("responses").size(); ++r) {
                const auto response = model.at("responses").at(r);
                const auto &beta = *field(fields, response.at("particular_coefficients"));
                const auto &fitted = *field(fields, response.at("represented"));
                const auto &residual = *field(fields, response.at("residual"));
                const auto &input = *field(fields, response.at("input"));
                const auto &normal = *field(fields, response.at("normal_equations"));
                check(fitted.output == metric.output && fitted.tensor.domain() == metric.tensor.domain() &&
                      fitted.tensor.support().count() == rows.size(), "fit lost original carrier or support");
                Rational quadrance{0};
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    const auto source = rows.at(i).integer().get_ui();
                    Rational value{0};
                    for (std::size_t j = 0; j < columns; ++j)
                        value += state(design.tensor, i * columns + j) * state(beta.tensor, j);
                    check(value == state(fitted.tensor, source), "native design action differs");
                    const auto error = state(residual.tensor, source);
                    check(value + error == state(input.tensor, source), "field reconstruction differs");
                    quadrance += error * error * state(metric.tensor, source);
                }
                check(quadrance == state(field(fields, response.at("residual_quadrance"))->tensor, 0),
                      "retained quadrance differs");
                for (std::size_t j = 0; j < columns; ++j) {
                    Rational pairing{0};
                    for (std::size_t i = 0; i < rows.size(); ++i) {
                        const auto source = rows.at(i).integer().get_ui();
                        pairing += state(design.tensor, i * columns + j) * state(metric.tensor, source) *
                                   state(residual.tensor, source);
                    }
                    check(pairing == 0 && state(normal.tensor, j) == 0,
                          "orthogonality or supported zero differs");
                }
            }
        }
    }
}
rfc::numeric::RationalVector values(const NativeRecord &owner, const NativeRecord &descriptor) {
    const auto &tensor = field(owner, descriptor)->tensor;
    rfc::numeric::RationalVector result;
    if (descriptor.has("coordinates")) {
        const auto &coordinates = field(owner, descriptor.at("coordinates"))->tensor;
        for (std::size_t i = 0; i < coordinates.domain().logical_size(); ++i)
            result.push_back(state(tensor, state(coordinates, i).get_num().get_ui()));
    } else {
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i) result.push_back(state(tensor, i));
    }
    return result;
}
rfc::numeric::RationalVector apply(const NativeRecord &owner, const NativeRecord &descriptor,
                                   const rfc::numeric::RationalVector &input) {
    const auto &tensor = field(owner, descriptor)->tensor;
    const auto shape = tensor.domain().extents();
    check(shape.size() == 2 && shape[1] == input.size(), "native matrix input shape differs");
    rfc::numeric::RationalVector output(shape[0]);
    for (std::size_t i = 0; i < output.size(); ++i)
        for (std::size_t j = 0; j < input.size(); ++j) output[i] += state(tensor, i * input.size() + j) * input[j];
    return output;
}
void audit_fields(const NativeRecord &owner, const NativeRecord &view) {
    if (view.kind() == NativeRecordKind::Object) {
        if (view.has("schema") && view.at("schema").text() == "gadget.scientific.tensor.view.v1") {
            const auto &tensor = *field(owner, view);
            const auto encoded = owner.encoded_tensor(tensor.key);
            check(encoded.has_value(), "computed field has no retained native encoding");
            rfc::storage::NativeAccessionRequest request; request.subject = tensor.key;
            const auto pin = owner.accession_directory()->pin(request);
            check(pin.source_namespace() == owner.source_namespace(),
                  "computed field pin has foreign namespace");
            check(pin.sequence() == owner.publication_history()->last_sequence(), "computed field pin lost cutoff");
            check(tensor.output.has_value(), "computed result lost output carrier");
            if (tensor.input) check(tensor.tensor.domain().rank() == 2, "map lost ordered axes");
            return;
        }
        if (view.has("schema") && view.at("schema").text() == "gadget.scientific.publication.view.v1") {
            const auto child = NativeRecord::decode(view.at("publication").bytes());
            audit_fields(child, child.at("view")); audit_fields(child, child.at("bindings"));
            return;
        }
        if (view.has("schema")) check(view.at("schema").text() != "gadget.scientific.tensor.stage.v1", "unpublished staged result escaped");
        for (const auto &name : view.names()) audit_fields(owner, view.at(name));
    } else if (view.kind() == NativeRecordKind::Array) {
        for (std::size_t i = 0; i < view.size(); ++i) audit_fields(owner, view.at(i));
    }
}
void verify_response(const NativeRecord &publication, const NativeRecord &original) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    check(owner.at("schema").text() == "gadget.scientific.fields.v1", "response lost native fields");
    check(rfc::storage::encode_rfcx(owner.endpoint()->source).bytes == rfc::storage::encode_rfcx(original.endpoint()->source).bytes,
          "response changed original complex");
    audit_fields(owner, view);
    check(field(owner, view.at("metric"))->role == rfc::storage::TemporalTensorRole::Metric,
          "response metric lost its native role");
    const auto preparations = view.at("preparations");
    for (std::size_t p = 0; p < preparations.size(); ++p) {
        const auto models = preparations.at(p).at("models");
        for (std::size_t m = 0; m < models.size(); ++m) {
            const auto model = models.at(m), design = model.at("design"), kernel = model.at("coefficient_kernel");
            const auto &a = *field(owner, design), &k = *field(owner, kernel), &f = *field(owner, model.at("independent_frame"));
            check(a.input == k.output && a.output == f.output && a.input != f.input && k.input != k.output,
                  "response coefficient axes merged");
            const auto n = a.tensor.domain().extents()[0], columns = a.tensor.domain().extents()[1], free = k.tensor.domain().extents()[1];
            check(columns == model.at("rank").integer() + model.at("coefficient_nullity").integer(), "response free family incomplete");
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < free; ++j) {
                    Rational sum{0};
                    for (std::size_t c = 0; c < columns; ++c) sum += state(a.tensor, i * columns + c) * state(k.tensor, c * free + j);
                    check(sum == 0, "response kernel changes fitted field");
                }
        }
    }
    for (std::size_t s = 0; s < view.at("states").size(); ++s) {
        const auto fits = view.at("states").at(s).at("fits");
        for (std::size_t i = 0; i < fits.size(); ++i) {
            const auto fit = fits.at(i);
            const auto input = values(owner, fit.at("objective_input")), represented = values(owner, fit.at("objective_represented")), residual = values(owner, fit.at("objective_residual"));
            check(rfc::numeric::add(represented, residual) == input, "response reconstruction differs");
            const auto weighted = apply(owner, view.at("metric"), residual);
            check(rfc::numeric::dot(residual, weighted) == values(owner, fit.at("residual_quadrance"))[0], "response metric quadrance differs");
            check(rfc::numeric::dot(represented, weighted) == 0, "response projection is not metric orthogonal");
            check(rfc::numeric::is_zero(values(owner, fit.at("normal_equations"))), "response normal equations failed");
            for (const auto *reading : {"observed", "represented", "residual"}) {
                const auto sectors = fit.at(reading).at("sectors");
                check(rfc::numeric::add(rfc::numeric::add(values(owner, sectors.at("gradient")), values(owner, sectors.at("curl"))), values(owner, sectors.at("harmonic"))) == values(owner, sectors.at("full")), "response sectors lost field");
            }
        }
    }
    if (view.has("physical_shape")) {
        const auto physical = view.at("physical_shape"), families = physical.at("families"), results = physical.at("results");
        for (std::size_t i = 0; i < results.size(); ++i) {
            const auto result = results.at(i);
            check(rfc::numeric::is_zero(values(owner, result.at("half_gradient_stationarity"))), "physical stationarity failed");
            for (const auto &q : values(owner, result.at("slack"))) check(q >= 0, "physical primal feasibility failed");
            for (const auto &q : values(owner, result.at("multipliers"))) check(q >= 0, "physical dual feasibility failed");
        }
        for (std::size_t i = 0; i < families.size(); ++i) {
            const auto family = families.at(i);
            check(field(owner, family.at("frame"))->input == field(owner, family.at("gram"))->input, "physical Gram owner differs");
            check(field(owner, family.at("source_normals"))->input == field(owner, family.at("frame"))->output, "physical source constraint owner differs");
        }
    }
}
void verify_ensemble(const NativeRecord &publication, const NativeRecord &original) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    check(rfc::storage::encode_rfcx(owner.endpoint()->source).bytes == rfc::storage::encode_rfcx(original.endpoint()->source).bytes, "ensemble changed original source");
    audit_fields(owner, view);
    check(field(owner, view.at("metric"))->role == rfc::storage::TemporalTensorRole::Metric,
          "ensemble metric lost its native role");
    for (std::size_t s = 0; s < view.at("samples").size(); ++s) {
        const auto sample = view.at("samples").at(s);
        const auto full = values(owner, sample.at("full")), gradient = values(owner, sample.at("gradient")), curl = values(owner, sample.at("curl")), harmonic = values(owner, sample.at("harmonic"));
        check(rfc::numeric::add(rfc::numeric::add(gradient, curl), harmonic) == full, "ensemble sectors lost field");
        const auto carrier = field(owner, sample.at("full"))->output;
        for (const auto *sector : {"gradient", "curl", "harmonic"}) check(field(owner, sample.at(sector))->output == carrier, "ensemble sector changed source owner");
    }
    for (std::size_t g = 0; g < view.at("groups").size(); ++g) {
        const auto group = view.at("groups").at(g);
        Rational variance{0};
        for (const auto *sector : {"gradient", "curl", "harmonic"}) variance += values(owner, group.at("ensembles").at(sector).at("variance"))[0];
        check(variance == values(owner, group.at("ensembles").at("full").at("variance"))[0], "ensemble finite variance budget differs");
    }
}
void verify_green(const std::shared_ptr<const ScientificStudy> &study) {
    const auto n = study->endpoint()->source.grade(rfc::core::Grade{0}).cell_count();
    auto query = Json::object(), coordinates = Json::object(), cells = Json::array(), interior = Json::array();
    query.set("operation", Json::string("green_response")); query.set("context", Json::string("all"));
    query.set("field", Json::string("all.signal")); query.set("lambda", Json::string("1"));
    coordinates.set("id", Json::string("81000"));
    for (std::size_t i = 0; i < n; ++i) {
        cells.push(Json::string(std::to_string(rfc::storage::stable_cell_key(study->endpoint()->source,
            {rfc::core::Grade{0}, static_cast<std::uint32_t>(i)}).key)));
        interior.push(Json::integer(static_cast<std::int64_t>(i)));
    }
    coordinates.set("cells", cells); query.set("coordinates", coordinates); query.set("interior", interior);
    const auto maximum = std::numeric_limits<std::size_t>::max();
    rfc::runtime::ExecutionHostPolicy policy;
    policy.name = "green-contract"; policy.limits = {maximum, maximum, maximum, maximum, maximum};
    auto resources = CpuResources::local();
    const auto bytes = resources->coordinator().available_declared_buffer_bytes();
    {
        ScientificSession session(study, resources);
        const auto result = session.execute_green_response(query, rfc::runtime::ExecutionLease(policy));
        check(result.at("result").at("interface_coordinates").size() == 0, "complete interior became a truncated partition");
        const auto field = result.at("result").at("field");
        check(field.size() == n, "Green truncated output coordinates");
        for (std::size_t i = 0; i < n; ++i)
            check(exact_rational(field.at(i).at("value").at("rational").text()) == Rational(std::to_string(i)), "Green line oracle differs");
        check(session.execute(query).at("result").at("factorization_reused").boolean(), "Green discarded source-bound preparation");
        if (n == 3) {
            for (std::size_t i = 1; i < 137; ++i) {
                auto next = Json::parse(query.dump()), space = Json::parse(coordinates.dump());
                space.set("id", Json::string(std::to_string(81000 + i))); next.set("coordinates", space);
                check(session.execute(next).at("result").at("retained_preparations").integer() == static_cast<std::int64_t>(i + 1),
                      "Green omitted a retained solve above prior count");
            }
        } else check(n > 128 && interior.size() > 128, "former Green dimension/interior limits untested");
        auto limited = policy; limited.limits.sparse_entries = 0;
        bool refused = false;
        try {(void)session.execute_green_response(query, rfc::runtime::ExecutionLease(limited));}
        catch (const std::length_error &) {refused = true;}
        check(refused, "cached Green ignored actual sparse-entry policy");
        auto cancel = std::make_shared<std::atomic_bool>(true); refused = false;
        try {(void)session.execute_green_response(query, rfc::runtime::ExecutionLease(policy, cancel));}
        catch (const rfc::runtime::ExecutionCancelled &) {refused = true;}
        check(refused, "Green cancellation published a result");
        auto foreign = Json::parse(query.dump()), space = Json::parse(coordinates.dump());
        auto reversed = Json::array(); for (std::size_t i = 0; i < n; ++i) reversed.push(cells.at(n - i - 1));
        space.set("cells", reversed); foreign.set("coordinates", space); refused = false;
        try {(void)session.execute(foreign);} catch (const std::invalid_argument &) {refused = true;}
        check(refused, "Green accepted a foreign coordinate order");
    }
    check(resources->coordinator().available_declared_buffer_bytes() == bytes, "Green retained resource charge survived owner destruction");
    auto hardware = rfc::runtime::detect_hardware();
    rfc::runtime::CpuResourcePolicy tight; tight.worker_limit = 1; tight.declared_buffer_byte_limit = 1;
    auto small = std::make_shared<CpuResources>(hardware, tight, [hardware] {return hardware;});
    ScientificSession session(study, small);
    bool refused = false;
    try {(void)session.execute_green_response(query, rfc::runtime::ExecutionLease(policy));}
    catch (const rfc::runtime::ExecutionLimitExceeded &error) {refused = error.resource() == "declared buffer bytes";}
    check(refused, "Green did not identify exhausted buffer resource");
}
Rational scalar(const NativeRecord &owner, const NativeRecord &descriptor) {
    if (descriptor.has("index")) return state(field(owner,descriptor)->tensor,descriptor.at("index").integer().get_ui());
    const auto coefficients = values(owner, descriptor);
    check(coefficients.size() == 1, "native scalar shape differs");
    return coefficients.front();
}
void verify_structural(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    audit_fields(owner, view);
    audit_fields(owner, owner.at("bindings"));
    if (view.at("status").text() == "planned") return;
    const auto &metric = *field(owner, view.at("metric"));
    check(metric.role == rfc::storage::TemporalTensorRole::Metric, "structural metric lost role");
    const auto states = view.at("states"), scenarios = view.at("scenarios");
    const auto blocks = view.at("exchange_blocks");
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const auto coordinates = blocks.at(b).at("coordinates");
        const auto assigned = values(owner, states.at(0).at("block_source_coordinates").at(b));
        check(coordinates.size() > 0 && coordinates.size() == assigned.size(),
              "exchange block coordinate witness lost original rows");
        for (std::size_t i = 0; i < coordinates.size(); ++i)
            check(Rational(coordinates.at(i).integer()) == assigned[i], "observed assignment differs from block coordinates");
    }
    bool missing = false, supported_zero = false;
    for (std::size_t c = 0; c < scenarios.size(); ++c) {
        const auto scenario = scenarios.at(c), records = scenario.at("statistics");
        NativeRecord statistics = owner.at("bindings").at(0);
        bool found = false;
        for (std::size_t b = 0; b < owner.at("bindings").size(); ++b) {
            const auto binding = owner.at("bindings").at(b);
            if (binding.has("kind") && binding.at("kind").text() == "structural_reference_statistics" &&
                binding.at("scenario").text() == scenario.at("id").text()) {
                statistics = binding.at("statistics"); found = true; break;
            }
        }
        check(found, "structural score family absent");
        const auto &scores = *field(owner, statistics);
        check(scores.input != scores.output, "reference states and hypotheses merged");
        const auto hypotheses = scenario.at("hypotheses");
        check(scores.tensor.domain().extents() == std::vector<std::size_t>{hypotheses.size(), states.size()},
              "structural score axes changed");
        for (std::size_t state_ = 0; state_ < states.size(); ++state_) {
            const auto action = states.at(state_).at("assignment");
            check(field(owner, action)->input == field(owner, action)->output,
                  "assignment changed original row carrier");
            const auto record = records.at(state_);
            for (std::size_t h = 0; h < hypotheses.size(); ++h) {
                const auto value = record.at("statistics").at(h);
                const auto index = h * states.size() + state_;
                check(scores.tensor.support().contains(index) == (value.kind() != NativeRecordKind::Null),
                      "structural missing score became supported");
                if (value.kind() == NativeRecordKind::Null) missing = true;
                else {
                    check(state(scores.tensor, index) == scalar(owner, value), "structural score owner differs");
                    supported_zero = supported_zero || scalar(owner, value) == 0;
                }
            }
            if (record.has("readings")) {
                const auto fields = record.at("readings").at("fields");
                for (std::size_t f = 0; f < fields.size(); ++f) {
                    const auto value = fields.at(f);
                    const auto original = values(owner, view.at("source_fields").at(value.at("source_field").text()));
                    const auto full = values(owner, value.at("full"));
                    check(apply(owner, action, original) == full, "assignment lost original full field");
                    auto reconstructed = rfc::numeric::add(values(owner, value.at("gradient")),
                        rfc::numeric::add(values(owner, value.at("curl")), values(owner, value.at("harmonic"))));
                    check(reconstructed == full, "structural sectors do not reconstruct");
                    check(field(owner, value.at("full"))->output == field(owner, view.at("source_fields").at(value.at("source_field").text()))->output,
                          "structural field lost original carrier");
                }
            }
        }
        for (std::size_t h = 0; h < hypotheses.size(); ++h) {
            const auto hypothesis = hypotheses.at(h);
            if (hypothesis.at("reference_tail_mass").kind() == NativeRecordKind::Null) continue;
            const auto threshold = scalar(owner, hypothesis.at("observed"));
            Rational tail{0}, ties{0};
            for (std::size_t state_ = 0; state_ < states.size(); ++state_) {
                const auto value = state(scores.tensor, h * states.size() + state_);
                const auto mass = scalar(owner, states.at(state_).at("mass"));
                if (value >= threshold) tail += mass;
                if (value == threshold) ties += mass;
            }
            check(tail == scalar(owner, hypothesis.at("reference_tail_mass")) &&
                  ties == scalar(owner, hypothesis.at("tie_mass")), "native finite tail differs from complete counting law");
        }
    }
    if (view.at("status").text() == "undefined_family_statistic") check(missing, "undefined law lost absent support");
    if (view.at("source_fields").has("all.zero")) check(supported_zero || missing, "zero/undefined fixture lost");
    if (view.has("scenario_sensitivity")) {
        const auto summaries = view.at("scenario_sensitivity");
        for (std::size_t h = 0; h < summaries.size(); ++h) {
            const auto summary = summaries.at(h);
            if (summary.at("status").text() != "completed") continue;
            Rational mean{0}, variance{0};
            for (std::size_t c = 0; c < scenarios.size(); ++c)
                mean += scalar(owner, scenarios.at(c).at("weight")) * scalar(owner, scenarios.at(c).at("hypotheses").at(h).at("observed"));
            for (std::size_t c = 0; c < scenarios.size(); ++c) {
                const Rational delta = scalar(owner, scenarios.at(c).at("hypotheses").at(h).at("observed")) - mean;
                variance += scalar(owner, scenarios.at(c).at("weight")) * delta * delta;
            }
            check(mean == scalar(owner, summary.at("weighted_mean")) && variance == scalar(owner, summary.at("weighted_variance")),
                  "native scenario moments differ");
        }
    }
}
void verify_sections(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    audit_fields(owner, view);
    audit_fields(owner, owner.at("bindings"));
    const auto stalks = owner.at("bindings").at(0).at("stalks");
    std::set<std::uint64_t> axes;
    for (std::size_t c = 0; c < stalks.size(); ++c) {
        const auto stalk = stalks.at(c);
        check(axes.insert(stalk.at("carrier").integer().get_ui()).second, "independent stalks share carriers");
        const auto &inclusion = *field(owner, stalk.at("direct_sum_inclusion"));
        check(inclusion.input == stalk.at("carrier").integer().get_ui(), "stalk inclusion lost actual domain");
        check(inclusion.input != inclusion.output, "stalk/direct sum identified by dimension");
        check(inclusion.tensor.support().count() == inclusion.tensor.domain().logical_size(), "stalk inclusion lost supported zero");
    }
    const auto index = rfc::storage::StableSourceIndex(owner.endpoint()->source);
    const auto carrier_index = rfc::storage::temporal_endpoint_detail::carriers(*owner.endpoint());
    const auto maps = view.at("maps");
    for (std::size_t m = 0; m < maps.size(); ++m) {
        const auto &map = *field(owner, maps.at(m).at("matrix"));
        check(map.owner.kind == rfc::storage::TemporalOwnerKind::Occurrence, "section map lost occurrence owner");
        const auto occurrence = map.owner.occurrence;
        check(occurrence.upper_grade.value == maps.at(m).at("grade").integer() &&
              occurrence.key.value == std::stoull(maps.at(m).at("key").text()), "section map flag identity differs");
        const auto slot = index.occurrence(occurrence).value();
        const auto &arena = owner.endpoint()->source.grade(occurrence.upper_grade);
        const auto upper = rfc::storage::stable_cell_key(owner.endpoint()->source,
            {occurrence.upper_grade, arena.upper_index[slot]});
        const auto lower = rfc::storage::stable_cell_key(owner.endpoint()->source,
            {{static_cast<std::uint16_t>(occurrence.upper_grade.value - 1u)}, arena.lower_index[slot]});
        const bool upward = view.at("construction").text() == "cellular_sheaf";
        check(carrier_index.at(*map.input)->cell == (upward ? lower : upper) &&
              carrier_index.at(*map.output)->cell == (upward ? upper : lower), "section map directions/carriers differ");
    }
    const auto sections = view.at("sections");
    for (std::size_t s = 0; s < sections.size(); ++s) {
        const auto section = sections.at(s);
        const auto input = values(owner, section.at("values"));
        if (section.has("glued")) {
            const auto glued = values(owner, section.at("glued")), residual = values(owner, section.at("residual"));
            check(rfc::numeric::add(glued, residual) == input && rfc::numeric::dot(glued, residual) == 0,
                  "native section reconstruction/orthogonality differs");
            check(apply(owner, section.at("coboundary"), input) == values(owner, section.at("coboundary_defect")),
                  "native section defect differs");
            check(rfc::numeric::is_zero(apply(owner, section.at("coboundary"), glued)), "glued section has nonzero defect");
            check(scalar(owner, section.at("residual_quadrance")) == rfc::numeric::dot(residual, residual), "native section quadrance differs");
        }
        if (section.has("relations")) {
            const auto branches = section.at("relations");
            for (std::size_t b = 0; b < branches.size(); ++b) {
                const auto branch = branches.at(b), legs = branch.at("legs");
                for (std::size_t l = 0; l < legs.size(); ++l) {
                    const auto leg = legs.at(l);
                    bool found = false;
                    for (std::size_t m = 0; m < maps.size(); ++m) {
                        const auto &map = *field(owner, maps.at(m).at("matrix"));
                        if (map.owner.occurrence.upper_grade.value != section.at("grade").integer() ||
                            index.occurrence(map.owner.occurrence).value() != leg.at("slot").integer()) continue;
                        check(apply(owner, maps.at(m).at("matrix"), values(owner, branch.at("upper"))) == values(owner, leg.at("value")),
                              "native relation leg differs from its actual ordered map");
                        found = true; break;
                    }
                    check(found, "native relation leg lost occurrence map");
                }
            }
        }
        if (section.at("input_origin").text() == "admitted source fields") {
            const auto sources = owner.at("bindings").at(1).at("sections");
            bool found = false;
            for (std::size_t b = 0; b < sources.size(); ++b) if (sources.at(b).at("id").text() == section.at("id").text()) {
                rfc::numeric::RationalVector assembled(input.size());
                const auto fields = sources.at(b).at("fields");
                for (std::size_t f = 0; f < fields.size(); ++f) {
                    const auto binding = fields.at(f);
                    check(field(owner, binding.at("realization"))->input == field(owner, binding.at("field"))->output,
                          "source section uses foreign field carrier");
                    assembled = rfc::numeric::add(assembled, apply(owner, binding.at("realization"), values(owner, binding.at("field"))));
                }
                check(assembled == input, "native source section operand differs"); found = true;
            }
            check(found, "source section lost original field bindings");
        }
    }
}
void verify_ablation(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    audit_fields(owner, view); audit_fields(owner, owner.at("bindings"));
    for (std::size_t a = 0; a < view.at("ablations").size(); ++a) {
        const auto variant = view.at("ablations").at(a);
        bool source_changed = false;
        for (std::size_t k = 0; k < variant.at("removed_cells").size(); ++k)
            source_changed = source_changed || variant.at("removed_cells").at(k).size() != 0;
        for (std::size_t g = 0; g < variant.at("grades").size(); ++g) {
            const auto grade = variant.at("grades").at(g);
            const auto &inclusion = *field(owner, grade.at("inclusion"));
            const auto &restriction = *field(owner, grade.at("cochain_restriction"));
            check(inclusion.input == restriction.output && inclusion.output == restriction.input,
                  "ablation inclusion/restriction lost ordered source axes");
            check((inclusion.input != inclusion.output) == source_changed, "actual source owner identity differs");
            const auto &boundary = *field(owner, grade.at("boundary"));
            check(boundary.input == inclusion.input, "retained boundary lost its actual source domain");
            if (grade.has("harmonic_representatives")) {
                const auto &harmonic = *field(owner, grade.at("harmonic_representatives"));
                const auto &included = *field(owner, grade.at("included_representatives"));
                check(harmonic.input == included.input && harmonic.output == inclusion.input &&
                      included.output == inclusion.output, "complete homology family source axes differ");
                const auto &kernel = *field(owner, grade.at("killed_class_coefficients"));
                check(kernel.output == harmonic.input && kernel.input != kernel.output,
                      "killed-class free axes merged with harmonic parameters");
            }
        }
        for (std::size_t f = 0; f < variant.at("fields").size(); ++f) {
            const auto record = variant.at("fields").at(f), sectors = record.at("sectors");
            const auto full = values(owner, sectors.at("full"));
            check(rfc::numeric::add(rfc::numeric::add(values(owner, sectors.at("gradient")),
                  values(owner, sectors.at("curl"))), values(owner, sectors.at("harmonic"))) == full,
                  "ablation full native sectors do not reconstruct");
            check(full == values(owner, record.at("restricted_input")), "retained source input lost coefficients");
            const auto quadrances = record.at("quadrances");
            check(scalar(owner, quadrances.at("gradient")) + scalar(owner, quadrances.at("curl")) +
                  scalar(owner, quadrances.at("harmonic")) == scalar(owner, quadrances.at("full")),
                  "ablation complete quadrance budget differs");
            const auto green = record.at("green");
            check(rfc::numeric::add(values(owner, green.at("structure_term")), values(owner, green.at("field_term"))) ==
                      values(owner, green.at("change")), "native Green cross-source terms do not reconstruct");
            check(rfc::numeric::sub(values(owner, green.at("value")), values(owner, green.at("restricted_original"))) ==
                      values(owner, green.at("change")), "Green transport lost actual retained coordinates");
        }
    }
    for (std::size_t b = 0; b < owner.at("bindings").size(); ++b) {
        const auto binding = owner.at("bindings").at(b);
        if (binding.at("kind").text() != "actual_ablation_source_axes") continue;
        const auto actual = NativeRecord::decode(binding.at("source").at("publication").bytes());
        const auto grades = binding.at("grades");
        for (std::size_t g = 0; g < grades.size(); ++g) {
            const auto grade = grades.at(g);
            check(grade.at("coordinate_carrier").integer() != grade.at("lattice_carrier").integer(),
                  "integral and realized coordinates identified by dimension");
            const auto &metric = *field(actual, actual.at("view").at("grades").at(g).at("metric"));
            check(metric.role == rfc::storage::TemporalTensorRole::Metric &&
                  metric.output == grade.at("source_carrier").integer().get_ui(), "actual source metric carrier missing");
        }
    }
    if (!view.has("original_sheaf")) return;
    auto request = Json::parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[]})");
    request.set("source_namespace", Json::string(std::to_string(publication.source_namespace())));
    request.set("cutoff", Json::integer(static_cast<std::int64_t>(publication.publication_history()->last_sequence())));
    std::vector<NativeRecord> selected_fields;
    for (const auto &path_text : {
            R"(["original_sheaf","sections",0,"glued"])",
            R"(["ablations",0,"sheaf","sections",0,"glued"])",
            R"(["ablations",0,"sheaf","sections",0,"projection_order_defect"])",
            R"(["ablations",0,"sheaf","sections",0,"projection_order_quadrance"])"}) {
        const auto path = Json::parse(path_text);
        const auto selected = select_scientific_result(view, admit_record_view(path));
        const auto encoded = selected.encoded_tensor(selected.at("tensor").integer().get_ui());
        check(encoded.has_value(), "archived or overlay claim lost actual tensor owner");
        check(encoded->tensor().support().contains(0), "claim fixture has no original supported coordinate");
        auto item = Json::object();
        item.set("path", path); item.set("coordinate", Json::integer(0));
        item.set("value_type", Json::string("rational")); item.set("comparison", Json::string("equal"));
        item.set("value", Json::string(encoded->tensor().value(0)->get_str()));
        request.at("items").push(item); selected_fields.push_back(selected);
    }
    const auto receipt = check_certificate_claims_record(publication, admit_record_view(request)).at("result");
    check(receipt.at("matched").integer() == selected_fields.size() && receipt.at("contradicted").integer() == 0,
          "native archived and overlay literal comparisons differ");
    for (std::size_t i = 0; i < selected_fields.size(); ++i) {
        const auto &selected = selected_fields[i];
        const auto encoded = selected.encoded_tensor(selected.at("tensor").integer().get_ui());
        const auto binding = receipt.at("checks").at(i).at("binding");
        check(binding.at("source_namespace").text() == std::to_string(encoded->publication().source_namespace()) &&
              binding.at("cutoff").integer() == encoded->publication().sequence(), "claim receipt changed actual selected owner");
    }
    check(selected_fields[0].source_namespace() != publication.source_namespace() &&
          selected_fields[1].source_namespace() != publication.source_namespace(), "child publications became enclosing fields");
    check(selected_fields[2].source_namespace() == publication.source_namespace(), "overlay additions lost enclosing field owner");
    auto foreign = Json::parse(request.dump());
    foreign.set("source_namespace", Json::string(std::to_string(selected_fields[0].source_namespace())));
    bool refused = false;
    try { (void)check_certificate_claims_record(publication, admit_record_view(foreign)); }
    catch (const std::invalid_argument &) { refused = true; }
    check(refused, "child namespace substituted for actual operation claim owner");
    for (const auto &bad_path : {R"(["ablations",-1])", R"(["original_sheaf","sections",999999,"glued"])",
                                 R"(["ablations",0,"sheaf","sections",0,"absent"])", R"(["ablations","0"])"}) {
        refused = false;
        try { (void)select_scientific_result(view, admit_record_view(Json::parse(bad_path))); }
        catch (const std::exception &) { refused = true; }
        check(refused, "invalid archived or overlay path accepted");
    }
}
void verify_direct(const NativeRecord &publication, const NativeRecord &original) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    audit_fields(owner, view); audit_fields(owner, owner.at("bindings"));
    const auto binding = owner.at("bindings").at(0);
    if (binding.at("kind").text() == "actual_reported_field_summary") {
        const auto &input = *field(owner, binding.at("source_field"));
        Rational sum{0}, squares{0}; std::size_t count{0};
        input.tensor.support().for_each_active([&](std::size_t i) {
            const auto x = *input.tensor.state(i); sum += x; squares += x*x; ++count;
        });
        check(view.at("reported").integer() == count &&
              view.at("missing").integer() == input.tensor.domain().logical_size()-count,
              "summary changed original supported versus missing coordinates");
        check(scalar(owner,view.at("sum")) == sum && scalar(owner,view.at("sum_of_squares")) == squares,
              "summary differs from original complete supported field");
        const auto &mean = field(owner,view.at("mean"))->tensor;
        check(mean.support().count() == (count ? 1u : 0u), "empty mean became supported zero");
        if (count) check(state(mean,0) == sum/Rational(count), "native mean differs");
        auto request = Json::parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[{"path":["mean"],"value_type":"rational","comparison":"equal","value":"0"}]})");
        request.set("source_namespace", Json::string(std::to_string(publication.source_namespace())));
        request.set("cutoff", Json::integer(static_cast<std::int64_t>(publication.publication_history()->last_sequence())));
        if (count) {
            request.at("items").at(0).set("value", Json::string(Rational(sum/Rational(count)).get_str()));
            check(check_certificate_claims_record(publication,admit_record_view(request)).at("result").at("matched").integer()==1,
                  "rational object native claim lost actual mean owner");
        } else {
            bool refused=false;
            try {(void)check_certificate_claims_record(publication,admit_record_view(request));}
            catch(const std::invalid_argument &e){refused=std::string(e.what()).find("missing")!=std::string::npos;}
            check(refused,"unsupported empty mean accepted as supported zero claim");
        }
        check(field(owner,binding.at("reported_values"))->output == field(owner,binding.at("weights"))->output,
              "summary rows and weights lost actual reported coordinate axis");
    } else {
        check(binding.at("kind").text() == "actual_direct_source_action", "direct action owner absent");
        const auto &input = *field(original,binding.at("input_tensor"));
        check(binding.at("input_carrier").integer() == *input.output, "original direct input carrier differs");
        check(NativeRecord::decode(binding.at("original_study_publication").bytes()).encode() == original.encode(),
              "direct action lost exact original study archive");
        const auto context_input = values(owner,binding.at("context_input"));
        const auto &map = field(owner,binding.at("input_realization"))->tensor;
        const auto shape = map.domain().extents();
        check(shape.size()==2 && shape[0]==context_input.size() && shape[1]==input.tensor.domain().logical_size(),
              "direct context realization shape differs");
        for (std::size_t i=0;i<shape[0];++i) {
            Rational value{0};
            for (std::size_t j=0;j<shape[1];++j) {
                const auto coefficient=state(map,i*shape[1]+j);
                if(coefficient!=0)value+=coefficient*state(input.tensor,j);
            }
            check(value==context_input[i], "context input differs from original tensor and actual lift");
        }
        const auto output=values(owner,binding.at("output"));
        check(output.size()==view.at("result").size(), "direct result lost complete coordinates");
        for(std::size_t i=0;i<output.size();++i) {
            const auto coordinate=view.at("result").at(i);
            check(coordinate.at("coordinate").integer()==i && scalar(owner,coordinate.at("value"))==output[i],
                  "direct native coordinate view changed actual output field");
        }
        std::size_t metrics{};
        for(const auto &t:owner.endpoint()->tensors)if(t.role==rfc::storage::TemporalTensorRole::Metric)++metrics;
        check(metrics==owner.endpoint()->source.grade_count(), "direct source lacks complete declared metrics");
    }
}
void verify_censored(const NativeRecord &publication, const NativeRecord &original, const NativeRecord &query) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    audit_fields(owner, view); audit_fields(owner, owner.at("bindings"));
    const auto binding = owner.at("bindings").at(0);
    check(binding.at("kind").text() == "actual_censored_endpoint_inputs", "censored actual input binding absent");
    check(std::ranges::equal(binding.at("original_study_publication").bytes(), original.encode()),
          "censored original publication changed");
    const auto source_axis = binding.at("source_carrier").integer().get_ui();
    const auto group_axis = binding.at("group_carrier").integer().get_ui();
    const auto contrast_axis = binding.at("contrast_carrier").integer().get_ui();
    check(source_axis != group_axis && group_axis != contrast_axis && source_axis != contrast_axis,
          "censored independently declared axes identified by dimension");
    const auto &lift = *field(owner, binding.at("input_realization"));
    const auto n = lift.tensor.domain().extents().at(0), ambient = lift.tensor.domain().extents().at(1);
    std::map<std::string, const rfc::storage::TemporalTensor *> inputs;
    for (std::size_t f = 0; f < binding.at("inputs").size(); ++f) {
        const auto input = binding.at("inputs").at(f);
        const auto &old = *field(original, input.at("tensor"));
        const auto &context = *field(owner, input.at("context_values"));
        check(old.output && input.at("carrier").integer() == *old.output && context.output == source_axis,
              "censored field carrier owner changed");
        rfc::storage::NativeAccessionRequest request; request.subject = old.key;
        const auto pin = original.accession_directory()->pin(request);
        check(input.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
              input.at("cutoff").integer() == pin.sequence(), "censored original pin/cutoff changed");
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t count{}; bool present{}; Rational expected{};
            for (std::size_t j = 0; j < ambient; ++j) if (state(lift.tensor, i * ambient + j) != 0) {
                ++count; check(state(lift.tensor, i * ambient + j) == 1, "censored lift not identity selection");
                present = old.tensor.support().contains(j);
                if (present) expected = *old.tensor.state(j);
            }
            check(count == 1 && context.tensor.support().contains(i) == present,
                  "censored original upper missing support became zero");
            if (present) check(state(context.tensor, i) == expected, "censored context field differs from original");
        }
        inputs.emplace(input.at("id").text(), &context);
    }
    const auto input_values = [&](const std::string &id) {
        rfc::numeric::RationalVector x(n);
        for (std::size_t i = 0; i < n; ++i) if (inputs.at(id)->tensor.support().contains(i)) x[i] = state(inputs.at(id)->tensor, i);
        return x;
    };
    const auto lower = input_values(query.at("lower").text()), upper = input_values(query.at("upper").text());
    const auto closed = input_values(query.at("lower_closed").text()), finite = input_values(query.at("upper_finite").text());
    rfc::numeric::RationalVector selected(n, Rational{1});
    for (std::size_t k = 0; k < query.at("selection").size(); ++k) {
        const auto x = input_values(query.at("selection").at(k).text());
        for (std::size_t i = 0; i < n; ++i) selected[i] *= x[i];
    }
    const auto groups = view.at("groups");
    const auto &readout = *field(owner, view.at("readout"));
    check(readout.input == source_axis && readout.output == group_axis, "censored readout axes differ");
    for (std::size_t g = 0; g < groups.size(); ++g) {
        auto membership = selected;
        if (g) {
            const auto x = input_values(query.at("strata").at(g-1).text());
            for (std::size_t i = 0; i < n; ++i) membership[i] *= x[i];
        }
        Rational mass{}; for (const auto &x : membership) mass += x;
        check(groups.at(g).at("defined").boolean() == (mass != 0) && scalar(owner, groups.at(g).at("size")) == mass,
              "censored group support/mass differs");
        for (std::size_t i = 0; i < n; ++i)
            check(state(readout.tensor, g*n+i) == (mass != 0 ? Rational(membership[i]/mass) : Rational{0}),
                  "censored complete original group map differs");
    }
    const Rational origin(query.at("origin").text());
    for (std::size_t t = 0; t < view.at("readings").size(); ++t) {
        const auto reading = view.at("readings").at(t); const Rational time(query.at("times").at(t).text());
        rfc::numeric::RationalVector lo(n), hi(n), rlo(n), rhi(n);
        std::size_t ambiguity{};
        for (std::size_t i = 0; i < n; ++i) if (selected[i] != 0) {
            lo[i] = finite[i] != 0 && upper[i] <= time ? 1 : 0;
            hi[i] = lower[i] < time || (lower[i] == time && closed[i] == 1) ? 1 : 0;
            rlo[i] = std::min(lower[i], time) - origin;
            rhi[i] = (finite[i] != 0 ? std::min(upper[i], time) : time) - origin;
            ambiguity += lo[i] != hi[i];
        }
        check(scalar(owner, reading.at("time")) == time && values(owner,reading.at("event_field_lower")) == lo &&
              values(owner,reading.at("event_field_upper")) == hi && values(owner,reading.at("restricted_time_lower_field")) == rlo &&
              values(owner,reading.at("restricted_time_upper_field")) == rhi, "censored native marginal fields differ");
        for (const auto &[name,x] : std::vector<std::pair<std::string,rfc::numeric::RationalVector>>{
                {"cumulative_event_fraction_lower",lo},{"cumulative_event_fraction_upper",hi},
                {"restricted_event_free_time_infimum",rlo},{"restricted_event_free_time_supremum",rhi}})
            check(values(owner,reading.at(name)) == apply(owner,view.at("readout"),x), "censored original aggregate differs");
        const auto family = reading.at("admissible_event_family"); const auto &frame = *field(owner,family.at("frame"));
        check(frame.input != frame.output && frame.output == source_axis && frame.tensor.domain().extents().at(1) == ambiguity,
              "censored full ambiguity frame/carrier differs");
        const auto sectors = family.at("sectors");
        for (std::size_t i = 0; i < n; ++i) {
            Rational reconstructed{};
            for (const auto *name : {"gradient","curl","harmonic"}) reconstructed += values(owner,sectors.at(name).at("offset")).at(i);
            check(reconstructed == lo[i], "censored sector offsets fail reconstruction");
            for (std::size_t j = 0; j < ambiguity; ++j) {
                Rational column{};
                for (const auto *name : {"gradient","curl","harmonic"}) column += state(field(owner,sectors.at(name).at("frame"))->tensor,i*ambiguity+j);
                check(column == state(frame.tensor,i*ambiguity+j), "censored full sector ambiguity reconstruction differs");
            }
        }
    }
    if (!view.has("joint_family")) return;
    const auto joint = view.at("joint_family"), states = joint.at("states");
    const auto metric = view.at("metric");
    const auto inner = [&](const auto &a, const auto &b) {
        const auto mb = apply(owner,metric,b); Rational result{};
        for (std::size_t i = 0; i < n; ++i) result += a[i]*mb[i];
        return result;
    };
    rfc::numeric::RationalVector weights;
    for (std::size_t s = 0; s < states.size(); ++s) {
        const auto record = states.at(s);
        weights.push_back(scalar(owner,record.at("weight")));
        check(weights.back() == Rational(query.at("joint_family").at("states").at(s).at("weight").text()),
              "joint finite law weight changed");
        const auto &times = field(owner,record.at("event_times"))->tensor;
        const auto &presence = field(owner,record.at("event_presence"))->tensor;
        const auto declared = query.at("joint_family").at("states").at(s).at("values");
        for (std::size_t e = 0; e < declared.size(); ++e) {
            const auto entry = declared.at(e); std::optional<std::size_t> coordinate;
            for (std::size_t i = 0; i < n; ++i)
                if (view.at("coordinate_carrier").at(i).at("id").text() == entry.at("entity").text()) coordinate = i;
            check(coordinate.has_value(),"joint original entity disappeared");
            const auto event = entry.at("event_time");
            check(times.support().contains(*coordinate) == (event.kind() != NativeRecordKind::Null),
                  "joint declared no event support changed");
            if (event.kind() != NativeRecordKind::Null)
                check(state(times,*coordinate) == Rational(event.text()), "joint exact original event time or supported zero changed");
        }
        for (std::size_t i = 0; i < n; ++i) {
            check(presence.support().contains(i) == (selected[i] != 0), "joint state original selection support differs");
            if (selected[i] != 0) check(times.support().contains(i) == (state(presence,i) == 1),
                                      "joint no event became supported zero time");
        }
    }
    constexpr const char *names[] = {"full","gradient","curl","harmonic"};
    for (const auto &[moment_key,sector_key] : std::vector<std::pair<std::string,std::string>>{
            {"event_moments","event_sectors"},{"restricted_time_moments","restricted_time_sectors"}}) {
        const auto moment = joint.at(moment_key), channels = moment.at("channels"); const auto count = channels.size();
        std::vector<std::vector<rfc::numeric::RationalVector>> samples(count);
        std::vector<rfc::numeric::RationalVector> means(count,rfc::numeric::RationalVector(n));
        for (std::size_t s = 0; s < states.size(); ++s)
            for (std::size_t c = 0; c < count; ++c) {
                auto x = values(owner,states.at(s).at("readings").at(c/4).at(sector_key).at(names[c%4]));
                samples[c].push_back(x); for (std::size_t i = 0; i < n; ++i) means[c][i] += weights[s]*x[i];
            }
        const auto &expected = field(owner,moment.at("expected_gram"))->tensor;
        const auto &mean_gram = field(owner,moment.at("mean_gram"))->tensor;
        const auto &centered = field(owner,moment.at("centered_gram"))->tensor;
        for (std::size_t a = 0; a < count; ++a) {
            check(values(owner,channels.at(a).at("mean_field")) == means[a], "joint complete channel mean differs");
            const auto &family = field(owner,moment.at("field_families").at(a).at("centered_fields"))->tensor;
            Rational variance{};
            for (std::size_t s = 0; s < states.size(); ++s) {
                auto delta = samples[a][s];
                for (std::size_t i = 0; i < n; ++i) {
                    delta[i] -= means[a][i]; check(state(family,s*n+i) == delta[i], "joint centered field differs");
                }
                variance += weights[s]*inner(delta,delta);
            }
            check(scalar(owner,moment.at("field_families").at(a).at("variance")) == variance,
                  "joint complete finite variance differs");
            for (std::size_t b = 0; b < count; ++b) {
                Rational second{}; for (std::size_t s = 0; s < states.size(); ++s) second += weights[s]*inner(samples[a][s],samples[b][s]);
                const auto first = inner(means[a],means[b]);
                check(state(expected,a*count+b) == second && state(mean_gram,a*count+b) == first && state(centered,a*count+b) == second-first,
                      "joint complete cross time/sector moments differ");
            }
        }
    }
}

void verify_censored_claim(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    auto request = Json::parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[{"path":["readings",0,"time"],"value_type":"rational","comparison":"equal","value":"0"}]})");
    request.set("source_namespace", Json::string(std::to_string(publication.source_namespace())));
    request.set("cutoff", Json::integer(static_cast<std::int64_t>(publication.publication_history()->last_sequence())));
    request.at("items").at(0).set("value", Json::string(scalar(owner,view.at("readings").at(0).at("time")).get_str()));
    check(check_certificate_claims_record(publication,admit_record_view(request)).at("result").at("matched").integer() == 1,
          "censored native scalar literal claim lost its actual field owner");
    if (!view.has("joint_family")) return;
    const auto states = view.at("joint_family").at("states");
    for (std::size_t s = 0; s < states.size(); ++s) {
        const auto &tensor = field(owner,states.at(s).at("event_times"))->tensor;
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i) {
            auto path = Json::array();
            path.push(Json::string("joint_family")); path.push(Json::string("states"));
            path.push(Json::integer(static_cast<std::int64_t>(s))); path.push(Json::string("event_times"));
            request.at("items").at(0).set("path",path);
            request.at("items").at(0).set("coordinate",Json::integer(static_cast<std::int64_t>(i)));
            if (tensor.support().contains(i)) {
                request.at("items").at(0).set("value",Json::string(tensor.state(i)->get_str()));
                check(check_certificate_claims_record(publication,admit_record_view(request)).at("result").at("matched").integer() == 1,
                      "joint original coordinate literal claim differs");
            } else {
                request.at("items").at(0).set("value",Json::string("0"));
                bool refused{};
                try {(void)check_certificate_claims_record(publication,admit_record_view(request));}
                catch (const std::invalid_argument &error) {refused = std::string(error.what()).find("missing") != std::string::npos;}
                check(refused,"joint no event/unselected support accepted as zero event time");
            }
        }
    }
}
void verify_censored_workflow(const std::shared_ptr<const ScientificStudy> &study,
                             const NativeRecord &query, const fs::path &output) {
    auto declaration = Json::parse(R"({"schema":"gadget.scientific.workflow.v1","id":"censored-families","source_digest":"","limits":{"depth":4,"nodes":3,"model_calls":0,"parallel_workers":2},"definitions":[],"entry":{"id":"two-families","type":"parallel","workers":2,"steps":[]}})");
    declaration.set("source_digest",Json::string(study->digest()));
    for (const auto *id : {"first","second"}) {
        auto node = Json::object(); node.set("id",Json::string(id)); node.set("type",Json::string("query"));
        node.set("query",render_record_view(query)); declaration.at("entry").at("steps").push(node);
    }
    const auto request = admit_record_view(declaration);
    const auto operation = execute_workflow_record(study,request);
    check(operation.at("result").at("native_execution").boolean() && operation.at("result").at("status").text() == "complete",
          "censored producer refused native workflow");
    const auto results = operation.at("result").at("result");
    check(results.size() == 2,"censored native workflow lost an independent family");
    const auto first = NativeRecord::decode(results.at(0).at("publication").bytes());
    const auto second = NativeRecord::decode(results.at(1).at("publication").bytes());
    check(first.source_namespace() != second.source_namespace(),"censored independent results merged by dimension");
    for (const auto &result : {first,second}) {
        check(NativeRecord::decode(result.at("request_publication").bytes()).encode() == request.encode(),
              "censored workflow leaf lost its original request owner");
        verify_censored(result,study->native_publication(),query); verify_censored_claim(result);
    }
    operation.write(output/"owning-workflow.rfcxs");
    const auto reopened = NativeRecord::read(output/"owning-workflow.rfcxs");
    check(NativeRecord::decode(reopened.at("result").at("result").at(0).at("publication").bytes()).encode() == first.encode(),
          "censored reopened workflow changed complete native publication");
}

void verify_replay(const NativeRecord &expected, const NativeRecord &actual) {
    const auto a = expected.endpoint(), b = actual.endpoint();
    check(a->tensors.size() == b->tensors.size(), "replay field family differs");
    for (std::size_t i = 0; i < a->tensors.size(); ++i)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(a->tensors[i], b->tensors[i]),
              "replay actual tensor differs");
    check(std::ranges::equal(expected.at("request_publication").bytes(), actual.at("request_publication").bytes()),
          "replay retained a result archive as the request owner");
    check(std::ranges::equal(expected.at("result").at("source_study_publication").bytes(),
                            actual.at("result").at("source_study_publication").bytes()),
          "replay changed or extended the original study archive");
    check(expected.raw_sources().size() == actual.raw_sources().size(), "replay retained extra byte owners");
}
void verify_workflow(const std::shared_ptr<const ScientificStudy> &study,
                     const NativeRecord &query, const fs::path &output, bool sections) {
    auto declaration = Json::parse(R"({"schema":"gadget.scientific.workflow.v1","id":"native-producers","source_digest":"","limits":{"depth":4,"nodes":3,"model_calls":0,"parallel_workers":2},"definitions":[],"entry":{"id":"two-owners","type":"parallel","workers":2,"steps":[]}})");
    declaration.set("source_digest", Json::string(study->digest()));
    for (const auto *id : {"first", "second"}) {
        auto leaf = Json::object();
        leaf.set("id", Json::string(id)); leaf.set("type", Json::string("query"));
        leaf.set("query", render_record_view(query));
        declaration.at("entry").at("steps").push(leaf);
    }
    const auto request = admit_record_view(declaration);
    const auto workflow = execute_workflow_record(study, request);
    const auto receipt = workflow.at("result");
    check(receipt.at("native_execution").boolean() && receipt.at("status").text() == "complete",
          "accepted producer refused native workflow execution");
    check(receipt.at("result").size() == 2, "parallel workflow omitted a native owner");
    const auto first = NativeRecord::decode(receipt.at("result").at(0).at("publication").bytes());
    const auto second = NativeRecord::decode(receipt.at("result").at(1).at("publication").bytes());
    check(first.source_namespace() != second.source_namespace(), "independent results merged by dimension");
    for (const auto &owner : {first, second}) {
        check(NativeRecord::decode(owner.at("request_publication").bytes()).encode() == request.encode(),
              "leaf lost complete original workflow request");
        check(NativeRecord::decode(owner.at("result").at("source_study_publication").bytes()).encode() ==
                  study->native_publication().encode(), "workflow leaf changed original study bytes");
        if (sections) verify_sections(owner); else verify_structural(owner);
    }
    request.write(output / "workflow-request.rfcxs");
    workflow.write(output / "owning-workflow.rfcxs");
    const auto reopened = NativeRecord::read(output / "owning-workflow.rfcxs");
    check(NativeRecord::decode(reopened.at("result").at("result").at(0).at("publication").bytes()).encode() ==
              first.encode(), "workflow reopening changed complete native fields");
    write_text(output / "workflow-display.json", workflow_record_view(reopened).dump());
    auto unsupported = Json::object();
    unsupported.set("operation", Json::string("coverage"));
    declaration.at("entry").at("steps").at(1).set("query", unsupported);
    bool refused = false;
    try { (void)execute_workflow_record(study, admit_record_view(declaration)); }
    catch (const std::invalid_argument &error) {
        refused = std::string(error.what()).find("native result publisher") != std::string::npos;
    }
    check(refused, "unconverted publisher passed strict admission");
}
}
int main(int argc, char **argv) try {
    if (argc != 5) throw std::invalid_argument("workspace, queries, output and mode required");
    const fs::path workspace_path(argv[1]), queries_path(argv[2]), output(argv[3]);
    fs::create_directories(output);
    const auto study = std::string(argv[4]) == "study" || workspace_path.filename() == "study.json"
        ? std::make_shared<const ScientificStudy>(workspace_path)
        : AssayWorkspace(workspace_path).study("measurements");
    if (std::string(argv[4]) == "green") {
        verify_green(study);
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
        return 0;
    }
    const auto original = study->native_publication();
    ScientificSession session(study);
    const auto queries = admit_record_view(Json::parse(read_text(queries_path)));
    check(!has_native_scientific_publisher("coverage") && !has_native_scientific_publisher("unknown"),
          "unconverted producer reported native availability");
    std::vector<NativeRecord> results;
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name);
        check(has_native_scientific_publisher(query.at("operation").text()),
              "implemented producer missing from native inventory");
        const auto result = session.execute_record(query);
        if (std::string(argv[4]) == "response") verify_response(result, original);
        else if (std::string(argv[4]) == "ensemble") verify_ensemble(result, original);
        else if (std::string(argv[4]) == "structural") verify_structural(result);
        else if (std::string(argv[4]) == "sections") verify_sections(result);
        else if (std::string(argv[4]) == "ablation") verify_ablation(result);
        else if (std::string(argv[4]) == "direct") verify_direct(result, original);
        else if (std::string(argv[4]) == "censored") {verify_censored(result, original, query); verify_censored_claim(result);}
        else if (std::string(argv[4]) != "finite") verify_model(result, original);
        else {
            const auto fields = result.at("result").at("fields");
            const auto binding = fields.at("bindings").at(0);
            const auto &statistics = *field(fields, binding.at("statistics"));
            check(statistics.input != statistics.output, "reference law axes merged");
            if (name == "undefined") check(statistics.tensor.support().count() < statistics.tensor.domain().logical_size(),
                                           "undefined statistic became zero");
        }
        const auto retained = NativeRecord::decode(result.at("result").at("source_study_publication").bytes());
        check(retained.encode() == original.encode(), "original study bytes changed");
        result.write(output / (name + ".rfcxs"));
        results.push_back(NativeRecord::read(output / (name + ".rfcxs")));
    }
    if (std::string(argv[4]) == "sections" || std::string(argv[4]) == "structural")
        verify_workflow(study, queries.at(queries.names().front()), output, std::string(argv[4]) == "sections");
    if (std::string(argv[4]) == "censored")
        verify_censored_workflow(study, queries.at(queries.names().front()), output);
    auto plan = Json::object(), steps = Json::array();
    plan.set("schema", Json::string("gadget.scientific.plan.v1"));
    plan.set("id", Json::string("native-models")); plan.set("question", Json::string("Retain declared model fields"));
    plan.set("source_digest", Json::string(study->digest()));
    for (const auto &name : queries.names()) {
        auto step = Json::object(); step.set("id", Json::string(name)); step.set("query", render_record_view(queries.at(name))); steps.push(step);
    }
    plan.set("steps", steps);
    const auto planned = session.execute_plan_record(admit_record_view(plan));
    planned.write(output / "plan.rfcxs");
    const auto moved = workspace_path.parent_path().string() + ".moved";
    fs::rename(workspace_path.parent_path(), moved);
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto replayed = replay_scientific_record(results[i]);
        verify_replay(results[i], replayed);
        replayed.write(output / ("replay-" + std::to_string(i) + ".rfcxs"));
        const auto repeated = replay_scientific_record(replayed);
        verify_replay(replayed, repeated);
    }
    const auto replayed_plan = replay_scientific_record(NativeRecord::read(output / "plan.rfcxs"));
    check(replayed_plan.at("result").at("steps").size() == queries.size(), "plan replay lost step");
    check(std::ranges::equal(planned.at("request_publication").bytes(), replayed_plan.at("request_publication").bytes()),
          "plan replay retained the enclosing result as request owner");
    check(planned.raw_sources().size() == replayed_plan.raw_sources().size(), "plan replay retained extra byte owners");
    replayed_plan.write(output / "replay-plan.rfcxs");
    auto changed_receipt = results.front().unbound_value();
    member(member(member(changed_receipt, "result"), "certificate"), "gadget_sources_sha256") =
        NativeRecordValue(std::string("different software receipt"));
    const auto relocated = results.front().derive(std::move(changed_receipt));
    verify_replay(results.front(), replay_scientific_record(relocated));
    auto changed_request = results.front().unbound_value();
    member(member(changed_request, "request"), "operation") = NativeRecordValue(std::string("unretained_operation"));
    bool request_refused = false;
    try {(void)replay_scientific_record(results.front().derive(std::move(changed_request)));}
    catch (const std::invalid_argument &error) {
        request_refused = std::string(error.what()).find("request owner") != std::string::npos;
    }
    check(request_refused, "query detached from original request owner was accepted");
    auto endpoint = *results.front().endpoint();
    const auto candidate = std::find_if(endpoint.tensors.rbegin(), endpoint.tensors.rend(), [](const auto &item) {
        return item.role != rfc::storage::TemporalTensorRole::Metric && item.tensor.domain().logical_size() > 0 && item.tensor.support().contains(0);
    });
    check(candidate != endpoint.tensors.rend(), "tamper fixture has no supported result");
    const auto tamper_index = static_cast<std::size_t>(std::distance(endpoint.tensors.begin(), candidate.base()) - 1);
    auto &tensor = endpoint.tensors.at(tamper_index).tensor;
    const auto size = tensor.domain().logical_size();
    check(size > 0 && tensor.support().contains(0), "tamper fixture has no supported field");
    auto entries = tensor.values().entries();
    std::vector<rfc::numeric::SparseRationalEntry> changed(entries.begin(), entries.end());
    const auto found = std::find_if(changed.begin(), changed.end(), [](const auto &entry) {return entry.index == 0;});
    if (found == changed.end()) changed.insert(changed.begin(), {0, Rational{1}});
    else found->value += 1;
    tensor = rfc::numeric::SparseRationalTensor(tensor.domain(), tensor.support(),
        rfc::numeric::SparseRationalVector(size, std::move(changed)));
    endpoint.declarations.reset();
    const auto forged = NativeRecord::create_from_endpoint(results.front().unbound_value(), new_native_record_namespace(),
                                                           results.front().raw_sources(), std::move(endpoint));
    forged.write(output / "changed-field.rfcxs");
    bool refused = false;
    try {(void)replay_scientific_record(forged);} catch (const std::invalid_argument &) {refused = true;}
    check(refused, "changed field replay accepted with unchanged source digest");
    auto absent_endpoint = *results.front().endpoint();
    auto &absent = absent_endpoint.tensors.at(tamper_index).tensor;
    auto stored = absent.values().entries();
    std::vector<rfc::numeric::SparseRationalEntry> retained;
    for (const auto &entry : stored) if (entry.index != 0) retained.push_back(entry);
    std::vector<std::size_t> support;
    absent.support().for_each_active([&](std::size_t i) {if (i != 0) support.push_back(i);});
    absent = rfc::numeric::SparseRationalTensor(absent.domain(),
        rfc::numeric::StructuralSupport::sparse(absent.domain().logical_size(), std::move(support)),
        rfc::numeric::SparseRationalVector(absent.domain().logical_size(), std::move(retained)));
    absent_endpoint.declarations.reset();
    const auto missing = NativeRecord::create_from_endpoint(results.front().unbound_value(), new_native_record_namespace(),
        results.front().raw_sources(), std::move(absent_endpoint));
    refused = false;
    try {(void)replay_scientific_record(missing);} catch (const std::invalid_argument &) {refused = true;}
    check(refused, "changed support replay accepted");
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    return 0;
} catch (const std::exception &error) {std::cerr << error.what() << '\n'; return 1;}
