#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <gadget/resources.hpp>
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
}
int main(int argc, char **argv) try {
    if (argc != 5) throw std::invalid_argument("workspace, queries, output and mode required");
    const fs::path workspace_path(argv[1]), queries_path(argv[2]), output(argv[3]);
    fs::create_directories(output);
    const auto study = std::string(argv[4]) == "study"
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
    std::vector<NativeRecord> results;
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name);
        const auto result = session.execute_record(query);
        if (std::string(argv[4]) == "response") verify_response(result, original);
        else if (std::string(argv[4]) == "ensemble") verify_ensemble(result, original);
        else if (std::string(argv[4]) == "structural") verify_structural(result);
        else if (std::string(argv[4]) == "sections") verify_sections(result);
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
