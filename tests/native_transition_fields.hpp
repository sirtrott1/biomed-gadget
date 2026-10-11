#pragma once
#include "native_scientific_fields.hpp"

namespace native_transition_test {
using namespace native_scientific_test;
inline NativeRecord previous(const NativeRecord &fields) {
    return NativeRecord::decode(fields.at("bindings").at(0).at("old_publication").at("publication").bytes());
}
inline Vector values(const NativeRecord &owner, const NativeRecord &descriptor) {
    return vector(*tensor(owner, descriptor));
}
inline void source(const NativeRecord &fields, const NativeRecord &original) {
    check(rfc::storage::encode_rfcx(fields.endpoint()->source).bytes ==
          rfc::storage::encode_rfcx(original.endpoint()->source).bytes, "original complex changed");
    for (std::size_t i = 0; i < original.endpoint()->tensors.size(); ++i)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(
            original.endpoint()->tensors[i], fields.endpoint()->tensors.at(i)), "original tensor changed");
    for (const auto &bytes : original.raw_sources()) {
        const auto retained = fields.raw_sources();
        const auto found = std::find_if(retained.begin(), retained.end(), [&](const auto &item) {return item.key() == bytes.key();});
        check(found != retained.end() && std::ranges::equal(found->bytes(), bytes.bytes()), "original byte owner changed");
    }
}
inline Rational pairing(const Vector &a, const Tensor &metric, const Vector &b) {
    const auto measured = apply_map(metric, b);
    check(a.size() == measured.size(), "pairing axis mismatch");
    Rational out;
    for (std::size_t i = 0; i < a.size(); ++i) out += a[i] * measured[i];
    return out;
}
inline void novelty(const NativeRecord &publication) {
    const auto fields = publication.at("result").at("fields"), old = previous(fields);
    const auto view = fields.at("view"), split = view.at("split"), binding = fields.at("bindings").at(0);
    const auto input = values(old, old.at("view").at("input"));
    const auto next = values(fields, binding.at("new_input"));
    const auto &j = *tensor(fields, binding.at("comparison")), &m = *tensor(fields, binding.at("metric"));
    const auto transported = values(fields, split.at("transported_old"));
    check(apply_map(j, input) == transported, "transported original field differs");
    const auto innovation = values(fields, split.at("innovation")), persistent = values(fields, split.at("persistent"));
    const auto complement = values(fields, split.at("complement"));
    const auto coordinates = values(old, old.at("view").at("persistent_coordinates"));
    check(tensor(old, old.at("view").at("persistent_coordinates"))->output ==
          tensor(old, old.at("view").at("input"))->output, "persistent coordinates lost actual old carrier");
    check(apply_map(j, coordinates) == persistent, "persistent image differs");
    const auto reconstructed = values(fields, split.at("reconstructed"));
    for (std::size_t i = 0; i < next.size(); ++i) {
        check(innovation[i] == next[i] - transported[i], "innovation differs");
        check(reconstructed[i] == persistent[i] + complement[i] && reconstructed[i] == innovation[i], "split reconstruction differs");
    }
    const auto &g = *tensor(fields, split.at("pairings"));
    const auto &pulled = *tensor(old, old.at("view").at("pulled_back_metric"));
    const auto measured = apply_map(m, complement);
    for (std::size_t i = 0; i < input.size(); ++i) {
        Rational zero;
        for (std::size_t k = 0; k < next.size(); ++k) zero += matrix(j, k, i) * measured[k];
        check(zero == 0, "complement not orthogonal to entire old image");
        for (std::size_t c = 0; c < input.size(); ++c) {
            Rational sum;
            for (const auto &a : j.tensor.values().entries()) if (a.index % input.size() == i)
                for (const auto &b : j.tensor.values().entries()) if (b.index % input.size() == c)
                    sum += a.value * matrix(m, a.index / input.size(), b.index / input.size()) * b.value;
            check(matrix(pulled, i, c) == sum, "pulled back metric differs");
        }
    }
    const std::vector<Vector> sectors{persistent, complement};
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 2; ++c)
            check(matrix(g, i, c) == pairing(sectors[i], m, sectors[c]), "split pairing differs");
    check(matrix(g, 0, 1) == 0, "split cross pairing differs");
    check(*tensor(fields, split.at("quadrance"))->tensor.state(0) == pairing(innovation, m, innovation), "innovation quadrance differs");
    if (view.has("causes")) {
        const auto causes = view.at("causes"), word = causes.at("word"), dictionary = causes.at("fields");
        const auto slots = values(fields, binding.at("cause_word"));
        std::vector<Vector> channels;
        Vector aggregate(next.size());
        for (std::size_t i = 0; i < word.size(); ++i) {
            const auto index = word.at(i).at("dictionary_index").integer().get_ui();
            check(slots[i] == Rational(std::to_string(index)), "ordered cause word differs");
            const auto descriptor = dictionary.at(index).at("values");
            check(descriptor.at("tensor").integer() == rfc::numeric::Integer(dictionary.at(index).at("field_key").text()), "cause field was copied");
            channels.push_back(values(fields, descriptor));
            for (std::size_t k = 0; k < next.size(); ++k) aggregate[k] += channels.back()[k];
        }
        check(aggregate == innovation && aggregate == values(fields, causes.at("aggregate")), "cause aggregate differs");
        const auto &gram = *tensor(fields, causes.at("pairings"));
        std::vector<Vector> measured_channels;
        for (const auto &channel : channels) measured_channels.push_back(apply_map(m, channel));
        for (std::size_t i = 0; i < channels.size(); ++i)
            for (std::size_t c = 0; c < channels.size(); ++c) {
                Rational sum;
                for (std::size_t k = 0; k < next.size(); ++k) sum += channels[i][k] * measured_channels[c][k];
                check(matrix(gram, i, c) == sum, "ordered cause cross pairing differs");
            }
    }
}
inline void response(const NativeRecord &publication) {
    const auto fields = publication.at("result").at("fields"), old = previous(fields);
    const auto view = fields.at("view"), old_view = old.at("view"), binding = fields.at("bindings").at(0);
    const auto aggregate = values(fields, view.at("aggregate"));
    Vector sum(aggregate.size());
    std::vector<Vector> channels;
    for (std::size_t i = 0; i < view.at("cause_order").size(); ++i) {
        const auto name = view.at("cause_order").at(i).text();
        channels.push_back(values(fields, view.at("causes").at(name)));
        for (std::size_t k = 0; k < sum.size(); ++k) sum[k] += channels.back()[k];
    }
    check(sum == aggregate, "ordered response reconstruction differs");
    const auto &j = *tensor(fields, binding.at("ambient_comparison"));
    const auto &old_l = *tensor(old, old_view.at("declared_maps").at("old_resolving"));
    const auto &new_l = *tensor(fields, view.at("declared_maps").at("new_resolving"));
    if (old_view.has("old_response")) {
        const auto old_g = values(old, old_view.at("old_response")), new_g = values(fields, view.at("new_response"));
        const auto old_lg = apply_map(old_l, old_g), new_lg = apply_map(new_l, new_g);
        const auto old_injected = apply_map(*tensor(old, old_view.at("declared_maps").at("old_realization")), values(old, old_view.at("input_values")));
        const auto new_injected = apply_map(*tensor(fields, view.at("declared_maps").at("new_realization")), values(fields, view.at("publications").at("new").at("input_values")));
        const auto lambda = Rational(publication.at("request").at("lambda").text());
        for (std::size_t i = 0; i < old_g.size(); ++i) check(old_g[i] + lambda * old_lg[i] == old_injected[i], "old Green solve differs");
        for (std::size_t i = 0; i < new_g.size(); ++i) check(new_g[i] + lambda * new_lg[i] == new_injected[i], "new Green solve differs");
        const auto transported = apply_map(j, old_g);
        for (std::size_t i = 0; i < aggregate.size(); ++i) check(aggregate[i] == new_g[i] - transported[i], "response delta differs");
        check(tensor(old, old_view.at("old_response"))->output == old_l.output, "old response lost actual old carrier");
    }
    if (view.has("pairing")) {
        const auto &metric = *tensor(fields, view.at("declared_maps").at("metric"));
        const auto &gram = *tensor(fields, view.at("pairing").at("gram"));
        for (std::size_t i = 0; i < channels.size(); ++i)
            for (std::size_t c = 0; c < channels.size(); ++c)
                check(matrix(gram, i, c) == pairing(channels[i], metric, channels[c]), "response cross pairing differs");
        check(*tensor(fields, view.at("pairing").at("aggregate_quadrance"))->tensor.state(0) == pairing(aggregate, metric, aggregate), "response quadrance differs");
    }
}
inline NativeRecord changed_field(const NativeRecord &record, std::uint64_t key, bool remove_support) {
    auto endpoint = *record.endpoint();
    auto found = std::find_if(endpoint.tensors.begin(), endpoint.tensors.end(), [&](const auto &item) {return item.key == key;});
    auto &changed = found->tensor;
    std::vector<rfc::numeric::SparseRationalEntry> entries;
    for (const auto &entry : changed.values().entries()) if (entry.index != 0) entries.push_back(entry);
    if (!remove_support) entries.push_back({0, *changed.state(0) + 1});
    std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) {return a.index < b.index;});
    std::vector<std::size_t> support;
    changed.support().for_each_active([&](auto i) {if (!remove_support || i) support.push_back(i);});
    changed = rfc::numeric::SparseRationalTensor(changed.domain(), rfc::numeric::StructuralSupport::sparse(changed.domain().logical_size(), support),
        rfc::numeric::SparseRationalVector(changed.domain().logical_size(), entries));
    endpoint.declarations.reset();
    return NativeRecord::create_from_endpoint(record.unbound_value(), new_native_record_namespace(), record.raw_sources(), std::move(endpoint));
}
inline void resize_span(NativeRecordValue &value, std::uint64_t source_key, std::size_t size) {
    if (auto *span = std::get_if<NativeRecordValue::ByteSpan>(&value.data)) {
        if (span->source == source_key) span->length = size;
    } else if (auto *object = std::get_if<NativeRecordValue::Object>(&value.data)) {
        for (auto &[name, item] : object->members) {(void)name; resize_span(item, source_key, size);}
    } else if (auto *array = std::get_if<NativeRecordValue::Array>(&value.data)) {
        for (auto &item : array->elements) resize_span(item, source_key, size);
    }
}
inline void run(const fs::path &root, bool is_novelty) {
    checks = 0;
    const auto study = std::make_shared<const ScientificStudy>(root / "study.json");
    ScientificSession session(study);
    const std::vector<std::string> requests = is_novelty ? std::vector<std::string>{"query"} : std::vector<std::string>{"response-query", "observed-query"};
    for (const auto &name : requests) {
        const auto query = admit_record_view(Json::parse(read_text(root / (name + ".json"))));
        const auto result = session.execute_record(query), fields = result.at("result").at("fields"), old = previous(fields);
        const ScientificStudy old_source(root / "p0/study.json"), new_source(root / "p1/study.json");
        source(old, old_source.native_publication()); source(fields, new_source.native_publication());
        if (is_novelty) novelty(result); else response(result);
        const auto folder = root / "native" / name; fs::create_directories(folder);
        result.write(folder / "result.rfcxs");
        fs::rename(root / "p0", root / "p0.offline"); fs::rename(root / "p1", root / "p1.offline");
        try {
            const auto replayed = replay_scientific_record(NativeRecord::read(folder / "result.rfcxs"));
            compare(result, replayed); compare(replayed, replay_scientific_record(replayed));
            source(previous(replayed.at("result").at("fields")), old_source.native_publication());
            if (is_novelty) novelty(replayed); else response(replayed);
        } catch (...) {fs::rename(root / "p0.offline", root / "p0"); fs::rename(root / "p1.offline", root / "p1"); throw;}
        fs::rename(root / "p0.offline", root / "p0"); fs::rename(root / "p1.offline", root / "p1");
        const auto descriptor = is_novelty ? fields.at("view").at("split").at("innovation") : fields.at("view").at("aggregate");
        const auto target = tensor(fields, descriptor);
        if (target->tensor.domain().logical_size()) {
            for (const bool remove_support : {false, true}) {
                const auto forged = changed_field(result, target->key, remove_support);
                if (!remove_support) forged.write(folder / "changed-field.rfcxs");
                bool refused = false;
                try {(void)replay_scientific_record(forged);} catch (const std::invalid_argument &) {refused = true;}
                check(refused, "replay accepted changed field or support");
            }
        }
        const auto old_input = old.at("view").at(is_novelty ? "input" : "input_values");
        const auto old_tensor = tensor(old, old_input);
        if (old_tensor->tensor.domain().logical_size()) {
            const auto span = fields.at("bindings").at(0).at("old_publication").at("publication");
            const auto key = span.byte_source().key();
            auto bytes = changed_field(old, old_tensor->key, false).archive_source(key);
            auto raw = result.raw_sources();
            for (auto &owner : raw) if (owner.key() == key) owner = bytes;
            auto value = result.unbound_value(); resize_span(value, key, bytes.bytes().size());
            auto endpoint = *result.endpoint(); endpoint.declarations.reset();
            const auto forged = NativeRecord::create_from_endpoint(std::move(value), new_native_record_namespace(), std::move(raw), std::move(endpoint));
            bool refused = false;
            try {(void)replay_scientific_record(forged);} catch (const std::invalid_argument &) {refused = true;}
            check(refused, "replay accepted a changed actual old field");
        }
        auto changed = result.unbound_value();
        auto &binding = std::get<NativeRecordValue::Array>(member(member(member(changed, "result"), "fields"), "bindings").data).elements.at(0);
        member(binding, "source_sequence") = NativeRecordValue(std::string("detached cutoff"));
        bool refused = false;
        try {(void)replay_scientific_record(result.derive(std::move(changed)));} catch (const std::invalid_argument &) {refused = true;}
        check(refused, "replay accepted detached new cutoff");
        const auto original_bytes = new_source.native_publication().raw_sources();
        if (!original_bytes.empty() && !original_bytes.front().bytes().empty()) {
            auto raw = result.raw_sources();
            for (auto &owner : raw) if (owner.key() == original_bytes.front().key()) {
                std::vector<std::byte> bytes(owner.bytes().begin(), owner.bytes().end());
                bytes.front() ^= std::byte{1}; owner = RawByteSource::own(owner.key(), std::move(bytes));
            }
            auto endpoint = *result.endpoint(); endpoint.declarations.reset();
            refused = false;
            try {
                const auto forged = NativeRecord::create_from_endpoint(result.unbound_value(), new_native_record_namespace(), std::move(raw), std::move(endpoint));
                (void)replay_scientific_record(forged);
            } catch (const std::invalid_argument &) {refused = true;}
            check(refused, "replay accepted changed original bytes");
        }
    }
    auto plan = Json::object(), steps = Json::array();
    plan.set("schema", Json::string("gadget.scientific.plan.v1"));
    plan.set("id", Json::string("native-transition")); plan.set("question", Json::string("Retain both endpoints"));
    plan.set("source_digest", Json::string(study->digest()));
    for (std::size_t i = 0; i < 2; ++i) {
        auto step = Json::object(); step.set("id", Json::string(std::to_string(i)));
        step.set("query", Json::parse(read_text(root / (requests[i % requests.size()] + ".json")))); steps.push(step);
    }
    plan.set("steps", steps);
    const auto planned = session.execute_plan_record(admit_record_view(plan));
    planned.write(root / "native/plan.rfcxs");
    const auto replayed = replay_scientific_record(NativeRecord::read(root / "native/plan.rfcxs"));
    check(replayed.at("result").at("steps").size() == 2, "plan replay lost ordered steps");
    check(planned.raw_sources().size() == replayed.raw_sources().size(), "plan replay accumulated owners");
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
}
} // namespace native_transition_test
