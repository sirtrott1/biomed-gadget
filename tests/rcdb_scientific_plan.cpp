#include <gadget/native_record.hpp>
#include <gadget/scientific.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace {
using R = gadget::NativeRecord;
using N = gadget::NativeRecordValue;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
N &member(N &value, const std::string &name) {
    auto &members = std::get<N::Object>(value.data).members;
    const auto found = std::find_if(members.begin(), members.end(), [&](const auto &entry) {
        return entry.first == name;
    });
    require(found != members.end(), "missing native plan member");
    return found->second;
}
void source_keys(const N &value, std::set<std::uint64_t> &keys) {
    if (const auto *span = std::get_if<N::ByteSpan>(&value.data))
        keys.insert(span->source);
    else if (const auto *object = std::get_if<N::Object>(&value.data))
        for (const auto &[name, child] : object->members) {
            (void)name;
            source_keys(child, keys);
        }
    else if (const auto *array = std::get_if<N::Array>(&value.data))
        for (const auto &child : array->elements)
            source_keys(child, keys);
}
R plan(const R &first, const R &second) {
    const auto original = first.at("result").at("source_study_publication").bytes();
    require(std::ranges::equal(original, second.at("result").at("source_study_publication").bytes()),
            "fixture operations require the same actual scientific source");
    auto study = std::make_shared<const gadget::ScientificStudy>(R::decode(original));
    std::vector<gadget::RawByteSource> sources;
    N::Array steps;
    std::uint64_t key = 1;
    for (const auto &[id, operation] : std::vector<std::pair<std::string, R>>{{"first", first}, {"second", second}}) {
        const auto request = operation.at("request");
        std::set<std::uint64_t> keys;
        source_keys(request.unbound_value(), keys);
        std::map<std::uint64_t, std::uint64_t> remap;
        for (const auto &source : request.raw_sources()) {
            if (!keys.contains(source.key()))
                continue;
            remap.emplace(source.key(), key);
            sources.push_back(source.rekey(key));
            key = rfc::machine::checked_size_add(key, 1);
        }
        require(remap.size() == keys.size(), "query declaration has an unbound original byte source");
        steps.elements.emplace_back(N::Object{{{"id", N(id)}, {"query", request.unbound_value(remap)}}});
    }
    const auto request = R::create(N(N::Object{{{"schema", N(std::string("gadget.scientific.plan.v1"))},
        {"id", N(std::string("stored-plan"))}, {"question", N(std::string("Replay complete native fields"))},
        {"source_digest", N(study->digest())}, {"steps", N(std::move(steps))}}}),
        gadget::new_native_record_namespace(), std::move(sources));
    gadget::ScientificSession session(std::move(study));
    return session.execute_plan_record(request);
}
}

int main(int argc, char **argv) {
    try {
        require(argc == 4 || argc == 5, "expected two source operations, fresh output directory and optional foreign operation");
        const auto first = R::read(argv[1]), second = R::read(argv[2]);
        const std::filesystem::path root(argv[3]);
        require(!std::filesystem::exists(root), "plan fixture directory exists");
        std::filesystem::create_directories(root);
        const auto forward = plan(first, second), reverse = plan(second, first);
        forward.write(root/"forward.rfcxs");
        reverse.write(root/"reverse.rfcxs");
        for (const auto &change : {"missing-step", "changed-order", "changed-id", "changed-source"}) {
            auto value = forward.unbound_value();
            auto &result = member(value, "result");
            auto &steps = std::get<N::Array>(member(result, "steps").data).elements;
            if (std::string(change) == "missing-step")
                steps.pop_back();
            else if (std::string(change) == "changed-order")
                std::swap(steps.front(), steps.back());
            else if (std::string(change) == "changed-id")
                member(steps.front(), "id") = N(std::string("unbound-step"));
            else
                member(result, "source_digest") = N(std::string("unbound-source"));
            forward.derive(std::move(value)).write(root/(std::string(change)+".rfcxs"));
        }
        if (argc == 5) {
            const auto foreign = R::read(argv[4]);
            require(!std::ranges::equal(first.at("result").at("source_study_publication").bytes(),
                                       foreign.at("result").at("source_study_publication").bytes()),
                    "mixed-source control must retain a different actual study");
            std::uint64_t key = 1;
            for (const auto &source : forward.raw_sources())
                key = std::max(key, rfc::machine::checked_size_add(source.key(), 1));
            auto archive = foreign.archive_source(key);
            auto value = forward.unbound_value();
            auto &steps = std::get<N::Array>(member(member(value, "result"), "steps").data).elements;
            member(steps.front(), "publication") = N(N::ByteSpan{key, 0, archive.bytes().size(), false});
            forward.derive(std::move(value), {std::move(archive)}).write(root/"mixed-source.rfcxs");
        }
        std::cout << "{\"status\":\"PASS\",\"ordered_steps\":2,\"native_controls\":4,\"mixed_source\":"
                  << (argc == 5 ? "true" : "false") << "}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
