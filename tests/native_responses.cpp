#include <atomic>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected response case");
        const fs::path root = argv[1];
        auto owner = std::make_shared<const ScientificStudy>(root / "study.json");
        ScientificSession session(owner);
        const auto response = Json::parse(read_text(root / "response-query.json")),
                   observation = Json::parse(read_text(root / "observed-query.json"));
        session.validate(observation);
        if (session.native_response_preparation_statistics()
                .at("green_pair_preparations")
                .integer())
            throw std::logic_error("admission prepared a Green solve");
        const auto first = session.execute(response), observed = session.execute(observation);
        const auto counts = [&] { return session.native_response_preparation_statistics(); };
        if (counts().at("green_pair_preparations").integer() != 1 ||
            counts().at("observation_preparations").integer() != 1)
            throw std::logic_error("observation failed to reuse its actual Green owner");
        if (session.execute(response).dump() != first.dump() ||
            session.execute(observation).dump() != observed.dump())
            throw std::logic_error("certificate changed on preparation reuse");
        auto equivalent = Json::parse(observation.dump());
        equivalent.set("lambda", Json::string("1/1"));
        (void)session.execute(equivalent);
        auto changed = Json::parse(observation.dump());
        changed.at("old").set("field", Json::string("u2"));
        changed.at("new").set("field", Json::string("u2"));
        (void)session.execute(changed);
        if (counts().at("green_pair_preparations").integer() != 1 ||
            counts().at("observation_preparations").integer() != 1 ||
            session.temporal_preparation_statistics().at("green_preparations").integer())
            throw std::logic_error("field change duplicated native or derived preparations");
        rfc::runtime::ExecutionHostPolicy low;
        const auto n = std::numeric_limits<std::size_t>::max();
        low.name = "explicit-current-lease";
        low.limits = {n, n, n, n, n};
        auto refusal = [&](const Json &q, const rfc::runtime::ExecutionLease &lease) {
            const auto before = counts().dump();
            bool refused = false;
            try {
                (void)session.execute_native_response(q, lease);
            } catch (const std::exception &) {
                refused = true;
            }
            if (!refused || counts().dump() != before)
                throw std::logic_error("current lease not rechecked on reuse");
        };
        auto cancellation = std::make_shared<std::atomic_bool>(true);
        refusal(observation, rfc::runtime::ExecutionLease(low, cancellation));
        low.limits.records = 0;
        refusal(observation, rfc::runtime::ExecutionLease(low));
        low.limits.records = n;
        low.limits.scalar_work = 0;
        refusal(observation, rfc::runtime::ExecutionLease(low));
        low.limits.scalar_work = n;
        low.limits.integer_bits = 1;
        refusal(observation, rfc::runtime::ExecutionLease(low));
        low.limits.integer_bits = n;
        const auto custom =
            session.execute_native_response(observation, rfc::runtime::ExecutionLease(low));
        if (custom.at("result").at("native_execution").at("execution_policy").text() != low.name ||
            custom.at("result").at("native_execution").at("green_preparation_policy").text() ==
                low.name)
            throw std::logic_error("preparation and execution policies conflated");
        // Immutable admitted sources and detached certificates survive source-file mutation.
        const auto original = read_text(root / "p0/study.rfcxs");
        write_text(root / "p0/study.rfcxs", "changed after admission");
        const auto retained = session.execute(observation);
        write_text(root / "p0/study.rfcxs", original);
        if (retained.dump() != observed.dump())
            throw std::logic_error("source owner borrowed mutable files");
        const auto original_certificate = first.dump();
        auto detached = first;
        detached.at("result")
            .at("publications")
            .at("old")
            .at("declaration")
            .set("id", Json::string("modified returned source metadata"));
        if (session.execute(response).dump() != original_certificate)
            throw std::logic_error("returned JSON aliased source or result cache");
        owner.reset();
        std::cout << counts().dump() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
