#include "native_scientific_fields.hpp"
#include <atomic>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2 && argc != 3)
            throw std::invalid_argument("expected Hodge case");
        const fs::path root = argv[1];
        if (argc == 3 && std::string(argv[2]) == "native") {native_scientific_test::run(root); return 0;}
        auto owner = std::make_shared<const ScientificStudy>(root / "study/study.json");
        ScientificSession session(owner);
        const auto query = Json::parse(read_text(root / "query.json"));
        auto stats = [&] { return session.native_hodge_preparation_statistics(); };
        session.validate(query);
        if (stats().at("hodge_frame_preparations").integer())
            throw std::logic_error("admission prepared native frames");
        const auto first = session.execute(query);
        if (stats().at("hodge_frame_preparations").integer() != 1 ||
            session.execute(query).dump() != first.dump())
            throw std::logic_error("frame preparation reuse changed certificate");
        auto serial = Json::parse(query.dump());
        serial.set("workers", Json::integer(1));
        if (session.execute(serial).at("result").at("fields").dump() !=
            first.at("result").at("fields").dump())
            throw std::logic_error("native CPU workers changed complete fields");
        auto changed = Json::parse(query.dump()), fields = Json::array();
        fields.push(Json::string("x2"));
        fields.push(Json::string("x"));
        fields.push(Json::string("x2"));
        changed.set("fields", fields);
        const auto ordered = session.execute(changed).at("result").at("fields");
        if (ordered.at(0).dump() != ordered.at(2).dump() ||
            ordered.at(1).at("sectors").dump() !=
                first.at("result").at("fields").at(0).at("sectors").dump() ||
            stats().at("hodge_frame_preparations").integer() != 1 || session.preparation_count())
            throw std::logic_error("ordered field change rebuilt native or derived Hodge owner");
        rfc::runtime::ExecutionHostPolicy policy;
        const auto n = std::numeric_limits<std::size_t>::max();
        policy.name = "current-Hodge-lease";
        policy.limits = {n, n, n, n, n};
        auto refuse = [&](const rfc::runtime::ExecutionLease &lease) {
            const auto before = stats().dump();
            bool refused = false;
            try {
                (void)session.execute_native_hodge(query, lease);
            } catch (const std::exception &) {
                refused = true;
            }
            if (!refused || before != stats().dump())
                throw std::logic_error("reused frame ignored current lease or mutated counters");
        };
        refuse(rfc::runtime::ExecutionLease(policy, std::make_shared<std::atomic_bool>(true)));
        policy.limits.records = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.records = n;
        policy.limits.scalar_work = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.scalar_work = n;
        policy.limits.rational_slots = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.rational_slots = n;
        policy.limits.integer_bits = 1;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.integer_bits = n;
        policy.backend = rfc::runtime::ExecutionBackend::ExactGpu;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.fallback = rfc::runtime::BackendFallback::AllowExactCpuFallback;
        const auto fallback =
            session.execute_native_hodge(query, rfc::runtime::ExecutionLease(policy));
        if (!fallback.at("result")
                 .at("fields")
                 .at(0)
                 .at("native_execution")
                 .at("fallback")
                 .boolean())
            throw std::logic_error("explicit exact CPU fallback omitted receipt");
        const auto bytes = read_text(root / "study/study.rfcxs");
        write_text(root / "study/study.rfcxs", "changed after admission");
        const auto retained = session.execute(query);
        write_text(root / "study/study.rfcxs", bytes);
        if (retained.dump() != first.dump())
            throw std::logic_error("Hodge source owner borrowed mutable files");
        auto detached = first;
        detached.at("result")
            .at("publication")
            .at("declaration")
            .set("id", Json::string("modified returned metadata"));
        if (session.execute(query).dump() != retained.dump())
            throw std::logic_error("Hodge result aliases immutable metadata");
        owner.reset();
        std::cout << stats().dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
