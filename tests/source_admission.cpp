#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <algorithm>
#include <iostream>
#include <thread>

namespace {
void json_admission(const gadget::fs::path &corpus_path) {
    using namespace gadget;
    const auto check = [](bool ok, const char *message) {
        if (!ok) throw std::logic_error(message);
    };
    const auto rejects = [&](auto action) {
        bool refused = false;
        try { action(); } catch (const std::invalid_argument &) { refused = true; }
        check(refused, "unrepresentable JSON operation admitted");
    };
    auto object = Json::parse(R"({"schema":"original","":"empty name"})");
    const std::string name("schema\0suffix", 13);
    rejects([&] { object.set(name, Json::string("replacement")); });
    rejects([&] { (void)object.has(name); });
    rejects([&] { (void)object.at(name); });
    check(object.at("schema").text() == "original" && object.at("").text() == "empty name",
          "invalid member access changed its C-string prefix or empty name");
    for (const auto *source : {R"("\ud800")", R"("\udc00")",
                              R"({"\ud800":1,"\ufffd":2})"}) {
        rejects([&] { (void)Json::parse(source); });
        const auto bytes = std::string(source);
        rejects([&] { (void)decode_record(std::as_bytes(std::span(bytes.data(), bytes.size()))); });
    }
    const auto distinct = Json::parse(R"({"\ud800\udc00":1,"\ufffd":2})");
    check(json_object_object_length(distinct.raw()) == 2 && distinct.at("\xf0\x90\x80\x80").integer() == 1 &&
              distinct.at("\xef\xbf\xbd").integer() == 2,
          "valid supplementary identity and literal replacement character collapsed");
    for (const auto *token : {"9223372036854775808", "18446744073709551615"}) {
        const auto value = Json::parse(token);
        bool refused = false;
        try { (void)value.integer(); } catch (const std::length_error &) { refused = true; }
        check(refused && value.dump_compact() == token, "signed access saturated or altered original unsigned value");
        const auto bytes = std::string(token);
        const auto native = decode_record(std::as_bytes(std::span(bytes.data(), bytes.size()))).native();
        check(native.integer().get_str() == token, "native unsigned integer value changed");
    }
    check(Json::parse("9223372036854775807").integer() == std::numeric_limits<std::int64_t>::max() &&
              Json::parse("-9223372036854775808").integer() == std::numeric_limits<std::int64_t>::min(),
          "exact signed integer boundaries refused");
    for (const auto *token : {"1e400", "-1e400", "1e-400", "-0e400"}) {
        const auto bytes = std::string(token);
        const auto native = decode_record(std::as_bytes(std::span(bytes.data(), bytes.size()))).native();
        check(Json::parse(bytes).dump_compact() == bytes && native.rational() == exact_decimal(bytes),
              "valid finite exact decimal rejected by its binary64 view");
    }
    const auto corpus = Json::parse(read_text(corpus_path));
    for (std::size_t i = 0; i < corpus.size(); ++i) {
        const auto item = corpus.at(i);
        const auto bytes = item.at("source").text();
        if (!item.at("admit").boolean()) {
            rejects([&] { (void)Json::parse(bytes); });
            rejects([&] { (void)decode_record(std::as_bytes(std::span(bytes.data(), bytes.size()))); });
            continue;
        }
        const auto expected = item.at("value").dump_compact();
        check(Json::parse(bytes).dump_compact() == expected,
              "complete JSON differs from independent decoded control oracle");
        const auto record = decode_record(std::as_bytes(std::span(bytes.data(), bytes.size())));
        check(record.dump_compact() == expected,
              "native admitted value differs from independent decoded control oracle");
        const auto originals = record.native().raw_sources();
        const auto original = std::find_if(originals.begin(), originals.end(),
                                           [](const auto &source) { return source.key() == 1; });
        check(original != originals.end() && original->bytes().size() == bytes.size() &&
                  std::ranges::equal(original->bytes(),
                                     std::as_bytes(std::span(bytes.data(), bytes.size()))),
              "JSON admission changed complete original source bytes");
    }
    using N = NativeRecordValue;
    const auto native = NativeRecord::create(N(N::Object{{{name, N(std::string("value\0bytes", 11))}}}),
                                              971031);
    check(native.names() == std::vector<std::string>{name} &&
              native.at(name).text() == std::string("value\0bytes", 11),
          "native names or values lost original NUL bytes");
    const auto destination = corpus_path.parent_path() / "native-nul-name.rfcxs";
    native.write(destination);
    const auto reopened = NativeRecord::read(destination);
    check(reopened.names() == native.names() && reopened.at(name).text() == native.at(name).text(),
          "native NUL name or value changed after archive reopen");
    rejects([&] { (void)render_record_view(reopened); });
    std::cout << "PASS: " << corpus.size()
              << " independent JSON corpus controls, complete original byte owners and native names\n";
}
} // namespace

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2 && argc != 3)
            throw std::invalid_argument("expected version-two study and optional JSON corpus");
        if (argc == 3)
            json_admission(argv[2]);
        const fs::path path = argv[1];
        auto owner = std::make_shared<const ScientificStudy>(path);
        const auto original = read_text(path.parent_path() / "table.txt");
        if (sha256_file(path.parent_path() / "table.txt") != sha256(original))
            throw std::logic_error("incremental identity differs from byte identity");
        const auto index = owner->source_index();
        auto request = Json::object();
        request.set("schema", Json::string("gadget.scientific.source.read.v1"));
        request.set("source_digest", index.at("source_digest"));
        request.set("publication", Json::string("self"));
        request.set("source_sequence", index.at("source_sequence"));
        request.set("stream", Json::parse(R"({"kind":"document","id":"table"})"));
        request.set("offset", Json::integer(0));
        request.set("limit", Json::integer(original.size()));
        const auto expected = owner->read_source(request).dump();
        write_text(path.parent_path() / "table.txt", "replacement after admission\n");
        std::exception_ptr failure[2];
        std::thread first([&] {
            try {
                if (owner->read_source(request).dump() != expected)
                    throw std::logic_error("admitted source followed external file mutation");
            } catch (...) {
                failure[0] = std::current_exception();
            }
        });
        std::thread second([&] {
            try {
                ScientificSession session(owner);
                if (session.read_source(request).dump() != expected)
                    throw std::logic_error("concurrent immutable source changed");
            } catch (...) {
                failure[1] = std::current_exception();
            }
        });
        first.join();
        second.join();
        write_text(path.parent_path() / "table.txt", original);
        for (const auto &error : failure)
            if (error)
                std::rethrow_exception(error);
        auto detached = owner->declaration();
        detached.at("cells").at(0).set("id", Json::string("changed"));
        if (owner->declaration().at("cells").at(0).at("id").text() == "changed")
            throw std::logic_error("detached declaration changed owner");
        std::cout << "PASS: streaming file identity, immutable admitted bytes, concurrent readers "
                     "and detached metadata\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
