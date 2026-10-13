#include <gadget/native_record.hpp>
#include <gadget/study.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <iostream>
#include <algorithm>
#include <barrier>
#include <exception>
#include <thread>

using namespace gadget;
using namespace rfc;
namespace {
void check(bool ok,const char *message) {if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F &&f) {bool threw=false;try{f();}catch(const std::exception &){threw=true;}check(threw,"invalid native record admitted");}
std::vector<storage::PinnedNativeAccession> concurrent_access(NativeRecord record) {
    constexpr std::size_t workers = 8, rounds = 64;
    const auto endpoint = record.endpoint();
    const auto name_space = record.source_namespace();
    std::barrier ready(static_cast<std::ptrdiff_t>(workers));
    std::vector<std::optional<storage::PinnedNativeAccession>> selected(
        workers * endpoint->tensors.size());
    std::vector<std::exception_ptr> errors(workers);
    std::vector<std::jthread> threads;
    for (std::size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back([&, worker] {
            ready.arrive_and_wait();
            try {
                for (std::size_t round = 0; round < rounds; ++round) {
                    const auto view = round % 2 ? record.at("cell") : record;
                    const auto current = view.historical_publication(endpoint->sequence);
                    check(current.accession_directory() == record.accession_directory() &&
                              current.endpoint() == endpoint && current.is_root(),
                          "current publication prepared a second directory or lost its owner");
                    storage::NativeAccessionRequest request;
                    request.kind = storage::NativeAccessionKind::Tensor;
                    request.subject = 999;
                    rejects([&] { (void)view.access(request); });
                    for (std::size_t i = 0; i < endpoint->tensors.size(); ++i) {
                        request.subject = endpoint->tensors[i].key;
                        auto pin = view.access(request);
                        check(pin.endpoint() == endpoint && pin.source_namespace() == name_space &&
                                  pin.sequence() == endpoint->sequence &&
                                  pin.tensor().get() == &endpoint->tensors[i].tensor,
                              "concurrent record access changed its actual tensor owner");
                        const auto tensor = pin.tensor();
                        for (std::size_t coordinate = 0; coordinate < tensor->domain().logical_size(); ++coordinate)
                            check(tensor->state(coordinate) == endpoint->tensors[i].tensor.state(coordinate),
                                  "concurrent record access changed a coefficient or support state");
                        selected[worker * endpoint->tensors.size() + i] = std::move(pin);
                    }
                }
            } catch (...) {
                errors[worker] = std::current_exception();
            }
        });
    threads.clear();
    for (const auto &error : errors)
        if (error) std::rethrow_exception(error);
    std::vector<storage::PinnedNativeAccession> pins;
    for (auto &pin : selected) {
        check(pin.has_value(), "concurrent record access omitted an original tensor");
        if (!pins.empty())
            check(pin->same_execution_authority(pins.front()),
                  "record copies or children prepared independent directories");
        pins.push_back(std::move(*pin));
    }
    return pins;
}
void equal_value(const NativeRecord &a, const NativeRecord &b) {
    check(a.kind() == b.kind(), "native value constructor changed");
    switch (a.kind()) {
    case NativeRecordKind::Null: break;
    case NativeRecordKind::Boolean: check(a.boolean() == b.boolean(), "Boolean changed"); break;
    case NativeRecordKind::Integer: check(a.integer() == b.integer(), "Integer changed"); break;
    case NativeRecordKind::Rational: check(a.rational() == b.rational(), "Rational changed"); break;
    case NativeRecordKind::ExactNumber:
        check(a.rational() == b.rational() && a.approximate_view() == b.approximate_view(), "exact number changed"); break;
    case NativeRecordKind::Decimal:
        check(a.rational() == b.rational() && a.text() == b.text(), "decimal changed"); break;
    case NativeRecordKind::Text:
    case NativeRecordKind::Bytes:
        check(std::ranges::equal(a.bytes(), b.bytes()), "native byte span changed"); break;
    case NativeRecordKind::Object:
        check(a.names() == b.names(), "native factor names or order changed");
        for (const auto &name : a.names()) equal_value(a.at(name), b.at(name));
        break;
    case NativeRecordKind::Array:
        check(a.size() == b.size(), "native cardinality changed");
        for (std::size_t i = 0; i < a.size(); ++i) equal_value(a.at(i), b.at(i));
        break;
    case NativeRecordKind::Tensor:
        check(a.tensor().domain() == b.tensor().domain(), "native tensor domain changed");
        for (std::size_t i = 0; i < a.tensor().domain().logical_size(); ++i)
            check(a.tensor().state(i) == b.tensor().state(i), "native tensor/support changed");
        break;
    }
}
void named_views() {
    using N = NativeRecordValue;
    N::Object input;
    std::vector<std::string> order;
    for (std::size_t i = 0; i < 256; ++i) {
        const auto name = "original_declared_parameter_" + std::to_string(255 - i);
        order.push_back(name);
        input.members.emplace_back(name, N(numeric::Integer(i)));
    }
    const std::string embedded("source\0parameter", 16);
    order.push_back(embedded);
    input.members.emplace_back(embedded, N(numeric::Integer(256)));
    const auto source = NativeRecord::create(N(std::move(input)), 847);
    check(source.names() == order, "native name selection reordered authored factors");
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto surrounded = "[" + order[i] + "]";
        const std::string_view selection(surrounded.data() + 1, order[i].size());
        check(source.has(selection) && source.at(selection).integer() == numeric::Integer(i),
              "native name selection ignored exact string view bounds");
        check(!source.at(selection).has(selection), "scalar claimed a product factor");
    }
    check(!source.has("source") && !source.has("original_declared_parameter_256"),
          "native name lookup accepted a prefix or absent factor");
    rejects([&] { (void)source.at("source"); });
    for (const bool bound : {false, true}) {
        const auto value = bound ? source.materialize_value() : source.unbound_value();
        const auto &members = std::get<N::Object>(value.data).members;
        check(members.size() == order.size(), "native object copy omitted authored factors");
        for (std::size_t i = 0; i < order.size(); ++i)
            check(members[i].first == order[i], "native object copy reordered original factors");
        equal_value(source, NativeRecord::create(value, bound ? 849 : 848));
    }
    const auto selected = source.at(order[127]);
    const auto retained = [&] {
        const auto decoded = NativeRecord::decode(source.encode());
        return decoded.at(order[127]);
    }();
    equal_value(selected, retained);
    check(storage::encode_rfcx(retained.endpoint()->source).bytes ==
              storage::encode_rfcx(source.endpoint()->source).bytes,
          "retained native member lost its immutable source owner");
}
}
int main(int argc,char **argv) {
 try {
    const auto directory=std::filesystem::path(argc>1?argv[1]:"native-record-contracts");
    std::filesystem::remove_all(directory);std::filesystem::create_directories(directory);
    named_views();
    core::RelationSignature signature;const auto relation=signature.add_type("authored-group");
    core::RelationalComplexBuilder builder(signature);builder.declare_grade(core::Grade{3});
    builder.add_grade_one(core::RelationKey{77},relation,{
        {11,true,false,numeric::Rational{0},{},core::OccurrenceKey{99}},
        {11,false,true,numeric::Rational{0},{},core::OccurrenceKey{100}}});
    storage::TemporalEndpoint endpoint(builder.freeze(),4);
    storage::TemporalTensor field;field.key=888;field.role=storage::TemporalTensorRole::OpaqueExact;
    field.interpretation="coefficient selection fixture";
    field.tensor=numeric::SparseRationalTensor(numeric::IndexDomain::linear(3),
        numeric::StructuralSupport::sparse(3,{0,2}),numeric::SparseRationalVector(3,{{2,numeric::Rational{-5,13}}}));
    endpoint.tensors.push_back(field);
    endpoint.carriers.push_back({881, storage::TemporalCarrierKind::SourceGrade, {0}, {}, 0,
                                  numeric::IndexDomain::linear(1)});
    endpoint.carriers.push_back({882, storage::TemporalCarrierKind::SourceGrade, {1}, {}, 0,
                                  numeric::IndexDomain::linear(1)});
    storage::TemporalTensor measured;
    measured.key = 891;
    measured.role = storage::TemporalTensorRole::Field;
    measured.output = 881;
    measured.tensor = numeric::SparseRationalTensor(
        numeric::IndexDomain::linear(1), numeric::StructuralSupport::full(1),
        numeric::SparseRationalVector(1, {{0, numeric::Rational{-7, 5}}}));
    endpoint.tensors.push_back(measured);
    measured.key = 892;
    measured.output = 882;
    measured.tensor = numeric::SparseRationalTensor::zero(
        numeric::IndexDomain::linear(1), numeric::StructuralSupport::full(1));
    endpoint.tensors.push_back(measured);
    const auto original=storage::encode_rfcx(endpoint.source).bytes;
    auto collided=endpoint;
    storage::TemporalDeclarationState retained_declarations;
    retained_declarations.namespace_key=810;
    const composite::DeclarationId retained_type{810,1,1};
    retained_declarations.registry=composite::DeclarationRegistry({composite::Declaration::integer(retained_type)});
    retained_declarations.values.push_back({0x4741444745540203ULL,retained_type,0,0,{},{},composite::DeclaredValue::integer(numeric::Integer{19})});
    collided.declarations=retained_declarations;
    storage::validate_temporal_endpoint(collided);
    rejects([&]{NativeRecord::create(NativeRecordValue{},810,{},&collided);});

    std::vector<std::byte> raw(200003,std::byte{0});
    for(std::size_t i=0;i<200000;++i)raw[i]=static_cast<std::byte>(i%256);
    using N=NativeRecordValue;
    N payload(N::Object{{{"name",N(std::string("tumor"))},
        {"evidence",N(N::ByteSpan{42,0,raw.size(),false})},
        {"value",N(numeric::Rational{"2/7"})}}});
    N cell(N::Owned{0,storage::StableCellKey{core::Grade{1},77},{},std::make_shared<const N>(std::move(payload))});
    N occurrence(N::Owned{0,{},storage::StableOccurrenceKey{core::Grade{1},core::OccurrenceKey{100}},std::make_shared<const N>(std::string("inactive occurrence retained"))});
    const auto record=NativeRecord::create(N(N::Object{{{"cell",std::move(cell)},{"occurrence",std::move(occurrence)}}}),811,{RawByteSource(42,raw)},&endpoint);
    check(record.publication_history().has_value() &&
          record.publication_history()->segment_count()==0 &&
          record.publication_history()->tail().pin()==record.endpoint(),
          "record directory did not retain its actual immutable publication");
    auto temporary_endpoint=endpoint;
    const auto consumed=NativeRecord::create_from_endpoint(N(numeric::Rational{"5/13"}),819,
                                                          {},std::move(temporary_endpoint));
    check(storage::encode_rfcx(consumed.endpoint()->source).bytes==original &&
          consumed.rational()==numeric::Rational{"5/13"} &&
          consumed.endpoint()==consumed.publication_history()->tail().pin(),
          "consumed endpoint changed its source or duplicated the publication");
    storage::NativeTemporalSeries origin(endpoint,811,native_data_limits());
    const auto origin_bytes=storage::encode_native_series(origin);
    const auto annotated=NativeRecord::annotate(N(N::Object{{{"raw",N(N::ByteSpan{42,0,raw.size(),false})},
        {"observation",N(numeric::Rational{"2/7"})}}}),origin,{RawByteSource(42,raw)});
    const auto annotation_history=*annotated.publication_history();
    check(annotation_history.segment_count()==1 && origin.segment_count()==0,
          "annotation mutated its original history");
    check(storage::encode_native_series(storage::NativeSeriesPrefix::capture(annotation_history,4).history())==origin_bytes,
          "annotation rewrote its original native history bytes");
    check(storage::NativeSeriesPrefix::capture(annotation_history,4).tail().pin()==origin.tail().pin(),
          "annotation failed to retain actual original ownership");
    const auto frame=storage::decode_temporal_segment(annotation_history.segment_bytes(0),native_data_limits().temporal);
    check(frame.source.grades.empty() && frame.tensors.empty() && frame.carriers.empty() && frame.declarations.has_value(),
          "annotation redundantly published its mathematical source");
    const auto origin_publication = origin.tail();
    for (const auto &field : origin_publication.endpoint().tensors)
        check(!field.tensor.values().materialized_entries(),
              "declaration annotation materialized an unchanged source field");
    const auto reference_frame = storage::make_temporal_segment(
        origin_publication.endpoint(), *annotated.endpoint(), {}, {},
        storage::RelationalCodeEncoding::Automatic, {}, frame.declarations);
    check(storage::encode_temporal_segment(reference_frame) ==
              storage::encode_temporal_segment(frame),
          "declared annotation differs from complete target comparison");
    check(annotated.endpoint()==annotation_history.tail().pin(),
          "record directory copied an actual native publication");
    annotated.write(directory/"annotated.bin");
    const auto annotated_restart=NativeRecord::read(directory/"annotated.bin");
    const auto& retained_field=annotated_restart.endpoint()->tensors.front().tensor.values();
    check(retained_field.retained_source() && !retained_field.materialized_entries(),
          "record loading expanded an unchanged field");
    const auto encoded=annotated_restart.encoded_tensor(888);
    check(encoded && encoded->publication().pin()==annotated_restart.endpoint() &&
          encoded->publication().source_namespace()==811 && encoded->encoding_sequence()==4,
          "encoded field lost its actual publication or retained base");
    storage::EncodedTensorReadStats read_stats;
    const auto selected=encoded->tensor().select(std::vector<std::size_t>{0,1,2},&read_stats);
    check(selected.state(0)==numeric::Rational{0} && !selected.state(1) &&
          selected.state(2)==numeric::Rational{-5,13} && read_stats.decoded_coefficients==1,
          "encoded field changed support or coefficient state");
    check(!retained_field.materialized_entries(),"selected operation expanded the publication field");
    auto changed_field_history = annotation_history;
    const auto changed_field_base = changed_field_history.tail();
    auto changed_field_endpoint = changed_field_base.endpoint();
    const auto &old_tensor = changed_field_endpoint.tensors.front().tensor;
    changed_field_endpoint.tensors.front().tensor = numeric::SparseRationalTensor(
        old_tensor.domain(), old_tensor.support(), numeric::SparseRationalVector(3, {{2, numeric::Rational{3,11}}}));
    ++changed_field_endpoint.sequence;
    changed_field_history.append(storage::make_temporal_segment(
        changed_field_base.endpoint(), changed_field_endpoint));
    const auto changed_field_record = NativeRecord::from_history(changed_field_history);
    check(!changed_field_record.encoded_tensor(888) &&
              changed_field_record.endpoint()->tensors.front().tensor.state(0) == numeric::Rational{0} &&
              !changed_field_record.endpoint()->tensors.front().tensor.state(1) &&
              changed_field_record.endpoint()->tensors.front().tensor.state(2) == numeric::Rational{3,11} &&
              annotated_restart.endpoint()->tensors.front().tensor.state(2) == numeric::Rational{-5,13},
          "changed field reused a stale encoding or changed its previous publication");
    changed_field_record.write(directory / "unencoded-field.bin");
    rejects([&]{(void)annotated_restart.encoded_tensor(999);});
    check(annotated_restart.encode()==annotated.encode() &&
          annotated_restart.at("observation").rational()==numeric::Rational{"2/7"},
          "ordered annotation changed after restart");
    const auto annotated_derived=annotated_restart.derive(N(numeric::Rational{"3/11"}));
    const auto historical=annotated_derived.historical_publication(annotated_restart.endpoint()->sequence);
    const auto historical_field=historical.encoded_tensor(888);
    check(historical_field && historical_field->publication().sequence()==historical.endpoint()->sequence &&
          historical_field->tensor().value(2)==numeric::Rational{-5,13},
          "historical encoded field used the later publication");
    check(annotated_derived.publication_history()->segment_count()==2 &&
          annotated_derived.parent()->at("observation").rational()==numeric::Rational{"2/7"},
          "derived record discarded ordered history or parent values");
    check(NativeRecord::decode(annotated_derived.encode()).encode()==annotated_derived.encode(),
          "derived ordered history changed after restart");
    check(storage::encode_rfcx(record.endpoint()->source).bytes==original,"complex changed by record admission");
    check(record.at("cell").at("value").rational()==numeric::Rational{"2/7"},"exact rational recoded");
    const auto owned=record.accession_directory()->cell({core::Grade{1},77}).declared_values();
    check(owned.size()==1,"description is not attached to its actual cell");
    check(record.declared_value(owned.front()->key).at("name").text()=="tumor","owned cell declaration unavailable");
    check(record.accession_directory()->occurrence({core::Grade{1},core::OccurrenceKey{100}}).declared_values().size()==1,"inactive occurrence description lost");
    record.write(directory/"data.bin");
    check(storage::read_file_exact(directory/"data.bin",native_data_limits().max_archive_bytes)==record.encode(),
          "streamed publication differs from the canonical complete archive");
    rejects([&]{record.write(directory/"data.bin");});
    check(storage::read_file_exact(directory/"data.bin",native_data_limits().max_archive_bytes)==record.encode(),
          "failed exclusive publication changed the previous file");
    const auto reopened=NativeRecord::read(directory/"data.bin");
    const auto current = reopened.at("cell").historical_publication(reopened.endpoint()->sequence);
    check(current.accession_directory() == reopened.accession_directory() &&
              current.source_namespace() == reopened.source_namespace() &&
              current.encode() == reopened.encode(),
          "current publication changed directory authority, original source or complete bytes");
    equal_value(current, reopened);
    rejects([&] { (void)reopened.historical_publication(reopened.endpoint()->sequence + 1); });
    const auto concurrent_pins = concurrent_access(NativeRecord::decode(reopened.encode()));
    check(concurrent_pins.size() == 24 &&
              concurrent_pins[0].tensor()->state(0) == numeric::Rational{0} &&
              !concurrent_pins[0].tensor()->state(1) &&
              concurrent_pins[0].tensor()->state(2) == numeric::Rational{-5,13} &&
              concurrent_pins[1].tensor()->state(0) == numeric::Rational{-7,5} &&
              concurrent_pins[2].tensor()->state(0) == numeric::Rational{0},
          "concurrent tensor pins lost exact coefficients, zero, absence or lifetime");
    const auto bytes=reopened.at("cell").at("evidence").bytes();
    check(std::equal(raw.begin(),raw.end(),bytes.begin(),bytes.end()),"original bytes, NUL or trailing zeros lost");
    check(reopened.encode()==record.encode(),"canonical native bytes changed on replay");
    const auto derived=reopened.derive(N(N::Object{{{"result",N(numeric::Rational{"3/11"})}}}));
    check(derived.parent()->at("cell").at("value").rational()==numeric::Rational{"2/7"},"derived publication discarded earlier source");
    check(storage::encode_rfcx(derived.endpoint()->source).bytes==original,"derived record replaced the complex");
    check(derived.accession_directory()->cell({core::Grade{1},77}).declared_values().size()==1,"derived record duplicated the original cell description");
    derived.write(directory/"derived.bin");
    const auto derived_reopened=NativeRecord::read(directory/"derived.bin");
    const auto parent_bytes=derived_reopened.parent()->at("cell").at("evidence").bytes();
    check(std::equal(raw.begin(),raw.end(),parent_bytes.begin(),parent_bytes.end()),
          "parent source bytes changed after restart");
    std::vector<std::byte> added{std::byte{255},std::byte{0},std::byte{128},std::byte{0}};
    const auto extended=reopened.derive(N(N::Object{{{"new",N(N::ByteSpan{43,0,added.size(),false})}}}),
                                        {RawByteSource(43,added)});
    added[0]=std::byte{0};
    check(extended.at("new").bytes().front()==std::byte{255},
          "caller can mutate an admitted byte owner");
    check(extended.parent()->at("cell").at("value").rational()==numeric::Rational{"2/7"},
          "owner reuse changes the parent value");
    const auto extended_reopened=NativeRecord::decode(extended.encode());
    check(std::equal(extended.at("new").bytes().begin(),extended.at("new").bytes().end(),
                     extended_reopened.at("new").bytes().begin(),extended_reopened.at("new").bytes().end()),
          "compiled byte owner differs from decoded native cargo");
    check(storage::encode_rfcx(extended_reopened.endpoint()->source).bytes==original,
          "owner reuse changes repeated/inactive/zero occurrences");
    std::vector<std::byte> transferred{std::byte{255},std::byte{0}};
    auto *prior_alias=transferred.data();
    const auto transferred_source=RawByteSource::own(44,std::move(transferred));
    const auto rebound=transferred_source.rekey(46);
    const auto frozen=NativeRecord::create(N(N::ByteSpan{44,0,2,false}),815,{transferred_source});
    prior_alias[0]=std::byte{0};
    check(frozen.bytes().front()==std::byte{255},
          "a caller's prior vector alias changed admitted provenance bytes");
    check(rebound.key()==46 && rebound.bytes().front()==std::byte{255} &&
          rebound.bytes().data()!=transferred_source.bytes().data(),
          "rekey did not freeze a publicly transferred owner");
    const auto immutable=RawByteSource(47,raw);
    const auto shared=immutable.rekey(48);
    check(shared.key()==48 && immutable.key()==47 && shared.bytes().data()==immutable.bytes().data(),
          "independent raw identities duplicated immutable byte storage");
    rejects([&]{(void)immutable.rekey(0);});
    const auto encoded_owner=record.archive_source(45);
    check(encoded_owner.key()==45 &&
          NativeRecord::decode(encoded_owner.bytes()).encode()==record.encode(),
          "fresh encoded owner lost its source identity or canonical publication");
    rejects([&]{(void)NativeRecord::create(N(N::ByteSpan{42,raw.size(),1,false}),812,{RawByteSource(42,raw)},&endpoint);});
    rejects([&]{(void)NativeRecord::create(N(N::Owned{0,storage::StableCellKey{core::Grade{1},78},{},std::make_shared<const N>(true)}),812,{},&endpoint);});
    // A lawful declaration payload can contain an invalid application reference
    // cycle. This does not restrict cycles in the retained scientific complex.
    auto leaf=std::make_shared<const N>(std::string("leaf"));
    auto nested=std::make_shared<const N>(N::Object{{{"next",N(N::Owned{0,storage::StableCellKey{core::Grade{1},77},{},leaf})}}});
    auto outer=std::make_shared<const N>(N::Object{{{"next",N(N::Owned{0,storage::StableCellKey{core::Grade{1},77},{},nested})}}});
    const auto finite=NativeRecord::create(N(N::Owned{0,storage::StableCellKey{core::Grade{1},77},{},outer}),813,{},&endpoint);
    const auto copied=NativeRecord::create(finite.materialize_value(),816,finite.raw_sources(),&endpoint);
    check(copied.accession_directory()->cell({core::Grade{1},77}).declared_values().size()==3,
          "materialized nested declarations lost their original owner bindings");
    const auto unbound=NativeRecord::create(finite.unbound_value(),817,finite.raw_sources());
    check(unbound.at("next").at("next").text()=="leaf",
          "unbound declaration lost nested payloads");
    const auto operation=retain_operation(endpoint,818,{},record,N(numeric::Rational{"3/11"}));
    check(operation.accession_directory()->cell({core::Grade{1},77}).declared_values().empty() &&
          operation.accession_directory()->occurrence({core::Grade{1},core::OccurrenceKey{100}}).declared_values().empty(),
          "foreign request bindings attached to equal keys in the result complex");
    check(operation.at("request").at("cell").at("value").rational()==numeric::Rational{"2/7"},
          "unbound operation request lost its exact payload");
    const auto original_request=NativeRecord::decode(operation.at("request_publication").bytes());
    check(original_request.source_namespace()==record.source_namespace() &&
          original_request.encode()==record.encode() &&
          original_request.accession_directory()->cell({core::Grade{1},77}).declared_values().size()==1 &&
          original_request.accession_directory()->occurrence({core::Grade{1},core::OccurrenceKey{100}}).declared_values().size()==1,
          "request archive lost its original namespace or owner bindings");
    const auto derived_operation=retain_operation(record,finite,N(numeric::Rational{"3/11"}));
    check(derived_operation.accession_directory()->cell({core::Grade{1},77}).declared_values().size()==1,
          "derived operation imported foreign request owners");
    auto cycle=*finite.endpoint();
    const auto rows=finite.accession_directory()->cell({core::Grade{1},77}).declared_values();
    const auto root=std::find_if(rows.begin(),rows.end(),[&](const auto &row){return row->type==finite.type();});
    check(root!=rows.end(),"owned root type is missing");
    std::vector<storage::TemporalDeclaredValue *> products;
    for(auto &row:cycle.declarations->values)if(row.cell && row.type==finite.type())products.push_back(&row);
    check(products.size()==2,"cycle fixture has wrong product family");
    auto &reference=std::get<composite::DeclaredValue::Product>(std::get<composite::DeclaredValue::Product>(products[0]->value.data).elements[0].data).elements;
    reference[0]=composite::DeclaredValue::integer(numeric::Integer(std::to_string(products[1]->key)));
    reference[1]=composite::DeclaredValue::integer(numeric::Integer(std::to_string(products[1]->type.name_space)));
    reference[2]=composite::DeclaredValue::integer(numeric::Integer(std::to_string(products[1]->type.key)));
    reference[3]=composite::DeclaredValue::integer(numeric::Integer(std::to_string(products[1]->type.version)));
    storage::validate_temporal_endpoint(cycle);
    rejects([&]{(void)NativeRecord::from_endpoint(cycle,813);});
    check(NativeRecord::decode(record.encode()).encode()==record.encode(),"immutable byte decode changed source");
    const auto booleans=NativeRecord::create(N(N::Array{{N(true),N(false),N(true)}}),814);
    const auto truth=NativeRecord::decode(booleans.encode());
    check(truth.size()==3 && truth.at(0).boolean() && !truth.at(1).boolean(),"fixed-width array declaration overflowed or lost values");
    {
        std::vector<std::byte> large_bytes(2 * 1024 * 1024);
        for (std::size_t i = 0; i < large_bytes.size(); ++i)
            large_bytes[i] = std::byte(i % 251);
        const auto bulk = NativeRecord::create(
            N(N::Object{{{"raw", N(N::ByteSpan{91, 0, large_bytes.size(), false})}}}), 771234,
            {RawByteSource(91, large_bytes)});
        const auto body = N(N::Object{{{"raw", N(N::ByteSpan{91, 0, large_bytes.size(), false})},
                                       {"message", N(std::string("Updated description"))}}});
        const auto reference = NativeRecord::annotate(body, bulk);
        const auto ref_history = *reference.publication_history();
        check(reference.format_version() == 3 &&
                  ref_history.segment_bytes(0).size() < large_bytes.size() / 8,
              "annotation repeated unchanged raw payloads");
        check(storage::encode_native_series(
                  storage::NativeSeriesPrefix::capture(ref_history, bulk.endpoint()->sequence)
                      .history()) == bulk.encode(),
              "historical byte references changed the original publication");
        const auto alias = reference.derive(
            N(N::Object{{{"raw", N(N::ByteSpan{91, 0, large_bytes.size(), false})},
                         {"alias", N(N::ByteSpan{92, 0, large_bytes.size(), false})}}}),
            {reference.at("raw").byte_source().rekey(92)});
        const auto aliases = NativeRecord::decode(alias.encode());
        check(aliases.at("raw").bytes().data() == aliases.at("alias").bytes().data() &&
                  std::equal(large_bytes.begin(), large_bytes.end(),
                             aliases.at("raw").bytes().begin()),
              "historical aliases duplicated or changed their actual byte owner");
        rejects([&] {
            (void)NativeRecord::from_endpoint(*reference.endpoint(), reference.source_namespace());
        });
        auto legacy_endpoint = *bulk.endpoint();
        constexpr std::uint64_t control = 0x4741444745540201ULL;
        for (auto &value : legacy_endpoint.declarations->values) {
            if (value.key == control)
                std::get<composite::DeclaredValue::Product>(value.value.data).elements[0] =
                    composite::DeclaredValue::integer(numeric::Integer{2});
            if (value.key == control + 3) {
                value.type = {control, 7, 1};
                for (auto &entry :
                     std::get<composite::DeclaredValue::Repeated>(value.value.data).elements) {
                    auto &payload =
                        std::get<composite::DeclaredValue::Product>(entry.data).elements[2];
                    payload =
                        std::get<composite::DeclaredValue::Alternative>(payload.data).element.at(0);
                }
            }
        }
        const auto legacy = NativeRecord::from_endpoint(legacy_endpoint, bulk.source_namespace());
        const auto legacy_bytes = legacy.encode();
        const auto upgraded = NativeRecord::annotate(body, legacy);
        check(legacy.format_version() == 2 && upgraded.format_version() == 3 &&
                  storage::encode_native_series(
                      storage::NativeSeriesPrefix::capture(*upgraded.publication_history(),
                                                           legacy.endpoint()->sequence)
                          .history()) == legacy_bytes &&
                  std::equal(large_bytes.begin(), large_bytes.end(),
                             NativeRecord::decode(upgraded.encode()).at("raw").bytes().begin()),
              "legacy source migration lost bytes or rewrote its original history");
        for (const auto field : {0, 1, 2}) {
            auto invalid = *reference.endpoint();
            for (auto &value : invalid.declarations->values)
                if (value.key == control + 3) {
                    auto &entry =
                        std::get<composite::DeclaredValue::Repeated>(value.value.data).elements[0];
                    auto &payload =
                        std::get<composite::DeclaredValue::Product>(entry.data).elements[2];
                    auto &address =
                        std::get<composite::DeclaredValue::Product>(
                            std::get<composite::DeclaredValue::Alternative>(payload.data)
                                .element[0]
                                .data)
                            .elements;
                    address[field] = composite::DeclaredValue::integer(
                        numeric::Integer(std::to_string(field == 0   ? bulk.source_namespace() + 1
                                                        : field == 1 ? invalid.sequence
                                                                     : 999999)));
                }
            auto malformed = *bulk.publication_history();
            const auto malformed_tail = malformed.tail();
            auto update = storage::make_temporal_declaration_update(
                malformed_tail.endpoint().declarations, invalid.declarations,
                composite::identity_declaration_maps(
                    malformed_tail.endpoint().declarations->registry));
            malformed.append_target(invalid, {}, {}, storage::RelationalCodeEncoding::Automatic, {},
                                    std::move(update));
            rejects([&] { (void)NativeRecord::from_history(malformed); });
        }
    }
    N::Array family;
    N::Object item;
    item.members.emplace_back("text", N(std::string("same\0text", 9)));
    item.members.emplace_back("q", N(numeric::Rational{"-17/101"}));
    item.members.emplace_back("raw", N(N::ByteSpan{42, 0, raw.size(), false}));
    item.members.emplace_back("tensor", N(field.tensor));
    family.elements.emplace_back(std::move(item));
    family.elements.emplace_back(numeric::Integer{numeric::Integer{1} << 3073});
    family.elements.emplace_back(N::Array{{N(true), N{}, N(N::ExactNumber{numeric::Rational{"2/7"}, false})}});
    const auto expanded = record.derive(N(family));
    std::size_t supplied{};
    const N streamed(N::PackedArray{family.elements.size(), [&](const auto &emit) {
        for (const auto &value : family.elements) { emit(value); ++supplied; }
    }});
    const auto compact = record.derive(streamed);
    check(supplied == family.elements.size(), "packed producer was replayed during admission");
    check(compact.format_version() == 4 && expanded.format_version() == 3, "packed format negotiation");
    equal_value(compact, expanded);
    check(storage::encode_rfcx(compact.endpoint()->source).bytes == original, "packing changed source relations");
    const auto compact_reopened = NativeRecord::decode(compact.encode());
    const auto encoding = compact_reopened.packed_encoding(0);
    check(encoding && encoding->bit_count && !encoding->bytes.empty(), "packed encoding was not exposed");
    const auto &declarations = compact_reopened.endpoint()->declarations.value();
    const composite::CompiledLayout packed_layout(declarations.registry, encoding->type,
        declarations.profile, declarations.resource_limits);
    const auto selected_bits = storage::encode_declared_bits(packed_layout, compact_reopened.at(0).value());
    check(selected_bits.bit_count == encoding->bit_count &&
              selected_bits.bytes.size() == encoding->bytes.size() &&
              std::equal(selected_bits.bytes.begin(), selected_bits.bytes.end(), encoding->bytes.begin()),
          "packed encoding does not bind its selected value");
    check(!expanded.packed_encoding(0), "ordinary array reported a packed encoding");
    rejects([&] { (void)compact_reopened.packed_encoding(compact_reopened.size()); });
    equal_value(compact_reopened, expanded);
    check(compact_reopened.encode() == compact.encode(), "packed canonical encoding changed");
    compact.write(directory / "packed.rfcxs");
    equal_value(NativeRecord::read(directory / "packed.rfcxs"), expanded);
    const auto escaped = [&] { const auto owner = NativeRecord::decode(compact.encode()); return owner.at(0).at("q"); }();
    check(escaped.rational() == numeric::Rational{"-17/101"}, "decoded child lost its payload owner");
    const auto nested_packed = record.derive(N(N::PackedArray{1, [&](const auto &emit) { emit(streamed); }}));
    equal_value(nested_packed.at(0), expanded);
    equal_value(NativeRecord::decode(nested_packed.encode()).at(0), expanded);
    const auto descendant = compact.derive(N(std::string("later")));
    check(descendant.parent().has_value(), "packed parent lost");
    equal_value(*descendant.parent(), expanded);
    check(NativeRecord::decode(descendant.encode()).parent()->at(0).at("q").rational() == numeric::Rational{"-17/101"},
          "packed ancestor did not survive publication");
    rejects([&] { (void)record.derive(N(N::PackedArray{1, [](const auto &) {}})); });
    rejects([&] { (void)record.derive(N(N::PackedArray{0, [](const auto &emit) { emit(N(true)); }})); });
    rejects([&] { (void)record.derive(N(N::PackedArray{1, [](const auto &) { throw std::runtime_error("producer failure"); }})); });
    const auto empty = record.derive(N(N::PackedArray{0, [](const auto &) {}}));
    check(empty.size() == 0 && NativeRecord::decode(empty.encode()).size() == 0, "empty packed array changed");
    for (const auto index : {std::size_t{0}, std::size_t{3}, std::size_t{6}}) {
        for (const auto &publication : {compact, descendant}) {
            auto bad = *publication.endpoint();
            for (auto &value : bad.declarations->values)
                if (value.type == compact.type()) {
                    auto &reference = std::get<composite::DeclaredValue::Product>(
                        std::get<composite::DeclaredValue::Repeated>(value.value.data).elements.at(0).data).elements;
                    reference.at(index) = composite::DeclaredValue::integer(numeric::Integer{"999999999999999999"});
                }
            rejects([&] { (void)NativeRecord::from_endpoint(bad, publication.source_namespace()); });
        }
    }
    std::cout
        << "native records preserve owned descriptions, original bytes and source relations\n";
    return 0;
 } catch (const std::exception &e) {
     std::cerr << e.what() << '\n';
     return 1;
 }
}
