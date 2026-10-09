#include <gadget/native_record.hpp>
#include <gadget/study.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <iostream>

using namespace gadget;
using namespace rfc;
namespace {
void check(bool ok,const char *message) {if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F &&f) {bool threw=false;try{f();}catch(const std::exception &){threw=true;}check(threw,"invalid native record admitted");}
}
int main(int argc,char **argv) {
 try {
    const auto directory=std::filesystem::path(argc>1?argv[1]:"native-record-contracts");
    std::filesystem::remove_all(directory);std::filesystem::create_directories(directory);
    core::RelationSignature signature;const auto relation=signature.add_type("authored-group");
    core::RelationalComplexBuilder builder(signature);builder.declare_grade(core::Grade{3});
    builder.add_grade_one(core::RelationKey{77},relation,{
        {11,true,false,numeric::Rational{0},{},core::OccurrenceKey{99}},
        {11,false,true,numeric::Rational{0},{},core::OccurrenceKey{100}}});
    storage::TemporalEndpoint endpoint(builder.freeze(),4);
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
    check(annotated.endpoint()==annotation_history.tail().pin(),
          "record directory copied an actual native publication");
    annotated.write(directory/"annotated.bin");
    const auto annotated_restart=NativeRecord::read(directory/"annotated.bin");
    check(annotated_restart.encode()==annotated.encode() &&
          annotated_restart.at("observation").rational()==numeric::Rational{"2/7"},
          "ordered annotation changed after restart");
    const auto annotated_derived=annotated_restart.derive(N(numeric::Rational{"3/11"}));
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
    std::cout<<"native records preserve owned descriptions, original bytes and source relations\n";
    return 0;
 }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
