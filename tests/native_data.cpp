#include <gadget/native_data.hpp>
#include <gadget/native_gluing.hpp>
#include <rfc/storage/native_atlas_codec.hpp>
#include <iostream>
#include <future>

namespace {
using namespace rfc;
using namespace gadget;
using Q = numeric::Rational;
using V = composite::DeclaredValue;
std::size_t checks = 0;
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
template<class F> void refuses(F operation, const char *message) {
    bool refused = false;
    try { operation(); } catch (const std::exception &) { refused = true; }
    check(refused, message);
}
composite::DeclarationId type(std::uint64_t key) { return {71, key, 2}; }
runtime::ExecutionLease lease() {
    runtime::ExecutionHostPolicy policy;
    const auto max = std::numeric_limits<std::size_t>::max();
    policy.name = "native-data-contract";
    policy.limits = {max, max, max, max, max};
    return runtime::ExecutionLease(policy);
}
runtime::CpuResourceCoordinator coordinator() {
    const auto hardware = runtime::detect_hardware();
    runtime::CpuResourcePolicy policy;
    policy.worker_limit = hardware.allowed_cpus.size();
    policy.tile_limit = std::numeric_limits<std::size_t>::max();
    policy.declared_buffer_byte_limit = machine::checked_size(*hardware.advisory_memory_headroom());
    return runtime::CpuResourceCoordinator(hardware,policy);
}
storage::TemporalEndpoint program_source() {
    core::RelationSignature signature;
    const auto producer = signature.aspects().add("producer", composite::AspectRole::RelationRole);
    const auto consumer = signature.aspects().add("consumer", composite::AspectRole::RelationRole);
    const auto relation = signature.add_type("authored-collective");
    core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(core::Grade{3});
    std::vector<core::BoundaryOccurrenceInput> occurrences;
    for (std::size_t i = 0; i < 130; ++i)
        occurrences.push_back({10 + i % 2, i % 5 != 0, i % 3 == 0, Q{0},
            {{producer, i % 2 != 0}, {consumer, i % 7 != 0}}, core::OccurrenceKey{1000 + i}});
    builder.add_grade_one(core::RelationKey{17}, relation, occurrences);
    storage::TemporalEndpoint endpoint(builder.freeze(), 4);
    storage::TemporalDeclarationState declarations;
    declarations.namespace_key = 71;
    declarations.registry = composite::DeclarationRegistry({
        composite::Declaration::boolean(type(1)), composite::Declaration::integer(type(2)),
        composite::Declaration::product(type(3), {{91,type(1),0},{12,type(1),0},{44,type(1),0},{61,type(1),0}}),
        composite::Declaration::product(type(4), {{700,type(2),0},{900,type(2),0},{500,type(2),0}})});
    declarations.source_layouts.push_back({16, type(3),
        {{91,"producer"},{12,"consumer"},{44,"exists"},{61,"orientation"}}});
    composite::BinaryProgram first(7);
    first.append({composite::BinaryOpcode::LoadAspect,0,0,0,0,0});
    first.append({composite::BinaryOpcode::LoadAspect,1,0,0,0,1});
    first.append({composite::BinaryOpcode::LoadAspect,2,0,0,0,2});
    first.append({composite::BinaryOpcode::LoadAspect,3,0,0,0,3});
    first.append({composite::BinaryOpcode::Xor,4,0,1});
    first.append({composite::BinaryOpcode::And,5,4,2});
    first.append({composite::BinaryOpcode::AndNot,6,5,3});
    declarations.programs.push_back({27,16,first,6});
    composite::BinaryProgram second(4);
    second.append({composite::BinaryOpcode::LoadAspect,0,0,0,0,0});
    second.append({composite::BinaryOpcode::LoadAspect,1,0,0,0,1});
    second.append({composite::BinaryOpcode::Not,2,1});
    second.append({composite::BinaryOpcode::Select,3,0,1,2});
    declarations.programs.push_back({28,16,second,3});
    declarations.values.push_back({33,type(4),0,0,storage::StableCellKey{core::Grade{0},10},{},
        V::product({V::integer(4),V::integer(1),V::integer(83)})});
    endpoint.declarations = declarations;
    storage::validate_temporal_endpoint(endpoint);
    return endpoint;
}
storage::TemporalEndpoint gluing_source(bool common) {
    core::RelationSignature signature;
    const auto edge = signature.add_type("edge");
    const auto multiway = signature.add_type("multiway");
    const auto filling = signature.add_type("filling");
    core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(core::Grade{2});
    builder.add_grade_one(core::RelationKey{100}, edge, {
        {10,true,true,Q{1}, {},core::OccurrenceKey{1}},
        {11,true,false,Q{1}, {},core::OccurrenceKey{2}},
        {10,true,false,Q{1,3}, {},core::OccurrenceKey{3}},
        {10,true,true,Q{1,3}, {},core::OccurrenceKey{4}},
        {11,false,true,Q{2}, {},core::OccurrenceKey{5}},
        {11,true,false,Q{0}, {},core::OccurrenceKey{6}}});
    if (!common) {
        builder.add_grade_one(core::RelationKey{200}, multiway, {
            {10,true,true,Q{1},{},core::OccurrenceKey{7}},
            {11,true,false,Q{1,2},{},core::OccurrenceKey{8}},
            {12,true,false,Q{1,2},{},core::OccurrenceKey{9}}});
        builder.add_higher(core::Grade{2},core::RelationKey{300},filling,{
            {1,true,false,Q{1},{},core::OccurrenceKey{10}},
            {1,true,true,Q{1},{},core::OccurrenceKey{11}}});
    }
    storage::TemporalEndpoint endpoint(builder.freeze());
    const auto n = endpoint.source.grade(core::Grade{1}).cell_count();
    endpoint.carriers.push_back({41,storage::TemporalCarrierKind::SourceGrade,core::Grade{1},{},0,
                                 numeric::IndexDomain::linear(n)});
    storage::TemporalTensor tensor;
    tensor.key=51; tensor.role=storage::TemporalTensorRole::Field; tensor.output=41;
    std::vector<numeric::SparseRationalEntry> values;
    if (!common) values.push_back({1,Q{4}});
    tensor.tensor=numeric::SparseRationalTensor(numeric::IndexDomain::linear(n),
        numeric::StructuralSupport::full(n),numeric::SparseRationalVector(n,values));
    endpoint.tensors.push_back(tensor);
    storage::validate_temporal_endpoint(endpoint);
    return endpoint;
}
void gluing(const std::filesystem::path &root) {
    const NativeStudy common(gluing_source(true),101), left(gluing_source(false),102), right(gluing_source(false),103);
    storage::NativeAtlasManifest manifest;
    manifest.key=71; manifest.revision=1;
    for (std::uint64_t side=2;side<=3;++side) {
        storage::NativeSourceEmbedding embedding;
        embedding.key=side; embedding.from=1; embedding.to=side;
        embedding.cells={{{core::Grade{0},10},{core::Grade{0},10}},
                         {{core::Grade{0},11},{core::Grade{0},11}},
                         {{core::Grade{1},100},{core::Grade{1},100}}};
        for(std::uint64_t i=1;i<=6;++i)
            embedding.occurrences.push_back({{core::Grade{1},core::OccurrenceKey{i}},
                                             {core::Grade{1},core::OccurrenceKey{i}}});
        manifest.embeddings.push_back(embedding);
    }
    manifest.pushouts.push_back({91,2,3,{},{}});
    storage::NativeAtlas atlas({common.atlas_source(1),left.atlas_source(2),right.atlas_source(3)},manifest);
    const inference::NativeGluingAvailability available{{101,0},{102,0},{103,0}};
    const auto plan=inference::PreparedNativeGluing::prepare({atlas,91,8},available,lease());
    const auto field=storage::NativeAccessionRequest::chain(51,core::Grade{1});
    const std::vector<NativeGluingSection> sections{{core::Grade{1},61,left.access(field),right.access(field)}};
    auto resources=coordinator();
    const auto result=execute_native_gluing(plan,sections,available,true,resources.reserve(4),lease());
    const auto serial=execute_native_gluing(plan,sections,available,true,resources.reserve(1),lease());
    check(result.sections.size()==1 && result.sections[0].size()==2,"native lineage lost");
    for(std::size_t i=0;i<2;++i)
        check(result.sections[0][i].chain().values()==serial.sections[0][i].chain().values(),"CPU worker parity");
    check(result.sections[0][0].plan().same_owner(plan),"native prepared owner lost");
    const auto endpoint=native_gluing_endpoint(result,{{0,0,71,81},{0,1,72,82}});
    check(storage::encode_rfcx(endpoint.source).bytes==storage::encode_rfcx(plan.constructed_source()->source).bytes,
          "binary export replaced the relational complex");
    const storage::NativeTemporalSeries series(endpoint,301,native_data_limits());
    storage::write_file_exact(root/"glued.dat",storage::encode_native_series(series));
    storage::write_file_exact(root/"inputs.rfcxa",storage::encode_native_atlas(plan.atlas()));
    const auto loaded=NativeStudy::read(root/"glued.dat");
    for(std::size_t i=0;i<2;++i) {
        const auto &chain=result.sections[0][i].chain();
        const auto pin=loaded.access(storage::NativeAccessionRequest::chain(71+i,chain.grade()));
        for(const auto &[coordinate,value]:chain.values()) {
            const auto key=storage::stable_cell_key(endpoint.source,{chain.grade(),static_cast<std::uint32_t>(coordinate)});
            check(pin.state(key.key)==value,"native binary supported coordinate differs");
        }
        check(pin.support_count()==chain.values().size(),"native binary support changed");
    }
    const auto restored=storage::read_native_atlas(root/"inputs.rfcxa");
    const auto replay=inference::PreparedNativeGluing::prepare({restored,91,8},available,lease());
    const auto replayed=execute_native_gluing(replay,{{core::Grade{1},61,
        restored.access(2,field),restored.access(3,field)}},available,true,resources.reserve(2),lease());
    for(std::size_t i=0;i<2;++i)
        check(replayed.sections[0][i].chain().values()==result.sections[0][i].chain().values(),"native atlas replay differs");
    refuses([&]{native_gluing_endpoint(result,{{0,0,71,81}});},"partial native lineage export accepted");
    refuses([&]{resources.reserve(0);},"zero CPU workers accepted");
    auto cancelled=std::make_shared<std::atomic_bool>(true);
    const auto current=lease();
    refuses([&]{execute_native_gluing(plan,sections,available,true,resources.reserve(1),runtime::ExecutionLease(current.policy(),cancelled));},
            "cancelled native gluing accepted");
}
} // namespace
int main(int argc,char **argv) {
    try {
        if(argc!=2)throw std::invalid_argument("expected output directory");
        const std::filesystem::path root(argv[1]); std::filesystem::create_directories(root);
        auto endpoint=program_source();
        std::vector<std::byte> original{std::byte{255},std::byte{0},std::byte{128},std::byte{10},std::byte{42},std::byte{13}};
        storage::write_file_exact(root/"raw.bin",original);
        const RawSpanDeclaration span{type(4),500,900,700};
        auto transferred=original;
        auto *prior_alias=transferred.data();
        const auto host_source=RawByteSource::own(83,std::move(transferred));
        const NativeStudy frozen_study(program_source(),701,{host_source});
        prior_alias[2]=std::byte{0};
        check(frozen_study.raw_span(33,span).bytes()[1]==std::byte{128},
              "prior caller alias changed a native study's provenance bytes");
        std::optional<PinnedRawSpan> retained;
        {
            NativeStudy study(endpoint,701,{RawByteSource::read(83,root/"raw.bin")});
            retained=study.raw_span(33,span);
            endpoint.declarations->values[0].value=V::product({V::integer(1),V::integer(0),V::integer(83)});
            check(retained->bytes().size()==4 && retained->bytes()[0]==std::byte{0},"mutable declaration alias");
            storage::write_file_exact(root/"raw.bin",std::vector<std::byte>{std::byte{99}});
            check(retained->bytes().size()==4 && retained->bytes()[1]==std::byte{128},"original raw bytes changed");
            for(std::uint64_t key:{27,28}) {
                const auto program=study.program(key,16,{core::Grade{1},core::OccurrenceKey{1000}});
                composite::BinaryExecutionWorkspace scratch;
                const auto bits=program.materialize(scratch);
                std::size_t expected_count=0;
                check(bits.size()==130,"program lost original occurrence domain");
                for(std::size_t i=0;i<130;++i) {
                    const bool producer=i%2!=0,consumer=i%7!=0;
                    const bool expected=key==27?((producer!=consumer) && i%5!=0 && i%3!=0):
                                                (producer?consumer:!consumer);
                    check(bits.get(i)==expected,"authored machine program differs");
                    if(expected)++expected_count;
                }
                check(program.count(scratch)==expected_count,"native packed count differs");
            }
            refuses([&]{study.program(999,16,{core::Grade{1},core::OccurrenceKey{1000}});},"unknown program accepted");
            refuses([&]{study.raw_span(33,{type(4),500,500,700});},"ambiguous raw span factors accepted");
            refuses([&]{study.raw_span(33,{type(3),500,900,700});},"foreign raw span type accepted");
        }
        check(retained->bytes()[3]==std::byte{42},"raw span lost native/source owners");
        endpoint=program_source();
        auto definitions=endpoint.declarations->registry.definitions();
        definitions.push_back(composite::Declaration::rational_tensor(type(5),{{1,257,{}}}));
        definitions.push_back(composite::Declaration::product(type(6),{{19,type(2),0},{7,type(1),0},{31,type(5),0}}));
        endpoint.declarations->registry=composite::DeclarationRegistry(definitions);
        const composite::CompiledLayout layout(endpoint.declarations->registry,type(6));
        const auto value=V::product({V::integer(numeric::Integer("-123456789012345678901234567890")),V::boolean(false),
            V::tensor(numeric::SparseRationalTensor(numeric::IndexDomain::linear(257),
                numeric::StructuralSupport::sparse(257,{0,1,256}),
                numeric::SparseRationalVector(257,{{0,Q{2,3}},{256,Q{-7,13}}})))});
        const auto encoded=storage::encode_declared_bits(layout,value);
        endpoint.declarations->values.push_back({34,type(4),0,0,{},
            storage::StableOccurrenceKey{core::Grade{1},core::OccurrenceKey{1000}},
            V::product({V::integer(numeric::Integer(std::to_string(encoded.bytes.size()))),V::integer(0),V::integer(84)})});
        NativeStudy typed(endpoint,701,{RawByteSource(84,encoded.bytes)});
        const auto decoded=typed.decode_raw_span(34,span,type(6),encoded.bit_count);
        check(composite::declared_value_equal(value,decoded.value),"custom raw composite payload differs");
        check(decoded.source.declaration().declared_owner()->occurrence.has_value(),"raw decoded value lost occurrence owner");
        const auto cached_layout = typed.compiled_layout(type(6));
        const auto alias = typed;
        check(cached_layout == typed.compiled_layout(type(6)) &&
                  cached_layout == alias.compiled_layout(type(6)),
              "study aliases recompiled the same immutable layout");
        check(&cached_layout->declarations() == &decoded.layout.declarations(),
              "raw decoding bypassed the study declaration owner");
        refuses([&] { (void)typed.compiled_layout(type(999)); },
                "unknown layout declaration accepted");
        NativeStudy concurrent(endpoint,701,{RawByteSource(84,encoded.bytes)});
        auto first_layout = std::async(std::launch::async,
            [&] { return concurrent.compiled_layout(type(6)); });
        auto second_layout = std::async(std::launch::async,
            [&] { return concurrent.compiled_layout(type(6)); });
        check(first_layout.get() == second_layout.get(),
              "concurrent callers published different compiled layout owners");
        auto first_decode=std::async(std::launch::async,[&]{return concurrent.decode_raw_span(34,span,type(6),encoded.bit_count);});
        auto second_decode=std::async(std::launch::async,[&]{return concurrent.decode_raw_span(34,span,type(6),encoded.bit_count);});
        const auto first_value=first_decode.get();
        const auto second_value=second_decode.get();
        check(&first_value.layout.declarations()==&second_value.layout.declarations(),
              "concurrent raw decoding copied the same endpoint registry");
        const auto reused=typed.decode_raw_span(34,span,type(6),encoded.bit_count);
        check(&decoded.layout.declarations()==&reused.layout.declarations(),
              "repeated raw decoding copied the same endpoint registry");
        check(&decoded.layout.declarations()!=&first_value.layout.declarations(),
              "independent study owners shared a declaration registry");
        auto limited_endpoint = endpoint;
        ++limited_endpoint.sequence;
        limited_endpoint.declarations->resource_limits.max_integer_bits = 64;
        NativeStudy limited(limited_endpoint,701,{RawByteSource(84,encoded.bytes)});
        const auto limited_layout = limited.compiled_layout(type(6));
        check(limited_layout != cached_layout && limited_layout->limits().max_integer_bits == 64,
              "another publication reused a broader cached declaration profile");
        refuses([&] { (void)limited.decode_raw_span(34,span,type(6),encoded.bit_count); },
                "raw decoding bypassed its actual source integer limit");
        check(limited_layout == limited.compiled_layout(type(6)),
              "failed value admission changed its immutable cached layout");
        const auto retained_layout = [&] {
            NativeStudy temporary(endpoint,701,{RawByteSource(84,encoded.bytes)});
            return temporary.compiled_layout(type(6));
        }();
        check(composite::declared_value_equal(
                  storage::decode_declared_bits(*retained_layout,encoded.bytes,encoded.bit_count), value),
              "compiled layout lost its registry after study destruction");
        const auto retained_decode=[&]{
            NativeStudy temporary(endpoint,702,{RawByteSource(84,encoded.bytes)});
            return temporary.decode_raw_span(34,span,type(6),encoded.bit_count);
        }();
        check(composite::declared_value_equal(retained_decode.value,value) &&
                  storage::encode_declared_bits(retained_decode.layout,retained_decode.value)==encoded,
              "decoded source lost its layout owner after study destruction");
        refuses([&]{typed.decode_raw_span(34,span,type(6),encoded.bit_count+1);},"invalid raw payload bit count accepted");
        endpoint=program_source();
        storage::NativeTemporalSeries series(endpoint,701,native_data_limits());
        const auto actual_endpoint=series.tail();
        NativeStudy adopted(actual_endpoint);
        check(adopted.endpoint()==actual_endpoint.pin(),"native study copied its actual canonical endpoint");
        const auto deferred=NativeStudy::from_publication(actual_endpoint,{RawByteSource(83,original)});
        check(deferred.endpoint()==actual_endpoint.pin(),"deferred study replaced its publication owner");
        auto updated_endpoint = endpoint;
        updated_endpoint.sequence = endpoint.sequence + 1;
        updated_endpoint.declarations->values.at(0).value = V::product(
            {V::integer(0), V::integer(0), V::integer(83)});
        auto changing_history = series;
        changing_history.append(storage::make_temporal_segment(endpoint, updated_endpoint));
        const auto retained_tail = changing_history.tail().pin();
        const auto latest = NativeStudy::from_history(changing_history, {RawByteSource(83,original)});
        const auto explicit_latest = NativeStudy::from_history(
            changing_history, {}, changing_history.last_sequence());
        check(latest.endpoint() == retained_tail && explicit_latest.endpoint() == retained_tail,
              "latest history selection replayed its retained endpoint");
        const auto earlier = NativeStudy::from_history(changing_history, {}, series.first_sequence());
        check(earlier.endpoint() == series.tail().pin(), "earlier base selection replaced its owner");
        changing_history = series;
        check(latest.endpoint() == retained_tail &&
              latest.raw_span(33,span).bytes().empty(), "history replacement changed a retained tail");
        const auto old_value = earlier.accession_directory()->pin_snapshot()->declarations->values.at(0).value;
        check(composite::declared_value_equal(old_value, endpoint.declarations->values.at(0).value),
              "earlier cutoff exposed a later declaration");
        refuses([&]{NativeStudy::from_history(series, {}, series.first_sequence() - 1);},
                "history selection admitted a cutoff before the source");
        refuses([&]{NativeStudy::from_history(series, {}, series.last_sequence() + 1);},
                "history selection admitted a future cutoff");
        storage::write_file_exact(root/"changed-history.payload",storage::encode_native_series(
            [&] {auto history = series; history.append(storage::make_temporal_segment(endpoint, updated_endpoint)); return history;}()));
        const auto loaded_latest = NativeStudy::read(root/"changed-history.payload");
        const auto loaded_earlier = NativeStudy::read(root/"changed-history.payload", {}, {}, series.first_sequence());
        check(loaded_latest.endpoint()->sequence == updated_endpoint.sequence &&
              composite::declared_value_equal(loaded_latest.endpoint()->declarations->values.at(0).value,
                                               updated_endpoint.declarations->values.at(0).value) &&
              composite::declared_value_equal(loaded_earlier.endpoint()->declarations->values.at(0).value,
                                               endpoint.declarations->values.at(0).value),
              "native file reading changed the selected publication");
        std::promise<void> release;
        const auto ready=release.get_future().share();
        std::vector<std::future<std::shared_ptr<const storage::NativeAccessionDirectory>>> directories;
        for(std::size_t i=0;i<4;++i)
            directories.push_back(std::async(std::launch::async,[deferred,ready] {
                ready.wait();return deferred.accession_directory();
            }));
        release.set_value();
        const auto directory=directories.front().get();
        for(std::size_t i=1;i<directories.size();++i)
            check(directories[i].get()==directory,"concurrent preparation created distinct accession directories");
        check(directory->pin_snapshot()==actual_endpoint.pin() &&
              deferred.raw_span(33,span).bytes()[2]==std::byte{10},
              "deferred accession changed its actual source or byte owner");
        auto restricted=native_accession_limits();restricted.max_index_entries=1;
        refuses([&]{NativeStudy(actual_endpoint,{},restricted);},
                "eager constructor deferred a caller's admission resource refusal");
        const auto bytes=storage::encode_native_series(series);
        storage::write_file_exact(root/"native.payload",bytes);
        check(detect_native_data(bytes)==NativeDataFormat::TemporalSeries,"series signature detection");
        const auto reopened=NativeStudy::read(root/"native.payload",{RawByteSource(83,original)});
        check(reopened.raw_span(33,span).bytes()[2]==std::byte{10},"binary declaration/raw byte restoration");
        check(reopened.endpoint()->source.grade(core::Grade{3}).cell_count()==0,"empty grade lost");
        const auto complex=storage::encode_rfcx(endpoint.source).bytes;
        storage::write_file_exact(root/"source.payload",complex);
        check(NativeStudy::read(root/"source.payload",{},702).endpoint()->source.grade(core::Grade{1}).incidence_count()==130,
              "RFCX original occurrence fibre lost");
        refuses([&]{NativeStudy::read(root/"source.payload");},"RFCX inferred a namespace");
        refuses([&]{NativeStudy::read(root/"native.payload",{},702);},"foreign namespace accepted");
        refuses([&]{NativeStudy::read(root/"native.payload",{},{},3);},"invalid cutoff accepted");
        auto limits=native_data_limits(); limits.max_archive_bytes=bytes.size()-1;
        refuses([&]{NativeStudy::read(root/"native.payload",{},{},{},limits);},"caller byte resource ignored");
        auto corrupt=bytes;corrupt.back()^=std::byte{1}; storage::write_file_exact(root/"corrupt.payload",corrupt);
        refuses([&]{NativeStudy::read(root/"corrupt.payload");},"corrupt native data accepted");
        NativeStudy no_raw(endpoint,701);
        refuses([&]{no_raw.raw_span(33,span);},"missing original bytes reconstructed");
        NativeStudy short_raw(endpoint,701,{RawByteSource(83,std::span<const std::byte>(original).first(2))});
        refuses([&]{short_raw.raw_span(33,span);},"out-of-source byte span accepted");
        gluing(root);
        std::cout<<checks<<" native data, machine program, raw provenance and binary gluing checks passed\n";
    } catch(const std::exception &error) {std::cerr<<error.what()<<'\n';return 1;}
}
