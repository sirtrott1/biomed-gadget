#include <gadget/homology_diagram.hpp>
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <rfc/inference/native_homology_diagram_codec.hpp>

namespace gadget {
namespace {
using N = NativeRecordValue;
using V = scientific_detail::ResultValue;
using Matrix = rfc::numeric::ExactMatrix;
using Product = rfc::inference::HomologyDiagramProduct;
using Provider = rfc::inference::PreparedNativeHomologyDiagram;
void require(bool ok, const char *message) { if (!ok) throw std::invalid_argument(message); }
void members(const NativeRecord &record, const std::set<std::string> &required) {
    require(record.kind() == NativeRecordKind::Object, "homology diagram requires an object");
    for (const auto &name : required) require(record.has(name), "homology diagram member absent");
    for (const auto &name : record.names()) require(required.contains(name), "unknown homology diagram member");
}
std::size_t count(const NativeRecord &record) {
    const auto n = record.integer();
    require(n >= 0 && n <= rfc::numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())), "homology diagram count is not machine addressable");
    return static_cast<std::size_t>(std::stoull(n.get_str()));
}
rfc::numeric::Rational rational(const NativeRecord &record) {
    if (record.kind() == NativeRecordKind::Integer) return rfc::numeric::Rational(record.integer());
    if (record.kind() != NativeRecordKind::Text) return record.rational();
    rfc::numeric::Rational q;
    require(q.set_str(record.text(), 10) == 0 && q.get_den() != 0, "invalid declared rational map coefficient");
    q.canonicalize(); return q;
}
Matrix map(const NativeRecord &record, rfc::numeric::ExactActionGuard &guard) {
    members(record, {"rows", "columns", "entries"});
    const auto rows = count(record.at("rows")), columns = count(record.at("columns"));
    guard.dimension(rows); guard.dimension(columns);
    const auto entries = record.at("entries");
    require(entries.kind() == NativeRecordKind::Array, "homology map entries require an array");
    guard.entries(entries.size()); Matrix result(rows, columns); std::set<std::pair<std::size_t,std::size_t>> seen;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        guard.step(); const auto entry = entries.at(i);
        require(entry.kind() == NativeRecordKind::Array && entry.size() == 3, "homology map entry requires row, column and coefficient");
        const auto row = count(entry.at(0)), col = count(entry.at(1));
        require(row < rows && col < columns && seen.emplace(row, col).second, "duplicate or out of range homology map entry");
        const auto q = rational(entry.at(2)); guard.value(q); result.set(row, col, q);
    }
    return result;
}
V integer(const rfc::numeric::Integer &n) { return V(N(n)); }
bool same_value(const NativeRecord &left, const NativeRecord &right,
    rfc::runtime::ExecutionLease lease) {
    std::vector<std::pair<NativeRecord,NativeRecord>> pending{{left,right}};
    while (!pending.empty()) {
        lease.check_cancelled(); const auto [a,b] = std::move(pending.back()); pending.pop_back();
        if (a.kind() != b.kind()) return false;
        switch (a.kind()) {
        case NativeRecordKind::Null: break;
        case NativeRecordKind::Boolean: if (a.boolean() != b.boolean()) return false; break;
        case NativeRecordKind::Integer: if (a.integer() != b.integer()) return false; break;
        case NativeRecordKind::Rational: if (a.rational() != b.rational()) return false; break;
        case NativeRecordKind::ExactNumber:
            if (a.rational() != b.rational() || a.approximate_view() != b.approximate_view()) return false;
            break;
        case NativeRecordKind::Decimal:
            if (a.rational() != b.rational()) return false;
            [[fallthrough]];
        case NativeRecordKind::Text:
        case NativeRecordKind::Bytes: {
            const auto x=a.bytes(), y=b.bytes(); if (!std::equal(x.begin(),x.end(),y.begin(),y.end())) return false;
            break;
        }
        case NativeRecordKind::Object: {
            const auto names=a.names(); if (names != b.names()) return false;
            for (const auto &name:names) pending.emplace_back(a.at(name),b.at(name));
            break;
        }
        case NativeRecordKind::Array:
            if (a.size() != b.size()) return false;
            for (std::size_t i=0;i<a.size();++i) pending.emplace_back(a.at(i),b.at(i));
            break;
        case NativeRecordKind::Tensor:
            if (!rfc::storage::temporal_endpoint_detail::tensor_equal(a.tensor(),b.tensor())) return false;
            break;
        }
    }
    return true;
}
V matrix(scientific_detail::ResultFields &fields, const Matrix &value) {
    const auto input = fields.axis(value.cols()), output = fields.axis(value.rows());
    auto result = fields.matrix(value, input, output);
    result.set("rows", V::count(value.rows())); result.set("columns", V::count(value.cols()));
    return result;
}
V matrix(scientific_detail::ResultFields &fields, const rfc::numeric::IntegerMatrix &value) {
    Matrix q(value.rows(), value.cols());
    for (const auto &[ij, n] : value.entries()) q.set(ij.first, ij.second, rfc::numeric::Rational(n));
    auto result = matrix(fields, q); result.set("coefficient_ring", V::string("Z")); return result;
}
V smith(scientific_detail::ResultFields &fields, const rfc::numeric::SmithNormalFormActions &value) {
    auto result = V::object(); result.set("diagonal", matrix(fields, value.diagonal)); result.set("rank", V::count(value.rank));
    auto factors = V::array(); for (const auto &n : value.invariant_factors) factors.push(integer(n)); result.set("invariant_factors", std::move(factors));
    for (const auto &[name, action] : std::array<std::pair<const char *,const rfc::numeric::IntegerUnimodularAction *>,2>{{{"left", &value.left}, {"right_inverse", &value.right_inverse}}}) {
        auto item = V::object(); item.set("dimension", V::count(action->dimension())); auto operations = V::array();
        for (const auto &op : action->operations()) { auto operation = V::object();
            operation.set("kind", V::count(static_cast<std::size_t>(op.kind))); operation.set("target", V::count(op.target)); operation.set("source", V::count(op.source));
            operation.set("multiple", integer(op.multiple)); operations.push(std::move(operation)); }
        item.set("operations", std::move(operations)); result.set(name, std::move(item));
    }
    return result;
}
NativeRecord publish(const NativeRecord &request, const std::vector<NativeRecord> &sources,
    const Provider &provider, std::size_t maximum, rfc::runtime::ExecutionLease lease) {
    scientific_detail::ResultFields fields(sources.front()); auto result = V::object();
    result.set("schema", V::string("gadget.homology.diagram.result.v1"));
    result.set("product", V::string(provider.declaration().product == Product::RationalZigzag ? "rational_zigzag" : "integral_modules"));
    result.set("grade", V::count(provider.declaration().grade.value)); result.set("request", fields.publication(request));
    auto vertices = V::array();
    for (std::size_t i = 0; i < sources.size(); ++i) { lease.check_cancelled(); auto vertex = V::object();
        vertex.set("source", fields.publication(sources[i])); vertex.set("source_namespace", V::count(sources[i].source_namespace()));
        vertex.set("sequence", V::count(sources[i].endpoint()->sequence)); auto carriers = V::array();
        for (auto c : provider.declaration().vertices[i].grade_carriers) carriers.push(V::count(c));
        vertex.set("grade_carriers", std::move(carriers));
        if (provider.declaration().product == Product::RationalZigzag) { const auto &h = provider.rational().homology()[i];
            vertex.set("cycles", matrix(fields, h.cycles)); vertex.set("boundaries", matrix(fields, h.boundaries));
            vertex.set("representatives", matrix(fields, h.representatives)); vertex.set("quotient_frame", matrix(fields, h.quotient_frame));
        } else { const auto &h = provider.integral().modules()[i];
            vertex.set("cycle_lattice", matrix(fields, h.cycle_basis)); vertex.set("boundary_coordinates", matrix(fields, h.boundary_coordinates));
            vertex.set("chain_smith", smith(fields, h.chain_smith)); vertex.set("boundary_smith", smith(fields, h.boundary_smith));
            vertex.set("free_rank", V::count(h.free_rank)); auto torsion = V::array(); for (const auto &n : h.torsion) torsion.push(integer(n)); vertex.set("torsion", std::move(torsion));
        }
        vertices.push(std::move(vertex));
    }
    result.set("vertices", std::move(vertices)); auto arrows = V::array();
    for (std::size_t i = 0; i < provider.declaration().arrows.size(); ++i) { lease.check_cancelled(); const auto &a = provider.declaration().arrows[i]; auto arrow = V::object();
        arrow.set("forward", V::boolean(a.forward)); auto maps = V::array(); for (const auto &m : a.maps) maps.push(matrix(fields, m)); arrow.set("chain_maps", std::move(maps));
        if (provider.declaration().product == Product::RationalZigzag) { const auto &h = provider.rational().induced()[i];
            arrow.set("cycle_coordinates", matrix(fields, h.cycle_coordinates)); arrow.set("boundary_coordinates", matrix(fields, h.boundary_coordinates)); arrow.set("homology_map", matrix(fields, h.homology_map));
        } else { const auto &h = provider.integral().induced()[i];
            arrow.set("cycle_map", matrix(fields, h.cycle_map)); arrow.set("smith_map", matrix(fields, h.smith_map)); arrow.set("relation_witness", matrix(fields, h.relation_witness)); }
        arrows.push(std::move(arrow));
    }
    result.set("arrows", std::move(arrows));
    if (provider.declaration().product == Product::RationalZigzag) { auto intervals = V::array();
        for (const auto &q : provider.rational().zigzag().intervals()) { auto interval = V::object(); interval.set("begin", V::count(q.begin)); interval.set("end", V::count(q.end));
            interval.set("multiplicity", V::count(q.generators.front().cols())); auto generators = V::array(); for (const auto &m : q.generators) generators.push(matrix(fields, m)); interval.set("generators", std::move(generators)); intervals.push(std::move(interval)); }
        result.set("intervals", std::move(intervals));
    }
    const auto key = fields.next_source_key();
    auto bytes = rfc::inference::encode_native_homology_diagram(provider, maximum, [&lease] {lease.check_cancelled();});
    result.set("provider", V(N(N::ByteSpan{key, 0, bytes.size(), false}))); fields.retain_sources({RawByteSource::own(key, std::move(bytes))});
    result.set("arithmetic", V::string(provider.declaration().product == Product::RationalZigzag ? "Q" : "Z"));
    return fields.finish(std::move(result));
}
} // namespace

namespace {
NativeRecord compute(const NativeRecord &request, const std::vector<NativeRecord> &sources,
    rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    members(request, {"schema", "product", "grade", "vertices", "arrows", "evidence"});
    require(request.at("schema").text() == "gadget.homology.diagram.v1", "homology diagram schema");
    require(!request.at("evidence").text().empty(), "homology diagram requires the map evidence declaration");
    const auto product = request.at("product").text(); require(product == "rational_zigzag" || product == "integral_modules", "homology diagram product");
    const auto grade = count(request.at("grade")); require(grade <= std::numeric_limits<std::uint16_t>::max(), "homology diagram grade representation");
    const auto vertices = request.at("vertices"), arrows = request.at("arrows");
    require(vertices.kind() == NativeRecordKind::Array && arrows.kind() == NativeRecordKind::Array && !sources.empty() && sources.size() == vertices.size(), "homology diagram source family");
    require(arrows.size() == sources.size()-1, "homology diagram consecutive arrow family");
    const auto &policy = lease.policy().limits;
    rfc::numeric::ExactActionGuard guard({std::numeric_limits<std::size_t>::max(), std::min(policy.rational_slots,policy.sparse_entries),policy.scalar_work,policy.integer_bits}, [&lease] {lease.check_cancelled();});
    rfc::inference::NativeHomologyDiagramRequest d; d.product = product == "rational_zigzag" ? Product::RationalZigzag : Product::IntegralModules; d.grade = rfc::core::Grade{static_cast<std::uint16_t>(grade)};
    for (std::size_t i = 0; i < vertices.size(); ++i) { members(vertices.at(i), {"source", "carriers"}); const auto ids = vertices.at(i).at("carriers"); require(ids.kind() == NativeRecordKind::Array, "homology diagram carrier family");
        require(!vertices.at(i).at("source").text().empty(), "homology diagram source name absent");
        std::vector<std::uint64_t> carriers; for (std::size_t k = 0; k < ids.size(); ++k) {guard.step();carriers.push_back(count(ids.at(k)));}
        d.vertices.push_back({sources[i].access(rfc::storage::NativeAccessionRequest::sparse_chain(rfc::core::Grade{0},{},{})),std::move(carriers)});
    }
    for (std::size_t i = 0; i < arrows.size(); ++i) { const auto a = arrows.at(i); members(a, {"forward", "maps"}); rfc::inference::NativeHomologyDiagramArrow arrow; arrow.forward = a.at("forward").boolean();
        const auto maps = a.at("maps"); require(maps.kind() == NativeRecordKind::Array, "homology diagram maps require an array");
        for (std::size_t k = 0; k < maps.size(); ++k) arrow.maps.push_back(map(maps.at(k),guard));
        d.arrows.push_back(std::move(arrow)); }
    const auto provider = Provider::prepare(std::move(d),lease);
    return publish(request,sources,provider,source_byte_limit(),lease);
}
} // namespace
NativeRecord homology_diagram_record(const NativeRecord &request, const std::vector<NativeRecord> &sources,
    std::shared_ptr<CpuResources> resources, std::optional<rfc::runtime::ExecutionLease> supplied) {
    const auto lease = supplied.value_or(scientific_detail::source_execution_lease("gadget-homology-diagram")); lease.check_cancelled();
    if (!resources) resources = CpuResources::local();
    auto worker = resources->reserve(1);
    return compute(request,sources,lease);
}
NativeRecord run_homology_diagram(const fs::path &specification, std::shared_ptr<CpuResources> resources) {
    const auto request = read_record(specification).native(); std::vector<NativeRecord> sources;
    const auto vertices = request.at("vertices"); require(vertices.kind() == NativeRecordKind::Array, "homology diagram vertices require an array");
    for (std::size_t i = 0; i < vertices.size(); ++i) sources.push_back(NativeRecord::read(specification.parent_path()/vertices.at(i).at("source").text()));
    return homology_diagram_record(request,sources,std::move(resources));
}
NativeRecord replay_homology_diagram(const NativeRecord &record, std::shared_ptr<CpuResources> resources,
    std::optional<rfc::runtime::ExecutionLease> supplied) {
    const auto lease = supplied.value_or(scientific_detail::source_execution_lease("gadget-homology-diagram-replay"));lease.check_cancelled();
    const auto view = record.at("view"); require(view.at("schema").text() == "gadget.homology.diagram.result.v1", "homology diagram result schema");
    if (!resources) resources = CpuResources::local();
    auto worker = resources->reserve(1);
    const auto request = NativeRecord::decode(view.at("request").at("publication").bytes()); std::vector<NativeRecord> sources;
    const auto vertices = view.at("vertices");
    for (std::size_t i = 0; i < vertices.size(); ++i) sources.push_back(NativeRecord::decode(vertices.at(i).at("source").at("publication").bytes()));
    const auto result = compute(request,sources,lease);
    const auto original = view.at("provider").bytes(), replay = result.at("view").at("provider").bytes();
    require(std::equal(original.begin(),original.end(),replay.begin(),replay.end()), "homology replay provider declaration or certificate differs");
    require(same_value(view,result.at("view"),lease), "homology replay typed result declaration differs");
    const auto before = record.endpoint(), after = result.endpoint();
    require(before->carriers == after->carriers && before->tensors.size() == after->tensors.size(), "homology replay result carriers differ");
    for (std::size_t i = 0; i < before->tensors.size(); ++i) { lease.check_cancelled(); const auto &a = before->tensors[i], &b = after->tensors[i];
        require(a.key == b.key && a.role == b.role && a.owner == b.owner && a.input == b.input && a.output == b.output && a.interpretation == b.interpretation &&
            rfc::storage::temporal_endpoint_detail::tensor_equal(a.tensor,b.tensor), "homology replay native result field differs"); }
    return result;
}
Json homology_diagram_view(const NativeRecord &record) { return scientific_detail::result_fields_view(record); }
} // namespace gadget
