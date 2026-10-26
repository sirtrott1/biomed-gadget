"""Retained native programs and results with explicit versioned owners."""
from dataclasses import dataclass, field, replace

from ._records import declaration, graph, natural, text


def check_link(link):
    """Check an address declaration without reading any scientific coefficients."""
    parent, record = link["dataset"], link["record"]
    for name in ("dataset", "store_id", "state_digest"):
        text(parent[name], "linked " + name)
    for name in ("id", "store_id", "state_digest"):
        text(record[name], "linked record " + name)
    natural(record["version"], "linked version", positive=True)
    natural(record["native_bytes"], "linked bytes", positive=True)
    if parent["id"] != "dataset" or parent["version"] != 1 or parent["store_id"] != record["store_id"]:
        raise ValueError("linked record and sealed dataset must share a store")


def record_link(store, reference, *, owners=()):
    """Admit an existing record while retaining its dataset manifest identity."""
    from copy import deepcopy
    matches = [owner for owner in owners if isinstance(owner, NativeOwner)
               and owner.link["dataset"] == store.reference and owner.link["record"] == reference]
    if len(matches) > 1:
        raise ValueError("one actual owner binding per linked source required")
    if matches:
        store.check_reference(reference)
        matches[0].value
        return deepcopy(matches[0].link)
    return _record_link(store, reference, store.read(reference))


def _record_link(store, reference, owner):
    from copy import deepcopy
    digest = reference["state_digest"]
    value = declaration(owner)
    if value["kind"] == "native" and value["native_type"] == "record":
        digest = owner.get_metadata(1, 0, "native_digest")
    if value["kind"] == "native" and value["native_type"] == "publication":
        digest = owner.get_metadata(1, 0, "native_value")["reference"]["state_digest"]
    link = {"dataset": dict(store.reference), "record": deepcopy(reference), "owner_digest": digest}
    check_link(link)
    return link


@dataclass(frozen=True)
class NativeOwner:
    """An explicitly retained actual owner, with its original store address.

    This is a caller owned lifetime, not a global cache or a process memory cap.
    Current record metadata and the actual source are rechecked before each use.
    """
    link: dict
    reference: object
    dependencies: tuple = ()
    _seal: bytes = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        from copy import deepcopy
        from rexgraph.value_codec import pack_value
        check_link(self.link)
        object.__setattr__(self, "link", deepcopy(self.link))
        object.__setattr__(self, "dependencies", tuple(deepcopy(self.dependencies)))
        object.__setattr__(self, "_seal", pack_value((self.link, self.dependencies)))

    @property
    def value(self):
        from rexgraph.value_codec import pack_value
        if pack_value((self.link, self.dependencies)) != self._seal:
            raise ValueError("retained native owner declaration changed")
        self.reference.check()
        if self.reference.state_digest != self.link["owner_digest"] or self.reference.source is None:
            raise ValueError("retained native owner differs from its original source")
        return self.reference.source


def native_owner(store, reference, *, dependencies=(), owners=()):
    """Restore one original carrier for reuse by explicitly related calculations."""
    from rexgraph import RexGraph
    from rexgraph.tensor_field import FieldSource
    from rexgraph.source_binding import CppSourceReference
    record = store.read(reference)
    link = _record_link(store, reference, record)
    value = declaration(record)
    owner = (restore_record(record, datasets=(store, *dependencies), limits=store.limits, owners=owners)
             if value["kind"] == "native" else record)
    if isinstance(owner, RexGraph):
        source = FieldSource(owner, reference["id"], reference["version"], link["owner_digest"])
    else:
        source = CppSourceReference(owner, reference["id"], reference["version"])
    result = NativeOwner(link, source, value["dependencies"] if value["kind"] == "native" else ())
    result.value
    return result


def resolve_link(link, stores, owners=()):
    """Restore a pinned owner through the caller's admitted dataset stores.

    Dataset, record and owner identities must match. Paths contained in stored
    declarations do not authorize opening additional stores.
    """
    check_link(link)
    store = stores[link["dataset"]["dataset"]]
    if store.reference != link["dataset"]:
        raise ValueError("linked dataset manifest differs from its retained version")
    matches = [owner for owner in owners if isinstance(owner, NativeOwner) and owner.link == link]
    if len(matches) > 1:
        raise ValueError("one actual owner binding per linked source required")
    if matches:
        store.check_reference(link["record"])
        for dependency in matches[0].dependencies:
            resolve_link(dependency, stores, owners)
        return matches[0].value
    owner = store.read(link["record"])
    value = declaration(owner)
    digest = link["record"]["state_digest"]
    if value["kind"] == "native" and value["native_type"] == "record":
        from rexgraph.object_identity import object_digest
        owner = restore_record(owner, datasets=tuple(stores.values()), limits=store.limits, owners=owners)
        digest = object_digest(owner)
    if value["kind"] == "native" and value["native_type"] == "publication":
        from rexgraph.cpp_model import restore_cpp_source
        from rexgraph.source_binding import source_from_record
        attachment = owner.get_metadata(1, 0, "native_value")
        ref = source_from_record(attachment["reference"], allow_cpp=True)
        size = attachment["archive"].size
        store.limits.admit(size, "decode_bytes")
        import numpy as np
        owner = restore_cpp_source(attachment, ref, max_bytes=store.limits.decode_bytes,
            max_records=size, max_dimension=size, max_logical_size=int(np.iinfo(np.uintp).max))
        digest = ref.state_digest
    if digest != link["owner_digest"]:
        raise ValueError("linked native owner identity differs")
    return owner


def _references(value):
    from rexgraph import RexGraph
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    if isinstance(value, (TensorField, SectionFamily)):
        return tuple(ref for ref in (value.source, *value.dependencies) if ref is not None)
    if isinstance(value, RexGraph):
        return ()
    from rcql.source_context import field_references
    return field_references(value)


def require_owner_links(references, links):
    """Require one original stored owner for every actual program dependency."""
    from rexgraph.source_binding import CppSourceReference
    for ref in references:
        ref.check()
        matches = [link for link in links if link["owner_digest"] == ref.state_digest and (
            isinstance(ref, CppSourceReference) or ref.record_id is None or
            (ref.record_id, ref.version) == (link["record"]["id"], link["record"]["version"]))]
        if len(matches) != 1:
            raise ValueError("each program dependency needs one original versioned owner link")


def retain(writer, name, value, *, dependencies=()):
    """Store actual provider values, with each contributing owner referenced once."""
    from rcql import Program
    from rexgraph import RexGraph
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    from rexgraph.source_binding import CppSourceReference
    from rexgraph.model_codec import MODEL_VALUES
    text(name, "native object name")
    if isinstance(value, (TensorField, SectionFamily)) and value.source is None:
        raise ValueError("retained calculations require a native source owner")
    links = writer._links(dependencies)
    for ref in _references(value):
        ref.check()
        matches = [link for link in links if link["owner_digest"] == ref.state_digest and (
            isinstance(ref, CppSourceReference) or
            (ref.record_id, ref.version) == (link["record"]["id"], link["record"]["version"]))]
        if len(matches) != 1:
            raise ValueError("each native source needs exactly one matching versioned owner link")
    if isinstance(value, Program):
        owner, kind = value.to_record(), "program"
    elif isinstance(value, RexGraph):
        from rexgraph.state import to_state
        from rexgraph.io.safetensors_bridge import state_to_safetensors_bytes
        from rexgraph.object_identity import object_digest
        digest = object_digest(value)
        raw = state_to_safetensors_bytes(to_state(value))
        writer.limits.admit(len(raw), "record_bytes")
        writer.limits.admit(len(raw), "decode_bytes")
        owner, kind = graph("native", ((writer.identity, name),)), "record"
        owner.attach_metadata(1, 0, "native_value", raw)
        owner.attach_metadata(1, 0, "native_digest", digest)
    elif isinstance(value, (TensorField, SectionFamily)):
        value.check_state()
        owner, kind = graph("native", ((writer.identity, name),)), type(value).__name__
        def detached(field):
            return replace(field, source=field.source.detached(),
                           dependencies=tuple(ref.detached() for ref in field.dependencies))
        if isinstance(value, TensorField):
            value = detached(value)
        else:
            value = replace(value, particular=detached(value.particular), rhs=detached(value.rhs))
        owner.attach_metadata(1, 0, "native_value", value)
    elif isinstance(value, MODEL_VALUES):
        value.check_state()
        owner, kind = graph("native", ((writer.identity, name),)), "model"
        owner.attach_metadata(1, 0, "native_value", value)
    elif type(value).__module__ == "rexgraph_cpp._bridge" and type(value).__name__ == "NativeSource":
        # The optional binding is imported only when a C++ publication is stored.
        from rexgraph_cpp import NativeSource
        if not isinstance(value, NativeSource):
            raise TypeError("native field, family, program, context or owning publication required")
        from rexgraph.cpp_model import cpp_source_record
        writer.limits.admit(len(value.archive), "decode_bytes")
        writer.limits.admit(len(value.archive), "record_bytes")
        owner, kind = graph("native", ((writer.identity, name),)), "publication"
        owner.attach_metadata(1, 0, "native_value", cpp_source_record(value, CppSourceReference(value)))
    else:
        from rcql import Result
        raw = Result((value,)).to_bytes(max_bytes=min(writer.limits.record_bytes, writer.limits.decode_bytes))
        owner, kind = graph("native", ((writer.identity, name),)), "value"
        owner.attach_metadata(1, 0, "native_value", raw)
    entry = writer._put("objects/" + name, owner,
        {"kind": "native", "native_type": kind, "dependencies": links})
    writer.objects.append(entry)
    return entry


def restore(store, reference, *, dependencies=(), owners=()):
    """Rebind a retained result to its original source and contributor versions."""
    return restore_record(store.read(reference), datasets=(store, *dependencies), limits=store.limits, owners=owners)


def restore_record(owner, *, datasets=(), limits=None, owners=()):
    """Restore an RCDB or RCBD native record using admitted dependency stores."""
    from rcql import Program
    from rexgraph.section_calculus import SectionFamily
    from ._records import Limits
    limits = limits or Limits()
    stores = {item.identity: item for item in datasets}
    owners = tuple(owners)
    metadata = declaration(owner)
    if metadata["kind"] != "native":
        raise TypeError("retained native object required")
    kind = metadata["native_type"]
    if kind in ("value", "model"):
        from rcql import Result
        from rcql.result_codec import bind_sources
        references = {}
        for link in metadata["dependencies"]:
            current = resolve_link(link, stores, owners)
            references[(link["owner_digest"], link["record"]["id"], link["record"]["version"])] = current
        value = owner.get_metadata(1, 0, "native_value")
        if kind == "value":
            limits.admit(len(value), "decode_bytes")
            result = Result.from_bytes(value, max_bytes=limits.decode_bytes)
            if len(result.values) != 1:
                raise ValueError("one complete retained RCQL value required")
            value = result.values[0]
        from rexgraph.source_binding import CppSourceReference
        for ref in _references(value):
            if isinstance(ref, CppSourceReference):
                matches = [link for link in metadata["dependencies"] if link["owner_digest"] == ref.state_digest]
                if len(matches) != 1:
                    raise ValueError("native model source needs one original owner")
                link = matches[0]
                references[(ref.state_digest, ref.record_id, ref.version)] = references[
                    (link["owner_digest"], link["record"]["id"], link["record"]["version"])]
        return bind_sources(value, references, strict=True)
    if kind == "record":
        from rexgraph.state import from_state
        from rexgraph.io.safetensors_bridge import state_from_safetensors_bytes
        from rexgraph.object_identity import object_digest
        for link in metadata["dependencies"]:
            resolve_link(link, stores, owners)
        raw = owner.get_metadata(1, 0, "native_value")
        limits.admit(len(raw), "decode_bytes")
        value = from_state(state_from_safetensors_bytes(raw))
        if object_digest(value) != owner.get_metadata(1, 0, "native_digest"):
            raise ValueError("retained native record identity differs")
        return value
    if kind == "publication":
        from rexgraph.cpp_model import restore_cpp_source
        from rexgraph.source_binding import source_from_record
        import numpy as np
        attachment = owner.get_metadata(1, 0, "native_value")
        ref = source_from_record(attachment["reference"], allow_cpp=True)
        size = attachment["archive"].size
        limits.admit(size, "decode_bytes")
        return restore_cpp_source(attachment, ref, max_bytes=limits.decode_bytes,
            max_records=size, max_dimension=size, max_logical_size=int(np.iinfo(np.uintp).max))
    if kind == "context":
        for link in metadata["dependencies"]:
            resolve_link(link, stores, owners)
        return owner
    if kind == "program":
        # Resolve every declared dependency before handing the program to RCQL.
        for link in metadata["dependencies"]:
            resolve_link(link, stores, owners)
        return Program.from_record(owner)
    value = owner.get_metadata(1, 0, "native_value")
    original_owners, references = {}, {}
    for link in metadata["dependencies"]:
        current = resolve_link(link, stores, owners)
        original_owners[link["owner_digest"]] = current
        references[(link["owner_digest"], link["record"]["id"], link["record"]["version"])] = current
    bound = []
    from rexgraph.source_binding import CppSourceReference
    for ref in _references(value):
        current = original_owners[ref.state_digest] if isinstance(ref, CppSourceReference) else references[
            (ref.state_digest, ref.record_id, ref.version)]
        bound.append(ref.bind(current))
    if isinstance(value, SectionFamily):
        return value.bind(bound[0], bound[1:])
    return replace(value, source=None if value.source is None else bound[0],
                   dependencies=tuple(bound[int(value.source is not None):]))


def field_program(reference, *, relation_family=False, source=None):
    """Build an RCQL program that retrieves the recorded field version.

    With relation_family=True, append the provider's complete exact family
    construction. Execution requires the original dataset and owner bindings.
    An explicit source alias uses its pinned native owner directly.
    """
    import json
    from rcql import OutputRef, Program, ProgramInput, ProgramStep, parse
    identifier = json.dumps(reference["id"])
    version = natural(reference["version"], "field version", positive=True)
    digest = json.dumps(reference["state_digest"])
    dependencies = ()
    selected = f"RCDB_VERSION($db,{identifier},{version})"
    if source is not None:
        if type(source) is not str or not source.isascii() or not source.isidentifier():
            raise ValueError("native program source must be an ASCII identifier")
        from rexgraph.tensor_field import FieldSource
        selected = "$" + source
        dependencies = (FieldSource(None, reference["id"], version, reference["state_digest"]),)
    field = ProgramStep("field", parse(f"FROM {selected} "
        f"LET owner=GADGET_RECORD_OWNER({digest}) RETURN GADGET_DATA_FIELD($columns)"))
    steps, outputs = (field,), (("field", OutputRef("field")),)
    if relation_family:
        family = ProgramStep("family", parse(f"FROM {selected} "
            "RETURN GADGET_RELATION_FAMILY($field)"), (("field", OutputRef("field")),))
        steps += (family,)
        outputs += (("family", OutputRef("family")),)
    return Program("gadget_native_field", steps, (ProgramInput("columns"),), outputs, dependencies=dependencies)
