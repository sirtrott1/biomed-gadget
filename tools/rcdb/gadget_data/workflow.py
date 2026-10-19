"""Retained native programs and results with explicit versioned owners."""
from dataclasses import replace

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


def record_link(store, reference):
    """Admit an existing record while retaining its dataset manifest identity."""
    from copy import deepcopy
    owner = store.read(reference)
    digest = reference["state_digest"]
    value = declaration(owner)
    if value["kind"] == "native" and value["native_type"] == "publication":
        digest = owner.get_metadata(1, 0, "native_value")["reference"]["state_digest"]
    link = {"dataset": dict(store.reference), "record": deepcopy(reference), "owner_digest": digest}
    check_link(link)
    return link


def resolve_link(link, stores):
    """Restore a pinned owner through the caller's admitted dataset stores.

    Dataset, record and owner identities must match. Paths contained in stored
    declarations do not authorize opening additional stores.
    """
    check_link(link)
    store = stores[link["dataset"]["dataset"]]
    if store.reference != link["dataset"]:
        raise ValueError("linked dataset manifest differs from its retained version")
    owner = store.read(link["record"])
    value = declaration(owner)
    digest = link["record"]["state_digest"]
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
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    if isinstance(value, (TensorField, SectionFamily)):
        return tuple(ref for ref in (value.source, *value.dependencies) if ref is not None)
    return ()


def retain(writer, name, value, *, dependencies=()):
    """Store actual provider values, with each contributing owner referenced once."""
    from rcql import Program
    from rexgraph import RexGraph
    from rexgraph.section_calculus import SectionFamily
    from rexgraph.tensor_field import TensorField
    from rexgraph.source_binding import CppSourceReference
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
        owner, kind = value.copy(), "context"
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
    else:
        # The optional binding is imported only when a C++ publication is stored.
        from rexgraph_cpp import NativeSource
        if not isinstance(value, NativeSource):
            raise TypeError("native field, family, program, context or owning publication required")
        from rexgraph.cpp_model import cpp_source_record
        writer.limits.admit(len(value.archive), "decode_bytes")
        writer.limits.admit(len(value.archive), "record_bytes")
        owner, kind = graph("native", ((writer.identity, name),)), "publication"
        owner.attach_metadata(1, 0, "native_value", cpp_source_record(value, CppSourceReference(value)))
    entry = writer._put("objects/" + name, owner,
        {"kind": "native", "native_type": kind, "dependencies": links})
    writer.objects.append(entry)
    return entry


def restore(store, reference, *, dependencies=()):
    """Rebind a retained result to its original source and contributor versions."""
    return restore_record(store.read(reference), datasets=(store, *dependencies), limits=store.limits)


def restore_record(owner, *, datasets=(), limits=None):
    """Restore an RCDB or RCBD native record using admitted dependency stores."""
    from rcql import Program
    from rexgraph.section_calculus import SectionFamily
    from ._records import Limits
    limits = limits or Limits()
    stores = {item.identity: item for item in datasets}
    metadata = declaration(owner)
    if metadata["kind"] != "native":
        raise TypeError("retained native object required")
    kind = metadata["native_type"]
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
            resolve_link(link, stores)
        return owner
    if kind == "program":
        # Resolve every declared dependency before handing the program to RCQL.
        for link in metadata["dependencies"]:
            resolve_link(link, stores)
        return Program.from_record(owner)
    value = owner.get_metadata(1, 0, "native_value")
    owners, references = {}, {}
    for link in metadata["dependencies"]:
        current = resolve_link(link, stores)
        owners[link["owner_digest"]] = current
        references[(link["owner_digest"], link["record"]["id"], link["record"]["version"])] = current
    bound = []
    from rexgraph.source_binding import CppSourceReference
    for ref in _references(value):
        current = owners[ref.state_digest] if isinstance(ref, CppSourceReference) else references[
            (ref.state_digest, ref.record_id, ref.version)]
        bound.append(ref.bind(current))
    if isinstance(value, SectionFamily):
        return value.bind(bound[0], bound[1:])
    return replace(value, source=None if value.source is None else bound[0],
                   dependencies=tuple(bound[int(value.source is not None):]))


def field_program(reference, *, relation_family=False):
    """Build an RCQL program that retrieves the recorded field version.

    With relation_family=True, append the provider's complete exact family
    construction. Execution requires the original dataset and owner bindings.
    """
    import json
    from rcql import OutputRef, Program, ProgramInput, ProgramStep, parse
    identifier = json.dumps(reference["id"])
    version = natural(reference["version"], "field version", positive=True)
    digest = json.dumps(reference["state_digest"])
    field = ProgramStep("field", parse(f"FROM RCDB_VERSION($db,{identifier},{version}) "
        f"LET owner=GADGET_RECORD_OWNER({digest}) RETURN GADGET_DATA_FIELD($columns)"))
    steps, outputs = (field,), (("field", OutputRef("field")),)
    if relation_family:
        family = ProgramStep("family", parse(f"FROM RCDB_VERSION($db,{identifier},{version}) "
            "RETURN GADGET_RELATION_FAMILY($field)"), (("field", OutputRef("field")),))
        steps += (family,)
        outputs += (("family", OutputRef("family")),)
    return Program("gadget_native_field", steps, (ProgramInput("columns"),), outputs)
