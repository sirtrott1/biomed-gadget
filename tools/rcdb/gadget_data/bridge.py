"""Transport selected native coefficients to GADGET without numerical JSON."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import hashlib
import subprocess


def export_gadget(selection, destination, *, importer, namespace, max_bytes, working_bytes=None, source_rows=4096):
    """Export one selected field as an owning native GADGET publication."""
    import rexgraph_cpp as cpp
    from rexgraph.column import ColumnDeclaration, declaration_of, exact_slot_coefficients
    from ._records import natural
    from .identity import runtime_identity
    natural(namespace, "native source namespace", positive=True)
    natural(max_bytes, "native export byte budget", positive=True)
    natural(source_rows, "native source row block", positive=True)
    if namespace >= 2**64 - 1:
        raise ValueError("native source namespace leaves no output identity")
    working_bytes = max_bytes if working_bytes is None else natural(working_bytes, "native export working bytes", positive=True)
    # This admission covers object slots during exact transport, independently
    # of the much smaller original typed matrix or encoded native archive.
    if selection.values.size * 1024 > working_bytes:
        raise MemoryError("native field export working admission exceeded; select less data or increase working_bytes")
    destination = Path(destination).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if selection.values.dtype.hasobject:
        bits = sum(int(value.numerator).bit_length() + int(value.denominator).bit_length()
                   for value, active in zip(selection.values.flat, selection.observed.flat) if active)
        if selection.values.size * 1024 + (bits + 1) // 2 > working_bytes:
            raise MemoryError("native field integer transport exceeds its working byte admission")
    source = selection.native()
    exact = source.get_metadata(1, 0, "exact_values")
    values = exact.values()
    relations = source.relations
    roles = tuple(dict.fromkeys(key[0] for key in relations.vertices.ids))
    role_positions = {role: i for i, role in enumerate(roles)}
    arities = relations.support_ptr[1:] - relations.support_ptr[:-1]
    types = ({"name": "original-measurement-row", "min_arity": int(arities.min()), "max_arity": int(arities.max())},)

    def component(start, length):
        begin, end = int(relations.support_ptr[start]), int(relations.support_ptr[start + length])
        declared = declaration_of(source)
        if declared is not None:
            declared = ColumnDeclaration(declared.head_slot[start:start + length], declared.share_num[begin:end], declared.share_den[begin:end])
        coefficients = exact_slot_coefficients(relations.support_ptr[start:start + length + 1] - begin,
                                               relations.support_idx[begin:end], declared)
        rows, used = [], set()
        for row in range(start, start + length):
            occurrences = []
            for slot in range(int(relations.support_ptr[row]), int(relations.support_ptr[row + 1])):
                vertex, value = int(relations.support_idx[slot]), coefficients[slot - begin]
                used.add(vertex + 1)
                occurrences.append({"key": slot + 1, "lower": vertex + 1, "exists": True,
                    "negative": value < 0, "magnitude": abs(value),
                    "aspects": (role_positions[relations.vertices.ids[vertex][0]],)})
            rows.append({"grade": 1, "key": row + 1, "type": 0, "occurrences": tuple(occurrences)})
        return cpp.NativeSource.from_declaration({
            "profile": "declared-role-source-v1", "namespace": namespace, "semantic_scale": 0,
            "max_grade": 1, "role_aspects": roles, "witnesses": tuple(sorted(used)),
            "types": types, "relations": tuple(rows), "carriers": (), "fields": (), "metrics": ()}, max_bytes=max_bytes).archive
    present_count = int(exact.presence.sum())
    bounds = (values.size, 2 * values.size + 2, 4 * values.size + 16, working_bytes * 8)
    tensor = cpp.encode_stalk_column_tensor(values.size,
        ((i, value) for i, value in enumerate(values.flat) if exact.presence.flat[i] and value), limits=bounds)
    support, support_kind = None, None
    if present_count != values.size:
        support_kind = "observed" if present_count <= values.size - present_count else "missing"
        support = cpp.encode_stalk_column_tensor(values.size,
            ((i, 1) for i, active in enumerate(exact.presence.flat) if bool(active) == (support_kind == "observed")), limits=bounds)
    if len(tensor) > max_bytes or support is not None and len(support) > max_bytes:
        raise MemoryError("native field tensor exceeds its encoded byte admission")
    metadata = {"schema": "gadget.native-data.field.v1", "dataset": selection.dataset,
        "exporter": runtime_identity(),
        "native_source": {"namespace": str(namespace),
            "values_sha256": hashlib.sha256(tensor).hexdigest(),
            "support_sha256": None if support is None else hashlib.sha256(support).hexdigest(),
            "support_kind": support_kind},
        "field": selection.field, "field_id": selection.specification["identity"],
        "shape": selection.values.shape, "rows": selection.rows, "columns": selection.columns,
        "axes": selection.specification["axes"], "row_labels": selection.row_labels,
        "column_labels": selection.column_labels,
        "vertex_map": tuple((i + 1, key) for i, key in enumerate(relations.vertices.ids)),
        "row_map": tuple((i + 1, row) for i, row in enumerate(selection.rows)),
        "parent": selection.parent, "dependencies": selection.dependencies,
        "quantity": selection.specification["quantity"], "units": selection.specification["units"],
        "number_rule": selection.specification["number_rule"], "observed_count": present_count}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="gadget-field-export-", dir=destination.parent) as temporary:
        temporary = Path(temporary)
        parts, offset = [], 0
        while offset < source.nE:
            length = min(source_rows, source.nE - offset)
            while True:
                try:
                    archive = component(offset, length)
                    break
                except ValueError:
                    if length == 1:
                        raise
                    length = max(1, length // 2)
            name = f"source-{len(parts)}.rfcxs"
            (temporary / name).write_bytes(archive)
            parts.append({"file": name, "offset": offset, "count": length, "sha256": hashlib.sha256(archive).hexdigest()})
            offset += length
        if len(parts) == 1:
            metadata["native_source"]["sha256"] = parts[0]["sha256"]
            source_input = ["--input", str(temporary / parts[0]["file"])]
        else:
            description = json.dumps({"schema": "gadget.native-data.components.v1", "parts": parts}) + "\n"
            metadata["native_source"]["components_sha256"] = hashlib.sha256(description.encode()).hexdigest()
            (temporary / "components.json").write_text(description)
            source_input = ["--components", str(temporary / "components.json")]
        (temporary / "values.tensor").write_bytes(tensor)
        if support is not None:
            (temporary / "support.tensor").write_bytes(support)
        # This is an identity and coordinate declaration, never field coefficients.
        try:
            control = json.dumps(metadata)
        except TypeError:
            from rexgraph.value_codec import pack_value
            identities = {key: metadata[key] for key in ("row_labels", "column_labels", "vertex_map")}
            metadata["identities_native"] = {"codec": "rexgraph.value", "payload_hex": pack_value(identities).hex()}
            metadata["row_labels"] = tuple({"native_identity": "row_labels", "index": i} for i in range(len(selection.rows)))
            metadata["column_labels"] = tuple({"native_identity": "column_labels", "index": i} for i in range(len(selection.columns)))
            metadata["vertex_map"] = tuple((i + 1, {"native_identity": "vertex_map", "index": i}) for i in range(source.nV))
            control = json.dumps(metadata)
        if len(control.encode()) > max_bytes:
            raise MemoryError("native coordinate declaration exceeds its byte admission")
        (temporary / "coordinates.json").write_text(control + "\n")
        command = [str(importer), *source_input,
            "--coordinates", str(temporary / "coordinates.json"), "--namespace", str(namespace + 1),
            "--values", str(temporary / "values.tensor"),
            "--max-bytes", str(max_bytes),
            "--out", str(destination)]
        if support is not None:
            command.extend(("--support", str(temporary / "support.tensor")))
        subprocess.run(command, check=True, stdout=subprocess.PIPE)
    return destination
