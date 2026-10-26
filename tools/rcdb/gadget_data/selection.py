"""Selected original coefficients, support and retained source maps."""
from dataclasses import dataclass, field
from math import prod

import numpy as np
import hashlib

from ._records import SCHEMA, declare, declaration, frozen, graph


def coordinates(values, extent, *, max_bytes=None):
    if max_bytes is not None:
        from itertools import islice
        limit = max_bytes // 64
        if hasattr(values, "__len__") and len(values) > limit:
            raise MemoryError("coordinate selection metadata exceeds its execution byte admission")
        result = tuple(islice(values, limit + 1))
        if len(result) > limit:
            raise MemoryError("coordinate selection metadata exceeds its execution byte admission")
    else:
        from itertools import islice
        if hasattr(values, "__len__") and len(values) > extent:
            raise ValueError("distinct coordinates cannot exceed the declared axis extent")
        result = tuple(islice(values, extent + 1))
    if (not result or any(type(i) is not int or not 0 <= i < extent for i in result)
            or len(set(result)) != len(result)):
        raise ValueError("distinct explicit original coordinates within the axis required")
    return result


def exact_slice(array, rows, columns):
    """Select native exact tensors, preserving absent slots and bigint escapes."""
    from rexgraph.exact_array import ExactArray
    row_positions = {row: i for i, row in enumerate(rows)}
    column_positions = {column: i for i, column in enumerate(columns)}
    width, selected_width = array.shape[1], len(columns)
    escapes = []
    for i, numerator, denominator in array.bigint:
        row, column = divmod(i, width)
        if row in row_positions and column in column_positions:
            escapes.append((row_positions[row] * selected_width + column_positions[column], numerator, denominator))
    index = np.ix_(rows, columns)
    return ExactArray(array.numerator[index], array.denominator[index], array.presence[index],
                      array.kind[index], tuple(sorted(escapes)))


def exact_take(array, positions):
    from rexgraph.exact_array import ExactArray
    selected = {int(position): i for i, position in enumerate(positions)}
    escapes = tuple(sorted((selected[i], n, d) for i, n, d in array.bigint if i in selected))
    return ExactArray(array.numerator[positions], array.denominator[positions], array.presence[positions],
                      array.kind[positions], escapes)


def block_values(owner, rows, columns):
    value = declaration(owner)
    if value["layout"] == "exact":
        array = exact_slice(owner.get_metadata(1, 0, "exact_values"), rows, columns)
        return array.values(), array.presence
    if value["layout"] != "dense":
        raise ValueError("sparse occurrence fields require occurrence selection")
    shape = tuple(value["shape"])
    raw = owner.get_metadata(1, 0, "values")
    packed = owner.get_metadata(1, 0, "observed")
    if (raw.shape != shape or packed.dtype != np.dtype("uint8") or
            len(packed) != (prod(shape) + 7) // 8):
        raise ValueError("native field coefficient or support extent differs")
    if prod(shape) % 8 and int(packed[-1]) >> (prod(shape) % 8):
        raise ValueError("native observed support has nonzero padding")
    present = np.empty((len(rows), len(columns)), dtype=bool)
    slots = 64 * 1024 // np.dtype(np.int64).itemsize
    for column_start in range(0, len(columns), slots):
        selected_columns = np.asarray(columns[column_start:column_start + slots], dtype=np.int64)
        row_count = max(1, slots // len(selected_columns))
        for row_start in range(0, len(rows), row_count):
            selected_rows = np.asarray(rows[row_start:row_start + row_count], dtype=np.int64)
            flat = selected_rows[:, None] * shape[1] + selected_columns
            present[row_start:row_start + len(selected_rows),
                    column_start:column_start + len(selected_columns)] = (
                        (packed[flat // 8] >> (flat % 8)) & 1).astype(bool)
    return raw[np.ix_(rows, columns)], present


def _check(owner, block, field):
    value = declaration(owner)
    for key in ("kind", "field_id", "layout", "offset", "shape", "axes", "number_rule"):
        if value[key] != block[key]:
            raise ValueError("native field block differs from its catalog declaration")
    if (value["field_id"] != field["identity"] or value["axes"] != field["axes"] or
            value["number_rule"] != field["number_rule"]):
        raise ValueError("foreign native field carrier or coefficient declaration")
    if value["layout"] == "dense":
        raw = owner.get_metadata(1, 0, "values")
        expected = field.get("storage_dtype", np.dtype(field["dtype"]).newbyteorder("<").str)
        if raw.dtype.str != expected or raw.dtype.kind not in "iuf":
            raise ValueError("native field scalar dtype differs from the declared source")
    elif value["layout"] == "exact":
        if tuple(owner.get_metadata(1, 0, "exact_values").shape) != tuple(value["shape"]):
            raise ValueError("native exact field shape differs from the declared source")


@dataclass(frozen=True)
class FieldSelection:
    dataset: str
    field: str
    specification: dict
    rows: tuple
    columns: tuple
    row_labels: tuple
    column_labels: tuple
    row_context: tuple
    column_context: tuple
    values: object
    observed: object
    dependencies: tuple
    parent: dict
    _seal: bytes = field(init=False, repr=False, compare=False)

    def _digest(self):
        from rexgraph.value_codec import pack_value
        digest = hashlib.sha256()
        digest.update(self.values.dtype.str.encode())
        if self.values.dtype.hasobject:
            for value in self.values.flat:
                digest.update(pack_value(value))
        else:
            digest.update(memoryview(self.values))
        digest.update(memoryview(self.observed))
        digest.update(pack_value((self.dataset, self.field, self.specification, self.rows, self.columns,
            self.row_labels, self.column_labels, self.row_context, self.column_context, self.dependencies, self.parent)))
        return digest.digest()

    def __post_init__(self):
        object.__setattr__(self, "_seal", self._digest())

    def check(self):
        if self._digest() != self._seal:
            raise ValueError("selected native data or its source mapping changed")

    def exact(self):
        """Convert only admitted selected coordinates through the native codec."""
        from rexgraph.exact_array import ExactArray
        from rexgraph.value import Absent
        self.check()
        values = np.asarray(self.values, dtype=object).copy()
        values[~self.observed] = Absent
        return ExactArray.from_values(values, rule=self.specification["number_rule"])

    def native(self):
        """Create an actual selected complex with original row identities."""
        owner = graph("selection/" + self.specification["identity"],
            tuple((self.dataset, self.specification["axes"][0], row) for row in self.rows),
            context=self.row_context)
        owner.attach_metadata(1, 0, "exact_values", self.exact())
        owner.attach_metadata(1, 0, "row_labels", self.row_labels)
        owner.attach_metadata(1, 0, "column_labels", self.column_labels)
        owner.attach_metadata(1, 0, "column_context", self.column_context)
        return declare(owner, {"schema": SCHEMA, "kind": "selection", "dataset": self.dataset,
            "field": self.field, "field_id": self.specification["identity"],
            "axes": self.specification["axes"], "rows": self.rows, "columns": self.columns,
            "shape": self.values.shape, "quantity": self.specification["quantity"],
            "units": self.specification["units"], "number_rule": self.specification["number_rule"],
            "dependencies": self.dependencies, "parent": self.parent})

    def tensor(self):
        """Return an exact field on the selected rows and declared column space.

        Missing coordinates require a declared missing data operation first.
        Publish and retrieve the selection to use it as a versioned query input.
        """
        from rexgraph.chain_map import CoordinateComplex
        from rexgraph.tensor_field import FieldSource, TensorField
        from rexgraph.type_accession import CoordinateSpace
        if not self.observed.all():
            raise ValueError("unobserved coordinates require a declared missing data operation")
        owner = self.native()
        value = owner.get_metadata(1, 0, "exact_values")
        references = tuple(FieldSource(None, ref["id"], ref["version"], ref["state_digest"])
                           for ref in self.dependencies)
        return TensorField(CoordinateComplex.from_rex(owner).spaces[1], value.values(),
            (CoordinateSpace(self.specification["axes"][1],
                tuple(f"{self.dataset}:{self.specification['axes'][1]}:{column}" for column in self.columns)),),
            FieldSource(owner), 1, "chain", dependencies=references)

    def materialize(self, output, *, name=None, limits=None):
        """Publish selected native coefficients with linked original dependencies."""
        from .store import DataWriter, context_links
        self.check()
        specification = self.specification
        with DataWriter(output, name=name or self.field, limits=limits,
                        parents=(dict(self.parent, dependencies=self.dependencies,
                            field_id=specification["identity"], rows=self.rows, columns=self.columns),)) as writer:
            row_axis, column_axis = specification["axes"]
            writer.axis(row_axis, self.row_labels, offset=0, context=self.row_context,
                        context_references=context_links(ref for ref in self.dependencies if ref.get("axis") == row_axis))
            writer.axis(column_axis, self.column_labels, offset=0, context=self.column_context,
                        context_references=context_links(ref for ref in self.dependencies if ref.get("axis") == column_axis))
            layout = "exact" if specification["layout"] == "exact" else "dense"
            writer.field(self.field, axes=specification["axes"], shape=self.values.shape, layout=layout,
                         number_rule=specification["number_rule"], quantity=specification["quantity"],
                         units=specification["units"])
            if layout == "exact":
                writer.exact(self.field, self.exact(), offset=(0, 0), dependencies=self.dependencies)
            else:
                writer.dense(self.field, self.values, self.observed, offset=(0, 0), dependencies=self.dependencies)
            return writer.finish()


@dataclass(frozen=True)
class OccurrenceSelection:
    dataset: str
    field: str
    specification: dict
    rows: tuple
    columns: tuple
    indices: object
    values: object
    original_entries: object
    dependencies: tuple
    parent: dict
    row_labels: tuple
    column_labels: tuple

    @property
    def implicit(self):
        return self.specification["implicit"]


def selection_item_bytes(specification):
    """Fixed selection slots, including separate observed and coverage masks."""
    return 96 if specification["layout"] == "exact" else np.dtype(specification["dtype"]).itemsize + 2


def select(store, name, *, rows, columns):
    specification = store.manifest["fields"][name]
    if specification["layout"] == "sparse_occurrences":
        raise ValueError("sparse occurrence arithmetic requires an explicit aggregation declaration")
    rows = coordinates(rows, specification["shape"][0], max_bytes=store.limits.selection_bytes)
    columns = coordinates(columns, specification["shape"][1], max_bytes=store.limits.selection_bytes)
    count = len(rows) * len(columns)
    # Exact selections retain fixed tensors and Python scalar objects. The
    # allocation admission counts slots conservatively, not total process RSS.
    item_bytes = selection_item_bytes(specification)
    store.limits.admit(count * item_bytes, "selection_bytes")
    dtype = object if specification["layout"] == "exact" else np.dtype(specification["dtype"])
    values = np.empty((len(rows), len(columns)), dtype=dtype)
    observed, covered = np.zeros(values.shape, bool), np.zeros(values.shape, bool)
    dependencies = []
    extra_integer_bytes = 0
    for block in store._blocks(name, rows, columns):
        ri = [(i, row - block["offset"][0]) for i, row in enumerate(rows)
              if block["offset"][0] <= row < block["offset"][0] + block["shape"][0]]
        ci = [(i, column - block["offset"][1]) for i, column in enumerate(columns)
              if block["offset"][1] <= column < block["offset"][1] + block["shape"][1]]
        if not ri or not ci:
            continue
        owner = store.read(block)
        _check(owner, block, specification)
        out = np.ix_([i for i, _ in ri], [i for i, _ in ci])
        if covered[out].any():
            raise ValueError("overlapping selected original field coordinates")
        selected, present = block_values(owner, [i for _, i in ri], [i for _, i in ci])
        if specification["layout"] == "exact":
            for value, active in zip(selected.flat, present.flat):
                if active:
                    numerator = int(value.numerator)
                    denominator = int(value.denominator)
                    extra_integer_bytes += max(0, (numerator.bit_length() + denominator.bit_length() + 7) // 8 - 16)
            store.limits.admit(count * item_bytes + 2 * extra_integer_bytes, "selection_bytes")
        values[out], observed[out], covered[out] = selected, present, True
        dependencies.append(block)
    if not covered.all():
        raise ValueError("selected original field coordinates are not covered")
    row_labels, row_context, row_refs = store.axis(specification["axes"][0], rows)
    column_labels, column_context, column_refs = store.axis(specification["axes"][1], columns)
    return FieldSelection(store.identity, name, {key: value for key, value in specification.items() if key != "blocks"}, rows, columns,
        row_labels, column_labels, row_context, column_context, frozen(values), frozen(observed),
        tuple(dependencies) + row_refs + column_refs, dict(store.reference))


def occurrences(store, name, *, rows, columns, blocks=None, include_axes=True, checked=False):
    specification = store.manifest["fields"][name]
    if specification["layout"] != "sparse_occurrences":
        raise ValueError("sparse occurrence field required")
    if not checked:
        rows = coordinates(rows, specification["shape"][0], max_bytes=store.limits.selection_bytes)
        columns = coordinates(columns, specification["shape"][1], max_bytes=store.limits.selection_bytes)
    indices, values, entries, dependencies = [], [], [], []
    total = 0
    for block in specification["blocks"] if blocks is None else blocks:
        if not all(selected is None or any(start <= i < start + count for i in selected)
                   for selected, start, count in zip((rows, columns), block["offset"], block["shape"], strict=True)):
            continue
        owner = store.read(block)
        _check(owner, block, specification)
        raw_indices = owner.get_metadata(1, 0, "indices")
        exact = block.get("value_layout") == "exact"
        raw_values = owner.get_metadata(1, 0, "exact_values" if exact else "values")
        if raw_indices.shape != (block["entry_count"], 2) or raw_values.size != block["entry_count"]:
            raise ValueError("native sparse occurrence extent differs")
        chosen = np.ones(len(raw_indices), dtype=bool)
        if rows is not None:
            chosen &= np.isin(raw_indices[:, 0], rows)
        if columns is not None:
            chosen &= np.isin(raw_indices[:, 1], columns)
        positions = np.flatnonzero(chosen)
        count = int(chosen.sum())
        total += count * (24 + (96 if exact else raw_values.dtype.itemsize))
        store.limits.admit(total, "selection_bytes")
        indices.append(raw_indices[chosen])
        values.append(exact_take(raw_values, positions).values() if exact else raw_values[chosen])
        if block.get("mapped_occurrences"):
            mapping = owner.get_metadata(1, 0, "original_entries")
            if mapping.shape != (block["entry_count"],) or mapping.dtype.kind not in "iu":
                raise ValueError("native original occurrence map differs")
            entries.append(mapping[positions])
        else:
            entries.append(positions + block["entry_offset"])
        dependencies.append(block)
    dtype = object if specification.get("value_layout") == "exact" else specification.get("dtype", "int64")
    if include_axes:
        row_labels, _, row_refs = store.axis(specification["axes"][0], rows)
        column_labels, _, column_refs = store.axis(specification["axes"][1], columns)
        dependencies.extend((*row_refs, *column_refs))
    else:
        row_labels, column_labels = (), ()
    return OccurrenceSelection(store.identity, name, {key: value for key, value in specification.items() if key != "blocks"},
        range(specification["shape"][0]) if rows is None else rows,
        range(specification["shape"][1]) if columns is None else columns,
        frozen(np.concatenate(indices) if indices else np.empty((0, 2), dtype=np.int64)),
        frozen(np.asarray(np.concatenate(values), dtype=dtype) if values else np.empty(0, dtype=dtype)),
        frozen(np.concatenate(entries) if entries else np.empty(0, dtype=np.int64)),
        tuple(dependencies), dict(store.reference), row_labels, column_labels)
