"""Generic original file intake into bounded native records."""
from contextlib import contextmanager
from itertools import islice
from pathlib import Path

import numpy as np
import mmap

from ._records import natural
from .store import DataWriter, Limits


def _axis(writer, name, extent, labels, *, block_size, source, selector):
    labels = iter(range(extent) if labels is None else labels)
    offset = 0
    while offset < extent:
        count = min(block_size, extent - offset)
        block = tuple(islice(labels, count))
        if len(block) != count:
            raise ValueError("axis labels do not cover the original extent")
        writer.axis(name, block, offset=offset,
            origins=({"source": source, "selector": dict(selector, axis=name, offset=offset, count=count)},))
        offset += count
    if tuple(islice(labels, 1)):
        raise ValueError("axis labels exceed the original extent")


@contextmanager
def _dense_file(path, *, reader, dataset):
    if reader == "npy":
        if dataset is not None:
            raise ValueError("npy intake does not select a container member")
        source = np.load(path, mmap_mode="r", allow_pickle=False)
        if not isinstance(source, np.memmap):
            if hasattr(source, "close"):
                source.close()
            raise ValueError("npy intake requires one typed NPY array")
        try:
            yield source
        finally:
            if isinstance(source, np.memmap):
                source._mmap.close()
    elif reader == "hdf5":
        import h5py
        if not dataset:
            raise ValueError("HDF5 intake requires an explicit dataset path")
        with h5py.File(path, "r") as container:
            # A single file content identity cannot own another file's values.
            current = container
            parts = [part for part in dataset.split("/") if part]
            for i, part in enumerate(parts):
                if not isinstance(current.get(part, getlink=True), h5py.HardLink):
                    raise ValueError("HDF5 selected dataset must use local hard links")
                current = current[part]
                if i + 1 < len(parts) and not isinstance(current, h5py.Group):
                    raise ValueError("HDF5 selected member path is not a group")
            if not isinstance(current, h5py.Dataset) or current.is_virtual or current.external:
                raise ValueError("HDF5 selected dataset must own its values in the original file")
            yield current
    else:
        raise ValueError("dense reader must be npy or hdf5")


def import_dense(source, output, *, reader, name, field="measurements", dataset=None,
                 axes=("rows", "columns"), row_labels=None, column_labels=None,
                 quantity, units, missing=None, tile=None, limits=None):
    """Map original typed matrices using bounded two dimensional selections.

    IEEE bits remain in their original dtype. A declared NaN marker becomes
    missing support with a zero placeholder, never an observed numerical zero.
    """
    limits = limits or Limits()
    if tile is not None:
        for count in tile:
            natural(count, "intake tile extent", positive=True)
    if tile is not None and len(tile) != 2 or missing not in (None, "nan"):
        raise ValueError("rank two tiles and an explicit supported missing rule required")
    source = Path(source).resolve(strict=True)
    with DataWriter(output, name=name, limits=limits) as writer:
        writer.source("original", source)
        with _dense_file(source, reader=reader, dataset=dataset) as matrix:
            if len(matrix.shape) != 2 or min(matrix.shape) <= 0 or matrix.dtype.kind not in "iuf":
                raise ValueError("nonempty original integer or IEEE matrix required")
            if reader == "hdf5" and matrix.chunks is not None:
                from math import prod
                limits.admit(prod(matrix.chunks) * matrix.dtype.itemsize, "decode_bytes")
            if tile is None:
                row_count = min(32, matrix.shape[0])
                target_bytes = min(limits.record_bytes, limits.decode_bytes) // 4
                column_count = max(1, min(matrix.shape[1], target_bytes // (row_count * (matrix.dtype.itemsize + 1))))
                tile = (row_count, column_count)
            selector = {"reader": reader, "dataset": dataset}
            _axis(writer, axes[0], matrix.shape[0], row_labels, block_size=256,
                  source="original", selector=selector)
            _axis(writer, axes[1], matrix.shape[1], column_labels, block_size=256,
                  source="original", selector=selector)
            writer.field(field, axes=axes, shape=tuple(map(int, matrix.shape)),
                number_rule="binary_exact" if matrix.dtype.kind == "f" else "exact",
                quantity=quantity, units=units)
            for row in range(0, matrix.shape[0], tile[0]):
                for column in range(0, matrix.shape[1], tile[1]):
                    shape = (min(tile[0], matrix.shape[0] - row), min(tile[1], matrix.shape[1] - column))
                    limits.admit(shape[0] * shape[1] * (matrix.dtype.itemsize + 1), "decode_bytes")
                    values = np.array(matrix[row:row + shape[0], column:column + shape[1]], copy=True)
                    observed = np.ones(values.shape, bool)
                    if values.dtype.kind == "f" and missing == "nan":
                        observed = ~np.isnan(values)
                        values[~observed] = 0
                    writer.dense(field, values, observed, offset=(row, column),
                        origins=({"source": "original", "selector": dict(selector,
                            rows=(row, row + shape[0]), columns=(column, column + shape[1]), missing=missing)},))
                if isinstance(matrix, np.memmap) and hasattr(matrix._mmap, "madvise"):
                    matrix._mmap.madvise(mmap.MADV_DONTNEED)
        return writer.finish()


def import_records(source, output, *, declaration, name, limits=None, grouping=None):
    """Use canonical declared readers and relation construction per native block."""
    from rexgraph import RexGraph
    from rexgraph.construction import construct
    from rexgraph.io.declaration import DatasetDeclaration
    from rexgraph.io.records import SourcePointer, TypedRecord
    if not isinstance(declaration, DatasetDeclaration):
        raise TypeError("a native DatasetDeclaration is required")
    limits = limits or Limits()
    if declaration.relations.shape == "long" and grouping != "contiguous":
        raise ValueError("long form relations require complete group construction before block publication")
    source = Path(source).resolve(strict=True)
    with DataWriter(output, name=name, limits=limits) as writer:
        original = writer.source("original", source)
        writer.source_declaration(declaration)
        def batches():
            if declaration.relations.shape != "long":
                for batch in declaration.read(source):
                    yield tuple(batch.records())
                return
            from ._records import UNSET
            from rexgraph.value_codec import pack_value
            current, pending, weight = UNSET, [], 0
            for batch in declaration.read(source):
                for record in batch.records():
                    identity = record.fields[declaration.relations.identity]
                    if pending and identity != current:
                        yield tuple(pending)
                        pending, weight = [], 0
                    current = identity
                    weight += 1024 + 2 * len(pack_value(dict(record.fields)))
                    limits.admit(weight, "decode_bytes")
                    pending.append(record)
            if pending:
                yield tuple(pending)

        for batch in batches():
            # Attach the immutable byte identity to every original row pointer.
            records = [TypedRecord(record.fields, SourcePointer(
                "gadget-source:original", bytes.fromhex(original["sha256"]), record.pointer.record,
                record.pointer.field, record.pointer.byte_range)) for record in batch]
            owner = RexGraph.from_relations(construct(records, declaration.relations, vertices=declaration.vertices))
            owner.attach_metadata(1, 0, "input_records", tuple(
                {"fields": dict(record.fields), "pointer": record.pointer.as_record()} for record in records))
            from ._records import UNSET
            writer.records_from(owner, declaration_digest=declaration.digest,
                origins=({"source": "original", "selector": {"reader": declaration.reader,
                    "records": tuple(record.pointer.record for record in records)}},),
                group=records[0].fields[declaration.relations.identity] if declaration.relations.shape == "long" else UNSET)
        return writer.finish()


def import_measurements(source, output, *, declaration, measurements, name,
                        quantity, units, field="measurements", label=None, context=(), limits=None,
                        index_labels=False, indexed_context=()):
    """Read declared table columns into exact blocks with original row pointers.

    The first bounded reader pass counts rows. The second publishes complete
    native relation records and exact measurement blocks. No whole table is held.
    """
    from rexgraph import RexGraph
    from rexgraph.construction import construct
    from rexgraph.io.declaration import DatasetDeclaration
    from rexgraph.io.records import SourcePointer, TypedRecord
    from rexgraph.value import Absent
    if not isinstance(declaration, DatasetDeclaration) or declaration.relations.shape == "long":
        raise ValueError("row or pair DatasetDeclaration required for measurement intake")
    measurements, context = tuple(measurements), tuple(context)
    fields = {entry.name: entry for entry in declaration.schema.fields}
    if (not measurements or len(set(measurements)) != len(measurements) or
            any(name not in fields or fields[name].kind != "number" for name in measurements)):
        raise ValueError("distinct declared numeric measurement columns required")
    if (label is not None and label not in fields or any(name not in fields for name in context)):
        raise ValueError("label and context must name original declared fields")
    source = Path(source).resolve(strict=True)
    with DataWriter(output, name=name, limits=limits) as writer:
        original = writer.source("original", source)
        writer.source_declaration(declaration)
        count = sum(len(batch) for batch in declaration.read(source))
        if not count:
            raise ValueError("measurement source has no original rows")
        writer.axis("measurements", measurements, offset=0, origins=({"source": "original",
                    "selector": {"reader": declaration.reader, "fields": measurements}},))
        writer.field(field, axes=("records", "measurements"), shape=(count, len(measurements)),
                     layout="exact", quantity=quantity, units=units)
        offset = 0
        for batch in declaration.read(source):
            records = [TypedRecord(record.fields, SourcePointer("gadget-source:original",
                bytes.fromhex(original["sha256"]), record.pointer.record, record.pointer.field,
                record.pointer.byte_range)) for record in batch.records()]
            if not records:
                continue
            owner = RexGraph.from_relations(construct(records, declaration.relations, vertices=declaration.vertices))
            owner.attach_metadata(1, 0, "input_records", tuple(
                {"fields": dict(record.fields), "pointer": record.pointer.as_record()} for record in records))
            provenance = ({"source": "original", "selector": {"reader": declaration.reader,
                           "records": tuple(record.pointer.record for record in records)}},)
            reference = writer.records_from(owner, declaration_digest=declaration.digest, origins=provenance)
            labels = tuple(record.pointer.record if label is None else record.fields[label] for record in records)
            if any(value is Absent for value in labels):
                raise ValueError("original row label is missing")
            context_rows = [{key: record.fields[key] for key in context if record.fields[key] is not Absent}
                            for record in records]
            writer.axis("records", labels, offset=offset, context=context_rows, origins=provenance,
                        index_labels=index_labels, indexed_context=indexed_context)
            values = [[record.fields[key] for key in measurements] for record in records]
            writer.exact(field, values, offset=(offset, 0), origins=provenance, dependencies=(reference,))
            offset += len(records)
        return writer.finish()


def import_coordinate(source, output, *, name, field="measurements", quantity, units,
                      implicit, block_entries=4096, axis_block=1024, limits=None):
    """Stream general Matrix Market coordinate occurrences without coalescing.

    Symmetry expansion is a separate mathematical operation. Only an explicitly
    general original matrix is admitted by this reader.
    """
    natural(block_entries, "occurrence block size", positive=True)
    natural(axis_block, "axis block size", positive=True)
    source = Path(source).resolve(strict=True)
    from rexgraph.exact_array import ExactArray
    from decimal import Decimal
    with DataWriter(output, name=name, limits=limits) as writer:
        writer.source("original", source)
        with source.open("r", encoding="ascii") as stream:
            header = stream.readline().strip().lower().split()
            if (len(header) != 5 or header[:3] != ["%%matrixmarket", "matrix", "coordinate"]
                    or header[3] not in ("integer", "real") or header[4] != "general"):
                raise ValueError("general integer or real Matrix Market coordinate source required")

            def lines():
                for line in stream:
                    if line.strip() and not line.lstrip().startswith("%"):
                        yield line.split()

            content = lines()
            dimensions = next(content, ())
            if len(dimensions) != 3:
                raise ValueError("coordinate source dimensions missing")
            rows, columns, count = map(int, dimensions)
            natural(rows, "row extent", positive=True); natural(columns, "column extent", positive=True)
            natural(count, "original occurrence count")
            _axis(writer, "rows", rows, None, block_size=axis_block, source="original", selector={"reader": "matrix_market"})
            _axis(writer, "columns", columns, None, block_size=axis_block, source="original", selector={"reader": "matrix_market"})
            writer.field(field, axes=("rows", "columns"), shape=(rows, columns), layout="sparse_occurrences",
                number_rule="exact" if header[3] == "integer" else "decimal_exact",
                quantity=quantity, units=units, implicit=implicit)
            offset = 0
            while offset < count:
                block = tuple(islice(content, min(block_entries, count - offset)))
                if not block or any(len(entry) != 3 for entry in block):
                    raise ValueError("original coordinate entries are truncated or malformed")
                indices = np.asarray([(int(entry[0]) - 1, int(entry[1]) - 1) for entry in block], dtype=np.int64)
                values = ExactArray.from_values([int(entry[2]) if header[3] == "integer" else Decimal(entry[2]) for entry in block],
                    rule="exact" if header[3] == "integer" else "decimal_exact")
                writer.sparse(field, indices, values, entry_offset=offset,
                    origins=({"source": "original", "selector": {"reader": "matrix_market",
                        "entries": (offset, offset + len(block)), "index_base": 1}},))
                offset += len(block)
            if next(content, None) is not None:
                raise ValueError("original coordinate entries exceed declared occurrence count")
        return writer.finish()
