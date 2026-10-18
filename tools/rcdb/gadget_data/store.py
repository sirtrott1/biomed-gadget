"""Versioned native blocks, source resolution and metadata exploration."""
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
import os
import stat
import uuid

import numpy as np

from ._records import (SCHEMA, Limits, declare, declaration, fingerprint, graph,
                       natural, text, index_term, UNSET, file_state)


def context_links(dependencies):
    """Retain shared context addresses once across selected original axis blocks."""
    from rexgraph.value_codec import pack_value
    result = {}
    for block in dependencies:
        for link in block.get("context_references", ()):
            result.setdefault(pack_value(link), link)
    return tuple(result.values())


def _open_store(path, *, create=False, backend="sqlite"):
    from rcdb import LocalStore, SQLStore
    path = Path(path)
    if create:
        if backend not in ("sqlite", "local"):
            raise ValueError("native data backend must be sqlite or local")
        path.mkdir(exist_ok=False)
        if backend == "local":
            return LocalStore(path)
        return SQLStore("sqlite:///" + str(path / "store.sqlite3"), native=True)
    if (path / "store.header").is_file():
        return LocalStore(path)
    if (path / "store.sqlite3").is_file():
        return SQLStore("sqlite:///" + str(path / "store.sqlite3"), native=True)
    raise ValueError("complete native data store required")


class SourceResolver:
    """Resolve original content by declared identity, independently of location."""
    def __init__(self, root=None, *, paths=None):
        self.root = None if root is None else Path(root).resolve(strict=True)
        self.paths = dict(paths or {})

    def resolve(self, reference, *, verify=True):
        identity = reference["id"]
        if identity in self.paths:
            path = Path(self.paths[identity]).resolve(strict=True)
        elif self.root is not None:
            locator = PurePosixPath(reference["locator"])
            if locator.is_absolute() or ".." in locator.parts:
                raise ValueError("source locator must be relative to the declared root")
            path = (self.root / str(locator)).resolve(strict=True)
            if not path.is_relative_to(self.root):
                raise ValueError("source locator leaves the declared root")
        else:
            raise KeyError("no original source location for " + identity)
        if verify and fingerprint(path) != {key: reference[key] for key in ("sha256", "bytes")}:
            raise ValueError("resolved original source differs from its retained content identity")
        return path

    @contextmanager
    def open(self, reference):
        """Pin and verify the actual opened original, including changes during use."""
        from ._records import file_state
        path = self.resolve(reference, verify=False)
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0))
        with os.fdopen(descriptor, "rb") as stream:
            if not stat.S_ISREG(os.fstat(stream.fileno()).st_mode):
                raise ValueError("original source must remain a regular file")
            before = file_state(os.fstat(stream.fileno()))
            # Verify the pinned descriptor rather than relying on the path check.
            import hashlib
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
            if (digest.hexdigest() != reference["sha256"] or
                    os.fstat(stream.fileno()).st_size != reference["bytes"] or
                    before != file_state(os.fstat(stream.fileno()))):
                raise ValueError("opened original source identity differs")
            stream.seek(0)
            yield stream
            if before != file_state(os.fstat(stream.fileno())):
                raise ValueError("original source changed during use")


class DataWriter:
    """Publish a complete dataset atomically through the canonical providers.

    A failed intake leaves no public dataset. Each native record is bounded;
    field values stay in typed native arrays and source containers stay external.
    """
    def __init__(self, output, *, name, limits=None, parents=(), backend="sqlite"):
        from .catalog import WriteIndex
        from rexgraph.io.publication import staged_publication
        self.output = Path(output).absolute()
        text(name, "dataset name")
        if os.path.lexists(self.output):
            raise FileExistsError(self.output)
        self.output.parent.mkdir(parents=True, exist_ok=True)
        self.limits = limits or Limits()
        self.identity, self.name = uuid.uuid4().hex, name
        self.parents = tuple(parents)
        self.sources, self._source_paths = {}, {}
        self._source_states = {}
        self.axes, self.fields = {}, {}
        self.axis_queries = {}
        self.records = WriteIndex(self, "records", "rows")
        self.declarations = WriteIndex(self, "declarations", "source")
        self.objects = WriteIndex(self, "objects", "native")
        self._closed = False
        self._staging = staged_publication(self.output, directory=True)
        self._path = self._staging.__enter__()
        try:
            (self._path / ".gitignore").write_text("*\n")
            self.store = _open_store(self._path / "data.rcdb", create=True, backend=backend)
        except BaseException:
            self._staging.__exit__(*__import__("sys").exc_info())
            self._closed = True
            raise

    def _put(self, identifier, owner, metadata, *, vocabulary=()):
        from rcdb.engine import encode_native_payload
        from rexgraph.object_identity import object_digest
        value = dict(metadata, schema=SCHEMA, dataset=self.identity)
        declare(owner, value)
        raw = encode_native_payload(owner, self.store.header.compression)
        self.limits.admit(len(raw), "record_bytes")
        self.limits.admit(len(raw), "decode_bytes")
        if self.store.get_record(identifier) is not None:
            raise ValueError("native data record already published")
        # Structural analytics are explicit query work, not an intake requirement.
        stored = self.store.put_prepared(identifier, raw,
            {"object_type": "RexGraph", "nV": int(owner.nV), "nE": int(owner.nE)},
            meta={"dataset": self.identity, "kind": value["kind"], "native_bytes": len(raw),
                  "vertex_labels": tuple(vocabulary)},
            tags=["gadget_data", value["kind"]])
        return {"id": identifier, "version": stored.version,
                "state_digest": object_digest(owner), "native_bytes": len(raw),
                "store_id": self.store.store_id}

    def source(self, identifier, path, *, locator=None, media_type=None):
        text(identifier, "source identity")
        if identifier in self.sources:
            raise ValueError("original source identity already registered")
        path = Path(path).resolve(strict=True)
        locator = path.name if locator is None else locator
        relative = PurePosixPath(text(locator, "source locator"))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("original source locator must be relative")
        reference = dict(fingerprint(path), id=identifier, locator=str(relative), media_type=media_type)
        self.sources[identifier] = reference
        self._source_paths[identifier] = path
        self._source_states[identifier] = file_state(path.stat())
        return dict(reference)

    def _origins(self, origins):
        from rexgraph.value_codec import pack_value, unpack_value
        origins = tuple(unpack_value(pack_value(origin)) for origin in origins)
        for origin in origins:
            if not isinstance(origin, dict) or origin.get("source") not in self.sources:
                raise ValueError("origin must name a registered immutable source")
            if not isinstance(origin.get("selector"), dict):
                raise ValueError("origin requires an explicit source selector")
        return origins

    def axis(self, name, labels, *, offset, context=None, origins=(), mapping=None,
             index_labels=False, indexed_context=(), context_references=()):
        from .catalog import WriteIndex
        text(name, "axis name"); natural(offset, "axis offset")
        labels = tuple(labels)
        indexed_context = tuple(indexed_context)
        if not labels:
            raise ValueError("axis block must be nonempty")
        if type(index_labels) is not bool or any(type(role) is not str or not role for role in indexed_context):
            raise ValueError("axis indexes require an explicit Boolean and named context roles")
        config = {"index_labels": bool(index_labels), "indexed_context": tuple(indexed_context)}
        if name in self.axis_queries and self.axis_queries[name] != config:
            raise ValueError("axis index declaration must be consistent across its blocks")
        self.axis_queries[name] = config
        if name not in self.axes:
            self.axes[name] = WriteIndex(self, "axis", name)
        blocks = self.axes[name]
        if offset != blocks.extent:
            raise ValueError("axis blocks must retain contiguous original coordinate order")
        owner = graph("axis/" + name,
                      tuple((self.identity, name, offset + i) for i in range(len(labels))), context=context)
        owner.attach_metadata(1, 0, "labels", labels)
        metadata = {"kind": "axis", "axis": name, "offset": offset, "count": len(labels),
                    "origins": self._origins(origins), "mapping": mapping,
                    "context_references": self._links(context_references),
                    "context_roles": tuple(dict.fromkeys(key for row in context or () for key in row))}
        vocabulary = {index_term("axis", self.identity, name)}
        if index_labels:
            vocabulary.update(index_term("label", value) for value in labels)
        for row in context or ():
            vocabulary.update(index_term("context", key, value) for key, value in row.items() if key in indexed_context)
        entry = self._put(f"axes/{name}/{offset}", owner, metadata, vocabulary=vocabulary)
        blocks.append(dict(entry, **metadata))
        return entry

    def _links(self, links):
        from copy import deepcopy
        from .workflow import check_link
        links = tuple(deepcopy(link) for link in links)
        parents = {parent["dataset"]: parent for parent in self.parents}
        for link in links:
            check_link(link)
            parent = link["dataset"]
            if parent["dataset"] == self.identity:
                if parent["store_id"] != self.store.store_id:
                    raise ValueError("local dependency has a foreign store")
                continue
            previous = parents.setdefault(parent["dataset"], parent)
            if any(previous[key] != parent[key] for key in ("dataset", "id", "version", "store_id", "state_digest")):
                raise ValueError("conflicting linked dataset versions")
        self.parents = tuple(parents.values())
        return links

    def native(self, name, value, *, dependencies=()):
        """Retain a native program, field, family or context with explicit owners."""
        from .workflow import retain
        return retain(self, name, value, dependencies=dependencies)

    def field(self, name, *, axes, shape, layout="dense", number_rule="exact",
              quantity, units, implicit=None):
        from .catalog import WriteIndex
        from rexgraph.value import NumberRule
        text(name, "field name")
        if name in self.fields:
            raise ValueError("field already declared")
        axes, shape = tuple(axes), tuple(shape)
        if len(axes) != 2 or len(set(axes)) != 2 or len(shape) != 2:
            raise ValueError("field requires two distinct named axes and their extents")
        for axis, size in zip(axes, shape, strict=True):
            text(axis, "axis name"); natural(size, "axis extent", positive=True)
        if layout not in ("dense", "exact", "sparse_occurrences"):
            raise ValueError("unknown native field layout")
        if layout == "sparse_occurrences" and implicit not in ("missing", "observed_zero"):
            raise ValueError("sparse field requires declared implicit support")
        if layout != "sparse_occurrences" and implicit is not None:
            raise ValueError("implicit support applies only to sparse occurrences")
        self.fields[name] = {"identity": uuid.uuid4().hex, "axes": axes, "shape": shape,
            "layout": layout, "number_rule": NumberRule(number_rule).value,
            "quantity": text(quantity, "quantity"), "units": text(units, "units"),
            "implicit": implicit, "blocks": WriteIndex(self, "field", name)}

    def dense(self, name, values, observed, *, offset, origins=(), dependencies=()):
        field = self.fields[name]
        if field["layout"] != "dense":
            raise ValueError("typed dense field declaration required")
        values, observed = np.asarray(values), np.asarray(observed)
        if (values.ndim != 2 or observed.shape != values.shape or observed.dtype != np.dtype(bool)
                or values.dtype.kind not in "iuf" or not values.size):
            raise ValueError("nonempty typed values and separate Boolean support required")
        if values.dtype.kind == "f" and not np.isfinite(values).all():
            raise ValueError("native typed coefficients must be finite; missing uses support")
        if values.dtype.kind == "f" and field["number_rule"] != "binary_exact":
            raise ValueError("typed IEEE fields require the binary_exact number rule")
        self._dtype(field, values.dtype)
        self.limits.admit(values.nbytes + observed.nbytes, "decode_bytes")
        owner = self._field_owner(name, offset, values.shape)
        owner.attach_metadata(1, 0, "values", np.array(values, copy=True, order="C"))
        owner.attach_metadata(1, 0, "observed", np.packbits(observed.reshape(-1), bitorder="little"))
        return self._field_put(name, owner, offset, values.shape, origins, dependencies=tuple(dependencies))

    def exact(self, name, values, *, offset, origins=(), dependencies=()):
        from rexgraph.exact_array import ExactArray
        field = self.fields[name]
        if field["layout"] != "exact":
            raise ValueError("exact field declaration required")
        if not isinstance(values, ExactArray):
            values = ExactArray.from_values(values, rule=field["number_rule"])
        if len(values.shape) != 2 or not values.size:
            raise ValueError("nonempty exact matrix required")
        self.limits.admit(sum(array.nbytes for array in
            (values.numerator, values.denominator, values.presence, values.kind)), "decode_bytes")
        owner = self._field_owner(name, offset, values.shape)
        owner.attach_metadata(1, 0, "exact_values", values)
        return self._field_put(name, owner, offset, values.shape, origins, dependencies=tuple(dependencies))

    def _field_owner(self, name, offset, shape):
        field = self.fields[name]
        offset = tuple(offset)
        if len(offset) != 2:
            raise ValueError("field offset requires both original coordinates")
        for start, count, extent in zip(offset, shape, field["shape"], strict=True):
            natural(start, "field offset"); natural(int(count), "block extent", positive=True)
            if start + count > extent:
                raise ValueError("field block outside declared original axes")
        return graph("field/" + field["identity"],
            tuple((self.identity, field["axes"][0], offset[0] + i) for i in range(shape[0])))

    def _field_put(self, name, owner, offset, shape, origins, **extra):
        field = self.fields[name]
        if field["layout"] != "sparse_occurrences":
            # Physical tiles cover complete row bands in original order. This
            # admits overlap and gap checks without comparing every old pair.
            previous = field["blocks"][-1] if field["blocks"] else None
            if previous is None:
                valid = tuple(offset) == (0, 0)
            elif offset[0] == previous["offset"][0]:
                valid = shape[0] == previous["shape"][0] and offset[1] == previous["offset"][1] + previous["shape"][1]
            else:
                valid = (offset[0] == previous["offset"][0] + previous["shape"][0] and offset[1] == 0 and
                         previous["offset"][1] + previous["shape"][1] == field["shape"][1])
            if not valid:
                raise ValueError("field tiles must cover contiguous original row bands in order")
        metadata = {"kind": "field", "field": name, "field_id": field["identity"],
                    "layout": field["layout"], "offset": tuple(offset), "shape": tuple(shape),
                    "axes": field["axes"], "number_rule": field["number_rule"],
                    "origins": self._origins(origins), **extra}
        entry = self._put(f"fields/{name}/{len(field['blocks'])}", owner, metadata)
        field["blocks"].append(dict(entry, **metadata))
        return entry

    def sparse(self, name, indices, values, *, entry_offset, origins=(), dependencies=(), original_entries=None):
        from rexgraph.exact_array import ExactArray
        field = self.fields[name]
        if field["layout"] != "sparse_occurrences":
            raise ValueError("sparse occurrence declaration required")
        natural(entry_offset, "original occurrence offset")
        indices = np.asarray(indices)
        exact = isinstance(values, ExactArray)
        if not exact:
            values = np.asarray(values)
        count = values.size
        if (indices.dtype.kind not in "iu" or indices.shape != (count, 2) or
                len(values.shape) != 1 or not count or
                np.any(indices < 0) or np.any(indices >= np.asarray(field["shape"]))):
            raise ValueError("original occurrence coordinates and typed values required")
        if exact and not values.presence.all():
            raise ValueError("stored sparse occurrences require observed coefficients")
        if not exact and values.dtype.kind not in "iuf":
            raise ValueError("sparse coefficients require typed values or an ExactArray")
        if not exact and values.dtype.kind == "f" and not np.isfinite(values).all():
            raise ValueError("finite sparse occurrence coefficients required")
        if not exact and values.dtype.kind == "f" and field["number_rule"] != "binary_exact":
            raise ValueError("typed IEEE fields require the binary_exact number rule")
        layout = "exact" if exact else "typed"
        if "value_layout" in field and field["value_layout"] != layout:
            raise ValueError("sparse field blocks must retain one coefficient layout")
        field["value_layout"] = layout
        if not exact:
            self._dtype(field, values.dtype)
        if entry_offset != field["blocks"].occurrences:
            raise ValueError("sparse blocks must retain contiguous original occurrence order")
        value_bytes = sum(a.nbytes for a in (values.numerator, values.denominator, values.presence, values.kind)) if exact else values.nbytes
        self.limits.admit(indices.nbytes + value_bytes, "decode_bytes")
        owner = graph("occurrences/" + field["identity"], ((self.identity, name, entry_offset),))
        owner.attach_metadata(1, 0, "indices", np.asarray(indices, dtype="<i8"))
        owner.attach_metadata(1, 0, "exact_values" if exact else "values", values if exact else np.array(values, copy=True))
        if original_entries is not None:
            original_entries = np.asarray(original_entries)
            if original_entries.shape != (count,) or original_entries.dtype.kind not in "iu" or np.any(original_entries < 0):
                raise ValueError("original occurrence mapping must retain every entry")
            owner.attach_metadata(1, 0, "original_entries", original_entries)
        low, high = tuple(map(int, indices.min(axis=0))), tuple(map(int, indices.max(axis=0) + 1))
        return self._field_put(name, owner, low, tuple(b - a for a, b in zip(low, high)), origins,
                               entry_offset=entry_offset, entry_count=count, value_layout=layout,
                               dependencies=tuple(dependencies), mapped_occurrences=original_entries is not None)

    @staticmethod
    def _dtype(field, dtype):
        if "dtype" in field and field["dtype"] != dtype.str:
            raise ValueError("field blocks must retain one original scalar dtype")
        field["dtype"] = dtype.str
        field["storage_dtype"] = dtype.newbyteorder("<").str

    def records_from(self, owner, *, declaration_digest, origins, group=UNSET):
        """Retain a complex constructed by the canonical declared record reader."""
        identifier = f"records/{len(self.records)}"
        if group is not UNSET:
            # This is a physical lookup bucket. Native identity equality decides
            # whether an existing group repeats; hash equality does not.
            base, suffix = f"records/group/{hash(group)}", 0
            identifier = base
            while self.store.get_record(identifier) is not None:
                previous = self.store.get(identifier)
                if declaration(previous)["group_identity"] == group:
                    raise ValueError("long form group reappears after another group; complete contiguous groups required")
                suffix += 1
                identifier = base + "/" + str(suffix)
        entry = self._put(identifier, owner,
            {"kind": "records", "declaration_digest": declaration_digest,
             "origins": self._origins(origins), "group_identity": None if group is UNSET else group})
        self.records.append(entry)
        return entry

    def source_declaration(self, value):
        """Keep a canonical DatasetDeclaration once through its existing codec."""
        from rcdb.codecs import DECLARATION_CODEC
        raw = value.to_bytes()
        self.limits.admit(len(raw), "record_bytes")
        self.limits.admit(len(raw), "decode_bytes")
        identifier = f"declarations/{len(self.declarations)}"
        row = self.store.put_record(identifier, value, codec=DECLARATION_CODEC,
            meta={"dataset": self.identity, "kind": "declaration", "native_bytes": len(raw)}, expected_version=0)
        ref = {"id": identifier, "version": row.version, "store_id": self.store.store_id,
               "state_digest": row.envelope.object_digest, "native_bytes": len(raw), "kind": "declaration"}
        self.declarations.append(ref)
        return ref

    def _pages(self, kind, name, entries, *, start_index=0):
        from rexgraph.value_codec import pack_value
        pages, pending, size = [], [], 0

        def flush():
            nonlocal pending, size
            if not pending:
                return
            index = start_index + len(pages)
            reference = self._put(f"catalog/{kind}/{name}/{index}",
                graph("catalog", ((kind, name, index),)), {"kind": "catalog", "entries": tuple(pending)})
            page = {"page": reference, "count": len(pending)}
            if kind == "axis":
                page["extent"] = (pending[0]["offset"], pending[-1]["offset"] + pending[-1]["count"])
            elif kind == "field":
                lower = tuple(min(block["offset"][i] for block in pending) for i in range(2))
                upper = tuple(max(block["offset"][i] + block["shape"][i] for block in pending) for i in range(2))
                page["bounds"] = (lower, upper)
            pages.append(page)
            pending, size = [], 0

        for entry in entries:
            count = len(pack_value(entry))
            if pending and (len(pending) == 128 or size + count > min(self.limits.record_bytes, self.limits.decode_bytes) // 2):
                flush()
            pending.append(entry); size += count
        flush()
        return pages

    def finish(self):
        from rexgraph.io.publication import publication_lock
        from .identity import runtime_identity
        if self._closed:
            raise ValueError("native dataset writer is closed")
        try:
            for name, field in self.fields.items():
                for axis, extent in zip(field["axes"], field["shape"], strict=True):
                    if axis not in self.axes or self.axes[axis].extent != extent:
                        raise ValueError("field axis does not cover the complete original extent")
                if field["layout"] != "sparse_occurrences":
                    covered = field["blocks"].area
                    if covered != field["shape"][0] * field["shape"][1]:
                        raise ValueError("field blocks leave uncovered original coordinates")
            for identity, reference in self.sources.items():
                path = self._source_paths[identity]
                if (file_state(path.stat()) != self._source_states[identity] or
                        fingerprint(path) != {k: reference[k] for k in ("sha256", "bytes")} or
                        file_state(path.stat()) != self._source_states[identity]):
                    raise ValueError("original source changed during native intake")
            axes = {name: blocks.finish() for name, blocks in self.axes.items()}
            fields = {name: dict(field, blocks=field["blocks"].finish(),
                                 block_count=len(field["blocks"])) for name, field in self.fields.items()}
            manifest = {"kind": "dataset", "name": self.name, "sources": self.sources,
                        "axes": axes, "fields": fields, "records": self.records.finish(),
                        "objects": self.objects.finish(),
                        "parents": self.parents, "declarations": self.declarations.finish(), "producer": runtime_identity(),
                        "axis_queries": self.axis_queries}
            result = self._put("dataset", graph("dataset", (self.identity,)), manifest)
            self.store.close()
            # The canonical directory publisher arbitrates the final fresh name.
            with publication_lock(self.output.parent):
                if os.path.lexists(self.output):
                    raise FileExistsError(self.output)
                self._staging.__exit__(None, None, None)
            self._closed = True
            return result
        except BaseException:
            self.close()
            raise

    def close(self):
        if not self._closed:
            self.store.close()
            error = RuntimeError("native data intake not committed")
            self._staging.__exit__(type(error), error, error.__traceback__)
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class DataStore:
    """Open only a sealed catalog; decode mathematical blocks on selection."""
    def __init__(self, path, *, limits=None, dataset=None):
        from collections import OrderedDict
        from .catalog import Entries
        self.path = Path(path).resolve(strict=True)
        self.limits = limits or Limits()
        self.store = _open_store(self.path / "data.rcdb")
        try:
            row = self.store.get_record("dataset")
            if row is None or row.version != 1:
                raise ValueError("sealed native dataset manifest missing")
            self.limits.admit(row.meta["native_bytes"], "decode_bytes")
            snapshot = self.store.read_record("dataset", version=1)
            from copy import deepcopy
            self.manifest = deepcopy(declaration(snapshot.value))
            if self.manifest["kind"] != "dataset" or dataset is not None and dataset != self.manifest["dataset"]:
                raise ValueError("dataset identity differs from the requested source")
            self.identity = self.manifest["dataset"]
            self._catalog = deepcopy(self.manifest)
            self._page_cache = OrderedDict()
            self._cache_bytes = 0
            self.manifest["axes"] = {name: Entries(self, pages) for name, pages in self.manifest["axes"].items()}
            self.manifest["fields"] = {name: dict(field, blocks=Entries(self, field["blocks"]))
                                       for name, field in self.manifest["fields"].items()}
            self.manifest["records"] = Entries(self, self.manifest["records"])
            self.manifest["declarations"] = Entries(self, self.manifest["declarations"])
            self.manifest["objects"] = Entries(self, self.manifest.get("objects", ()))
            self.reference = {"id": "dataset", "version": 1, "state_digest": snapshot.state_digest,
                              "store_id": self.store.store_id, "dataset": self.identity}
        except BaseException:
            self.store.close()
            raise

    def read(self, reference):
        """Validate address, version and state before returning a native owner."""
        if reference["store_id"] != self.store.store_id:
            raise ValueError("foreign native store reference")
        row = next((row for row in self.store.history(reference["id"])
                    if row.id == reference["id"] and row.version == reference["version"]), None)
        if (row is None or row.meta.get("native_bytes") != reference["native_bytes"] or
                row.envelope is None or row.envelope.object_digest != reference["state_digest"]):
            raise ValueError("native record address or retained envelope differs")
        self.limits.admit(max(reference["native_bytes"], row.envelope.payload_size), "decode_bytes")
        snapshot = self.store.read_record(reference["id"], version=reference["version"])
        if snapshot is None or snapshot.state_digest != reference["state_digest"]:
            raise ValueError("missing or substituted native record version")
        if row.meta.get("dataset") != self.identity:
            raise ValueError("foreign native dataset block")
        if reference.get("kind") != "declaration" and declaration(snapshot.value)["dataset"] != self.identity:
            raise ValueError("foreign native dataset owner")
        return snapshot.value

    def catalog(self):
        """Return administrative metadata without opening numerical blocks."""
        from copy import deepcopy
        value = deepcopy(self._catalog)
        value["reference"] = dict(self.reference)
        return value

    def axis(self, name, coordinates):
        from copy import deepcopy
        from .selection import coordinates as checked
        blocks = self.manifest["axes"][name]
        extent = blocks.extent
        selected = checked(coordinates, extent, max_bytes=self.limits.selection_bytes)
        labels, contexts, dependencies = {}, {}, []
        metadata_bytes = 0
        for block in blocks.axis(selected):
            local = [i for i in selected if block["offset"] <= i < block["offset"] + block["count"]]
            if not local:
                continue
            actual, context_rows = self._axis_metadata(block)
            metadata_bytes = self._admit_metadata(block, metadata_bytes)
            dependencies.append(block)
            for i in local:
                metadata_bytes = self._admit_metadata(
                    (actual[i - block["offset"]], context_rows[i - block["offset"]]), metadata_bytes)
                labels[i] = deepcopy(actual[i - block["offset"]])
                contexts[i] = deepcopy(context_rows[i - block["offset"]])
        if len(labels) != len(selected):
            raise ValueError("native axis has missing original coordinates")
        return (tuple(labels[i] for i in selected), tuple(contexts[i] for i in selected), tuple(dependencies))

    def _admit_metadata(self, value, retained):
        from rexgraph.value_codec import pack_value
        retained += len(pack_value(value)) + 64
        self.limits.admit(retained, "selection_bytes")
        return retained

    def select(self, field, *, rows, columns):
        from .selection import select
        return select(self, field, rows=rows, columns=columns)

    def native(self, reference, *, dependencies=()):
        """Restore a retained native value using explicitly supplied dataset stores."""
        from .workflow import restore
        return restore(self, reference, dependencies=dependencies)

    def contexts(self, axis_reference, *, dependencies=()):
        """Resolve shared context owners without inferring joins or merging assertions."""
        from .workflow import resolve_link
        metadata = declaration(self.read(axis_reference))
        if metadata["kind"] != "axis":
            raise TypeError("native axis record required")
        stores = {store.identity: store for store in (self, *dependencies)}
        for link in metadata.get("context_references", ()):
            yield resolve_link(link, stores)

    def find_axis(self, name, *, label=UNSET, context=None):
        """Stream original axis matches using declared RCDB vocabulary where available."""
        from rexgraph.value_codec import pack_value
        context = dict(context or {})
        config = self.manifest.get("axis_queries", {}).get(name, {})
        terms = []
        if label is not UNSET and config.get("index_labels"):
            terms.append(index_term("label", label))
        terms.extend(index_term("context", key, value) for key, value in context.items()
                     if key in config.get("indexed_context", ()))
        candidates = None
        if terms:
            cap = max(1, self.limits.catalog_bytes // 2048)
            found = self.store.query(limit=cap, labels_all=[index_term("axis", self.identity, name), *terms])
            if len(found) < cap:
                candidates = {row.id for row in found}
        for block in self.manifest["axes"][name]:
            if candidates is not None and block["id"] not in candidates:
                current = self.store.get_record(block["id"])
                if current is not None and current.version == block["version"]:
                    continue
            labels, contexts = self._axis_metadata(block)
            for i, (actual, row) in enumerate(zip(labels, contexts, strict=True)):
                if label is not UNSET and pack_value(actual) != pack_value(label):
                    continue
                if any(key not in row or pack_value(row[key]) != pack_value(value) for key, value in context.items()):
                    continue
                from copy import deepcopy
                yield {"coordinate": block["offset"] + i, "label": deepcopy(actual),
                       "context": deepcopy(row), "source": dict(block)}

    def _blocks(self, name, rows, columns):
        yield from self.manifest["fields"][name]["blocks"].field(rows, columns)

    def _page(self, reference):
        def load():
            owner = self.read(reference)
            value = declaration(owner)
            if value["kind"] != "catalog":
                raise ValueError("native catalog page kind differs")
            return tuple(value["entries"])
        return self._metadata(reference, "page", load)

    def _metadata(self, reference, kind, load):
        while self._cache_bytes > self.limits.catalog_bytes:
            _, (_, old_weight) = self._page_cache.popitem(last=False)
            self._cache_bytes -= old_weight
        key = (kind, reference["id"], reference["version"], reference["state_digest"],
               reference.get("axis"), reference.get("offset"), reference.get("count"))
        if key not in self._page_cache:
            value, weight = load(), reference["native_bytes"]
            if weight > self.limits.catalog_bytes:
                return value
            while self._cache_bytes + weight > self.limits.catalog_bytes:
                _, (_, old_weight) = self._page_cache.popitem(last=False)
                self._cache_bytes -= old_weight
            self._page_cache[key] = (value, weight)
            self._cache_bytes += weight
        self._page_cache.move_to_end(key)
        return self._page_cache[key][0]

    def _axis_metadata(self, reference):
        def load():
            owner = self.read(reference)
            value = declaration(owner)
            if value["kind"] != "axis":
                raise ValueError("native axis record kind differs")
            if any(value[key] != reference[key] for key in ("axis", "offset", "count")):
                raise ValueError("native axis record coordinates differ")
            labels = tuple(owner.get_metadata(1, 0, "labels"))
            if len(labels) != reference["count"] or owner.nE != len(labels):
                raise ValueError("native axis record extent differs")
            rows = []
            for row in range(len(labels)):
                start, end = owner.relations.support_ptr[row:row + 2]
                context = {}
                for vertex in owner.relations.support_idx[int(start):int(end)]:
                    key = owner.relations.vertices.ids[int(vertex)]
                    if key[0].startswith("context/"):
                        context[key[0][8:]] = key[2]
                rows.append(context)
            return labels, tuple(rows)
        return self._metadata(reference, "axis", load)

    def clear_metadata_cache(self):
        self._page_cache.clear()
        self._cache_bytes = 0

    def occurrences(self, field, *, rows, columns):
        from .selection import occurrences
        return occurrences(self, field, rows=rows, columns=columns)

    def iter_occurrences(self, field, *, rows=None, columns=None):
        """Stream stored occurrences without expanding a sparse logical domain."""
        from .selection import coordinates, occurrences
        specification = self.manifest["fields"][field]
        if specification["layout"] != "sparse_occurrences":
            raise ValueError("sparse occurrence field required")
        rows = None if rows is None else coordinates(rows, specification["shape"][0], max_bytes=self.limits.selection_bytes)
        columns = None if columns is None else coordinates(columns, specification["shape"][1], max_bytes=self.limits.selection_bytes)
        for block in specification["blocks"].field(rows, columns):
            yield occurrences(self, field, rows=rows, columns=columns, blocks=(block,), include_axes=False, checked=True)

    def provenance(self, field, *, rows, columns):
        """Map selected coordinates to retained records without decoding values."""
        from copy import deepcopy
        from .selection import coordinates
        specification = self.manifest["fields"][field]
        rows = coordinates(rows, specification["shape"][0], max_bytes=self.limits.selection_bytes)
        columns = coordinates(columns, specification["shape"][1], max_bytes=self.limits.selection_bytes)
        result = {"dataset": self.reference, "field": specification["identity"],
            "rows": rows, "columns": columns, "axes": {}, "field_records": [],
            "originals": self.manifest["sources"], "parents": self.manifest["parents"]}
        retained = self._admit_metadata(result, 0)

        def bounded(entries):
            nonlocal retained
            output = []
            for entry in entries:
                retained = self._admit_metadata(entry, retained)
                output.append(entry)
            return output

        if specification["layout"] == "sparse_occurrences":
            result["field_records"] = bounded(specification["blocks"].field(rows, columns))
        else:
            result["field_records"] = bounded(self._blocks(field, rows, columns))
        for axis, selected in zip(specification["axes"], (rows, columns), strict=True):
            result["axes"][axis] = bounded(self.manifest["axes"][axis].axis(selected))
        return deepcopy(result)

    def iter_select(self, field, *, rows=None, columns=None):
        """Visit physical native field windows without building a whole matrix."""
        specification = self.manifest["fields"][field]
        if specification["layout"] == "sparse_occurrences":
            raise ValueError("use occurrence selection for sparse fields")
        from .selection import coordinates
        rows = None if rows is None else coordinates(rows, specification["shape"][0], max_bytes=self.limits.selection_bytes)
        columns = None if columns is None else coordinates(columns, specification["shape"][1], max_bytes=self.limits.selection_bytes)
        for block in specification["blocks"]:
            selected = []
            for requested, start, count in zip((rows, columns), block["offset"], block["shape"], strict=True):
                selected.append(tuple(range(start, start + count)) if requested is None else
                                tuple(i for i in requested if start <= i < start + count))
            if selected[0] and selected[1]:
                yield self.select(field, rows=selected[0], columns=selected[1])

    def materialize(self, field, output, *, rows, columns, tile=(256, 256), name=None):
        """Stream selected native tiles and preserve original axis maps in records."""
        from .selection import coordinates
        specification = self.manifest["fields"][field]
        if specification["layout"] == "sparse_occurrences":
            return self._materialize_occurrences(field, output, rows=rows, columns=columns, tile=tile, name=name)
        rows = coordinates(rows, specification["shape"][0], max_bytes=self.limits.selection_bytes)
        columns = coordinates(columns, specification["shape"][1], max_bytes=self.limits.selection_bytes)
        if len(tile) != 2:
            raise ValueError("materialization tile requires two extents")
        for count in tile:
            natural(count, "materialization tile extent", positive=True)
        with DataWriter(output, name=name or field, limits=self.limits,
                        parents=(dict(self.reference, field_id=specification["identity"]),)) as writer:
            for axis, selected, width in zip(specification["axes"], (rows, columns), tile, strict=True):
                for offset in range(0, len(selected), width):
                    original = selected[offset:offset + width]
                    labels, context, dependencies = self.axis(axis, original)
                    writer.axis(axis, labels, offset=offset, context=context,
                        context_references=context_links(dependencies),
                        mapping={"parent": self.reference, "axis": axis,
                                 "coordinates": original, "dependencies": dependencies})
            writer.field(field, axes=specification["axes"], shape=(len(rows), len(columns)),
                layout=specification["layout"], number_rule=specification["number_rule"],
                quantity=specification["quantity"], units=specification["units"])
            for row in range(0, len(rows), tile[0]):
                for column in range(0, len(columns), tile[1]):
                    selected = self.select(field, rows=rows[row:row + tile[0]], columns=columns[column:column + tile[1]])
                    if specification["layout"] == "dense":
                        writer.dense(field, selected.values, selected.observed, offset=(row, column),
                                     dependencies=selected.dependencies)
                    else:
                        writer.exact(field, selected.exact(), offset=(row, column), dependencies=selected.dependencies)
            return writer.finish()

    def _materialize_occurrences(self, field, output, *, rows, columns, tile, name):
        from .selection import coordinates, occurrences
        from rexgraph.exact_array import ExactArray
        specification = self.manifest["fields"][field]
        rows = coordinates(rows, specification["shape"][0], max_bytes=self.limits.selection_bytes)
        columns = coordinates(columns, specification["shape"][1], max_bytes=self.limits.selection_bytes)
        if len(tile) != 2:
            raise ValueError("materialization requires two physical axis tile extents")
        for count in tile:
            natural(count, "axis tile extent", positive=True)
        row_map, column_map = {value: i for i, value in enumerate(rows)}, {value: i for i, value in enumerate(columns)}
        parent = dict(self.reference, field_id=specification["identity"])
        with DataWriter(output, name=name or field, limits=self.limits, parents=(parent,)) as writer:
            for axis, selected, width in zip(specification["axes"], (rows, columns), tile, strict=True):
                for offset in range(0, len(selected), width):
                    original = selected[offset:offset + width]
                    labels, context, dependencies = self.axis(axis, original)
                    writer.axis(axis, labels, offset=offset, context=context,
                        context_references=context_links(dependencies),
                        mapping={"parent": self.reference, "axis": axis, "coordinates": original, "dependencies": dependencies})
            writer.field(field, axes=specification["axes"], shape=(len(rows), len(columns)), layout="sparse_occurrences",
                number_rule=specification["number_rule"], quantity=specification["quantity"], units=specification["units"],
                implicit=specification["implicit"])
            writer.fields[field]["original_entry_owner"] = specification.get("original_entry_owner", parent)
            offset = 0
            for block in specification["blocks"]:
                selected = occurrences(self, field, rows=rows, columns=columns, blocks=(block,), include_axes=False, checked=True)
                if not selected.values.size:
                    continue
                local = np.asarray([(row_map[int(row)], column_map[int(column)]) for row, column in selected.indices], dtype=np.int64)
                values = (ExactArray.from_values(selected.values, rule=specification["number_rule"])
                          if specification.get("value_layout") == "exact" else selected.values)
                writer.sparse(field, local, values, entry_offset=offset, dependencies=selected.dependencies,
                              original_entries=selected.original_entries)
                offset += selected.values.size
            return writer.finish()

    def verify_sources(self, resolver):
        return {identity: resolver.resolve(ref) for identity, ref in self.manifest["sources"].items()}

    def verify(self):
        """Check all native record versions one bounded owner at a time."""
        self.clear_metadata_cache()
        count = 1
        for blocks in self.manifest["axes"].values():
            for block in blocks:
                self.read(block)
                count += 1
        for field in self.manifest["fields"].values():
            for block in field["blocks"]:
                from .selection import _check
                _check(self.read(block), block, field)
                count += 1
        for block in self.manifest["records"]:
            self.read(block)
            count += 1
        for block in self.manifest["declarations"]:
            self.read(block)
            count += 1
        for block in self.manifest["objects"]:
            self.read(block)
            count += 1
        return count

    def close(self):
        self.clear_metadata_cache()
        self.store.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
