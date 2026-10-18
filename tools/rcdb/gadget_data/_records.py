"""Application declarations carried by canonical native records."""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import os
import stat

import numpy as np

SCHEMA = "gadget.native-data.v1"
UNSET = object()


def index_term(*parts):
    from rexgraph.value_codec import pack_value
    return "gadgetdata" + hashlib.sha256(pack_value(parts)).hexdigest()


def natural(value, name, *, positive=False):
    if type(value) is not int or value < int(positive):
        raise ValueError(name + " requires a " + ("positive" if positive else "nonnegative") + " integer")
    return value


def text(value, name):
    if type(value) is not str or not value:
        raise ValueError(name + " requires nonempty text")
    return value


def file_state(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def fingerprint(path):
    """Hash one regular original file without retaining its contents."""
    digest = hashlib.sha256()
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0))
    with os.fdopen(descriptor, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("original source must be a regular file")
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        if file_state(before) != file_state(os.fstat(stream.fileno())):
            raise ValueError("original source changed while hashing")
        if file_state(before) != file_state(Path(path).stat()):
            raise ValueError("original source path changed while hashing")
    return {"sha256": digest.hexdigest(), "bytes": before.st_size}


def graph(kind, keys, *, context=None):
    """Independent original row uses with explicitly typed participants."""
    from rexgraph import Relations, RexGraph, VertexTable
    keys = tuple(keys)
    if not keys or len(set(keys)) != len(keys):
        raise ValueError("distinct nonempty row occurrence identities required")
    vertices = [("row", key) for key in keys] + [("frame", kind)]
    positions = {key: i for i, key in enumerate(vertices)}
    supports = [[i, len(keys)] for i in range(len(keys))]
    if context is not None:
        if len(context) != len(keys):
            raise ValueError("context must retain every original row")
        for support, row in zip(supports, context, strict=True):
            for role, value in row.items():
                text(role, "context role")
                if value is None:
                    continue
                # Type remains part of identity: false, zero and text differ.
                key = ("context/" + role, type(value).__name__, value)
                if key not in positions:
                    positions[key] = len(vertices)
                    vertices.append(key)
                support.append(positions[key])
    return RexGraph.from_relations(Relations.from_supports(
        supports, vertices=VertexTable(tuple(vertices)), heads=[0] * len(keys),
        relation_ids=keys, relation_types=[kind] * len(keys)))


def declare(owner, value):
    owner.attach_metadata(1, 0, "gadget_data", value)
    return owner


def declaration(owner):
    from rexgraph import RexGraph
    if not isinstance(owner, RexGraph):
        raise TypeError("native data requires an actual relational complex")
    value = owner.get_metadata(1, 0, "gadget_data")
    if not isinstance(value, dict) or value.get("schema") != SCHEMA:
        raise ValueError("native data declaration missing")
    return value


def frozen(array):
    result = np.array(array, copy=True, order="C")
    if not result.dtype.hasobject:
        return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class Limits:
    """Execution admissions, adjustable independently of scientific dimensions.

    record_bytes bounds a serialized native block. decode_bytes bounds the
    declared uncompressed block before RCDB reconstruction. selection_bytes
    bounds selected typed arrays. These are not total process RSS guarantees.
    """
    record_bytes: int = 16 * 1024 * 1024
    decode_bytes: int = 64 * 1024 * 1024
    selection_bytes: int = 64 * 1024 * 1024
    catalog_bytes: int = 4 * 1024 * 1024

    def __post_init__(self):
        for name in ("record_bytes", "decode_bytes", "selection_bytes", "catalog_bytes"):
            natural(getattr(self, name), name, positive=True)

    def admit(self, count, resource):
        if count > getattr(self, resource):
            raise MemoryError(resource + " admission exceeded; select less data or change the execution budget")
