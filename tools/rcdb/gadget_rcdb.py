"""RCDB storage and RCQL retrieval of canonical GADGET publications."""
from pathlib import Path
import filecmp
import os
import subprocess
from tempfile import TemporaryDirectory
import zlib

_MAGIC = b"RFCXS\x01\x00\x00"
_CODEC_NAME = "gadget.native-publication"
_OPERATOR_NAME = "GADGET_PUBLICATION_BYTES"


def validate_publication_file(path, *, gadget):
    """Check canonical native bytes without a complete Python archive buffer."""
    return _validate_publication_file(path, gadget=gadget)


def _validate_publication_file(path, *, gadget, expected=None):
    source = Path(path).resolve(strict=True)
    binary = str(Path(gadget).resolve(strict=True))
    with source.open("rb") as stream:
        if stream.read(len(_MAGIC)) != _MAGIC:
            raise ValueError("GADGET publication requires RFCXS bytes")
    with TemporaryDirectory(prefix="gadget-native-admission-") as folder:
        canonical = Path(folder)/"canonical.rfcxs"
        process = subprocess.run([binary, "native-record", "--input", str(source),
                                  "--out", str(canonical)], capture_output=True)
        if process.returncode:
            raise ValueError("GADGET native admission failed: " +
                             process.stderr.decode("utf-8", errors="replace").strip())
        if expected is None:
            unchanged = filecmp.cmp(source, canonical, shallow=False)
        else:
            with canonical.open("rb") as stream:
                unchanged = all(stream.read(1024*1024) == expected[offset:offset+1024*1024]
                                for offset in range(0, len(expected), 1024*1024))
                unchanged = unchanged and not stream.read(1)
        if process.stdout or not unchanged:
            raise ValueError("GADGET publication does not have canonical native bytes")
    return source


def _read_publication_file(path, *, gadget):
    raw = Path(path).read_bytes()
    _validate_publication_file(path, gadget=gadget, expected=raw)
    return raw


def compress_publication(source, destination, *, gadget, codec, level=None):
    """Write an RFCXS file in RCDB's existing lossless compressed blob frame."""
    from rcdb.core import compress_blob
    if codec not in ("zstd", "zlib"):
        raise ValueError("compressed publication requires zstd or zlib")
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    raw = _read_publication_file(source, gadget=gadget)
    encoded = compress_blob(raw, level=level, codec=codec)
    with destination.open("xb") as stream:
        stream.write(encoded)
    return destination


def read_compressed_publication(path, *, gadget, max_encoded_bytes, max_output_bytes):
    """Admit an RCDB blob under explicit encoded and decoded byte budgets."""
    from rcdb.core import decompress_blob
    if type(max_encoded_bytes) is not int or max_encoded_bytes <= 0:
        raise ValueError("encoded byte budget must be a positive integer")
    if type(max_output_bytes) is not int or max_output_bytes <= 0:
        raise ValueError("decoded byte budget must be a positive integer")
    with Path(path).open("rb") as stream:
        size = os.fstat(stream.fileno()).st_size
        if size > max_encoded_bytes:
            raise ValueError("compressed publication exceeds its encoded byte budget")
        encoded = stream.read(size + 1)
    if len(encoded) != size:
        raise ValueError("compressed publication changed size while reading")
    try:
        raw = decompress_blob(encoded, max_output_bytes=max_output_bytes)
    except zlib.error as exc:
        raise ValueError("compressed publication has invalid zlib bytes") from exc
    return validate_publication(raw, gadget=gadget)


def install_compressed_format(*, gadget, max_encoded_bytes, max_output_bytes):
    """Register RCDB-framed RFCXS files with RexGraph's format registry."""
    from rexgraph.io import register_format
    binary = str(Path(gadget).resolve(strict=True))

    def load(path):
        return read_compressed_publication(path, gadget=binary,
                                            max_encoded_bytes=max_encoded_bytes,
                                            max_output_bytes=max_output_bytes)

    return register_format("gadget-native-compressed", load=load,
                           extensions=(".rcblob",))


def validate_publication(raw, *, gadget):
    """Validate and re-encode an RFCXS publication with GADGET's native reader."""
    if type(raw) is not bytes or not raw.startswith(_MAGIC):
        raise ValueError("GADGET publication requires RFCXS bytes")
    with TemporaryDirectory(prefix="gadget-rcdb-") as folder:
        source = Path(folder)/"input.rfcxs"
        source.write_bytes(raw)
        _validate_publication_file(source, gadget=gadget, expected=raw)
    return raw


def install_codec(*, gadget, max_bytes, replace=False):
    """Install an RFCXS provider under the caller's declared RCDB byte budget."""
    from rcdb import CodecRef, RecordCodec, register_record_codec
    binary = str(Path(gadget).resolve(strict=True))

    def encode(raw):
        if type(raw) is not bytes or not raw.startswith(_MAGIC):
            raise ValueError("GADGET publication requires RFCXS bytes")
        return raw

    def decode(raw):
        return validate_publication(raw, gadget=binary)

    reference = CodecRef(_CODEC_NAME, 1)
    return register_record_codec(RecordCodec(reference, "GadgetNativePublication", encode,
                                             decode, max_bytes=max_bytes), replace=replace)


def put_publication(store, record_id, path, *, expected_version, **selectors):
    """Publish unchanged native bytes through RCDB's versioned record engine."""
    from rcdb import CodecRef
    from rcdb.codecs import record_codec
    reference = CodecRef(_CODEC_NAME, 1)
    provider = record_codec(reference)
    path = Path(path)
    if path.stat().st_size > provider.max_bytes:
        raise ValueError("publication exceeds the caller's declared RCDB byte budget")
    return store.put_record(record_id, path.read_bytes(), codec=reference,
                            expected_version=expected_version, **selectors)


def read_publication(store, record_id, *, version=None, as_of=None, valid_at=None):
    """Read one selected RFCXS payload without reconstructing a Python graph."""
    snapshot = store.read_record(record_id, version=version, as_of=as_of, valid_at=valid_at)
    if snapshot is None:
        raise KeyError(f"GADGET publication {record_id!r} is not present")
    envelope = snapshot.record.envelope
    if envelope is None or envelope.codec != _CODEC_NAME or envelope.codec_version != 1:
        raise ValueError("selected RCDB record does not use the GADGET publication codec")
    if type(snapshot.value) is not bytes or not snapshot.value.startswith(_MAGIC):
        raise ValueError("selected RCDB publication has no native RFCXS payload")
    return snapshot


def install_rcql_operator():
    """Register typed current/versioned native-byte retrieval on an RCDB source."""
    from rcql import Domain, Exactness, OperatorSignature, RCType, TypePattern, ValueKind
    from rcql.execution_trace import record_method
    from rcql.operators import register as register_operator
    from rcql.signatures import catalogued, lookup, register as register_signature
    signature = OperatorSignature(
        name=_OPERATOR_NAME, source_kind=ValueKind.RCDB_STORE,
        inputs=(TypePattern("record_id", literal=str),
                TypePattern("version", literal=int, optional=True)),
        result=RCType("GadgetPublicationBytes", kind=ValueKind.BYTES,
                      domain=Domain.BYTES, exactness=Exactness.STRUCTURAL),
        implementation_key="gadget.native-publication.bytes.v1",
        requires=frozenset({"identity", "records"}), source_methods=frozenset({"read_record"}))
    if _OPERATOR_NAME in catalogued():
        if lookup(_OPERATOR_NAME) != signature:
            raise ValueError("RCQL operator name already has a different declaration")
        return _OPERATOR_NAME

    def retrieve(source, record_id, version=None):
        snapshot = read_publication(source, record_id, version=version)
        record_method("gadget-native-publication-read", record_id=snapshot.record.id,
                      record_version=snapshot.record.version, state_digest=snapshot.state_digest)
        return snapshot.value

    register_signature(signature)
    register_operator(_OPERATOR_NAME)(retrieve)
    return _OPERATOR_NAME


def execute_publication(workspace, request, *, gadget):
    """Execute an assay request with binary input and output publications."""
    binary = str(Path(gadget).resolve(strict=True))
    for raw in (workspace, request):
        if type(raw) is not bytes or not raw.startswith(_MAGIC):
            raise ValueError("GADGET publication requires RFCXS bytes")
    with TemporaryDirectory(prefix="gadget-native-query-") as folder:
        workspace_path = Path(folder)/"workspace.rfcxs"
        request_path = Path(folder)/"request.rfcxs"
        output = Path(folder)/"result.rfcxs"
        workspace_path.write_bytes(workspace)
        request_path.write_bytes(request)
        _validate_publication_file(workspace_path, gadget=binary, expected=workspace)
        _validate_publication_file(request_path, gadget=binary, expected=request)
        process = subprocess.run([binary, "assay-query", "--workspace", str(workspace_path),
                                  "--request", str(request_path), "--out-record", str(output)],
                                 capture_output=True)
        if process.returncode:
            raise ValueError("GADGET native query failed: " +
                             process.stderr.decode("utf-8", errors="replace").strip())
        if process.stdout:
            raise ValueError("GADGET binary query produced an unexpected display response")
        return _read_publication_file(output, gadget=binary)


def install_rcql_executor(*, gadget):
    """Register binary assay execution on explicitly selected RCDB versions."""
    import hashlib
    from rcql import Domain, Exactness, OperatorSignature, RCType, TypePattern, ValueKind
    from rcql.types import Effect
    from rcql.execution_trace import record_method
    from rcql.operators import register as register_operator
    from rcql.signatures import register as register_signature
    binary = str(Path(gadget).resolve(strict=True))
    name = "GADGET_ASSAY_QUERY"
    signature = OperatorSignature(
        name=name, source_kind=ValueKind.RCDB_STORE,
        inputs=(TypePattern("workspace_id", literal=str), TypePattern("workspace_version", literal=int),
                TypePattern("request_id", literal=str), TypePattern("request_version", literal=int)),
        result=RCType("GadgetOperationBytes", kind=ValueKind.BYTES,
                      domain=Domain.BYTES, exactness=Exactness.STRUCTURAL),
        implementation_key="gadget.native-assay-query.v1",
        requires=frozenset({"identity", "records", "files", "gadget_execute"}),
        effects=frozenset({Effect.READ, Effect.FILESYSTEM}),
        source_methods=frozenset({"read_record"}))

    def execute(source, workspace_id, workspace_version, request_id, request_version):
        workspace = read_publication(source, workspace_id, version=workspace_version)
        request = read_publication(source, request_id, version=request_version)
        result = execute_publication(workspace.value, request.value, gadget=binary)
        record_method("gadget-native-assay-query", workspace_id=workspace.record.id,
                      workspace_version=workspace.record.version, request_id=request.record.id,
                      request_version=request.record.version, result_sha256=hashlib.sha256(result).hexdigest())
        return result

    register_signature(signature)
    register_operator(name)(execute)
    return name


def install_format(*, gadget):
    """Register RFCXS admission with RexGraph's existing file-format registry."""
    from rexgraph.io import register_format
    binary = str(Path(gadget).resolve(strict=True))

    def load(path):
        return _read_publication_file(path, gadget=binary)

    return register_format("gadget-native", load=load, extensions=(".rfcxs",))


def select_publication_field(publication, request, *, selector, gadget):
    """Select an exact original field owner with the native C++ field provider."""
    binary = str(Path(selector).resolve(strict=True))
    for raw in (publication, request):
        if type(raw) is not bytes or not raw.startswith(_MAGIC):
            raise ValueError("GADGET field selection requires RFCXS bytes")
    with TemporaryDirectory(prefix="gadget-native-field-") as folder:
        source = Path(folder) / "source.rfcxs"
        control = Path(folder) / "request.rfcxs"
        output = Path(folder) / "selection.rfcxs"
        source.write_bytes(publication)
        control.write_bytes(request)
        _validate_publication_file(source, gadget=gadget, expected=publication)
        _validate_publication_file(control, gadget=gadget, expected=request)
        process = subprocess.run([binary, "--input", str(source), "--request", str(control),
                                  "--out", str(output)], capture_output=True)
        if process.returncode:
            raise ValueError("GADGET native field selection failed: " +
                             process.stderr.decode("utf-8", errors="replace").strip())
        if process.stdout:
            raise ValueError("GADGET native field selection produced an unexpected display")
        return _read_publication_file(output, gadget=gadget)


def install_rcql_field_selector(*, selector, gadget):
    """Register original field selection at explicit versions, returning a native owner."""
    import hashlib
    from rcql import Domain, Exactness, OperatorSignature, RCType, TypePattern, ValueKind
    from rcql.types import Effect
    from rcql.execution_trace import record_method
    from rcql.operators import register as register_operator
    from rcql.signatures import catalogued, lookup, register as register_signature
    binary = str(Path(selector).resolve(strict=True))
    native_reader = str(Path(gadget).resolve(strict=True))
    name = "GADGET_NATIVE_FIELD"
    signature = OperatorSignature(
        name=name, source_kind=ValueKind.RCDB_STORE,
        inputs=(TypePattern("publication_id", literal=str),
                TypePattern("publication_version", literal=int),
                TypePattern("request_id", literal=str), TypePattern("request_version", literal=int)),
        result=RCType("GadgetNativeFieldSelectionBytes", kind=ValueKind.BYTES,
                      domain=Domain.BYTES, exactness=Exactness.STRUCTURAL),
        implementation_key="gadget.native-field-selection.v1",
        requires=frozenset({"identity", "records", "files", "gadget_field_select"}),
        effects=frozenset({Effect.READ, Effect.FILESYSTEM}),
        source_methods=frozenset({"read_record"}))
    if name in catalogued():
        if lookup(name) != signature:
            raise ValueError("RCQL operator name already has a different declaration")
        return name

    def select(source, publication_id, publication_version, request_id, request_version):
        publication = read_publication(source, publication_id, version=publication_version)
        request = read_publication(source, request_id, version=request_version)
        result = select_publication_field(publication.value, request.value,
                                          selector=binary, gadget=native_reader)
        record_method("gadget-native-field-selection", publication_id=publication.record.id,
                      publication_version=publication.record.version,
                      source_state_digest=publication.state_digest,
                      request_id=request.record.id, request_version=request.record.version,
                      request_state_digest=request.state_digest,
                      result_sha256=hashlib.sha256(result).hexdigest())
        return result

    register_signature(signature)
    register_operator(name)(select)
    return name
