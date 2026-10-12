"""Lossless RCDB compression, native admission and file-budget contracts."""
from contextlib import closing
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

from rcdb import LocalStore
from rcdb.core import decompress_blob
from rcdb.envelope import RecordEnvelope
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS

spec = importlib.util.spec_from_file_location(
    "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools"/"rcdb"/"gadget_rcdb.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
binary, publication, root = (Path(p).resolve() for p in sys.argv[1:4])
raw = publication.read_bytes()
codec = adapter.install_codec(gadget=binary, max_bytes=len(raw))


def rejects(action):
    try:
        action()
    except (ValueError, FileExistsError):
        return
    raise AssertionError("Expected refusal")


if len(sys.argv) > 4 and sys.argv[4] == "--reopen":
    for name in ("zstd", "zlib"):
        with closing(LocalStore(root/(name+".rcdb"))) as store:
            assert store.header.compression.name == name
            assert adapter.read_publication(store, "coverage", version=1).value == raw
    print("Fresh-process compressed RCDB replay PASS")
    sys.exit(0)

root.mkdir(parents=True, exist_ok=True)
from rexgraph import io
adapter.install_format(gadget=binary)
adapter.install_compressed_format(gadget=binary, max_encoded_bytes=len(raw),
                                  max_output_bytes=len(raw))

replacement_control = root / "replacement.json"
replacement_control.write_text('{"purpose":"canonical replacement control"}')
replacement_archive = root / "replacement.rfcxs"
subprocess.run([str(binary), "native-record", "--input", str(replacement_control),
                "--out", str(replacement_archive)], check=True)
replacement = replacement_archive.read_bytes()
source_copy = root / "changing-source.rfcxs"
original_run = subprocess.run
replacements = []


def replace_source(command, **kwargs):
    if command[1] == "native-record":
        path = Path(command[command.index("--input") + 1])
        if path == source_copy:
            path.write_bytes(replacement)
            replacements.append(path)
    return original_run(command, **kwargs)


for action in [lambda: io.load(source_copy),
               lambda: adapter.compress_publication(source_copy, root / "changed.rcblob",
                                                     gadget=binary, codec="zlib")]:
    source_copy.write_bytes(raw)
    subprocess.run = replace_source
    try:
        rejects(action)
    finally:
        subprocess.run = original_run
assert len(replacements) == 2 and not (root / "changed.rcblob").exists()
assert publication.read_bytes() == raw
for name in ("zstd", "zlib"):
    path = root/(name+".rfcxs.rcblob")
    assert adapter.compress_publication(publication, path, gadget=binary, codec=name) == path
    encoded = path.read_bytes()
    assert len(encoded) < len(raw)
    assert io.load(path) == raw
    budget = dict(gadget=binary, max_encoded_bytes=len(encoded), max_output_bytes=len(raw))
    assert adapter.read_compressed_publication(path, **budget) == raw
    assert adapter.read_compressed_publication(
        path, **{**budget, "max_encoded_bytes": 2**80}) == raw
    rejects(lambda: adapter.compress_publication(publication, path, gadget=binary, codec=name))
    assert path.read_bytes() == encoded
    rejects(lambda: adapter.read_compressed_publication(
        path, **{**budget, "max_encoded_bytes": len(encoded)-1}))
    rejects(lambda: adapter.read_compressed_publication(
        path, **{**budget, "max_output_bytes": len(raw)-1}))
    for value in (0, -1, True, "100"):
        rejects(lambda: adapter.read_compressed_publication(
            path, **{**budget, "max_output_bytes": value}))
        rejects(lambda: adapter.read_compressed_publication(
            path, **{**budget, "max_encoded_bytes": value}))
    for bad in (encoded[:-1], encoded+b"trailing", encoded[:-1]+bytes([encoded[-1]^1])):
        damaged = root/(name+"-damaged.rcblob")
        damaged.write_bytes(bad)
        rejects(lambda: adapter.read_compressed_publication(
            damaged, **{**budget, "max_encoded_bytes": len(bad)}))
    with closing(LocalStore(root/(name+".rcdb"), compression=BlobCodecSpec(name=name),
                            record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
        adapter.put_publication(store, "coverage", publication, expected_version=0)
        assert adapter.read_publication(store, "coverage", version=1).value == raw
    blobs = list((root/(name+".rcdb")/"blobs").iterdir())
    assert blobs
    encoded_size = sum(p.stat().st_size for p in blobs)
    for blob in blobs:
        envelope, payload = RecordEnvelope.from_bytes(blob.read_bytes())
        assert envelope.codec == codec.name
        assert decompress_blob(payload, max_output_bytes=len(raw)) == raw
    assert encoded_size < len(raw)

rejects(lambda: adapter.compress_publication(publication, root/"invalid.rcblob",
                                            gadget=binary, codec="unknown"))
assert not (root/"invalid.rcblob").exists()
process = subprocess.run([sys.executable, __file__, str(binary), str(publication), str(root),
                          "--reopen"], capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)})
assert process.returncode == 0, process.stderr
print(process.stdout.strip())
print("RCDB compression, unchanged native bytes, format registry and malformed/budget refusals PASS")
