"""Lossless RCDB compression, native admission and file-budget contracts."""
from contextlib import closing
import importlib.util
import os
from pathlib import Path
import signal
import subprocess
import sys
from unittest.mock import patch

from rcdb import LocalStore
from rcdb.core import decompress_blob
from rcdb.envelope import RecordEnvelope
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS

spec = importlib.util.spec_from_file_location(
    "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools"/"rcdb"/"gadget_rcdb.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
binary, publication, root = (Path(p).resolve() for p in sys.argv[1:4])
publisher = (Path(sys.argv[sys.argv.index("--publisher") + 1]).resolve()
             if "--publisher" in sys.argv else None)
raw = publication.read_bytes()
codec = adapter.install_codec(gadget=binary, max_bytes=len(raw))


def rejects(action):
    try:
        action()
    except (ValueError, FileExistsError):
        return
    raise AssertionError("Expected refusal")


if len(sys.argv) > 4 and sys.argv[4] == "--nonregular":
    fifo = root / "input.fifo"
    rejects(lambda: adapter.validate_publication_file(fifo, gadget=binary))
    from rexgraph import io
    adapter.install_format(gadget=binary)
    rejects(lambda: io.load(fifo, format="gadget-native"))
    rejects(lambda: adapter.read_compressed_publication(
        fifo, gadget=binary, max_encoded_bytes=len(raw), max_output_bytes=len(raw)))
    with closing(LocalStore(root / "fifo.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
        rejects(lambda: adapter.put_publication(store, "fifo", fifo, expected_version=0))
        assert store.read_record("fifo") is None
    print("Nonregular native validation, file loading, compression and RCDB input refusals PASS")
    sys.exit(0)


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
        if path == source_copy or path.name == "source.rfcxs":
            source_copy.write_bytes(replacement)
            replacements.append(source_copy)
    return original_run(command, **kwargs)


for action in [lambda: io.load(source_copy),
               lambda: adapter.compress_publication(source_copy, root / "changed.rcblob",
                                                     gadget=binary, codec="zlib", publisher=publisher)]:
    source_copy.write_bytes(raw)
    subprocess.run = replace_source
    try:
        rejects(action)
    finally:
        subprocess.run = original_run
assert len(replacements) == 2 and not (root / "changed.rcblob").exists()
assert publication.read_bytes() == raw
rejects(lambda: io.load(root, format="gadget-native"))
rejects(lambda: adapter.validate_publication_file(root, gadget=binary))
rejects(lambda: adapter.read_compressed_publication(
    root, gadget=binary, max_encoded_bytes=len(raw), max_output_bytes=len(raw)))
if hasattr(os, "mkfifo"):
    os.mkfifo(root / "input.fifo")
    process = subprocess.run([sys.executable, __file__, str(binary), str(publication),
                              str(root), "--nonregular"], capture_output=True, text=True,
                             timeout=10, env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)})
    assert process.returncode == 0, process.stderr
    print(process.stdout.strip())
for name in ("zstd", "zlib"):
    path = root/(name+".rfcxs.rcblob")
    assert adapter.compress_publication(publication, path, gadget=binary, codec=name,
                                        publisher=publisher) == path
    encoded = path.read_bytes()
    assert len(encoded) < len(raw)
    assert io.load(path) == raw
    budget = dict(gadget=binary, max_encoded_bytes=len(encoded), max_output_bytes=len(raw))
    assert adapter.read_compressed_publication(path, **budget) == raw
    assert adapter.read_compressed_publication(
        path, **{**budget, "max_encoded_bytes": 2**80}) == raw
    changing = root / (name + "-changing.rcblob")
    original_fstat = os.fstat
    for mutation in ("grow", "truncate"):
        changing.write_bytes(encoded)
        identity = changing.stat()
        altered = []

        def checked_fstat(descriptor):
            info = original_fstat(descriptor)
            if (info.st_dev, info.st_ino) == (identity.st_dev, identity.st_ino) and not altered:
                with changing.open("r+b") as stream:
                    stream.truncate(len(encoded) + 8 * 1024 * 1024 if mutation == "grow"
                                    else len(encoded) - 1)
                altered.append(mutation)
            return info

        with patch.object(os, "fstat", checked_fstat):
            rejects(lambda: adapter.read_compressed_publication(changing, **budget))
        assert altered == [mutation]
        assert path.read_bytes() == encoded and publication.read_bytes() == raw
    rejects(lambda: adapter.compress_publication(publication, path, gadget=binary, codec=name,
                                                publisher=publisher))
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

    if os.name == "posix":
        import resource
        failed = root / (name + "-failed-publication.rcblob")
        actual_publisher = (binary.with_name("gadget_publish_file") if publisher is None else publisher)
        ordinary_run = subprocess.run

        def limited():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            cap = min(1024, len(encoded) - 1)
            resource.setrlimit(resource.RLIMIT_FSIZE, (cap, cap))

        def write_failure(command, **kwargs):
            if Path(command[0]).resolve() == actual_publisher.resolve() and Path(command[2]) == failed:
                kwargs["preexec_fn"] = limited
            return ordinary_run(command, **kwargs)

        with patch.object(subprocess, "run", write_failure):
            try:
                adapter.compress_publication(publication, failed, gadget=binary, codec=name,
                                             publisher=publisher)
            except subprocess.CalledProcessError as error:
                assert error.returncode and any(message in error.stderr for message in
                    (b"file write failed", b"staged output flush failed", b"staged output close failed"))
            else:
                raise AssertionError("Expected actual compressed output file size refusal")
        assert not failed.exists() and not list(root.glob(".rfc-stage-*"))
        adapter.compress_publication(publication, failed, gadget=binary, codec=name, publisher=publisher)
        assert adapter.read_compressed_publication(failed, **budget) == raw
        assert publication.read_bytes() == raw

    altered_output = root / (name + "-altered-publication.rcblob")
    ordinary_run = subprocess.run
    actual_publisher = (binary.with_name("gadget_publish_file") if publisher is None else publisher)

    def alter_publisher_input(command, **kwargs):
        if Path(command[0]).resolve() == actual_publisher.resolve() and Path(command[2]) == altered_output:
            body = kwargs["input"]
            kwargs["input"] = body[:-1] + bytes([body[-1] ^ 1])
        return ordinary_run(command, **kwargs)

    with patch.object(subprocess, "run", alter_publisher_input):
        try:
            adapter.compress_publication(publication, altered_output, gadget=binary, codec=name,
                                         publisher=publisher)
        except ValueError as error:
            assert "differs from the complete supplied bytes" in str(error)
        else:
            raise AssertionError("Changed publisher input escaped final byte comparison")
    assert altered_output.read_bytes() != encoded and publication.read_bytes() == raw

missing = root / "missing-publisher.rcblob"
try:
    adapter.compress_publication(publication, missing, gadget=binary, codec="zlib",
                                 publisher=root / "missing-publisher")
except FileNotFoundError:
    pass
else:
    raise AssertionError("Missing native publisher used a direct-write fallback")
assert not missing.exists()
empty_output = root / "empty-publisher.rcblob"
if Path("/bin/true").is_file():
    try:
        adapter.compress_publication(publication, empty_output, gadget=binary, codec="zlib",
                                     publisher="/bin/true")
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("A successful process with missing output was accepted")
    assert not empty_output.exists()
existing = root / "existing-publisher.rcblob"
existing.write_bytes(b"original destination")
try:
    adapter.compress_publication(publication, existing, gadget=binary, codec="zlib",
                                 publisher=root / "missing-publisher")
except FileExistsError:
    pass
else:
    raise AssertionError("Existing destination was changed before publisher admission")
assert existing.read_bytes() == b"original destination"

rejects(lambda: adapter.compress_publication(publication, root/"invalid.rcblob",
                                            gadget=binary, codec="unknown", publisher=publisher))
assert not (root/"invalid.rcblob").exists()
process = subprocess.run([sys.executable, __file__, str(binary), str(publication), str(root),
                          "--reopen"], capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)})
assert process.returncode == 0, process.stderr
print(process.stdout.strip())
print("RCDB compression, unchanged native bytes, format registry and malformed/budget refusals PASS")
