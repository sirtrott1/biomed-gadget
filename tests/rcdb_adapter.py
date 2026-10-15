"""Native RCDB publication, RCQL retrieval and restart contracts."""
import importlib.util
from contextlib import closing
from unittest.mock import patch
import os
from pathlib import Path
import subprocess
import sys

from rcdb import LocalStore, VALUE_CODEC, VersionConflictError
from rcdb.header import DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result

module_path = Path(__file__).resolve().parents[1]/"tools"/"rcdb"/"gadget_rcdb.py"
spec = importlib.util.spec_from_file_location("gadget_rcdb", module_path)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
binary, publication, root = map(lambda x: Path(x).resolve(), sys.argv[1:4])
raw = publication.read_bytes()
codec = adapter.install_codec(gadget=binary, max_bytes=len(raw))


def rejects(action):
    try:
        action()
    except (ValueError, KeyError, TypeError, PermissionError, VersionConflictError):
        return
    raise AssertionError("Expected refusal")


if len(sys.argv) > 4 and sys.argv[4] == "--reopen":
    with closing(LocalStore(root/"native.rcdb")) as store:
        assert adapter.read_publication(store, "coverage", version=1).value == raw
        assert adapter.read_publication(store, "coverage").record.version == 2
    print("Fresh-process native RCDB replay PASS")
    sys.exit(0)

root.mkdir(parents=True, exist_ok=True)
assert adapter.validate_publication(raw, gadget=binary) == raw
from rexgraph import io
adapter.install_format(gadget=binary)
assert io.load(publication) == raw
rejects(lambda: adapter.validate_publication(b'{"schema":"pretend"}', gadget=binary))
rejects(lambda: adapter.validate_publication(raw[:-1], gadget=binary))
rejects(lambda: adapter.validate_publication(raw[:-1]+bytes([raw[-1]^1]), gadget=binary))
with closing(LocalStore(root/"native.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
    rejects(lambda: store.put_record("truncated", raw[:-1], codec=codec, expected_version=0))
    rejects(lambda: store.put_record("corrupt", raw[:-1]+bytes([raw[-1]^1]),
                                    codec=codec, expected_version=0))
    assert store.read_record("truncated") is None and store.read_record("corrupt") is None
    first = adapter.put_publication(store, "coverage", publication, expected_version=0, tx_time=10.0)
    assert first.version == 1
    second = adapter.put_publication(store, "coverage", publication, expected_version=1, tx_time=20.0)
    assert second.version == 2
    rejects(lambda: adapter.put_publication(store, "coverage", publication, expected_version=0))
    changing = root / "changing-publication.rfcxs"
    original_put = store.put_record
    attempts = []

    def track_put(*args, **kwargs):
        attempts.append(args[0])
        return original_put(*args, **kwargs)

    store.put_record = track_put
    original_fstat, original_stat = os.fstat, Path.stat
    for mutation in ("grow", "truncate"):
        changing.write_bytes(raw)
        identity = changing.stat()
        altered = []

        def change_after_size_check(info):
            if (info.st_dev, info.st_ino) == (identity.st_dev, identity.st_ino) and not altered:
                with changing.open("r+b") as stream:
                    stream.truncate(len(raw) + 8 * 1024 * 1024 if mutation == "grow" else len(raw) - 1)
                altered.append(mutation)
            return info

        def checked_fstat(descriptor):
            return change_after_size_check(original_fstat(descriptor))

        def checked_stat(path, **kwargs):
            return change_after_size_check(original_stat(path, **kwargs))

        with patch.object(os, "fstat", checked_fstat), patch.object(Path, "stat", checked_stat):
            rejects(lambda: adapter.put_publication(store, "coverage", changing,
                                                     expected_version=2))
        assert altered == [mutation]
        assert attempts == [], "changed files must refuse before store publication"
        assert adapter.read_publication(store, "coverage").record.version == 2
        assert adapter.read_publication(store, "coverage", version=1).value == raw
    changing.write_bytes(raw)
    identity = changing.stat()
    observed = []

    def edit_during_capture(descriptor):
        info = original_fstat(descriptor)
        if (info.st_dev, info.st_ino) == (identity.st_dev, identity.st_ino):
            observed.append(info)
            if len(observed) == 2:
                with changing.open("r+b") as stream:
                    stream.seek(len(raw)-1)
                    stream.write(bytes([raw[-1] ^ 1]))
        return info

    with patch.object(os, "fstat", edit_during_capture):
        rejects(lambda: adapter.put_publication(store, "coverage", changing,
                                                 expected_version=2))
    assert len(observed) >= 3 and attempts == []
    assert adapter.read_publication(store, "coverage", version=1).value == raw
    rejects(lambda: adapter.put_publication(store, "coverage", root, expected_version=2))
    assert attempts == []
    store.put_record = original_put
    assert adapter.read_publication(store, "coverage", version=1).value == raw
    assert adapter.read_publication(store, "coverage", as_of=15.0).record.version == 1
    assert adapter.read_publication(store, "coverage").record.version == 2
    store.put_record("other", raw, codec=VALUE_CODEC, expected_version=0)
    rejects(lambda: adapter.read_publication(store, "other"))
    rejects(lambda: adapter.read_publication(store, "absent"))
    assert adapter.install_rcql_operator() == adapter.install_rcql_operator()
    declaration = query(source("db"), call("GADGET_PUBLICATION_BYTES", "coverage", 1))
    result = run(declaration, sources={"db": BoundSource(store, SourcePolicy.allow("identity", "records"))})
    assert result.values == (raw,)
    assert result.execution[0]["methods"][0]["record_version"] == 1
    assert result.execution[0]["methods"][0]["method"] == "gadget-native-publication-read"
    replay = unpack_result(pack_result(result))
    assert replay.values == (raw,)
    assert replay.execution == result.execution
    rejects(lambda: run(declaration, sources={"db": BoundSource(store, SourcePolicy.allow("records"))}))
    current = run(query(source("db"), call("GADGET_PUBLICATION_BYTES", "coverage")), sources={"db": store})
    assert current.execution[0]["methods"][0]["record_version"] == 2

process = subprocess.run([sys.executable, __file__, str(binary), str(publication), str(root), "--reopen"],
                         capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)})
assert process.returncode == 0, process.stderr
print(process.stdout.strip())
print("RCDB unchanged bytes, canonical admission, versions, refusals and typed RCQL receipts PASS")
