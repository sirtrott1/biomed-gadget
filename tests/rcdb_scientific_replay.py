"""Explicit stored versions replay complete native scientific source and fields."""
from contextlib import closing
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

from rcdb import LocalStore, RecordCodec
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result

binary, verifier, original, replacement, root = (Path(p).resolve() for p in sys.argv[1:6])
compression = sys.argv[6] if len(sys.argv) > 6 else "none"
spec = importlib.util.spec_from_file_location(
    "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
raw, different = original.read_bytes(), replacement.read_bytes()
assert raw != different
root.mkdir(parents=True)


def refuses(action):
    try:
        action()
    except (ValueError, KeyError, PermissionError, TypeError):
        return
    raise AssertionError("Expected scientific replay refusal")


def inventory(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def verify(expected, replayed):
    path = root/("replay-"+str(len(comparisons))+".rfcxs")
    path.write_bytes(replayed)
    process = subprocess.run([str(verifier), str(expected), str(path)], check=True,
                             capture_output=True, text=True)
    receipt = json.loads(process.stdout)
    assert receipt["status"] == "PASS" and receipt["original_request_and_study_bytes_equal"]
    comparisons.append(receipt)
    return path


comparisons = []
for invalid in (b"invalid", bytearray(raw), raw[:-1]):
    refuses(lambda: adapter.replay_publication(invalid, gadget=binary))
control = root/"control.json"
control.write_text('{"purpose":"native bytes without scientific replay sources"}')
unsupported = root/"control.rfcxs"
subprocess.run([str(binary), "native-record", "--input", str(control), "--out", str(unsupported)], check=True)
refuses(lambda: adapter.replay_publication(unsupported.read_bytes(), gadget=binary))
verify(original, adapter.replay_publication(raw, gadget=binary))
ordinary_run = subprocess.run
replacements = []


def replace_answer(command, **keywords):
    if len(command) > 1 and command[1] == "native-record":
        path = Path(command[command.index("--input") + 1])
        if path.name == "replay.rfcxs":
            path.write_bytes(different)
            replacements.append(str(path))
    return ordinary_run(command, **keywords)


with patch.object(subprocess, "run", replace_answer):
    refuses(lambda: adapter.replay_publication(raw, gadget=binary))
assert len(replacements) == 1


def display_answer(command, **keywords):
    process = ordinary_run(command, **keywords)
    if len(command) > 1 and command[1] == "replay-scientific" and process.returncode == 0:
        process.stdout = b"unexpected display"
    return process


with patch.object(subprocess, "run", display_answer):
    refuses(lambda: adapter.replay_publication(raw, gadget=binary))

codec = adapter.install_codec(gadget=binary, max_bytes=max(len(raw), len(different),
                              RecordCodec.__dataclass_fields__["max_bytes"].default))
name = adapter.install_rcql_scientific_replay(gadget=binary)
assert adapter.install_rcql_scientific_replay(gadget=binary) == name
refuses(lambda: adapter.install_rcql_scientific_replay(gadget=verifier))
record_id = "operation@2/β?one#%"
store_path = root/"scientific.rcdb"
with closing(LocalStore(store_path, compression=BlobCodecSpec(name=compression),
                        record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
    assert adapter.put_publication(store, record_id, original, expected_version=0).version == 1
    assert adapter.put_publication(store, record_id, replacement, expected_version=1).version == 2
    assert adapter.read_publication(store, record_id).value == different
    before = inventory(store_path)
    restricted = BoundSource(store, SourcePolicy.allow("identity", "records", "files"))
    declaration = query(source("db"), call(name, record_id, 1))
    with patch.object(adapter, "replay_publication", side_effect=AssertionError("unadmitted replay executed")):
        refuses(lambda: run(declaration, sources={"db": restricted}))
    admitted = BoundSource(store, SourcePolicy.allow("identity", "records", "files", "gadget_replay"))
    received = []
    replay = adapter.replay_publication

    def observed(value, **keywords):
        received.append(value)
        return replay(value, **keywords)

    with patch.object(adapter, "replay_publication", observed):
        result = run(declaration, sources={"db": admitted})
    assert received == [raw]
    payload = result.values[0]
    assert type(payload) is bytes
    method = result.execution[0]["methods"][0]
    assert method["method"] == "gadget-native-scientific-replay"
    assert method["record_id"] == record_id and method["record_version"] == 1
    assert method["source_state_digest"] == adapter.read_publication(store, record_id, version=1).state_digest
    assert method["result_sha256"] == hashlib.sha256(payload).hexdigest()
    restored = unpack_result(pack_result(result))
    assert restored.values == result.values and restored.execution == result.execution
    verify(original, payload)
    result2 = run(query(source("db"), call(name, record_id, 2)), sources={"db": admitted})
    verify(replacement, result2.values[0])
    for invalid in (0, True, 3):
        refuses(lambda: run(query(source("db"), call(name, record_id, invalid)), sources={"db": admitted}))
    refuses(lambda: run(query(source("db"), call(name, record_id)), sources={"db": admitted}))
    assert inventory(store_path) == before
    assert adapter.read_publication(store, record_id, version=1).value == raw
    output = verify(original, adapter.replay_publication(payload, gadget=binary))
    adapter.put_publication(store, "replayed", output, expected_version=0)
    store.checkpoint(max_frames=1)
with closing(LocalStore(store_path)) as store:
    assert adapter.read_publication(store, record_id, version=1).value == raw
    assert adapter.read_publication(store, record_id, version=2).value == different
    selected = adapter.read_publication(store, "replayed", version=1)
    verify(original, selected.value)
assert original.read_bytes() == raw and replacement.read_bytes() == different
receipt = {"status": "PASS", "compression": compression, "source_bytes": len(raw),
           "source_sha256": hashlib.sha256(raw).hexdigest(), "replacement_sha256": hashlib.sha256(different).hexdigest(),
           "complete_native_comparisons": comparisons, "selected_version_before_replay": 1,
           "replay_repeated_without_extra_original_owners": True, "stored_originals_unchanged": True,
           "capability_refusal_before_replay": True, "result_packet_and_version_digest_receipt": True,
           "scope": "Native provider replay and complete stored tensor/source/request ownership; no Python numerical result import"}
(root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
print("Versioned native scientific replay, capabilities, complete fields/support/source bytes and restart PASS")
