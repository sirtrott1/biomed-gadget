"""Native publication copies, conditional writers and local recovery contracts."""
from contextlib import closing
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

from rcdb import CopyReceipt, LocalStore, VersionConflictError, copy_record
from rcdb.engine import ChangeCursor
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from rcdb.journal import LocalJournal, TornJournalError, _block
from rcdb.retention import RetentionPolicy

binary, publication, root = (Path(p).resolve() for p in sys.argv[1:4])
raw = publication.read_bytes()
spec = importlib.util.spec_from_file_location(
    "gadget_rcdb", Path(__file__).resolve().parents[1] / "tools/rcdb/gadget_rcdb.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
codec = adapter.install_codec(gadget=binary, max_bytes=len(raw))
record_id = "original/\u03b1@1"


def inventory(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob("*") if p.is_file()}


def refuses(action, exception=ValueError):
    try:
        action()
    except exception:
        return
    raise AssertionError("Expected provider refusal")


def selected(store):
    return adapter.read_publication(store, record_id, version=1)


def checked_copy(source, destination, **options):
    snapshot = selected(source)
    return copy_record(source, destination, snapshot.record,
                       expected_digest=snapshot.state_digest, preserve_signature=True,
                       return_receipt=True, **options)


if len(sys.argv) > 4 and sys.argv[4] in ("--fault", "--writer"):
    mode, destination = sys.argv[4], Path(sys.argv[5])
    with closing(LocalStore(root / "source-none.rcdb")) as source:
        with closing(LocalStore(destination)) as target:
            if mode == "--writer":
                expected = ChangeCursor.from_record(json.loads(sys.argv[6]))
                print("READY", flush=True)
                assert sys.stdin.read(1) == "g"
                try:
                    receipt = checked_copy(source, target, expected_cursor=expected, tx_time=30.)
                except VersionConflictError:
                    print("CONFLICT", flush=True)
                else:
                    print("COPIED " + str(receipt.destination_version), flush=True)
                sys.exit(0)
            stage = sys.argv[6]
            if stage == "after-blob":
                write = target._write_blob

                def cut(digest, blob):
                    write(digest, blob)
                    os._exit(41)

                target._write_blob = cut
            elif stage == "after-journal":
                publish = target._journal.publish

                def cut(frame):
                    publish(frame)
                    os._exit(42)

                target._journal.publish = cut
            elif stage == "torn-tail":
                def cut(frame):
                    encoded = _block(frame)
                    with target._journal.path.open("ab") as stream:
                        stream.write(encoded[:len(encoded) // 2])
                        stream.flush()
                        os.fsync(stream.fileno())
                    os._exit(43)

                target._journal.publish = cut
            else:
                raise AssertionError("Unknown process cut")
            checked_copy(source, target, expected_cursor=target.change_cursor, tx_time=30.)
    raise AssertionError("Process cut did not execute")

if len(sys.argv) > 4 and sys.argv[4] == "--reopen":
    for name in json.loads((root / "destinations.json").read_text()):
        with closing(LocalStore(root / name)) as store:
            assert selected(store).value == raw
            assert adapter.read_publication(store, "keep", version=1).value == raw
            assert len(store.history(record_id)) == 1
            assert store.change_cursor.sequence == 2
    print("Fresh-process native copy checkpoint replay PASS")
    sys.exit(0)

root.mkdir(parents=True)
control = root / "different.json"
control.write_text('{"purpose":"distinct current native publication"}')
replacement = root / "different.rfcxs"
subprocess.run([str(binary), "native-record", "--input", str(control),
                "--out", str(replacement)], check=True)
different = replacement.read_bytes()
assert different != raw and len(different) <= len(raw)
destinations = []
original_profile = len(sys.argv) > 4 and sys.argv[4] == "--original"
if not original_profile and os.name != "posix":
    raise RuntimeError("Local process publication qualification requires POSIX")
source_codecs = ("none",) if original_profile else ("none", "zstd", "zlib")


def opened(path, compression="none"):
    return LocalStore(path, compression=BlobCodecSpec(name=compression),
                      record_codecs=(*DEFAULT_RECORD_CODECS, codec))


def seed_target(path, compression="zstd"):
    with closing(opened(path, compression)) as store:
        adapter.put_publication(store, "keep", publication, expected_version=0, tx_time=5.)
        store.checkpoint(max_frames=1)
        return store.change_cursor


for source_codec in source_codecs:
    path = root / ("source-" + source_codec + ".rcdb")
    with closing(opened(path, source_codec)) as source:
        adapter.put_publication(source, record_id, publication, expected_version=0, tx_time=10.,
                                valid_from=3., valid_to=12., meta={"source": "original"},
                                tags=["native", "original"])
        adapter.put_publication(source, record_id, replacement, expected_version=1, tx_time=20.)
        assert adapter.read_publication(source, record_id).value == different
        source_before = inventory(path)
        snapshot = selected(source)
        for destination_codec in ("none", "zstd", "zlib"):
            name = source_codec + "-to-" + destination_codec + ".rcdb"
            destination = root / name
            cursor = seed_target(destination, destination_codec)
            with closing(LocalStore(destination)) as target:
                receipt = checked_copy(source, target, expected_cursor=cursor, tx_time=30.,
                                       expected_policy_digest=target.transfer_policy_digest())
                actual = selected(target)
                assert actual.value == raw and actual.value != different
                assert actual.record.meta == snapshot.record.meta
                assert actual.record.signature == snapshot.record.signature
                assert (actual.record.valid_from, actual.record.valid_to) == (3., 12.)
                assert receipt == CopyReceipt.from_bytes(receipt.to_bytes())
                assert (receipt.source_store_id, receipt.source_record_id, receipt.source_version,
                        receipt.source_digest) == (source.store_id, record_id, 1, snapshot.state_digest)
                assert (receipt.destination_store_id, receipt.destination_record_id,
                        receipt.destination_version, receipt.destination_digest) == (
                            target.store_id, record_id, 1, actual.state_digest)
                assert receipt.source_store_id != receipt.destination_store_id
                before = inventory(destination)
                refuses(lambda: checked_copy(source, target, expected_cursor=cursor), VersionConflictError)
                refuses(lambda: copy_record(source, target, snapshot.record, expected_digest="0" * 64))
                refuses(lambda: checked_copy(source, target, expected_policy_digest="0" * 64))
                refuses(lambda: checked_copy(source, target, expected_version=1))
                refuses(lambda: checked_copy(source, target, governed=True))
                for damaged in (receipt.to_bytes()[:-1], receipt.to_bytes() + b"extra"):
                    refuses(lambda: CopyReceipt.from_bytes(damaged))
                assert inventory(destination) == before
                target.checkpoint(max_frames=1)
                assert selected(target).value == raw and target.change_cursor.sequence == 2
                (root / (name + ".receipt")).write_bytes(receipt.to_bytes())
            destinations.append(name)
        assert inventory(path) == source_before

(root / "destinations.json").write_text(json.dumps(destinations))
command = [sys.executable, __file__, str(binary), str(publication), str(root)]
process = subprocess.run([*command, "--reopen"], capture_output=True, text=True)
assert process.returncode == 0, process.stderr
print(process.stdout.strip(), flush=True)

if not original_profile:
    with closing(LocalStore(root / "source-none.rcdb")) as source:
        missing = root / "missing-codec.rcdb"
        with closing(LocalStore(missing)) as target:
            before = inventory(missing)
            refuses(lambda: checked_copy(source, target))
            assert inventory(missing) == before

        destination = root / "exception.rcdb"
        cursor = seed_target(destination)
        with closing(LocalStore(destination)) as target:
            before = inventory(destination)
            with patch.object(target, "_write_blob", side_effect=OSError("controlled write failure")):
                refuses(lambda: checked_copy(source, target, expected_cursor=cursor), OSError)
            assert inventory(destination) == before
            assert checked_copy(source, target, expected_cursor=cursor).destination_version == 1

    for stage, exit_code in (("after-blob", 41), ("after-journal", 42), ("torn-tail", 43)):
        destination = root / (stage + ".rcdb")
        cursor = seed_target(destination)
        journal_before = (destination / "records.journal").read_bytes()
        source_before = inventory(root / "source-none.rcdb")
        process = subprocess.run([*command, "--fault", str(destination), stage],
                                 capture_output=True, text=True)
        assert process.returncode == exit_code, (stage, process.stdout, process.stderr)
        assert inventory(root / "source-none.rcdb") == source_before
        if stage == "torn-tail":
            before = inventory(destination)
            refuses(lambda: LocalStore(destination), TornJournalError)
            assert inventory(destination) == before
            journal = LocalJournal(destination / "records.journal", store_id=cursor.store_id)
            status = journal.inspect()
            assert status.torn and status.valid_end == len(journal_before)
            (root / "torn-tail.bin").write_bytes((destination / "records.journal").read_bytes())
            journal.repair_torn_tail()
            assert (destination / "records.journal").read_bytes() == journal_before
        with closing(LocalStore(root / "source-none.rcdb")) as source:
            with closing(LocalStore(destination)) as target:
                assert adapter.read_publication(target, "keep").value == raw
                if stage == "after-journal":
                    assert selected(target).value == raw
                    refuses(lambda: checked_copy(source, target, expected_cursor=cursor), VersionConflictError)
                else:
                    assert target.read_record(record_id) is None and target.change_cursor == cursor
                    plan = target.plan_retention(RetentionPolicy(grace_seconds=0))
                    assert plan.candidates
                    retained_before = (destination / "records.journal").read_bytes()
                    result = target.apply_retention(plan)
                    assert result["history_retained"] and result["deleted_objects"] == len(plan.candidates)
                    assert (destination / "records.journal").read_bytes() == retained_before
                    assert checked_copy(source, target, expected_cursor=cursor).destination_version == 1
                target.checkpoint(max_frames=1)
                assert selected(target).value == raw and len(target.history(record_id)) == 1
        print(stage, "native copy recovery PASS", flush=True)

    for damage in ("checkpoint", "segment"):
        destination = root / ("damaged-" + damage + ".rcdb")
        shutil.copytree(root / destinations[0], destination)
        path = (destination / "replay.checkpoint" if damage == "checkpoint"
                else next((destination / "replay").iterdir()))
        path.write_bytes(path.read_bytes()[:-1])
        before = inventory(destination)
        refuses(lambda: LocalStore(destination))
        assert inventory(destination) == before

    damaged_source = root / "damaged-source.rcdb"
    shutil.copytree(root / "source-none.rcdb", damaged_source)
    destination = root / "corrupt-source-target.rcdb"
    seed_target(destination)
    with closing(LocalStore(damaged_source)) as source:
        blob = damaged_source / "blobs" / source.changes()[0].mutation.blob_digest
        encoded = blob.read_bytes()
        blob.write_bytes(encoded[:-1] + bytes([encoded[-1] ^ 1]))
        with closing(LocalStore(destination)) as target:
            source_before, before = inventory(damaged_source), inventory(destination)
            record = source.history(record_id)[0]
            refuses(lambda: copy_record(source, target, record, return_receipt=True))
            assert inventory(damaged_source) == source_before and inventory(destination) == before

    destination = root / "writers.rcdb"
    cursor = seed_target(destination)
    children = [subprocess.Popen([*command, "--writer", str(destination),
                                  json.dumps(cursor.as_record())], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) for _ in range(4)]
    try:
        for child in children:
            assert child.stdout.readline().strip() == "READY"
        for child in children:
            child.stdin.write("g")
            child.stdin.flush()
        results = []
        for child in children:
            output, error = child.communicate(timeout=30)
            assert child.returncode == 0, error
            results.append(output.strip())
        assert sorted(results) == ["CONFLICT"] * 3 + ["COPIED 1"]
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
                child.communicate()
    with closing(LocalStore(destination)) as target:
        assert selected(target).value == raw and len(target.history(record_id)) == 1
        assert target.change_cursor.sequence == 2
    print("Four conditional process writers preserve one complete native copy PASS", flush=True)

assert publication.read_bytes() == raw
print("Native copy identities, exact old payload, compression, checkpoints and source bytes PASS", flush=True)
