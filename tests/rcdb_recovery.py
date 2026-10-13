"""Native record history and measured local journal/checkpoint write failures."""
import argparse
from contextlib import closing, contextmanager
import errno
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

from rcdb import LocalStore, PublicationUncertainError, VersionConflictError
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from rcdb.journal import LocalJournal, TornJournalError
import rcdb.localstore as localstore


def refuses(action, exception=ValueError):
    try:
        action()
    except exception as error:
        return error
    raise AssertionError("Expected provider refusal")


def inventory(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


@contextmanager
def file_size_limit(size):
    import resource
    import signal
    previous = resource.getrlimit(resource.RLIMIT_FSIZE)
    handler = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    resource.setrlimit(resource.RLIMIT_FSIZE, (size, previous[1]))
    try:
        yield
    finally:
        resource.setrlimit(resource.RLIMIT_FSIZE, previous)
        signal.signal(signal.SIGXFSZ, handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gadget", type=Path)
    parser.add_argument("publication", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--compression", choices=("none", "zstd", "zlib"), default="none")
    parser.add_argument("--original", action="store_true")
    parser.add_argument("--child", choices=("journal-empty", "journal-partial", "checkpoint", "reopen"))
    args = parser.parse_args()
    binary, publication, root = (path.resolve() for path in
                                 (args.gadget, args.publication, args.root))
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    raw = publication.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    codec = adapter.install_codec(gadget=binary, max_bytes=len(raw))
    record_id = "literal@2/β?one#%"
    meta = {"source": "unchanged publication", "exact": Fraction(1, 7),
            "integer": 2**80, "text": "raw\x00metadata"}
    store_path = root/"history.rcdb"

    def publish(store, path, expected, when):
        return adapter.put_publication(store, record_id, path, expected_version=expected,
                                       tx_time=when, valid_from=0., valid_to=100., meta=meta)

    def check_history(store, different, current=3):
        values = {1: raw, 2: different, 3: raw, 4: different}
        versions = list(range(1, current + 1))
        assert [record.version for record in store.history(record_id)] == versions
        for version in versions:
            snapshot = adapter.read_publication(store, record_id, version=version)
            assert snapshot.value == values[version] and snapshot.record.id == record_id
            assert all(snapshot.record.meta[key] == value for key, value in meta.items())
        assert adapter.read_publication(store, record_id).record.version == current
        assert adapter.read_publication(store, record_id).value == values[current]
        assert adapter.read_publication(store, record_id, as_of=15.).value == raw
        assert adapter.read_publication(store, record_id, as_of=25.).value == different
        refuses(lambda: adapter.read_publication(store, record_id, as_of=35.), KeyError)
        assert adapter.read_publication(store, record_id, as_of=45.).value == raw
        assert store.read_record(record_id, version=current + 1) is None
        assert store.change_cursor.sequence == current + 1

    replacement = root/"different.rfcxs"
    if args.child:
        different = replacement.read_bytes()
        with closing(LocalStore(store_path)) as store:
            if args.child == "reopen":
                check_history(store, different)
                print("Fresh process native tombstone/history/checkpoint replay PASS", flush=True)
                return
            before = (store_path/"records.journal").read_bytes()
            checkpoint = (store_path/"replay.checkpoint").read_bytes()
            if args.child == "checkpoint":
                assert publish(store, replacement, 3, 50.).version == 4
                original_stage = localstore.staged_publication
                attempts = []

                @contextmanager
                def limited_stage(path, *positional, **keywords):
                    with original_stage(path, *positional, **keywords) as staged:
                        if Path(path) == store_path/"replay.checkpoint":
                            attempts.append(str(path))
                            with file_size_limit(1):
                                yield staged
                        else:
                            yield staged

                with patch.object(localstore, "staged_publication", limited_stage):
                    failure = refuses(lambda: store.checkpoint(max_frames=1), OSError)
                assert failure.errno == errno.EFBIG and len(attempts) == 1
                assert (store_path/"replay.checkpoint").read_bytes() == checkpoint
                check_history(store, different, current=4)
                outcome = {"phase": args.child, "errno": failure.errno,
                           "previous_checkpoint_unchanged": True, "current_version": 4}
            else:
                original_publish = store._journal.publish
                extra = 0 if args.child == "journal-empty" else 32
                attempts = []

                def limited_publish(frame):
                    attempts.append(frame.sequence)
                    with file_size_limit(len(before) + extra):
                        return original_publish(frame)

                with patch.object(store._journal, "publish", limited_publish):
                    failure = refuses(lambda: publish(store, replacement, 3, 50.),
                                      (OSError, PublicationUncertainError))
                assert len(attempts) == 1
                chain, error = [], failure
                while error is not None:
                    chain.append({"type": type(error).__name__, "message": str(error),
                                  "errno": getattr(error, "errno", None)})
                    error = error.__cause__ or error.__context__
                assert any(error["errno"] == errno.EFBIG for error in chain)
                after = (store_path/"records.journal").read_bytes()
                assert after.startswith(before) and len(after) <= len(before) + extra
                assert (store_path/"replay.checkpoint").read_bytes() == checkpoint
                journal = LocalJournal(store_path/"records.journal", store_id=store.store_id)
                status = journal.inspect()
                if status.torn:
                    refuses(lambda: adapter.read_publication(store, record_id),
                            (TornJournalError, PublicationUncertainError))
                    stored = inventory(store_path)
                    refuses(lambda: LocalStore(store_path), TornJournalError)
                    assert inventory(store_path) == stored
                    assert status.valid_end == len(before) and status.last_sequence == 4
                    repair = journal.repair_torn_tail()
                    assert repair == status and (store_path/"records.journal").read_bytes() == before
                else:
                    assert after == before
                outcome = {"phase": args.child, "failure_chain": chain, "partial_bytes": len(after)-len(before),
                           "explicit_torn_tail_repair": status.torn, "committed_prefix_unchanged": True}
        with closing(LocalStore(store_path)) as store:
            if args.child != "checkpoint":
                check_history(store, different)
                assert publish(store, replacement, 3, 50.).version == 4
            check_history(store, different, current=4)
            store.checkpoint(max_frames=1)
        with closing(LocalStore(store_path)) as store:
            check_history(store, different, current=4)
        assert not list(store_path.rglob(".*.tmp-*"))
        assert publication.read_bytes() == raw
        (root/(args.child+".json")).write_text(json.dumps(outcome, indent=2)+"\n")
        print(args.child, "actual file size failure and complete native retry PASS", flush=True)
        return

    root.mkdir(parents=True)
    control = root/"different.json"
    control.write_text('{"purpose":"distinct native historical control"}')
    subprocess.run([str(binary), "native-record", "--input", str(control), "--out", str(replacement)],
                   check=True)
    different = replacement.read_bytes()
    assert different != raw and len(different) <= len(raw)
    with closing(LocalStore(store_path, compression=BlobCodecSpec(name=args.compression),
                            record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
        assert publish(store, publication, 0, 10.).version == 1
        pinned = adapter.read_publication(store, record_id, version=1)
        assert publish(store, replacement, 1, 20.).version == 2
        refuses(lambda: store.delete(record_id, expected_version=1, tx_time=30.), VersionConflictError)
        assert store.delete(record_id, expected_version=2, tx_time=30.)
        assert store.read_record(record_id) is None
        assert adapter.read_publication(store, record_id, version=1).value == raw
        deleted_prefix = (store_path/"records.journal").read_bytes()
        assert publish(store, publication, 0, 40.).version == 3
        assert pinned.value == raw and pinned.record.version == 1 and pinned.record.tx_to is None
        check_history(store, different)
        store.checkpoint(max_frames=1)
        check_history(store, different)
        committed = (store_path/"records.journal").read_bytes()
        identity = store.store_id

    command = [sys.executable, "-B", __file__, str(binary), str(publication), str(root),
               "--compression", args.compression]
    process = subprocess.run([*command, "--child", "reopen"], capture_output=True, text=True, timeout=120)
    assert process.returncode == 0, process.stderr
    print(process.stdout.strip(), flush=True)
    receipt = {"source_bytes": len(raw), "source_sha256": source_hash, "compression": args.compression,
               "native_versions": [1, 2, 3], "tombstone_gap": True, "pinned_old_snapshot_unchanged": True,
               "complete_native_bytes_and_metadata": True, "fresh_process_checkpoint_replay": True}

    if not args.original:
        if os.name != "posix":
            raise RuntimeError("Measured file size failure contracts require POSIX")
        journal_path = root/"corpus.journal"
        journal = LocalJournal(journal_path, store_id=identity)
        tail = committed[len(deleted_prefix):]
        for size in range(1, len(tail)):
            journal_path.write_bytes(deleted_prefix+tail[:size])
            status = journal.inspect()
            assert status.torn and status.valid_end == len(deleted_prefix) and status.last_sequence == 3
            assert journal.repair_torn_tail() == status
            assert journal_path.read_bytes() == deleted_prefix
        corrupt_store = root/"corrupt-complete.rcdb"
        shutil.copytree(store_path, corrupt_store)
        for offset in (0, 8, 16, len(tail)//2, len(tail)-1):
            damaged = bytearray(committed)
            damaged[len(deleted_prefix)+offset] ^= 1
            journal_path.write_bytes(damaged)
            error = refuses(journal.repair_torn_tail)
            assert not isinstance(error, TornJournalError) and journal_path.read_bytes() == damaged
            (corrupt_store/"records.journal").write_bytes(damaged)
            before = inventory(corrupt_store)
            refuses(lambda: LocalStore(corrupt_store))
            assert inventory(corrupt_store) == before
        receipt.update(exhaustive_incomplete_tail_lengths=len(tail)-1, damaged_complete_frames_refused=5)
        for phase in ("journal-empty", "journal-partial", "checkpoint"):
            destination = root/phase
            destination.mkdir()
            shutil.copytree(store_path, destination/"history.rcdb")
            shutil.copyfile(replacement, destination/"different.rfcxs")
            process = subprocess.run([sys.executable, "-B", __file__, str(binary), str(publication),
                                      str(destination), "--compression", args.compression, "--child", phase],
                                     capture_output=True, text=True, timeout=120)
            assert process.returncode == 0, (phase, process.stdout, process.stderr)
            print(process.stdout.strip(), flush=True)
            receipt[phase] = json.loads((destination/(phase+".json")).read_text())

    assert publication.read_bytes() == raw
    receipt.update(status="PASS", complete_original_source_unchanged=True,
                   storage_bytes={str(path.relative_to(store_path)): path.stat().st_size
                                  for path in store_path.rglob("*") if path.is_file()})
    (root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("Native history, checked journal recovery and complete source bytes PASS", flush=True)


if __name__ == "__main__":
    main()
