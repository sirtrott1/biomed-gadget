"""Complete native publications retain versions across installed local backends."""
import argparse
from contextlib import closing
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

from rcdb import (BlobCodecSpec, FileStore, LocalStore, MemoryStore, NativeObjectStore,
                  SQLStore, VersionConflictError, copy_record)
from rcdb.header import DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gadget", type=Path)
    parser.add_argument("publication", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--compression", choices=("none", "zlib", "zstd"), default="zstd")
    parser.add_argument("--backend", action="append",
                        choices=("memory", "local", "sql", "object-file", "object-memory"))
    parser.add_argument("--result-max-bytes", type=int)
    parser.add_argument("--reopen", action="store_true")
    args = parser.parse_args()
    if args.result_max_bytes is not None and args.result_max_bytes <= 0:
        parser.error("portable result byte budget must be positive")
    binary, publication, root = (p.resolve() for p in
                                 (args.gadget, args.publication, args.root))
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    raw = publication.read_bytes()
    control = root/"different.rfcxs"
    if not args.reopen:
        root.mkdir(parents=True, exist_ok=False)
        request = root/"different.json"
        request.write_text('{"purpose":"distinct current record version"}')
        subprocess.run([str(binary), "native-record", "--input", str(request),
                        "--out", str(control)], check=True, capture_output=True)
    different = control.read_bytes()
    assert different != raw
    reference = adapter.install_codec(gadget=binary, max_bytes=max(len(raw), len(different)))
    adapter.install_rcql_operator()
    kinds = args.backend or ["memory", "local", "sql", "object-file", "object-memory"]
    record_id = "publication@1"
    metadata = {"exact": Fraction(2, 7), "description": "native source"}

    def opened(kind, *, create=False):
        options = {"compression": BlobCodecSpec(args.compression),
                   "record_codecs": (*DEFAULT_RECORD_CODECS, reference)} if create else {}
        path = root/kind
        if kind == "memory":
            return MemoryStore(**options)
        if kind == "local":
            return LocalStore(path, **options)
        if kind == "sql":
            return SQLStore(f"sqlite:///{path}.sqlite", **options)
        prefix = "file://" if kind == "object-file" else "memory://"
        return NativeObjectStore(prefix+str(path), **options)

    def check(store):
        first = adapter.read_publication(store, record_id, version=1)
        assert first.value == raw and first.record.meta == metadata
        assert first.state_digest == first.record.envelope.object_digest
        current = adapter.read_publication(store, record_id)
        assert current.record.version == 3 and current.value == different
        historical = adapter.read_publication(store, record_id, as_of=15.0, valid_at=4.0)
        assert historical.value == raw and historical.record.version == 1
        result = run(query(source("db"), call("GADGET_PUBLICATION_BYTES", record_id, 1)),
                     sources={"db": BoundSource(store, SourcePolicy.allow("identity", "records"))})
        assert result.values == (raw,)
        method = result.execution[0]["methods"][0]
        assert (method["record_id"], method["record_version"], method["state_digest"]) == (
            record_id, 1, first.state_digest)
        limits = {} if args.result_max_bytes is None else {"max_bytes": args.result_max_bytes}
        portable = unpack_result(pack_result(result, **limits), **limits)
        assert portable.values == result.values and portable.execution == result.execution
        return first

    if args.reopen:
        receipts = json.loads((root/"receipt.json").read_text())
        for kind in kinds:
            if kind in ("memory", "object-memory"):
                continue
            with closing(opened(kind)) as store:
                check(store)
                before = receipts["backends"][kind]
                assert store.store_id == before["store_id"]
                assert store.state_digest() == before["state_digest"]
                assert store.header.compression.name == args.compression
        print("Fresh-process complete native backend histories and typed query replay PASS")
        return

    receipts = {}
    for kind in kinds:
        with closing(opened(kind, create=True)) as store:
            first = adapter.put_publication(store, record_id, publication, expected_version=0,
                                             meta=metadata, tx_time=10.0, valid_from=1.0, valid_to=9.0)
            assert first.version == 1
            store.put_record(record_id, different, codec=reference, expected_version=1, tx_time=20.0)
            before = store.state_digest()
            try:
                adapter.put_publication(store, record_id, publication, expected_version=0)
            except VersionConflictError:
                pass
            else:
                raise AssertionError("A stale writer must refuse")
            assert store.state_digest() == before
            assert store.delete(record_id, tx_time=30.0)
            assert store.read_record(record_id) is None
            assert adapter.read_publication(store, record_id, version=1).value == raw
            store.put_record(record_id, different, codec=reference, expected_version=0, tx_time=40.0)
            historical = check(store)
            with closing(MemoryStore(record_codecs=(*DEFAULT_RECORD_CODECS, reference))) as target:
                receipt = copy_record(store, target, historical.record,
                                      expected_cursor=target.change_cursor,
                                      expected_digest=historical.state_digest,
                                      return_receipt=True, destination_id="copied@1")
                copied = adapter.read_publication(target, "copied@1", version=1)
                assert copied.value == raw and copied.record.meta == metadata
                assert (receipt.source_store_id, receipt.source_record_id, receipt.source_version,
                        receipt.source_digest) == (store.store_id, record_id, 1, historical.state_digest)
                assert receipt.destination_digest == copied.state_digest
            before = store.state_digest()
            if kind in ("local", "object-file", "object-memory"):
                store.checkpoint()
            else:
                try:
                    store.checkpoint()
                except NotImplementedError:
                    pass
                else:
                    raise AssertionError("Undeclared checkpoint capability")
            assert store.state_digest() == before
            receipts[kind] = {"store_id": store.store_id, "state_digest": before,
                              "historical_digest": historical.state_digest}
        print(kind, args.compression, "native bytes/versions/tombstone/copy/query and checkpoint capability PASS", flush=True)
    with closing(FileStore(str(root/"legacy"), read_only=False)) as legacy:
        try:
            legacy.put_record(record_id, raw, codec=reference, expected_version=0)
        except ValueError:
            pass
        else:
            raise AssertionError("Legacy storage must refuse unsupported generic native records")
        assert legacy.list() == []
    result = {"status": "PASS", "compression": args.compression, "backends": receipts,
              "source_bytes": len(raw), "source_sha256": hashlib.sha256(raw).hexdigest(),
              "legacy_generic_record_refused": True,
              "scope": "Installed local providers; memory backends persist only in this process"}
    (root/"receipt.json").write_text(json.dumps(result, indent=2)+"\n")
    subprocess.run([sys.executable, "-I", "-B", str(Path(__file__).resolve()),
                    *sys.argv[1:], "--reopen"], check=True)
    assert publication.read_bytes() == raw


if __name__ == "__main__":
    main()
