"""Complete ordered native plans survive stored versions, copies and replay."""
import argparse
from contextlib import closing
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

from rcdb import LocalStore, RecordCodec, copy_record
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS

COLD_REPLAY = """
from contextlib import closing
import hashlib, importlib.util, json, sys
from pathlib import Path
from rcdb import LocalStore
spec = importlib.util.spec_from_file_location('gadget_rcdb', sys.argv[1])
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
adapter.install_codec(gadget=sys.argv[2], max_bytes=int(sys.argv[8]))
with closing(LocalStore(sys.argv[3])) as store:
    captured = adapter.read_publication(store, sys.argv[4], version=int(sys.argv[5]))
    assert store.store_id == sys.argv[7]
    result = adapter.replay_publication(captured.value, gadget=sys.argv[2])
    with Path(sys.argv[6]).open('xb') as output:
        output.write(result)
    print(json.dumps({'store_id': store.store_id, 'version': captured.record.version,
                     'source_sha256': hashlib.sha256(captured.value).hexdigest()}))
"""


def inventory(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gadget", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("verifier", type=Path)
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--foreign", type=Path)
    parser.add_argument("--copies", action="store_true")
    args = parser.parse_args()
    binary, fixture, verifier, first, second, root = (path.resolve() for path in
        (args.gadget, args.fixture, args.verifier, args.first, args.second, args.root))
    root.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    inputs = [first, second, *([args.foreign.resolve()] if args.foreign else [])]
    before_inputs = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    plans = root/"plans"
    command = [fixture, first, second, plans, *([args.foreign.resolve()] if args.foreign else [])]
    process = subprocess.run(list(map(str, command)), capture_output=True, text=True)
    assert process.returncode == 0, (process.stdout, process.stderr)
    declaration = json.loads(process.stdout)
    assert declaration["status"] == "PASS" and declaration["ordered_steps"] == 2
    forward, reverse = plans/"forward.rfcxs", plans/"reverse.rfcxs"
    originals = {str(path): path.read_bytes() for path in (forward, reverse)}
    assert originals[str(forward)] != originals[str(reverse)]
    baseline = inventory(plans)
    max_bytes = max(RecordCodec.__dataclass_fields__["max_bytes"].default,
                    *(path.stat().st_size for path in plans.glob("*.rfcxs")))

    def cold_replay(store, record_id, version, identity, expected, output):
        before = inventory(store)
        process = subprocess.run([sys.executable, "-I", "-B", "-c", COLD_REPLAY,
                                  str(spec.origin), str(binary), str(store), record_id, str(version),
                                  str(output), identity, str(max_bytes)], check=True,
                                 capture_output=True, text=True, timeout=300)
        selected = json.loads(process.stdout)
        assert selected["source_sha256"] == hashlib.sha256(expected.read_bytes()).hexdigest()
        assert selected["version"] == version
        comparison = subprocess.run([str(verifier), str(expected), str(output)], check=True,
                                    capture_output=True, text=True)
        fields = json.loads(comparison.stdout)
        assert fields["status"] == "PASS" and fields["plans"] == 1 and fields["ordered_steps"] == 2
        assert inventory(store) == before
        return {"selected": selected, "comparison": fields}

    profiles = []
    cold_profiles = []
    for codec in ("none", "zstd", "zlib"):
        output = root/("versioned-"+codec)
        process = subprocess.run([sys.executable, "-B", str(Path(__file__).with_name("rcdb_scientific_replay.py")),
                                  str(binary), str(verifier), str(forward), str(reverse), str(output), codec],
                                 capture_output=True, text=True, timeout=300)
        assert process.returncode == 0, (codec, process.stdout, process.stderr)
        receipt = json.loads((output/"receipt.json").read_text())
        assert receipt["status"] == "PASS"
        assert all(check["plans"] == 1 and check["ordered_steps"] == 2
                   for check in receipt["complete_native_comparisons"])
        profiles.append(receipt)
        store_path = output/"scientific.rcdb"
        with closing(LocalStore(store_path)) as store:
            identity = store.store_id
        for version, expected in ((1, forward), (2, reverse)):
            cold_profiles.append({"codec": codec, "version": version,
                                  **cold_replay(store_path, "operation@2/β?one#%", version, identity,
                                                expected, root/(codec+"-cold-"+str(version)+".rfcxs"))})
        print(codec, "complete ordered plan versions/repeated replay/restart PASS", flush=True)
    codec = adapter.install_codec(gadget=binary, max_bytes=max_bytes)
    source_path = root/"refusals.rcdb"
    refusals = []
    with closing(LocalStore(source_path, record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
        adapter.put_publication(store, "valid", forward, expected_version=0)
        for case in ("missing-step", "changed-order", "changed-id", "changed-source", "mixed-source"):
            path = plans/(case+".rfcxs")
            if not path.exists():
                continue
            adapter.put_publication(store, case, path, expected_version=0)
            before = inventory(source_path)
            raw = adapter.read_publication(store, case, version=1).value
            try:
                adapter.replay_publication(raw, gadget=binary)
            except ValueError as error:
                assert "scientific replay" in str(error)
                refusals.append({"case": case, "message": str(error)})
            else:
                raise AssertionError("Changed native plan passed semantic replay")
            assert inventory(source_path) == before
            assert adapter.read_publication(store, "valid", version=1).value == originals[str(forward)]
    copies = []
    if args.copies:
        for source_codec in ("none", "zstd", "zlib"):
            with closing(LocalStore(root/("versioned-"+source_codec)/"scientific.rcdb")) as source:
                selected = adapter.read_publication(source, "operation@2/β?one#%", version=1)
                source_before = inventory(Path(source.root))
                for destination_codec in ("none", "zstd", "zlib"):
                    destination = root/(source_codec+"-to-"+destination_codec+".rcdb")
                    with closing(LocalStore(destination, compression=BlobCodecSpec(name=destination_codec),
                                            record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as target:
                        copied = copy_record(source, target, selected.record,
                                             expected_digest=selected.state_digest,
                                             expected_cursor=target.change_cursor, preserve_signature=True,
                                             return_receipt=True)
                        captured = adapter.read_publication(target, selected.record.id, version=1)
                        assert captured.value == originals[str(forward)]
                        assert copied.source_digest == selected.state_digest
                        assert copied.source_store_id == source.store_id and copied.destination_store_id == target.store_id
                        target.checkpoint(max_frames=1)
                        identity = target.store_id
                    with closing(LocalStore(destination)) as target:
                        assert target.store_id == identity
                        assert adapter.read_publication(target, selected.record.id, version=1).value == captured.value
                    reopened = cold_replay(destination, selected.record.id, 1, identity, forward,
                                           root/(source_codec+"-to-"+destination_codec+"-replay.rfcxs"))
                    copies.append({"source_codec": source_codec, "destination_codec": destination_codec,
                                   "native_original_bytes_equal": True, **reopened})
                assert inventory(Path(source.root)) == source_before
    assert inventory(plans) == baseline
    assert all(hashlib.sha256(path.read_bytes()).hexdigest() == before_inputs[str(path)] for path in inputs)
    report = {"status": "PASS", "ordered_steps": 2, "input_sha256": before_inputs,
              "plan_bytes": {path.name: path.stat().st_size for path in (forward, reverse)},
              "plan_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (forward, reverse)},
              "profiles": profiles, "cold_profiles": cold_profiles, "refusals": refusals, "copies": copies,
              "whole_original_inputs_unchanged": True,
              "scope": "Existing native ScientificSession plan and replay providers; complete ordered child fields and original request/study owners through canonical record copies"}
    (root/"receipt.json").write_text(json.dumps(report, indent=2)+"\n")
    print("Complete native plans, ordered child owners, semantic refusals and canonical copies PASS", flush=True)


if __name__ == "__main__":
    main()
