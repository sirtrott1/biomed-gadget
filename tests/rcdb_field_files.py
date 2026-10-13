"""Native field selection retains complete file owners without archive buffers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tracemalloc
from unittest.mock import patch


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gadget", "selector", "publisher", "verifier", "fixture", "source", "root"):
        parser.add_argument(name, type=Path)
    parser.add_argument("tensor", type=int)
    parser.add_argument("--whole", action="store_true")
    args = parser.parse_args()
    binary, selector, publisher, verifier, fixture, source, root = (
        getattr(args, name).resolve() for name in
        ("gadget", "selector", "publisher", "verifier", "fixture", "source", "root"))
    root.mkdir(parents=True)
    before = digest(source)
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    requests = root/"requests"
    process = subprocess.run([str(fixture), str(source), str(args.tensor), str(requests)],
                             capture_output=True, text=True)
    assert process.returncode == 0, (process.stdout, process.stderr)
    declaration = json.loads(process.stdout)
    inputs = {str(path): digest(path) for path in requests.iterdir()}
    budgets = dict(selector=selector, gadget=binary, publisher=publisher,
                   max_source_bytes=source.stat().st_size,
                   max_request_bytes=max(path.stat().st_size for path in requests.iterdir()),
                   max_output_bytes=source.stat().st_size+16*1024*1024)
    comparisons, refusals = [], []

    def refuses(name, action, error=ValueError):
        try:
            action()
        except error as caught:
            refusals.append({"case": name, "message": str(caught)})
        else:
            raise AssertionError("Expected file field selection refusal: "+name)

    for name in ("full", "window"):
        request = requests/(name+".rfcxs")
        output = root/(name+"-selection.rfcxs")
        tracemalloc.start()
        assert adapter.select_publication_field_file(source, request, output, **budgets) == output
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        if args.whole:
            assert peak < 8*1024*1024, peak
        process = subprocess.run([str(verifier), str(source), str(output), str(request)],
                                 capture_output=True, text=True)
        assert process.returncode == 0, (process.stdout, process.stderr)
        original = digest(output)
        refuses("existing-"+name, lambda: adapter.select_publication_field_file(
            source, request, output, **budgets), FileExistsError)
        assert digest(output) == original
        comparisons.append({"case": name, "output_bytes": output.stat().st_size,
                            "output_sha256": original, "complete_native_verification": process.stdout.strip(),
                            "python_peak_traced_allocation_bytes": peak})
    full = requests/"full.rfcxs"
    for budget in ("max_source_bytes", "max_request_bytes", "max_output_bytes"):
        for value in (0, -1, True, "100"):
            refuses(budget+"-"+str(value), lambda: adapter.select_publication_field_file(
                source, full, root/"invalid-budget.rfcxs", **{**budgets, budget: value}))
    assert not (root/"invalid-budget.rfcxs").exists()
    for budget, size in (("max_source_bytes", source.stat().st_size),
                         ("max_request_bytes", full.stat().st_size)):
        output = root/(budget+"-exhausted.rfcxs")
        refuses(budget, lambda: adapter.select_publication_field_file(
            source, full, output, **{**budgets, budget: size-1}))
        assert not output.exists()
    if not args.whole:
        output = root/"output-budget-exhausted.rfcxs"
        refuses("max_output_bytes", lambda: adapter.select_publication_field_file(
            source, full, output, **{**budgets, "max_output_bytes": 1}))
        assert not output.exists()
        ordinary = subprocess.run

        def display(command, **keywords):
            answer = ordinary(command, **keywords)
            if Path(command[0]).resolve() == selector and answer.returncode == 0:
                answer.stdout = b"unexpected display"
            return answer

        with patch.object(subprocess, "run", display):
            refuses("unexpected-display", lambda: adapter.select_publication_field_file(
                source, full, root/"display.rfcxs", **budgets))
        assert not (root/"display.rfcxs").exists()
        changed = []

        def mutate_request(command, **keywords):
            if Path(command[0]).resolve() == selector:
                captured = Path(command[command.index("--request")+1])
                captured.write_bytes((requests/"window.rfcxs").read_bytes())
                changed.append(str(captured))
            return ordinary(command, **keywords)

        with patch.object(subprocess, "run", mutate_request):
            refuses("changed-captured-request", lambda: adapter.select_publication_field_file(
                source, full, root/"changed-request.rfcxs", **budgets))
        assert len(changed) == 1 and not (root/"changed-request.rfcxs").exists()
        winner = root/"winner.rfcxs"

        def select():
            try:
                adapter.select_publication_field_file(source, full, winner, **budgets)
            except FileExistsError:
                return False
            return True

        with ThreadPoolExecutor(max_workers=2) as workers:
            outcomes = list(workers.map(lambda _: select(), range(2)))
        assert sum(outcomes) == 1
        subprocess.run([str(verifier), str(source), str(winner), str(full)], check=True)
        failed = root/"write-exhausted.rfcxs"

        def limit_output():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))

        def write_failure(command, **keywords):
            if Path(command[0]).resolve() == publisher and Path(command[2]) == failed:
                keywords["preexec_fn"] = limit_output
            return ordinary(command, **keywords)

        with patch.object(subprocess, "run", write_failure):
            refuses("real-output-byte-exhaustion", lambda: adapter.select_publication_field_file(
                source, full, failed, **budgets), subprocess.CalledProcessError)
        assert not failed.exists() and not list(root.glob(".rfc-stage-*"))
        adapter.select_publication_field_file(source, full, failed, **budgets)
        subprocess.run([str(verifier), str(source), str(failed), str(full)], check=True)
    assert digest(source) == before and all(digest(Path(path)) == sha for path, sha in inputs.items())
    receipt = {"status": "PASS", "declaration": declaration, "source_bytes": source.stat().st_size,
               "source_sha256": before, "comparisons": comparisons, "refusals": refusals,
               "all_original_inputs_unchanged": True,
               "scope": "Complete native field/source/cutoff/carrier/support verification using existing provider; bounded Python file flow, no source partition or numerical rendering import"}
    (root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("Complete source/native window field file selection and original byte ownership PASS", flush=True)


if __name__ == "__main__":
    main()
