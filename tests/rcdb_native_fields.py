"""Exact source/cutoff/carrier and complete native field selection contracts."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import resource as resource_limits
import shutil
import signal
import subprocess
import sys
import threading

selector, fixture, verifier, folder = map(Path, sys.argv[1:])
shutil.rmtree(folder, ignore_errors=True)
folder.mkdir(parents=True)


def run(*command, succeeds=True):
    process = subprocess.run(list(map(str, command)), capture_output=True, text=True)
    assert (process.returncode == 0) == succeeds, (command, process.stdout, process.stderr)
    return process


run(fixture, folder / "original")
source = folder / "original/data.bin"
original = source.read_bytes()
base = dict(schema="gadget.rcdb.field.query.v1", source_namespace="811", cutoff=4, tensor=888,
            input_carrier=None, output_carrier=None, interpretation="coefficient selection fixture")
execution = dict(name="native-field-contract", records=100, rational_slots=100,
                 sparse_entries=100, scalar_work=1000, integer_bits=1024)
requests = [base, {**base, "window": dict(lower=[0], upper=[1]), "execution": execution},
            {**base, "window": dict(lower=[1], upper=[2]), "execution": execution},
            {**base, "window": dict(lower=[2], upper=[3]), "execution": execution},
            {**base, "window": dict(lower=[1], upper=[3]), "execution": execution},
            {**base, "window": dict(lower=[0], upper=[3]),
             "execution": {**execution, "sparse_entries": 0}},
            {**base, "tensor": 891, "output_carrier": 881, "interpretation": ""}]
for i, request in enumerate(requests):
    control, output = folder / f"request-{i}.json", folder / f"selection-{i}.rfcxs"
    control.write_text(json.dumps(request))
    run(selector, "--input", source, "--request", control, "--out", output)
    run(verifier, source, output, control)
    before = output.read_bytes()
    run(selector, "--input", source, "--request", control, "--out", output, succeeds=False)
    assert output.read_bytes() == before

bad = [{**base, "source_namespace": "812"}, {**base, "cutoff": 3}, {**base, "cutoff": 5},
       {**base, "tensor": 999}, {**base, "input_carrier": 881}, {**base, "output_carrier": 881},
       {**base, "interpretation": "claimed foreign units"}, {**base, "tensor": 891,
       "output_carrier": 882, "interpretation": ""}, {**base, "cutoff": 4.0},
       {**base, "ignored": True}, {**base, "window": dict(lower=[0], upper=[3])},
       {**base, "window": dict(lower=[0], upper=[4]), "execution": execution},
       {**base, "window": dict(lower=[2], upper=[1]), "execution": execution},
       {**base, "window": dict(lower=[], upper=[]), "execution": execution}]
for resource in ["records", "rational_slots", "scalar_work", "integer_bits"]:
    bad.append({**base, "window": dict(lower=[0], upper=[3]),
                "execution": {**execution, resource: 0}})
for i, request in enumerate(bad):
    control, output = folder / f"bad-{i}.json", folder / f"refused-{i}.rfcxs"
    control.write_text(json.dumps(request))
    refused = run(selector, "--input", source, "--request", control, "--out", output, succeeds=False)
    assert not output.exists() and refused.stderr.strip()

# Native CreateNew publication elects one winner without replacing its bytes.
control = folder / "race-request.json"
control.write_text(json.dumps(base))
for attempt in range(3):
    output = folder / f"race-{attempt}.rfcxs"
    ready = threading.Barrier(2)

    def publish():
        ready.wait()
        return subprocess.run([str(selector), "--input", str(source), "--request", str(control),
                               "--out", str(output)], capture_output=True, text=True)

    with ThreadPoolExecutor(max_workers=2) as writers:
        futures = [writers.submit(publish) for _ in range(2)]
        outcomes = [future.result() for future in futures]
    assert sum(process.returncode == 0 for process in outcomes) == 1
    assert all(process.returncode == 0 or process.stderr.strip() for process in outcomes)
    run(verifier, source, output, control)
    published = output.read_bytes()
    run(selector, "--input", source, "--request", control, "--out", output, succeeds=False)
    assert output.read_bytes() == published

for name, target in [("existing-link", source), ("dangling-link", folder / "never-created")]:
    output = folder / name
    output.symlink_to(target)
    run(selector, "--input", source, "--request", control, "--out", output, succeeds=False)
    assert output.is_symlink() and output.readlink() == target
    assert source.read_bytes() == original
    if name == "dangling-link":
        assert not target.exists()
output = folder / "source-hard-link"
output.hardlink_to(source)
run(selector, "--input", source, "--request", control, "--out", output, succeeds=False)
assert output.read_bytes() == original

def exhaust_file_bytes():
    signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    resource_limits.setrlimit(resource_limits.RLIMIT_FSIZE, (64, 64))

output = folder / "failed-publication.rfcxs"
staging_before = set(folder.glob(".rfc-stage-*"))
failed = subprocess.run([str(selector), "--input", str(source), "--request", str(control),
                         "--out", str(output)], capture_output=True, text=True,
                        preexec_fn=exhaust_file_bytes)
assert failed.returncode != 0 and failed.stderr.strip() and not output.exists()
assert set(folder.glob(".rfc-stage-*")) == staging_before
run(selector, "--input", source, "--request", control, "--out", output)
run(verifier, source, output, control)
assert source.read_bytes() == original
print("Native field selection, original owners, support, exact carriers, writer races, aliases and failed publication PASS")
