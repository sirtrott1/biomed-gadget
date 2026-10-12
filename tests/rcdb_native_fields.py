"""Exact source/cutoff/carrier and complete native field selection contracts."""
from pathlib import Path
import json
import shutil
import subprocess
import sys

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
assert source.read_bytes() == original
print("Native field selection, original owners, supported zero/missing, exact carriers and resource refusals PASS")
