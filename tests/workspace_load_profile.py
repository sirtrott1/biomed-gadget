"""Check full workspace admission, exact rewriting and profiling refusals."""
from pathlib import Path
import json
import shutil
import subprocess
import sys

profile, fixture, directory = map(Path, sys.argv[1:])
shutil.rmtree(directory, ignore_errors=True)
directory.mkdir(parents=True)


def run(*args, succeeds=True):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert (result.returncode == 0) == succeeds, (args, result.stdout, result.stderr)
    return result


run(fixture, directory / "fixture")
source = directory / "fixture/derived/workspace.rfcxs"
output = directory / "rewritten.rfcxs"
result = run(profile, "--input", source, "--rewrite", output)
phases, values, bindings = [], {}, {}
for line in result.stdout.splitlines():
    columns = line.split("\t")
    if columns[0] == "phase":
        assert len(columns) == 5 and float(columns[2]) >= 0, columns
        if columns[3] != "unavailable":
            assert int(columns[4]) >= int(columns[3]) > 0, columns
        phases.append(columns[1])
    elif columns[0] == "value":
        assert len(columns) == 3 and columns[1] not in values, columns
        values[columns[1]] = int(columns[2])
    else:
        assert columns[0] == "binding" and len(columns) == 3, columns
        assert columns[1] not in bindings and columns[2].startswith("sha256:"), columns
        bindings[columns[1]] = columns[2]
assert phases[0] == "workspace_record" and phases[-1] == "total"
assert set(["workspace_study_record", "study_prefix", "study_integrity",
            "study_accessions", "study_documents", "study_fields", "study_complete",
            "workspace_complete", "admission", "rewrite", "verify_bytes", "release"]) <= set(phases)
assert values["workspace.assays"] == 1 and values["workspace.partitions"] == 1
assert values["process.cpu_ticks_per_second"] > 0
assert values["workspace.study_endpoint_owners"] == 1
assert values["workspace.raw_owner_bytes"] > 0
assert values["workspace.source_references"] >= values["workspace.unique_byte_owners"] > 0
if sys.platform == "linux":
    assert len([name for name in values if name.startswith("phase.cpu.")]) == len(phases)
assert values["workspace.verified_bytes"] == source.stat().st_size
assert source.read_bytes() == output.read_bytes()
again = run(profile, "--input", output)
assert bindings["workspace.digest"] in again.stdout
for args in [[], ["--input"], ["--unknown", source],
             ["--input", source, "--input", source],
             ["--input", source, "--rewrite", output],
             ["--input", directory / "fixture/import.json"]]:
    run(profile, *args, succeeds=False)
assert source.read_bytes() == output.read_bytes()

request = directory / "query.json"
request.write_text(json.dumps(dict(
    schema="gadget.assays.query.v1", operation="chromatin_contrast",
    samples=[dict(id="base", resource="control", scale="1", scale_provenance="fixture"),
             dict(id="treated", resource="treated", scale="1", scale_provenance="fixture")],
    contrasts=[dict(id="change", treated="treated", control="base")],
    strand_policy="separate", thresholds=["0"], absence_policy="measured_zero_within_union")))
publication = directory / "query.rfcxs"
answer = run(profile, "--input", directory / "fixture/input-away/workspace.rfcxs",
             "--request", request, "--result", publication)
for name in ["request_admission", "execute_record", "account_result", "publication",
             "release_result", "release_request"]:
    assert "phase\t" + name + "\t" in answer.stdout
assert publication.read_bytes().startswith(b"RFCXS\x01\0\0")
canonical = directory / "canonical-query.rfcxs"
run(profile.parent / "gadget", "native-record", "--input", publication, "--out", canonical)
assert canonical.read_bytes() == publication.read_bytes()
for args in [["--request", request], ["--result", directory / "no-request.rfcxs"],
             ["--request", request, "--result", publication],
             ["--request", request, "--rewrite", directory / "same.rfcxs",
              "--result", directory / "same.rfcxs"]]:
    run(profile, "--input", source, *args, succeeds=False)
bad_request = directory / "bad-query.json"
bad_request.write_text(json.dumps({"schema": "gadget.assays.query.v1", "operation": "unknown"}))
run(profile, "--input", source, "--request", bad_request,
    "--result", directory / "refused.rfcxs", succeeds=False)
assert not (directory / "refused.rfcxs").exists()
assert not (directory / "no-request.rfcxs").exists()
assert not (directory / "same.rfcxs").exists()
print("Workspace load profile contracts passed")
