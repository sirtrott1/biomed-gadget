"""Check full workspace admission, exact rewriting and profiling refusals."""
from pathlib import Path
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
print("Workspace load profile contracts passed")
