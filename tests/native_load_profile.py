"""Check complete native rewrites, selected cutoffs and profiling refusals."""
import shutil
import subprocess
import sys
from pathlib import Path


def run(command, succeeds=True):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True)
    assert (result.returncode == 0) == succeeds, (command, result.stdout, result.stderr)
    return result


def observations(output):
    phases, values = {}, {}
    for line in output.splitlines():
        columns = line.split("\t")
        if columns[0] == "phase":
            assert len(columns) == 5 and columns[1] not in phases, columns
            assert float(columns[2]) >= 0
            if columns[3] != "unavailable":
                assert int(columns[4]) >= int(columns[3]) > 0
            phases[columns[1]] = columns[2:]
        else:
            assert columns[0] == "value" and len(columns) == 3, columns
            assert columns[1] not in values, columns
            values[columns[1]] = int(columns[2])
    return phases, values


def main():
    profile, fixture, records, directory = sys.argv[1:]
    root = Path(directory)
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    data = root / "native data"
    raw = root / "native records"
    run([fixture, data])
    run([records, raw])
    for label, extra in [("latest", []), ("previous", ["--cutoff", "4"])]:
        source = data / "changed-history.payload"
        rewritten = root / (label + ".rfcxs")
        command = [profile, "--input", source, "--accessions", "--rewrite", rewritten, *extra]
        phases, values = observations(run(command).stdout)
        assert set(phases) == {"decode", "select", "account", "accessions",
                               "accessions_reused", "rewrite", "verify_bytes"}
        assert values["source.grades"] == 4 and values["source.cells"] == 3
        assert values["source.occurrences"] == 130
        assert values["history.components"] == 3
        assert values["history.verified_bytes"] == source.stat().st_size
        assert source.read_bytes() == rewritten.read_bytes()
        run(command, succeeds=False)
    source = raw / "data.bin"
    phases, values = observations(run([profile, "--input", source, "--record", "--accessions"]).stdout)
    assert "record" in phases and values["record.unique_byte_owners"] > 0
    assert values["record.raw_owner_bytes"] > 0
    assert values["record.source_references"] >= values["record.unique_byte_owners"]
    _, values = observations(run([profile, "--input", data / "glued.dat", "--accessions"]).stdout)
    assert values["fields.tensors"] == 2 and values["fields.logical_coordinates"] == 7
    assert values["fields.supported_coordinates"] == 7 and values["fields.nonzero_entries"] == 6
    assert values["fields.support_capacity_bytes"] == 0
    assert values["fields.entry_capacity_bytes"] == 0 and values["fields.limb_capacity_bytes"] == 0
    assert values["fields.retained_value_fields"] == 2 and values["fields.materialized_value_fields"] == 0
    for extra in [["--cutoff", "3"], ["--cutoff", "6"], ["--cutoff", "-1"],
                  ["--cutoff", "18446744073709551616"], ["--record"],
                  ["--accessions", "--accessions"], ["--unknown"]]:
        run([profile, "--input", data / "changed-history.payload", *extra], succeeds=False)
    run([profile, "--input", data / "corrupt.payload"], succeeds=False)
    run([profile, "--input", data / "source.payload"], succeeds=False)
    run([profile, "--input", data / "native.payload", "--rewrite", data / "native.payload"], succeeds=False)
    print("Native load profile contracts passed")


if __name__ == "__main__":
    main()
