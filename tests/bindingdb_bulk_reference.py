"""Check every original TSV cell and affinity coordinate with independent readers.

This verifies the format reader only. Native evidence publication, interpretation,
ownership and replay remain separate acceptance checks.
"""
import argparse
import collections
import csv
from decimal import Decimal
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reader", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "qualification output already exists"
    source = args.source.read_bytes()
    source_sha = hashlib.sha256(source).hexdigest()
    assert source_sha == args.source_sha256, "original source checksum differs"
    lines = source.splitlines(keepends=True)
    original = csv.reader(io.StringIO(source.decode("utf-8"), newline=""),
                          delimiter="\t", quoting=csv.QUOTE_NONE)
    header = next(original)
    process = subprocess.Popen([str(args.reader.resolve()), str(args.source.resolve())],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        def take():
            line = process.stdout.readline()
            assert line, "format reader ended before complete source verification"
            return json.loads(line)

        reported_header = take()
        assert reported_header["header"] == header
        offset, rows, cells, observed = len(lines[0]), 0, 0, 0
        endpoints, relations = collections.Counter(), collections.Counter()
        affinity_columns = {i: name.removesuffix(" (nM)") for i, name in enumerate(header)
                            if name in {"Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"}}
        for ordinal, reference in enumerate(original):
            actual = take()
            assert actual["ordinal"] == ordinal and actual["line"] == ordinal + 2
            assert actual["byte_offset"] == offset
            assert actual["byte_length"] == len(lines[ordinal + 1])
            assert actual["cells"] == reference, f"source cells changed at row {ordinal + 2}"
            assert source[offset:offset + actual["byte_length"]] == lines[ordinal + 1]
            assert [value["column"] for value in actual["affinities"]] == list(affinity_columns), (
                f"declared affinity coordinates changed at row {ordinal + 2}")
            for value in actual["affinities"]:
                column = value["column"]
                assert value["endpoint"] == affinity_columns[column]
                assert value["reported"] == reference[column]
                token = reference[column].strip(" ")
                if not token:
                    assert value["number"] == value["relation"] == ""
                    continue
                match = re.fullmatch(r"(<=|>=|<|>|=|~)? *(.*?) *", token)
                assert match
                relation, number = match.group(1) or "=", match.group(2)
                assert Decimal(number).is_finite()
                assert value["number"] == number and value["relation"] == relation
                endpoints[value["endpoint"]] += 1
                relations[relation] += 1
                observed += 1
            offset += actual["byte_length"]
            rows += 1
            cells += len(reference)
        summary = take()
        assert summary == {"rows": rows, "reported_affinities": observed}
        assert offset == len(source)
        assert process.stdout.read() == ""
        stderr = process.stderr.read()
        assert process.wait() == 0, stderr
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
    result = {
        "schema": "gadget.bindingdb.format.reference.v1",
        "status": "PASS",
        "scope": "Complete original source format, occurrence, byte interval and lexical affinity verification; native admission is a separate operation",
        "source_sha256": source_sha,
        "reader_sha256": hashlib.sha256(args.reader.read_bytes()).hexdigest(),
        "source_bytes": len(source), "rows": rows, "columns": len(header),
        "original_cells_checked": cells, "reported_affinities": observed,
        "endpoints": dict(endpoints), "relations": dict(relations),
    }
    with args.out.open("x") as output:
        output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
