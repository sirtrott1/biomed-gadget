#!/usr/bin/env python3
"""Independently verify a source format reader against original PubChem captures.

Uses Python's JSON decoder (Decimal and raw numeric lexemes) and csv reader,
not the C++ normalization helpers. This is not native publication qualification.
"""
import argparse
import csv
from decimal import Decimal
import hashlib
import io
import json
from pathlib import Path
import subprocess


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(data):
    return "sha256:" + hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reader", required=True, type=Path)
    parser.add_argument("--assay", required=True, type=Path)
    parser.add_argument("--concise", required=True, type=Path)
    parser.add_argument("--aid", required=True)
    parser.add_argument("--assay-sha256", required=True)
    parser.add_argument("--concise-sha256", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    checks = 0

    def check(condition, message):
        nonlocal checks
        require(condition, message)
        checks += 1

    results = {}
    for mode, path, expected in [("assay", args.assay, args.assay_sha256),
                                 ("concise", args.concise, args.concise_sha256)]:
        data = path.read_bytes()
        check(sha(data) == expected, f"{mode}: input fingerprint differs")
        result = subprocess.run([str(args.reader.resolve()), mode, str(path.resolve()), args.aid],
                                text=True, capture_output=True, check=True)
        (args.out / f"{mode}-reader.json").write_text(result.stdout)
        output = json.loads(result.stdout, parse_float=Decimal)
        check(output["source_sha256"] == expected, f"{mode}: reader fingerprint differs")
        check(output["source_bytes"].encode() == data, f"{mode}: original bytes changed")
        results[mode] = output

    source_text = args.assay.read_bytes().decode()
    original = json.loads(source_text, parse_float=Decimal)
    # Each JSON numeric token is retained by an independent standard decoder.
    lexemes = json.loads(source_text, parse_int=lambda s: {"number": s},
                         parse_float=lambda s: {"number": s})
    output = results["assay"]
    check(output["original"] == original, "complete original assay document changed")
    submit = original["PC_AssaySubmit"]
    raw_submit = lexemes["PC_AssaySubmit"]
    description = submit["assay"]["descr"]
    projection = output["projection"]
    check(projection["aid"] == args.aid, "assay AID changed")
    check(projection["description"]["original"] == description, "full description changed")
    check(projection["description"]["source_pointer"] == "/PC_AssaySubmit/assay/descr",
          "description pointer changed")
    definitions = description.get("results", [])
    check(len(projection["definitions"]) == len(definitions), "descriptor count differs")
    for i, (actual, definition) in enumerate(zip(projection["definitions"], definitions)):
        check(actual["original"] == definition, f"descriptor {i}: metadata/units changed")
        check(actual["source_pointer"] == f"/PC_AssaySubmit/assay/descr/results/{i}",
              f"descriptor {i}: original position changed")
        check(actual["tid"] == str(definition["tid"]), f"descriptor {i}: TID differs")
        check(actual["value_member"] == {1: "fval", 2: "ival", 3: "bval", 4: "sval"}[definition["type"]],
              f"descriptor {i}: type differs")
    rows = submit.get("data", [])
    check(len(projection["substances"]) == len(rows), "substance occurrences lost")
    reported = absent = numeric = 0
    for i, (actual, row) in enumerate(zip(projection["substances"], rows)):
        check(actual["original"] == row, f"substance {i}: original row changed")
        check(actual["source_pointer"] == f"/PC_AssaySubmit/data/{i}", f"substance {i}: pointer changed")
        check(actual["sid"] == str(row["sid"]), f"substance {i}: SID changed")
        check(actual["outcome_reported"] == ("outcome" in row), f"substance {i}: outcome missingness changed")
        check(actual.get("outcome") == row.get("outcome"), f"substance {i}: outcome changed")
        check(len(actual["readouts"]) == len(definitions), f"substance {i}: declared missing fields dropped")
        values = {value["tid"]: (j, value) for j, value in enumerate(row.get("data", []))}
        for d, (readout, definition) in enumerate(zip(actual["readouts"], definitions)):
            tid = definition["tid"]
            check(readout["tid"] == str(tid), f"substance {i}: readout identity changed")
            check(readout["definition_pointer"] == f"/PC_AssaySubmit/assay/descr/results/{d}",
                  f"substance {i}, TID {tid}: descriptor binding changed")
            if tid not in values:
                check(readout["support"] == "absent" and "number_token" not in readout,
                      f"substance {i}, TID {tid}: absent value fabricated")
                absent += 1
                continue
            j, value = values[tid]
            check(readout["support"] == "reported", f"substance {i}, TID {tid}: reported value lost")
            check(readout["original"] == value, f"substance {i}, TID {tid}: typed value changed")
            check(readout["source_pointer"] == f"/PC_AssaySubmit/data/{i}/data/{j}",
                  f"substance {i}, TID {tid}: source occurrence changed")
            member = {1: "fval", 2: "ival", 3: "bval", 4: "sval"}[definition["type"]]
            check(readout["value_member"] == member, f"substance {i}, TID {tid}: type changed")
            if definition["type"] in (1, 2):
                token = raw_submit["data"][i]["data"][j]["value"][member]["number"]
                check(readout["number_token"] == token, f"substance {i}, TID {tid}: original numeric spelling changed")
                numeric += 1
            check("comparison" not in readout, f"substance {i}, TID {tid}: modifier pairing inferred")
            reported += 1
    check([r["sid"] for r in projection["revocations"]] == [str(s) for s in submit.get("revoke", [])],
          "revocation declarations changed")

    table = list(csv.reader(io.StringIO(args.concise.read_bytes().decode(), newline="")))
    header, source_rows = table[0], table[1:]
    concise = results["concise"]
    check(concise["original"]["header"] == header, "concise original header changed")
    check([r["cells"] for r in concise["original"]["rows"]] == source_rows, "concise source cells changed")
    expected_links = []
    for record, cells in enumerate(source_rows, 2):
        row = dict(zip(header, cells))
        check(row["AID"] == args.aid, "concise response contains another assay")
        pairs = []
        if row["SID"] and row["CID"]:
            pairs.append((f"pubchem:sid:{row['SID']}", f"pubchem:cid:{row['CID']}", "reported_substance_compound_link"))
        if row.get("Target GeneID"):
            pairs.append((f"pubchem:aid:{row['AID']}", f"ncbi:gene:{row['Target GeneID']}", "reported_assay_target"))
        if row.get("Target Accession"):
            pairs.append((f"pubchem:aid:{row['AID']}", f"pubchem:target_accession:{row['Target Accession']}", "reported_assay_target_accession"))
        for left, right, meaning in pairs:
            expected_links.append({"source_record": record, "cells": cells,
                                   "left": left, "right": right, "relationship": meaning})
    check(concise["projection"]["links"] == expected_links, "concise explicit identity occurrences differ")
    report = {"scope": "Independent original-source format verification only; native admission is a separate operation",
              "checks": checks, "aid": args.aid, "substance_occurrences": len(rows),
              "definitions": len(definitions), "reported_readouts": reported,
              "numeric_readouts": numeric, "absent_readout_coordinates": absent,
              "concise_rows": len(source_rows), "explicit_links": len(expected_links),
              "assay_sha256": args.assay_sha256, "concise_sha256": args.concise_sha256}
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
