"""Stream-check every original TSV cell and affinity coordinate independently.

Independent checks cover lexical round-trip, occurrences, byte intervals and
affinity tokens. Structural refusals are enforced by the C++ reader only.
Native evidence publication, interpretation, ownership and replay remain
separate acceptance checks.
"""
import argparse
import collections
import csv
from decimal import Decimal
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import threading


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("limits must be positive integers")
    return number


def sha256(value):
    if re.fullmatch(r"[0-9a-fA-F]{64}", value) is None:
        raise argparse.ArgumentTypeError("expected a 64-digit SHA-256 checksum")
    return value.lower()


def timeout_seconds(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("timeout must be a finite positive number of seconds")
    return number


def publish_result(path, result):
    # Same-filesystem hard-link publication is atomic and refuses an existing
    # destination. A failed write never exposes a partial PASS file.
    pending = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix="." + path.name + ".", delete=False) as output:
            pending = Path(output.name)
            json.dump(result, output, indent=2)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.link(pending, path)
    finally:
        if pending is not None:
            pending.unlink(missing_ok=True)


def main():
    require(sys.version_info >= (3, 11), "BindingDB qualification requires Python 3.11 or newer")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reader", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-sha256", type=sha256, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--minimum-rows", type=positive,
                        help="require at least this many observations (default: allow header-only input)")
    parser.add_argument("--timeout-seconds", type=timeout_seconds, default=3600.0,
                        help="reader wall-clock deadline, in seconds (default: 3600)")
    for name in ("source-bytes", "rows", "columns", "field-bytes"):
        parser.add_argument("--" + name, type=positive)
    args = parser.parse_args()
    require(not args.out.exists(), "qualification output already exists")
    size = args.source.stat().st_size
    limits = {
        "source_bytes": args.source_bytes or max(1, size),
        "rows": args.rows or max(1, size),
        "columns": args.columns or 4096,
        "field_bytes": args.field_bytes or 1024 * 1024,
        "allow_unterminated_final_line": False,
        "expected_source_bytes": size,
    }
    require(size <= limits["source_bytes"], "source exceeds declared byte budget")
    source_sha = digest(args.source)
    require(source_sha == args.source_sha256, "original source checksum differs")
    reader_sha = digest(args.reader)
    command = [str(args.reader.resolve()), str(args.source.resolve())]
    for name, value in limits.items():
        if name not in {"allow_unterminated_final_line", "expected_source_bytes"}:
            command += ["--" + name.replace("_", "-"), str(value)]
    csv.field_size_limit(limits["field_bytes"])
    # One source line and one reported JSON row at a time. A malformed long line
    # is bounded by the declared field/column policy and source byte budget.
    line_limit = min(limits["source_bytes"],
                     limits["columns"] * (limits["field_bytes"] + 1) + 1)
    checked_sha = hashlib.sha256()
    endpoints, relations = collections.Counter(), collections.Counter()
    offset = rows = cells = observed = 0
    with args.source.open("rb") as source, tempfile.TemporaryFile(mode="w+t") as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        timed_out = threading.Event()

        def expire():
            if process.poll() is None:
                timed_out.set()
                process.kill()

        deadline = threading.Timer(args.timeout_seconds, expire)
        deadline.daemon = True
        deadline.start()
        try:
            def failed(message, wait=False):
                if not wait and process.poll() is None:
                    process.kill()
                process.wait()
                errors.seek(0)
                detail = errors.read().strip()
                if timed_out.is_set():
                    detail = f"reader timed out after {args.timeout_seconds:g} seconds; " + detail
                return ValueError(message + (": " + detail if detail else
                                             f" (reader exit status {process.returncode})"))

            def take(required, context):
                # Decode one JSON record at a time, so an invalid source byte in
                # a later row cannot break decoding of an earlier header record.
                line = process.stdout.readline()
                if not line:
                    raise failed("format reader ended before complete source verification at " + context,
                                 wait=True)
                try:
                    record = json.loads(line.decode("utf-8"))
                except ValueError as error:
                    raise failed(f"invalid reader JSON at {context}: {error}") from error
                if not isinstance(record, dict) or not required <= record.keys():
                    raise failed("incomplete reader record at " + context)
                return record

            def original_line():
                raw = source.readline(line_limit + 1)
                if not raw:
                    return None, None
                require(len(raw) <= line_limit, "source line exceeds declared budgets")
                require(raw.endswith(b"\n"), "BindingDB source must end with a newline")
                checked_sha.update(raw)
                content = raw[:-1].removesuffix(b"\r")
                require(b"\r" not in content and b"\0" not in content,
                        "source contains an unsupported control byte")
                try:
                    decoded = content.decode("utf-8")
                except UnicodeDecodeError as error:
                    raise ValueError(f"source is not UTF-8 at source line {rows + (1 if offset == 0 else 2)}, "
                                     f"byte offset {offset + error.start}; encoding must be checked before qualification") from error
                values = next(csv.reader(io.StringIO(decoded, newline=""),
                                         delimiter="\t", quoting=csv.QUOTE_NONE))
                require(len(values) <= limits["columns"], "source exceeds declared column budget")
                return raw, values

            first, header = original_line()
            require(first is not None, "source is empty")
            reported_header = take({"header", "limits"}, "source line 1, byte offset 0")
            require(reported_header["header"] == header, "original header changed")
            require(reported_header["limits"] == limits, "reader limits differ from declared limits")
            offset = len(first)
            affinity_columns = {i: name.removesuffix(" (nM)") for i, name in enumerate(header)
                                if name in {"Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"}}
            while True:
                raw, reference = original_line()
                if raw is None:
                    break
                require(rows < limits["rows"], "source exceeds declared row budget")
                context = f" at source line {rows + 2}"
                actual = take({"ordinal", "line", "byte_offset", "byte_length", "cells", "affinities"},
                              f"source line {rows + 2}, byte offset {offset}")
                require(actual["ordinal"] == rows and actual["line"] == rows + 2,
                        "row occurrence changed" + context)
                require(actual["byte_offset"] == offset, "byte offset changed" + context)
                require(actual["byte_length"] == len(raw), "byte length changed" + context)
                require(actual["cells"] == reference, "source cells changed" + context)
                require([value["column"] for value in actual["affinities"]] == list(affinity_columns),
                        "declared affinity coordinates changed" + context)
                for value in actual["affinities"]:
                    column = value["column"]
                    require(value["endpoint"] == affinity_columns[column], "endpoint changed" + context)
                    require(value["reported"] == reference[column], "reported token changed" + context)
                    token = reference[column].strip(" ")
                    if not token:
                        require(value["number"] == value["relation"] == "",
                                "missing affinity acquired a value" + context)
                        continue
                    match = re.fullmatch(r"(<=|>=|<|>|=|~)? *((?:[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)) *", token)
                    require(match is not None, "invalid lexical affinity" + context)
                    relation, number = match.group(1) or "=", match.group(2)
                    require(Decimal(number).is_finite(), "non-finite affinity" + context)
                    require(value["number"] == number and value["relation"] == relation,
                            "affinity relation or exact number changed" + context)
                    endpoints[value["endpoint"]] += 1
                    relations[relation] += 1
                    observed += 1
                offset += len(raw)
                rows += 1
                cells += len(reference)
            summary = take({"rows", "reported_affinities", "ends_with_newline", "bytes_consumed"},
                           f"source summary after {rows} rows, byte offset {offset}")
            require(summary == {"rows": rows, "reported_affinities": observed,
                                "ends_with_newline": True, "bytes_consumed": size}, "source summary changed")
            require(rows >= (args.minimum_rows or 0), "source has fewer observations than --minimum-rows")
            require(offset == size, "source byte count changed")
            require(checked_sha.hexdigest() == source_sha, "source changed during verification")
            if process.stdout.read(1) != b"":
                raise failed("format reader emitted extra output")
            status = process.wait()
            if status != 0:
                raise failed(f"format reader failed after source summary at byte offset {offset}",
                             wait=True)
        finally:
            deadline.cancel()
            if process.poll() is None:
                process.kill()
            process.wait()
            deadline.join()
            process.stdout.close()
    require(digest(args.reader) == reader_sha, "reader changed during verification")
    result = {
        "schema": "gadget.bindingdb.format.reference.v2",
        "status": "PASS",
        "scope": "Lexical round-trip, occurrence, byte interval and affinity-token verification; structural refusals by reader only; native adapter qualification pending",
        "source_sha256": source_sha, "reader_sha256": reader_sha,
        "source_bytes": size, "rows": rows, "columns": len(header),
        "limits": limits, "ends_with_newline": True,
        "bytes_consumed": offset, "source_encoding": "utf-8",
        "minimum_rows": args.minimum_rows or 0, "reader_timeout_seconds": args.timeout_seconds,
        "original_cells_checked": cells, "reported_affinities": observed,
        "endpoints": dict(endpoints), "relations": dict(relations),
    }
    publish_result(args.out, result)
    print(json.dumps(result))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, csv.Error) as error:
        raise SystemExit(str(error)) from None
