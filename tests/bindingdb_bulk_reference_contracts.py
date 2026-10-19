"""Exercise the real format inspector and independent qualification command."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

READER = Path(sys.argv.pop(1)).resolve()
REFERENCE = Path(__file__).with_name("bindingdb_bulk_reference.py")
HEADER = ["BindingDB Reactant_set_id", "BindingDB MonomerID", "Target Name",
          "Target Source Organism According to Curator or DataSource",
          "EC50 (nM)", "Ki (nM)", "IC50 (nM)", "Kd (nM)", "Original extra"]
ROW = ["reaction-a", "compound-a", "construct", "", ">1e309", "0", " <= 0.50 ", "", '"quotes"']


def line(values, ending="\n"):
    return ("\t".join(values) + ending).encode()


class ReferenceContracts(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="bindingdb-reference-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "source.tsv"
        self.source.write_bytes(line(HEADER, "\r\n") + line(ROW))
        self.serial = 0

    def qualify(self, optimized=False, reader=READER, extra=(), checksum=None, output=None):
        self.serial += 1
        output = output or self.root / f"qualification-{self.serial}.json"
        with self.source.open("rb") as source:
            checksum = checksum or hashlib.file_digest(source, "sha256").hexdigest()
        command = [sys.executable] + (["-O"] if optimized else []) + [
            str(REFERENCE), "--reader", str(reader), "--source", str(self.source),
            "--source-sha256", checksum, "--out", str(output), *extra]
        return subprocess.run(command, capture_output=True, text=True, timeout=90), output

    def test_exact_values_limits_and_source_identity(self):
        for optimized in (False, True):
            with self.subTest(optimized=optimized):
                size = self.source.stat().st_size
                field = max(len(value.encode()) for value in HEADER + ROW)
                result, output = self.qualify(optimized, extra=[
                    "--source-bytes", str(size), "--rows", "1", "--columns", "9",
                    "--field-bytes", str(field)])
                self.assertEqual(result.returncode, 0, result.stderr)
                data = json.loads(output.read_text())
                self.assertEqual(data["status"], "PASS")
                self.assertEqual(data["limits"], {"source_bytes": size, "rows": 1,
                    "columns": 9, "field_bytes": field, "allow_unterminated_final_line": False,
                    "expected_source_bytes": size})
                self.assertEqual(data["bytes_consumed"], size)
                self.assertEqual(data["rows"], 1)
                self.assertEqual(data["reported_affinities"], 3)
                self.assertEqual(data["endpoints"], {"EC50": 1, "Ki": 1, "IC50": 1})
                self.assertTrue(data["ends_with_newline"])

    def test_header_only(self):
        self.source.write_bytes(line(HEADER))
        for optimized in (False, True):
            result, output = self.qualify(optimized)
            self.assertEqual(result.returncode, 0, result.stderr)
            data = json.loads(output.read_text())
            self.assertEqual(data["rows"], 0)
            self.assertEqual(data["reported_affinities"], 0)
            self.assertEqual(data["columns"], len(HEADER))
        result, output = self.qualify(True, extra=["--minimum-rows", "1"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("fewer observations", result.stderr)
        self.assertFalse(output.exists())

    def test_streams_beyond_old_byte_and_row_caps(self):
        with self.source.open("wb") as output:
            output.write(line(HEADER))
            for i in range(100001):
                row = ROW.copy()
                row[0], row[-1] = f"reaction-{i}", "source evidence " + "x" * 128
                output.write(line(row))
        self.assertGreater(self.source.stat().st_size, 16 * 1024 * 1024)
        result, output = self.qualify(optimized=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        data = json.loads(output.read_text())
        self.assertEqual(data["rows"], 100001)
        self.assertEqual(data["reported_affinities"], 300003)
        self.assertEqual(data["limits"]["source_bytes"], self.source.stat().st_size)

    def test_refuses_missing_newline_empty_and_bad_checksum(self):
        for source in (line(HEADER) + line(ROW).rstrip(b"\n"), b""):
            self.source.write_bytes(source)
            for optimized in (False, True):
                result, output = self.qualify(optimized)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())
        self.source.write_bytes(line(HEADER) + line(ROW))
        result, output = self.qualify(True, checksum="0" * 64)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum differs", result.stderr)
        self.assertFalse(output.exists())

    def test_limit_refusals_never_qualify(self):
        for option, value in (("--source-bytes", self.source.stat().st_size - 1),
                              ("--columns", 8), ("--field-bytes", 4), ("--rows", 0)):
            result, output = self.qualify(True, extra=[option, str(value)])
            self.assertNotEqual(result.returncode, 0, option)
            self.assertFalse(output.exists(), option)
        self.source.write_bytes(line(HEADER) + line(ROW) * 2)
        result, output = self.qualify(True, extra=["--rows", "1"])
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(output.exists())

    def test_reader_refusal_reports_stderr_and_location(self):
        row = ROW.copy(); row[5] = "bad"
        self.source.write_bytes(line(HEADER) + line(row))
        for optimized in (False, True):
            result, output = self.qualify(optimized)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("BindingDB affinity requires a complete decimal token", result.stderr)
            self.assertIn("source line 2", result.stderr)
            self.assertIn(f"byte offset {len(line(HEADER))}", result.stderr)
            self.assertFalse(output.exists())

    def test_non_utf8_source_has_an_explicit_encoding_refusal(self):
        self.source.write_bytes(line(HEADER) + line(ROW).replace(b"construct", b"caf\xe9"))
        result, output = self.qualify(True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source is not UTF-8 at source line 2", result.stderr)
        self.assertIn("byte offset", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertFalse(output.exists())

    def test_checksum_is_case_insensitive(self):
        checksum = hashlib.sha256(self.source.read_bytes()).hexdigest()
        result, output = self.qualify(True, checksum=checksum.upper())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(output.read_text())["source_sha256"], checksum)

    def test_hung_reader_is_stopped_with_a_timeout_reason(self):
        fake = self.root / "hung-reader"
        fake.write_text(f"#!{sys.executable}\nimport time\ntime.sleep(60)\n")
        fake.chmod(0o755)
        result, output = self.qualify(True, reader=fake, extra=["--timeout-seconds", "0.2"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("reader timed out after 0.2 seconds", result.stderr)
        self.assertFalse(output.exists())

    def test_timeout_after_summary_keeps_the_timeout_reason(self):
        baseline = subprocess.run([str(READER), str(self.source)], capture_output=True,
                                  text=True, check=True)
        for close_stdout in (False, True):
            with self.subTest(close_stdout=close_stdout):
                fake = self.root / "reader-stalled-after-summary"
                fake.write_text(
                    f"#!{sys.executable}\nimport os, sys, time\n"
                    f"sys.stdout.write({baseline.stdout!r})\nsys.stdout.flush()\n"
                    "sys.stderr.write('stalled after writing summary\\n')\nsys.stderr.flush()\n"
                    + ("os.close(sys.stdout.fileno())\n" if close_stdout else "")
                    + "time.sleep(60)\n")
                fake.chmod(0o755)
                result, output = self.qualify(True, reader=fake, extra=["--timeout-seconds", "1"])
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("format reader failed after source summary", result.stderr)
                self.assertIn("reader timed out after 1 seconds", result.stderr)
                self.assertIn("stalled after writing summary", result.stderr)
                self.assertFalse(output.exists())

    def test_invalid_incomplete_and_early_reader_output_preserve_reasons(self):
        for report, reason in (("{", "line 1 column 2 (char 1)"),
                               ('{"unexpected":1}\n', "incomplete reader record"),
                               ("", "format reader ended before complete source verification")):
            with self.subTest(report=report):
                fake = self.root / "invalid-reader-report"
                fake.write_text(
                    f"#!{sys.executable}\nimport sys\n"
                    "sys.stderr.write('deliberate reader diagnostic\\n')\nsys.stderr.flush()\n"
                    f"sys.stdout.write({report!r})\nsys.stdout.flush()\n")
                fake.chmod(0o755)
                result, output = self.qualify(True, reader=fake)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(reason, result.stderr)
                self.assertIn("deliberate reader diagnostic", result.stderr)
                self.assertIn("source line 1, byte offset 0", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertFalse(output.exists())

    def test_partial_result_write_never_publishes_a_pass_file(self):
        spec = importlib.util.spec_from_file_location("bindingdb_reference", REFERENCE)
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        output = self.root / "result.json"

        def fail_write(value, file, **kwargs):
            file.write('{"status": "PASS",')
            raise OSError("injected full disk")

        with mock.patch.object(reference.json, "dump", side_effect=fail_write):
            with self.assertRaisesRegex(OSError, "injected full disk"):
                reference.publish_result(output, {"status": "PASS"})
        self.assertFalse(output.exists())
        self.assertEqual(list(self.root.glob(".result.json.*")), [])
        output.write_text("retained evidence")
        with self.assertRaises(FileExistsError):
            reference.publish_result(output, {"status": "PASS"})
        self.assertEqual(output.read_text(), "retained evidence")

    def test_reader_cannot_qualify_tokens_the_oracle_refuses(self):
        baseline = subprocess.run([str(READER), str(self.source)], capture_output=True,
                                  text=True, check=True)
        for token, reason in (("NaN", "invalid lexical affinity"), ("1,000", "invalid lexical affinity")):
            row = ROW.copy(); row[5] = token
            self.source.write_bytes(line(HEADER, "\r\n") + line(row))
            records = [json.loads(item) for item in baseline.stdout.splitlines()]
            records[0]["limits"]["source_bytes"] = self.source.stat().st_size
            records[0]["limits"]["rows"] = self.source.stat().st_size
            records[0]["limits"]["expected_source_bytes"] = self.source.stat().st_size
            records[1]["cells"] = row
            records[1]["byte_length"] = len(line(row))
            records[1]["affinities"][1]["reported"] = token
            report = "\n".join(json.dumps(record) for record in records) + "\n"
            fake = self.root / "overpermissive-reader"
            fake.write_text(f"#!{sys.executable}\nimport sys\nsys.stdout.write({report!r})\n")
            fake.chmod(0o755)
            result, output = self.qualify(True, reader=fake)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(reason, result.stderr)
            self.assertFalse(output.exists())

    def test_optimized_verifier_detects_corrupted_reports(self):
        baseline = subprocess.run([str(READER), str(self.source)], capture_output=True,
                                  text=True, check=True)
        records = [json.loads(item) for item in baseline.stdout.splitlines()]
        for mutation in ("offset", "length", "endpoint", "column", "relation", "number", "cells", "summary", "limits"):
            with self.subTest(mutation=mutation):
                changed = json.loads(json.dumps(records))
                row, affinity = changed[1], changed[1]["affinities"][0]
                if mutation == "offset": row["byte_offset"] += 1
                elif mutation == "length": row["byte_length"] -= 1
                elif mutation == "endpoint": affinity["endpoint"] = "Kd"
                elif mutation == "column": affinity["column"] += 1
                elif mutation == "relation": affinity["relation"] = "="
                elif mutation == "number": affinity["number"] = "1"
                elif mutation == "cells": row["cells"][0] = "wrong occurrence"
                elif mutation == "summary": changed[-1]["rows"] += 1
                elif mutation == "limits": changed[0]["limits"]["rows"] -= 1
                fake = self.root / "incorrect-reader"
                report = "\n".join(json.dumps(record) for record in changed) + "\n"
                fake.write_text(f"#!{sys.executable}\nimport sys\nsys.stdout.write({report!r})\n")
                fake.chmod(0o755)
                result, output = self.qualify(True, reader=fake)
                self.assertNotEqual(result.returncode, 0, mutation)
                self.assertFalse(output.exists(), mutation)

    def test_existing_qualification_is_preserved(self):
        output = self.root / "retained.json"
        output.write_text("retained evidence\n")
        result, _ = self.qualify(True, output=output)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output.read_text(), "retained evidence\n")

    def test_inspector_option_and_io_errors(self):
        for extra in (["--rows", "0"], ["--rows", "-1"], ["--rows"],
                      ["--unknown", "1"], ["--rows", "1", "--rows", "2"]):
            result = subprocess.run([str(READER), str(self.source), *extra], capture_output=True)
            self.assertNotEqual(result.returncode, 0, extra)
        result = subprocess.run([str(READER), str(self.root / "absent.tsv")], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        # The inspector ignores SIGPIPE: this must be its checked ostream error,
        # rather than termination by a signal (macOS also has no /dev/full).
        import os
        read_end, write_end = os.pipe()
        os.close(read_end)
        try:
            result = subprocess.run([str(READER), str(self.source)], stdout=write_end,
                                    stderr=subprocess.PIPE)
        finally:
            os.close(write_end)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stderr)


if __name__ == "__main__":
    if sys.version_info < (3, 11):
        raise SystemExit("BindingDB qualification tests require Python 3.11 or newer")
    unittest.main()
