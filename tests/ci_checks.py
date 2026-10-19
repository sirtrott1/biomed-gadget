"""Check that CI retains native qualification for every executable input change."""
import importlib.util
import json
import sys
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "ci_checks", Path(__file__).resolve().parents[1] / "tools/ci_checks.py")
checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checks)


class CIContracts(unittest.TestCase):
    def test_missing_dependencies_fail_before_native_work(self):
        with patch.object(checks.importlib, "import_module", side_effect=ImportError("unavailable")):
            with self.assertRaisesRegex(ValueError, "Missing test dependencies: .*matplotlib"):
                checks.preflight()

    def test_documentation_requires_a_completed_compiler_qualification(self):
        import json
        runs = {"workflow_runs": [{"id": 1, "event": "push", "head_sha": "original", "html_url": "passed-run"}]}
        success = {"jobs": [{"name": "contracts (g++)", "conclusion": "success"},
                            {"name": "contracts (clang++)", "conclusion": "success"}]}
        for changes, jobs, expected in (
                (["README.md"], success, "passed-run"),
                ([], success, "passed-run"),
                (["CMakeLists.txt", "README.md"], success, None),
                (None, success, None),
                (["README.md"], {"jobs": success["jobs"][:1]}, None),
                (["README.md"], {"jobs": [{"name": "contracts (g++)", "conclusion": None}]}, None)):
            with self.subTest(changes=changes, jobs=jobs):
                with patch.object(checks, "changed_paths", return_value=changes), \
                     patch.object(checks.subprocess, "check_output",
                                  side_effect=[json.dumps(runs), json.dumps(jobs)]):
                    self.assertEqual(checks.matching_qualification("HEAD", "fixture/repo"), expected)
        with patch.object(checks.subprocess, "check_output", side_effect=OSError("unavailable")):
            self.assertIsNone(checks.matching_qualification("HEAD", "fixture/repo"))

    def test_prior_pr_head_is_not_a_tested_merge_qualification(self):
        for event in ("pull_request", "pull_request_target", None):
            runs = {"workflow_runs": [{"id": 1, "event": event,
                                      "head_sha": "untested-feature-head", "html_url": "pr-run"}]}
            with self.subTest(event=event), \
                 patch.object(checks, "changed_paths", return_value=["README.md"]) as changed, \
                 patch.object(checks.subprocess, "check_output", return_value=json.dumps(runs)) as api:
                self.assertIsNone(checks.matching_qualification("current", "fixture/repo"))
                changed.assert_not_called()
                self.assertEqual(api.call_count, 1)

    def test_documentation_can_avoid_native_rebuild(self):
        self.assertFalse(checks.needs_contracts([
            "README.md", "tools/rcdb/README.md", "examples/scientific/README.md"]))

    def test_every_other_input_requires_native_contracts(self):
        for path in ("src/statistics/query.cpp", "include/gadget/native_data.hpp", "CMakeLists.txt",
                     "CMakePresets.json", "cmake/RexGraph.lock.cmake", "schemas/query.json",
                     "tests/requirements.txt", "tools/rcdb/pyproject.toml",
                     "tests/fixtures/source.rfcxs", "tools/reader.py",
                     ".github/workflows/contracts.yml", "LICENSE", "new/provider.py"):
            with self.subTest(path=path):
                self.assertTrue(checks.needs_contracts(["README.md", path]))
        self.assertTrue(checks.needs_contracts(None))
        self.assertTrue(checks.needs_contracts([]))

    def test_git_changes_and_unavailable_history(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def git(*args):
                return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
            git("init", "--quiet")
            git("config", "user.name", "CI contract fixture")
            git("config", "user.email", "fixture@example.invalid")
            (root / "README.md").write_text("Original\n")
            git("add", "README.md")
            git("commit", "--quiet", "-m", "Original fixture")
            base = git("rev-parse", "HEAD")
            (root / "README.md").write_text("Updated\n")
            git("commit", "--quiet", "-am", "Update fixture")
            self.assertFalse(checks.needs_contracts(checks.changed_paths(base, "HEAD", root)))
            (root / "source.cpp").write_text("int value;\n")
            git("add", "source.cpp")
            git("commit", "--quiet", "-m", "Add executable input")
            self.assertTrue(checks.needs_contracts(checks.changed_paths(base, "HEAD", root)))
            self.assertIsNone(checks.changed_paths("0" * 40, "HEAD", root))
            self.assertIsNone(checks.changed_paths("f" * 40, "HEAD", root))

    def test_native_acceptance_cannot_be_replaced_by_input_checks(self):
        self.assertEqual(checks.qualification_status("success", "true", "success", ""),
                         "PASS_CURRENT_NATIVE_CONTRACTS")
        self.assertEqual(checks.qualification_status("success", "false", "skipped", "passed-run"),
                         "PASS_PREVIOUS_EQUIVALENT_NATIVE_CONTRACTS")
        for inputs, required, native, matching in (
                ("success", "true", "failure", "passed-run"),
                ("success", "true", "skipped", "passed-run"),
                ("success", "true", "cancelled", ""),
                ("success", "false", "skipped", ""),
                ("success", "false", "failure", "passed-run"),
                ("failure", "false", "skipped", "passed-run"),
                ("cancelled", "true", "success", ""),
                ("success", "", "skipped", "")):
            with self.subTest(inputs=inputs, required=required, native=native):
                with self.assertRaises(ValueError):
                    checks.qualification_status(inputs, required, native, matching)

    def test_pr_selection_uses_the_tested_merge_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def git(*args):
                return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
            git("init", "--quiet")
            git("config", "user.name", "PR selection fixture")
            git("config", "user.email", "fixture@example.invalid")
            (root / "README.md").write_text("Original\n")
            (root / "source.cpp").write_text("int original;\n")
            git("add", ".")
            git("commit", "--quiet", "-m", "Original")
            original = git("rev-parse", "HEAD")
            git("checkout", "--quiet", "-b", "feature-docs")
            (root / "README.md").write_text("PR documentation\n")
            git("commit", "--quiet", "-am", "PR documentation")
            feature = git("rev-parse", "HEAD")
            git("checkout", "--quiet", "-b", "develop", original)
            (root / "source.cpp").write_text("int newer_base;\n")
            git("commit", "--quiet", "-am", "New base executable input")
            base = git("rev-parse", "HEAD")
            git("merge", "--quiet", "--no-ff", feature, "-m", "Tested PR merge")
            merge = git("rev-parse", "HEAD")
            self.assertEqual(checks.changed_paths(base, merge, root), ["README.md"])
            self.assertFalse(checks.needs_contracts(checks.changed_paths(base, merge, root)))
            # A feature-head qualification cannot qualify newer code inherited from develop.
            self.assertTrue(checks.needs_contracts(checks.changed_paths(feature, merge, root)))
            for before in ("", "0" * 40, "f" * 40):
                self.assertTrue(checks.needs_contracts(checks.changed_paths(before, merge, root)))

    def test_report_cli_marks_missing_native_sdk_unqualified(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "summary.md"
            result = subprocess.run([
                sys.executable, str(checks.ROOT / "tools/ci_checks.py"), "report",
                "--checks-result", "success", "--required", "true",
                "--contracts-result", "failure", "--summary", str(summary)],
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(json.loads(result.stdout)["status"], "UNQUALIFIED")
            self.assertIn("native acceptance is unqualified", summary.read_text())

    def test_document_links_and_code_examples(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "README.md"
            path.write_text("[missing](missing.md)\n")
            with self.assertRaises(ValueError):
                checks.check_document(path)
            (path.parent / "present.md").write_text("Present\n")
            path.write_text("[local](present.md#anchor)\n[web](https://example.com/)\n"
                            "```python\nvalue = mapping[key](argument)\n```\n")
            self.assertEqual(checks.check_document(path), 1)


if __name__ == "__main__":
    unittest.main()
