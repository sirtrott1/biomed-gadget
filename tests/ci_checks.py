"""Check that CI retains native qualification for every executable input change."""
import importlib.util
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
        runs = {"workflow_runs": [{"id": 1, "head_sha": "original", "html_url": "passed-run"}]}
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

    def test_documentation_can_avoid_native_rebuild(self):
        self.assertFalse(checks.needs_contracts([
            "README.md", "tools/rcdb/README.md", "examples/scientific/README.md"]))

    def test_every_other_input_requires_native_contracts(self):
        for path in ("src/query.cpp", "include/gadget/native_data.hpp", "CMakeLists.txt",
                     "CMakePresets.json", "cmake/RexGraph.lock.cmake", "schemas/query.json",
                     "tests/requirements.txt", "tools/rcdb/requirements.lock",
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
