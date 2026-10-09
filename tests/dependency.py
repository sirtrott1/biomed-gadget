import importlib.util
import hashlib
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("dependency", Path(__file__).resolve().parents[1] / "tools/check_dependency.py")
dependency = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dependency)


class BuildInputContracts(unittest.TestCase):
    def test_consumed_input_changes_are_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "include/rfc").mkdir(parents=True)
            used, unused = root / "include/rfc/used.hpp", root / "include/rfc/unused.hpp"
            used.write_text("original used input")
            unused.write_text("original unused input")
            inventory = {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (used, unused)}
            consumed = {"include/rfc/used.hpp"}
            dependency.validate_consumed(root, inventory, consumed)
            unused.write_text("ongoing independent work")
            dependency.validate_consumed(root, inventory, consumed)
            used.write_text("changed mathematical implementation")
            with self.assertRaises(SystemExit):
                dependency.validate_consumed(root, inventory, consumed)
            with self.assertRaises(SystemExit):
                dependency.validate_consumed(root, inventory, {"include/rfc/unrecorded.hpp"})
            with self.assertRaises(SystemExit):
                dependency.validate_consumed(root, inventory, set())

    def test_missing_and_stale_compiler_records_reject(self):
        root = Path("/canonical").resolve()
        objects = {"CMakeFiles/gadget_core.dir/src/io.cpp.o"}
        output = "CMakeFiles/gadget_core.dir/src/io.cpp.o: #deps 1, deps mtime 1 (VALID)\n    /canonical/include/rfc/core.hpp\n"
        self.assertEqual(dependency.ninja_headers(output, root, objects), {"include/rfc/core.hpp"})
        for invalid in ("", output.replace("VALID", "STALE")):
            with self.assertRaises(SystemExit):
                dependency.ninja_headers(invalid, root, objects)


if __name__ == "__main__":
    unittest.main()
