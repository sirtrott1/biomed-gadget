#!/usr/bin/env python3
"""Check native SDK agreement, inventory damage and a mismatched provider receipt."""
import argparse
import copy
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_info", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "check_native_stack", root / "tools/rcdb/check_native_stack.py")
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    build = json.loads(args.build_info.read_text())
    report = checker.check_stack(build)
    assert report["status"] == "PASS"

    def refuses(candidate, message):
        try:
            checker.check_stack(candidate)
        except ValueError as error:
            assert message in str(error), str(error)
        else:
            raise AssertionError("incompatible native stack was accepted")

    candidate = copy.deepcopy(build)
    candidate.pop("rexgraph_header_files")
    refuses(candidate, "no SDK header inventory")
    candidate = copy.deepcopy(build)
    name = next(iter(candidate["rexgraph_header_files"]))
    candidate["rexgraph_header_files"][name] = "0" * 64
    refuses(candidate, "differs from its recorded fingerprint")
    entries = "".join(f"{name}:{digest}\n"
                      for name, digest in sorted(candidate["rexgraph_header_files"].items()))
    candidate["rexgraph_headers_sha256"] = hashlib.sha256(entries.encode()).hexdigest()
    refuses(candidate, "Python native SDK differs from GADGET")

    module = importlib.import_module("rexgraph_cpp._model_storage")
    original = module.model_storage_build_info
    receipt = dict(original(), engine_headers="0" * 64)
    module.model_storage_build_info = lambda: receipt
    try:
        refuses(build, "rexgraph_cpp._model_storage.model_storage_build_info")
    finally:
        module.model_storage_build_info = original
    print(f"PASS: {len(report['native_providers'])} native providers and four mismatch controls")


if __name__ == "__main__":
    main()
