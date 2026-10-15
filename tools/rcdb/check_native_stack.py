#!/usr/bin/env python3
"""Compare every installed Python native provider with a GADGET SDK inventory."""
import argparse
import hashlib
import importlib
import importlib.machinery
import importlib.metadata
import json
from pathlib import Path


def expected_engine(build):
    inventory = build.get("rexgraph_header_files")
    if not isinstance(inventory, dict) or not inventory:
        raise ValueError("GADGET build information has no SDK header inventory")
    for name, digest in inventory.items():
        if (not name.startswith("include/rfc/") or not name.endswith(".hpp")
                or not isinstance(digest, str) or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)):
            raise ValueError("GADGET SDK header inventory is invalid")
    entries = "".join(f"{name}:{digest}\n" for name, digest in sorted(inventory.items()))
    if hashlib.sha256(entries.encode()).hexdigest() != build.get("rexgraph_headers_sha256"):
        raise ValueError("GADGET SDK header inventory differs from its recorded fingerprint")
    binding_entries = "".join(f"{name.removeprefix('include/')}:{digest}\n"
                              for name, digest in sorted(inventory.items()))
    return hashlib.sha256(binding_entries.encode()).hexdigest()


def check_stack(build):
    expected = expected_engine(build)
    binding = importlib.import_module("rexgraph_cpp")
    directory = Path(binding.__file__).resolve().parent
    providers = {}
    for path in sorted(directory.iterdir()):
        suffix = next((s for s in importlib.machinery.EXTENSION_SUFFIXES
                       if path.name.endswith(s)), None)
        if suffix is None:
            continue
        name = "rexgraph_cpp." + path.name.removesuffix(suffix)
        module = importlib.import_module(name)
        if Path(module.__file__).resolve() != path.resolve():
            raise ValueError(f"{name}: imported provider differs from the installed extension")
        receipts = {key: value() for key, value in vars(module).items()
                    if callable(value) and (key == "build_info" or key.endswith("_build_info"))}
        if not receipts:
            raise ValueError(f"{name}: native provider has no SDK build receipt")
        for key, receipt in receipts.items():
            digests = [receipt[k] for k in ("engine_digest", "engine_headers") if k in receipt]
            if not digests or any(digest != expected for digest in digests):
                raise ValueError(f"{name}.{key}: Python native SDK differs from GADGET; "
                                 f"expected {expected}, received {digests}")
        providers[name] = {"path": str(path), "receipts": receipts}
    if not providers or "rexgraph_cpp._bridge" not in providers:
        raise ValueError("installed Python native binding has no compiled bridge")
    projects = {}
    for package, distribution in (("rexgraph", "rexgraph"),
                                  ("rexgraph_cpp", "rexgraph-cpp-bindings"),
                                  ("rcdb", "rexgraph-rcdb"), ("rcql", "rexgraph-rcql"),
                                  ("agent", "rexgraph-agent"), ("system", "rexgraph-system")):
        module = importlib.import_module(package)
        projects[package] = {"path": str(Path(module.__file__).resolve()),
                             "version": importlib.metadata.version(distribution)}
    return {"status": "PASS", "gadget_sources_sha256": build["gadget_sources_sha256"],
            "rexgraph_headers_sha256": build["rexgraph_headers_sha256"],
            "python_native_engine_digest": expected, "projects": projects,
            "native_providers": providers}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-info", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        report = check_stack(json.loads(args.build_info.read_text()))
    except (ValueError, ImportError, OSError) as error:
        parser.exit(1, f"native stack check failed: {error}\n")
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(report['native_providers'])} Python native providers match GADGET's SDK")


if __name__ == "__main__":
    main()
