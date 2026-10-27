#!/usr/bin/env python3
"""Compare every installed Python native provider with a GADGET SDK inventory."""
import argparse
import hashlib
import importlib
import importlib.machinery
import importlib.metadata
import importlib.util
import json
from pathlib import Path


def file_digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def project_inventory(package, distribution, origin):
    """Identify installed package contents independently of their version label."""
    project = importlib.metadata.distribution(distribution)
    if Path(project.locate_file(package + "/__init__.py")).resolve() != Path(origin).resolve():
        raise ValueError(f"{distribution}: imported package differs from its installed file inventory")
    files = {}
    for entry in project.files or ():
        name = entry.as_posix()
        if (name.startswith(package + "/") and not name.endswith(".pyc")
                and "__pycache__" not in entry.parts):
            files[name] = file_digest(Path(project.locate_file(entry)))
    if not files:
        raise ValueError(f"{distribution}: installed package file inventory is unavailable")
    entries = "".join(f"{name}:{digest}\n" for name, digest in sorted(files.items()))
    return {"version": project.version, "files": len(files),
            "files_sha256": hashlib.sha256(entries.encode()).hexdigest()}


def require_profile(report, expected):
    """Compare exact installed contents and receipts without binding local paths."""
    for key in ("gadget_sources_sha256", "rexgraph_headers_sha256", "python_native_engine_digest"):
        if report.get(key) != expected.get(key):
            raise ValueError(f"receiving profile differs: {key}")
    for group in ("projects", "native_providers", "auxiliary_providers"):
        if set(report[group]) != set(expected[group]):
            raise ValueError(f"receiving profile differs: {group}")
        for name, actual in report[group].items():
            keys = ("version", "files", "files_sha256") if group == "projects" else ("sha256", "receipts")
            if any(actual.get(key) != expected[group][name].get(key) for key in keys):
                raise ValueError(f"receiving profile differs: {name}")


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


def check_stack(build, *, services=False):
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
            if not isinstance(receipt, dict):
                raise ValueError(f"{name}.{key}: native provider receipt must be a mapping")
            digests = [receipt[k] for k in ("engine_digest", "engine_headers") if k in receipt]
            engine = receipt.get("engine")
            if (isinstance(engine, str) and len(engine) == 64
                    and all(c in "0123456789abcdef" for c in engine)):
                digests.append(engine)
            if not digests or any(digest != expected for digest in digests):
                raise ValueError(f"{name}.{key}: Python native SDK differs from GADGET; "
                                 f"expected {expected}, received {digests}")
        providers[name] = {"path": str(path), "sha256": file_digest(path), "receipts": receipts}
    if not providers or "rexgraph_cpp._bridge" not in providers:
        raise ValueError("installed Python native binding has no compiled bridge")
    projects = {}
    required = [("rexgraph", "rexgraph"), ("rexgraph_cpp", "rexgraph-cpp-bindings"),
                ("rcdb", "rexgraph-rcdb"), ("rcql", "rexgraph-rcql")]
    if services:
        required.extend((("agent", "rexgraph-agent"), ("system", "rexgraph-system")))
    auxiliary = {}
    if importlib.util.find_spec("rexgraph_temporal") is not None:
        module = importlib.import_module("rexgraph_temporal._native")
        receipt = module.build_info()
        digest = receipt.get("temporal_digest", "")
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("temporal provider has no implementation fingerprint")
        path = Path(module.__file__).resolve()
        auxiliary[module.__name__] = {"path": str(path), "sha256": file_digest(path),
            "receipts": {"build_info": receipt}, "sdk_agreement": "separate receiving profile required"}
        required.append(("rexgraph_temporal", "rexgraph-temporal-bindings"))
    for package, distribution in required:
        module = importlib.import_module(package)
        projects[package] = {"path": str(Path(module.__file__).resolve()),
                             **project_inventory(package, distribution, module.__file__)}
    return {"status": "PASS", "gadget_sources_sha256": build["gadget_sources_sha256"],
            "rexgraph_headers_sha256": build["rexgraph_headers_sha256"],
            "python_native_engine_digest": expected, "projects": projects,
            "native_providers": providers, "auxiliary_providers": auxiliary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-info", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--compare", type=Path,
                        help="require exact installed contents and receipts from a previous report")
    parser.add_argument("--operator", action="append", default=[],
                        help="also inspect a named RCQL contract without executing it")
    parser.add_argument("--services", action="store_true",
                        help="also require the optional Agent and System installations")
    args = parser.parse_args()
    try:
        report = check_stack(json.loads(args.build_info.read_text()), services=args.services)
        if args.compare:
            require_profile(report, json.loads(args.compare.read_text()))
        if args.operator:
            from gadget_data.capabilities import operation_capabilities
            report["operations"] = operation_capabilities(args.operator)
    except (ValueError, ImportError, OSError) as error:
        parser.exit(1, f"native stack check failed: {error}\n")
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(report['native_providers'])} Python native providers match GADGET's SDK")


if __name__ == "__main__":
    main()
