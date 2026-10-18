#!/usr/bin/env python3
"""Verify recorded build inputs, with an explicit live development profile."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def fingerprint(root, paths):
    entries = "".join(f"{p.relative_to(root).as_posix()}:{hashlib.sha256(p.read_bytes()).hexdigest()}\n"
                      for p in sorted(paths))
    return hashlib.sha256(entries.encode()).hexdigest()


def inventory_fingerprint(inventory):
    entries = "".join(f"{name}:{value}\n" for name, value in sorted(inventory.items()))
    return hashlib.sha256(entries.encode()).hexdigest()


def adapter_sources(root, suffix):
    for path in root.glob("tools/rcdb/**/*." + suffix):
        directory = path.relative_to(root / "tools/rcdb").parts[0]
        if directory not in ("build", "dist") and not directory.endswith(".egg-info"):
            yield path


def validate_consumed(native, inventory, consumed):
    if not consumed:
        raise SystemExit("no canonical compiler dependencies found")
    for name in sorted(consumed):
        if name not in inventory:
            raise SystemExit(f"{name}: compiler dependency absent from recorded inventory; rebuild")
        path = native / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != inventory[name]:
            raise SystemExit(f"{name}: consumed canonical header differs from recorded build; rebuild")


def ninja_headers(output, native, objects):
    consumed, seen = set(), set()
    selected = False
    for line in output.splitlines():
        if not line.startswith("    "):
            name = line.split(": #deps ", 1)[0]
            selected = name in objects
            if selected:
                if not line.endswith("(VALID)"):
                    raise SystemExit(f"{name}: compiler dependency record is stale; rebuild")
                seen.add(name)
        elif selected:
            path = Path(line.strip()).resolve()
            if path.is_relative_to(native / "include/rfc") and path.suffix == ".hpp":
                consumed.add(path.relative_to(native).as_posix())
    if seen != objects:
        raise SystemExit(f"missing compiler dependency records: {sorted(objects - seen)}; rebuild")
    return consumed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=Path("build/gadget"))
    source = Path(__file__).resolve().parents[2] / "rexgraphcpp"
    if not (source / "CMakeLists.txt").is_file():
        source = source.parent / "rexgraph-cpp/rexgraphcpp"
    parser.add_argument("--rexgraph", type=Path, default=source)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--consumed-only", action="store_true",
                        help="Ninja profile: verify every consumed header; report unrelated live-tree changes")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    native = args.rexgraph.resolve()
    expected = json.loads((args.build / "build-info.json").read_text())
    backend = expected.get("rexgraph_backend_source_files", {})
    if backend:
        if inventory_fingerprint(backend) != expected.get("rexgraph_backend_sources_sha256"):
            raise SystemExit("canonical compiled backend inventory is inconsistent")
        validate_consumed(native, backend, set(backend))
    actual = dict(rexgraph_headers_sha256=fingerprint(native, native.glob("include/rfc/**/*.hpp")),
                  gadget_sources_sha256=fingerprint(root, [root / "CMakeLists.txt",
                      *root.glob("include/gadget/**/*.hpp"), *root.glob("src/**/*.cpp"), *root.glob("src/**/*.hpp"),
                      *root.glob("apps/**/*.cpp"), *root.glob("examples/**/*.cpp"),
                      *adapter_sources(root, "cpp"), *adapter_sources(root, "py"),
                      *root.glob("examples/**/*.hpp"), *root.glob("schemas/**/*.json"),
                      *root.glob("cmake/**/*.cmake"),
                      *root.glob("cmake/**/*.cmake.in")]))
    if actual["gadget_sources_sha256"] != expected["gadget_sources_sha256"]:
        raise SystemExit("GADGET runtime sources differ from compiled build; rebuild before acceptance")
    same_tree = actual["rexgraph_headers_sha256"] == expected["rexgraph_headers_sha256"]
    consumed = set()
    if args.consumed_only:
        inventory = expected.get("rexgraph_header_files")
        if not inventory or inventory_fingerprint(inventory) != expected["rexgraph_headers_sha256"]:
            raise SystemExit("missing or inconsistent canonical header inventory; rebuild")
        objects = set(expected.get("runtime_objects", [])) or {f"CMakeFiles/gadget_core.dir/{p.relative_to(root).as_posix()}.o"
                   for p in root.glob("src/**/*.cpp")
                   if p.name != "model.cpp" and (p.name != "store.cpp" or expected.get("durable_store_enabled", False))}
        objects.add("CMakeFiles/gadget.dir/apps/main.cpp.o")
        deps = subprocess.run(["ninja", "-C", str(args.build.resolve()), "-t", "deps"],
                              check=True, capture_output=True, text=True).stdout
        consumed = ninja_headers(deps, native, objects)
        validate_consumed(native, inventory, consumed)
    elif not same_tree:
        raise SystemExit("rexgraph_headers_sha256: current tree differs from compiled build; rebuild or verify --consumed-only")
    if args.results:
        analytical = json.loads((args.results / "manifests/analysis.json").read_text())
        for key in actual:
            if analytical[key] != expected[key]:
                raise SystemExit(f"{key}: analytical output belongs to a different build")
    report = dict(status="PASS", mode="consumed_headers" if args.consumed_only else "complete_tree",
                  full_native_tree_matches=same_tree,
                  compiled_rexgraph_headers_sha256=expected["rexgraph_headers_sha256"],
                  current_rexgraph_headers_sha256=actual["rexgraph_headers_sha256"],
                  gadget_sources_sha256=actual["gadget_sources_sha256"],
                  consumed_headers=sorted(consumed))
    report["rexgraph_backend_sources_sha256"] = expected.get("rexgraph_backend_sources_sha256", "disabled")
    if backend:
        shader = args.build / "rexgraph/generated/vulkan/packed_binary.spv"
        report["compiled_shader_sha256"] = hashlib.sha256(shader.read_bytes()).hexdigest()
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: GADGET inputs and {'all consumed canonical headers' if args.consumed_only else 'complete canonical tree'} match recorded build")
    if not same_tree:
        print("The live tree has unrelated edits; the certificate identifies the recorded build revision.")


if __name__ == "__main__":
    main()
