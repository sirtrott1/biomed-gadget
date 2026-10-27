#!/usr/bin/env python3
"""Check CI inputs before native compilation."""
import argparse
import importlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


def build_jobs(requested, cpus, available_bytes):
    if requested < 1 or cpus < 1:
        raise ValueError("Compiler counts must be positive integers")
    # Leave room for the linker, build tools and the operating system.
    capacity = (available_bytes - 1024**3) // (3 * 1024**3)
    if capacity < 1:
        raise ValueError("Native compilation requires at least 4 GiB of available memory")
    return min(requested, cpus, capacity)


def linux_build_resources():
    cpus = len(os.sched_getaffinity(0))
    memory = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available = int(memory["MemAvailable"].split()[0]) * 1024
    # Hosted and container runners can expose more host resources than admitted.
    relative = next((line.split(":", 2)[2] for line in Path("/proc/self/cgroup").read_text().splitlines()
                     if line.startswith("0::")), "/")
    base = Path("/sys/fs/cgroup")
    group = (base / relative.lstrip("/")).resolve()
    if not group.is_relative_to(base):
        group = base
    for directory in (group, *group.parents):
        if directory != base and base not in directory.parents:
            break
        limit = directory / "memory.max"
        current = directory / "memory.current"
        if limit.is_file() and current.is_file():
            value = limit.read_text().strip()
            if value != "max":
                available = min(available, max(0, int(value) - int(current.read_text())))
        quota = directory / "cpu.max"
        if quota.is_file():
            value, period = quota.read_text().split()
            if value != "max":
                cpus = min(cpus, max(1, int(value) // int(period)))
    return cpus, available


def needs_contracts(paths):
    if not paths:
        return True
    for name in paths:
        path = Path(name)
        if name == "README.md":
            continue
        if path.name == "README.md" and path.parts[0] in ("examples", "tools"):
            continue
        return True
    return False


def changed_paths(base, head, root=ROOT):
    if not base or not base.strip("0"):
        return None
    result = subprocess.run(
        ["git", "-C", str(root), "diff", "--no-ext-diff", "--no-renames",
         "--name-only", "-z", base, head, "--"], capture_output=True)
    if result.returncode:
        return None
    return [os.fsdecode(name) for name in result.stdout.split(b"\0") if name]


def matching_qualification(head, repository, root=ROOT):
    if not repository:
        return None
    def api(path):
        return json.loads(subprocess.check_output(["gh", "api", path], stderr=subprocess.PIPE))
    try:
        runs = api(f"repos/{repository}/actions/workflows/contracts.yml/runs?status=success&per_page=20")
        for run in runs["workflow_runs"]:
            # PR API head_sha is the feature head, while checkout tests a merge
            # snapshot. Only these events identify the actual tested checkout.
            if run.get("event") not in ("push", "workflow_dispatch"):
                continue
            paths = changed_paths(run["head_sha"], head, root)
            if paths is None or (paths and needs_contracts(paths)):
                continue
            jobs = api(f"repos/{repository}/actions/runs/{run['id']}/jobs?per_page=100")
            passed = {job["name"] for job in jobs["jobs"] if job["conclusion"] == "success"}
            if {"contracts (g++)", "contracts (clang++)", "Python bindings"} <= passed:
                return run["html_url"]
    except (OSError, subprocess.CalledProcessError, ValueError, KeyError, TypeError):
        return None
    return None


def qualification_status(checks_result, required, contracts_result, matching, bindings_result="skipped"):
    if checks_result != "success":
        raise ValueError("Input/documentation checks did not succeed; native acceptance is unqualified")
    if required == "true":
        if contracts_result != "success":
            raise ValueError("Required native compiler contracts did not succeed; native acceptance is unqualified")
        if bindings_result != "success":
            raise ValueError("Required Python binding contracts did not succeed; native acceptance is unqualified")
        return "PASS_CURRENT_NATIVE_CONTRACTS"
    if required == "false" and contracts_result == "skipped" and bindings_result == "skipped" and matching:
        return "PASS_PREVIOUS_EQUIVALENT_NATIVE_CONTRACTS"
    raise ValueError("No successful native compiler qualification was established")


def check_document(path):
    text = re.sub(r"(?ms)^```[^\n]*\n.*?^```[ \t]*$", "", path.read_text())
    checked = 0
    for href in re.findall(r"(?<!!)\[[^\]\n]+\]\(([^()\n]+)\)", text):
        url = urlsplit(href.strip("<>"))
        if url.scheme or url.netloc or not url.path:
            continue
        target = path.parent / unquote(url.path)
        if not target.exists():
            raise ValueError(f"{path}: missing link target {href}")
        checked += 1
    return checked


def preflight():
    versions, failures = {}, []
    for name in ("numpy", "scipy", "sympy", "jsonschema", "matplotlib"):
        try:
            importlib.import_module(name)
            versions[name] = importlib.metadata.version(name)
        except ImportError as error:
            failures.append(f"{name}: {error}")
    if failures:
        raise ValueError("Missing test dependencies: " + "; ".join(failures))
    import sympy
    matrix = sympy.Matrix([[2, 1], [1, 1]])
    if matrix.inv() * matrix != sympy.eye(2):
        raise ValueError("SymPy exact inverse check failed")
    from matplotlib.figure import Figure
    output = io.BytesIO()
    Figure().savefig(output, format="png")
    if not output.getvalue().startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("Matplotlib PNG rendering check failed")
    return versions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--base", default="")
    plan.add_argument("--head", default="HEAD")
    plan.add_argument("--output", type=Path)
    report = commands.add_parser("report")
    report.add_argument("--checks-result", required=True)
    report.add_argument("--required", required=True)
    report.add_argument("--contracts-result", required=True)
    report.add_argument("--bindings-result", default="skipped")
    report.add_argument("--matching-qualification", default="")
    report.add_argument("--test-preset", choices=("cpu", "contracts"), default="cpu")
    report.add_argument("--summary", type=Path)
    commands.add_parser("docs")
    commands.add_parser("preflight")
    resources = commands.add_parser("resources")
    resources.add_argument("--requested", type=int, default=3)
    resources.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.command == "plan":
        paths = changed_paths(args.base, args.head)
        required = needs_contracts(paths)
        qualified = None
        if not required:
            qualified = matching_qualification(args.head, os.environ.get("GITHUB_REPOSITORY"))
            required = qualified is None
        if args.output:
            with args.output.open("a") as output:
                output.write(f"contracts={str(required).lower()}\n")
                output.write(f"matching_qualification={qualified or ''}\n")
        print(json.dumps({"contracts": required, "base": args.base, "head": args.head,
                          "changed_paths": paths, "matching_qualification": qualified}))
    elif args.command == "report":
        try:
            status = qualification_status(args.checks_result, args.required,
                                          args.contracts_result, args.matching_qualification,
                                          args.bindings_result)
            detail = (f"Both native compilers passed CTest profile {args.test_preset} and Python bindings passed for this snapshot."
                      if args.required == "true" else
                      "Native compilation was reused from an executable-input-equivalent run: "
                      + args.matching_qualification)
        except ValueError as error:
            status, detail = "UNQUALIFIED", str(error)
        if args.summary:
            with args.summary.open("a") as output:
                output.write(f"Native qualification: {status}\n\n{detail}\n")
                output.write("\nPassing input/documentation checks alone does not qualify native behavior.\n")
        print(json.dumps({"status": status, "detail": detail}))
        if status == "UNQUALIFIED":
            raise SystemExit(1)
    elif args.command == "docs":
        paths = (ROOT / "README.md", *sorted((ROOT / "examples").rglob("README.md")),
                 *sorted((ROOT / "tools").rglob("README.md")))
        count = sum(check_document(path) for path in paths)
        print(f"PASS: {count} local README link paths")
    elif args.command == "resources":
        if sys.platform != "linux":
            raise ValueError("Compiler resource observation requires Linux")
        cpus, available = linux_build_resources()
        jobs = build_jobs(args.requested, cpus, available)
        if args.output:
            with args.output.open("a") as output:
                output.write(f"CMAKE_BUILD_PARALLEL_LEVEL={jobs}\n")
        print(json.dumps({"compiler_jobs": jobs, "available_cpus": cpus,
                          "available_memory_bytes": available}))
    else:
        print(json.dumps({"status": "PASS", "test_dependencies": preflight()}))


if __name__ == "__main__":
    main()
