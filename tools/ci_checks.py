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
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


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
            paths = changed_paths(run["head_sha"], head, root)
            if paths is None or (paths and needs_contracts(paths)):
                continue
            jobs = api(f"repos/{repository}/actions/runs/{run['id']}/jobs?per_page=100")
            passed = {job["name"] for job in jobs["jobs"] if job["conclusion"] == "success"}
            if {"contracts (g++)", "contracts (clang++)"} <= passed:
                return run["html_url"]
    except (OSError, subprocess.CalledProcessError, ValueError, KeyError, TypeError):
        return None
    return None


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
    commands.add_parser("docs")
    commands.add_parser("preflight")
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
        print(json.dumps({"contracts": required, "changed_paths": paths,
                          "matching_qualification": qualified}))
    elif args.command == "docs":
        paths = (ROOT / "README.md", *sorted((ROOT / "examples").rglob("README.md")),
                 *sorted((ROOT / "tools").rglob("README.md")))
        count = sum(check_document(path) for path in paths)
        print(f"PASS: {count} local README link paths")
    else:
        print(json.dumps({"status": "PASS", "test_dependencies": preflight()}))


if __name__ == "__main__":
    main()
