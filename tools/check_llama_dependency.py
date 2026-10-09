#!/usr/bin/env python3
"""Verify the local llama.cpp source inventory recorded by the consumer build."""
import argparse
import hashlib
import json
from pathlib import Path

SUFFIXES = {".h", ".hpp", ".c", ".cpp", ".cc", ".inc", ".cmake", ".in", ".cu", ".mm", ".metal", ".glsl", ".comp"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=Path("build/presets/cpu"))
    parser.add_argument("--llama", type=Path, default=Path.home() / "llama.cpp")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    info = json.loads((args.build / "build-info.json").read_text())
    if not info.get("local_model_enabled"):
        raise SystemExit("local model adapter was not compiled")
    root = args.llama.resolve()
    files = [root / "CMakeLists.txt"]
    directories = ["include", "src", "common", "ggml", "vendor", "cmake"]
    if info.get("multimodal_enabled"):
        directories.append("tools/mtmd")
    for directory in directories:
        files.extend(p for p in (root / directory).rglob("*")
                     if p.is_file() and (p.suffix in SUFFIXES or p.name == "CMakeLists.txt"))
    inventory = {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    if inventory != info["llama_source_files"]:
        raise SystemExit("local llama.cpp inputs differ from configured build; rebuild before acceptance")
    fingerprint = hashlib.sha256("".join(f"{name}:{value}\n" for name, value in sorted(inventory.items())).encode()).hexdigest()
    if fingerprint != info["llama_sources_sha256"]:
        raise SystemExit("llama.cpp source inventory fingerprint mismatch")
    report = dict(status="PASS", llama_sources_sha256=fingerprint, source_files=len(files),
                  profile="direct C++ local GGUF inference; source identity, not binary/toolchain reproducibility")
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(files)} local llama.cpp inputs match configured source fingerprint")


if __name__ == "__main__":
    main()
