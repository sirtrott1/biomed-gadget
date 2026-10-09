#!/usr/bin/env python3
"""Fetch explicitly selected public GGUFs; C++ runtime itself stays offline."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import urllib.request

MODELS = [
    ("mistralai/Ministral-3-3B-Instruct-2512-GGUF", "eb599d408350ea2bb60452cb86be7c7b2fc28227", "Ministral-3-3B-Instruct-2512-Q4_K_M.gguf", 2147023008, "9ed150d4367e68df0ac8e1540f6ddc65b42d0ee26378329d1ecbca60f93fc5f8"),
    ("ggml-org/SmolVLM2-500M-Video-Instruct-GGUF", "ccd7aae53bcb1997355c2f094959e72b3642ce17", "SmolVLM2-500M-Video-Instruct-Q8_0.gguf", 436808704, "6f67b8036b2469fcd71728702720c6b51aebd759b78137a8120733b4d66438bc"),
    ("ggml-org/SmolVLM2-500M-Video-Instruct-GGUF", "ccd7aae53bcb1997355c2f094959e72b3642ce17", "mmproj-SmolVLM2-500M-Video-Instruct-Q8_0.gguf", 108785184, "921dc7e259f308e5b027111fa185efcbf33db13f6e35749ddf7f5cdb60ef520b"),
]
MINISTRAL_PROJECTOR = ("mistralai/Ministral-3-3B-Instruct-2512-GGUF", "eb599d408350ea2bb60452cb86be7c7b2fc28227", "Ministral-3-3B-Instruct-2512-BF16-mmproj.gguf", 841501856, "70ca26e41ecea9da997e2887462ad68b848fbf64da011ff4107bbe7353519ee0")

def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("data/models"))
    parser.add_argument("--ministral-vision", action="store_true", help="also fetch the official 3B vision projector")
    args = parser.parse_args()
    receipts = []
    for repo, revision, name, size, expected in MODELS + ([MINISTRAL_PROJECTOR] if args.ministral_vision else []):
        folder = args.out / repo.replace("/", "--")
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / name
        if path.exists():
            if not path.is_file() or path.stat().st_size != size or digest(path) != expected:
                raise SystemExit(f"existing file differs from pinned model: {path}")
        else:
            temp = folder / (name + ".download")
            url = f"https://huggingface.co/{repo}/resolve/{revision}/{name}"
            try:
                with temp.open("xb") as dst, urllib.request.urlopen(url, timeout=120) as src:
                    remaining = size
                    while block := src.read(4*1024*1024):
                        remaining -= len(block)
                        if remaining < 0: raise RuntimeError("download exceeds pinned file size")
                        dst.write(block)
                if temp.stat().st_size != size or digest(temp) != expected:
                    raise RuntimeError("downloaded model failed pinned size/SHA-256")
                # Never replace an existing user model, including concurrent creation.
                os.link(temp, path)
            finally:
                temp.unlink(missing_ok=True)
        receipts.append(dict(repository=repo, revision=revision, filename=name, bytes=size, sha256=expected, path=str(path.resolve())))
        print(f"verified {path}", flush=True)
    (args.out / "model-receipts.json").write_text(json.dumps(receipts, indent=2)+"\n")

if __name__ == "__main__": main()
