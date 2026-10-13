"""Whole original archives retain identical canonical mapped compression frames."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tracemalloc


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gadget", "publisher", "source", "baseline", "root"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    binary, publisher, source, baseline, root = (getattr(args, name).resolve() for name in
                                                ("gadget", "publisher", "source", "baseline", "root"))
    root.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    original, size = digest(source), source.stat().st_size
    checks = []
    for codec in ("zstd", "zlib"):
        expected = baseline/(codec+".rcblob")
        output = root/(codec+".rcblob")
        tracemalloc.start()
        assert adapter.compress_publication(source, output, gadget=binary, publisher=publisher,
                                            codec=codec) == output
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        assert digest(output) == digest(expected) and output.stat().st_size == expected.stat().st_size
        returned = adapter.read_compressed_publication(output, gadget=binary,
                                                       max_encoded_bytes=output.stat().st_size,
                                                       max_output_bytes=size)
        assert len(returned) == size and hashlib.sha256(returned).hexdigest() == original
        del returned
        checks.append({"codec": codec, "encoded_bytes": output.stat().st_size,
                       "encoded_sha256": digest(output), "canonical_baseline_frame_equal": True,
                       "complete_original_decompressed_bytes_equal": True,
                       "python_peak_traced_allocation_bytes": peak})
        print(codec, "complete baseline frame/decompressed native bytes PASS", flush=True)
    assert digest(source) == original
    receipt = {"status": "PASS", "source_bytes": size, "source_sha256": original,
               "checks": checks, "whole_original_source_unchanged": True,
               "scope": "Existing RCDB compressor consumes admitted immutable mapped bytes; encoded output, mapped address space/resident pages and separate native reader memory still apply. No new codec or biological recalculation."}
    (root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")


if __name__ == "__main__":
    main()
