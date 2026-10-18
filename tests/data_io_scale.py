"""Bounded native intake and complete original comparison under allocation admission."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time


def run(root, *, child=False):
    if not child:
        import numpy as np
        root.mkdir(parents=True, exist_ok=False)
        original = np.lib.format.open_memmap(root / "original.npy", mode="w+", dtype="int64", shape=(32768, 1024))
        columns = np.arange(1024, dtype=np.int64)[None, :]
        for row in range(0, len(original), 128):
            original[row:row + 128] = np.arange(row, row + 128, dtype=np.int64)[:, None] * 17 + columns - 23
        original.flush()
        del original
        process = subprocess.run([sys.executable, "-B", __file__, str(root), "--child"], capture_output=True, text=True)
        (root / "child.stdout").write_text(process.stdout)
        (root / "child.stderr").write_text(process.stderr)
        if process.returncode:
            raise RuntimeError("native intake resource child failed: " + process.stderr)
        report = json.loads(process.stdout)
        (root / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report))
        return
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "rcdb"))
    import numpy as np
    import rcdb
    from gadget_data import DataStore, Limits
    from gadget_data.readers import import_dense
    allocation_bytes = 192 * 1024**2
    baseline = next(int(line.split()[1]) * 1024 for line in Path("/proc/self/status").read_text().splitlines()
                    if line.startswith("VmData:"))
    resource.setrlimit(resource.RLIMIT_DATA, (baseline + allocation_bytes, baseline + allocation_bytes))
    started = time.monotonic()
    limits = Limits(record_bytes=2 * 1024**2, decode_bytes=8 * 1024**2, selection_bytes=4 * 1024**2)
    import_dense(root / "original.npy", root / "native", reader="npy", name="complete",
                 quantity="original_count", units="count", limits=limits)
    converted = time.monotonic()
    coordinates = 0
    with DataStore(root / "native", limits=limits) as store:
        for selected in store.iter_select("measurements"):
            expected = np.asarray(selected.rows, dtype=np.int64)[:, None] * 17 + np.asarray(selected.columns, dtype=np.int64) - 23
            if not np.array_equal(selected.values, expected) or not selected.observed.all():
                raise ValueError("complete native original coefficient or support comparison failed")
            coordinates += selected.values.size
        store.materialize("measurements", root / "selected", rows=tuple(range(0, 32768, 1009)),
                          columns=(1023, 0, 17), tile=(8, 2))
    with DataStore(root / "selected", limits=limits) as selected:
        selected.verify()
    print(json.dumps({"status": "PASS", "coordinates": coordinates,
        "original_bytes": (root / "original.npy").stat().st_size, "additional_allocation_bytes": allocation_bytes,
        "baseline_data_bytes": baseline, "allocation_limit_bytes": baseline + allocation_bytes,
        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "intake_seconds": converted - started, "complete_seconds": time.monotonic() - started,
        "scope": "Complete original integer and support comparison, native selected materialization and cold reader; one process allocation admission"}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("root", type=Path); parser.add_argument("--child", action="store_true")
    args = parser.parse_args()
    run(args.root.resolve(), child=args.child)
