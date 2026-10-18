"""Compare complete original axes while lookup leaves scientific arrays unopened."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import resource
import time
from unittest.mock import patch

import numpy as np

from gadget_data import AxisIndex, DataStore, DataWriter, Limits, build_axis_index


def execute(root):
    started = time.monotonic()
    root.mkdir(parents=True, exist_ok=False)
    rows, columns, width = 32768, 64, 256
    paths = []
    for number in range(2):
        original = root / (str(number) + ".npy")
        values = np.lib.format.open_memmap(original, mode="w+", dtype=np.int64, shape=(rows, columns))
        for offset in range(0, rows, width):
            values[offset:offset + width] = np.arange(offset * columns, (offset + width) * columns).reshape(width, columns)
        values.flush()
        del values
        table = root / (str(number) + ".tsv")
        with table.open("w") as stream:
            for i in range(rows):
                stream.write(str(i) + "\t" + ("TP53" if i % 1024 == 0 else "gene" + str(i)) + "\n")
        values = np.load(original, mmap_mode="r")
        target = root / ("data" + str(number))
        with DataWriter(target, name="original" + str(number)) as writer:
            writer.source("counts", original)
            writer.source("axis", table)
            writer.axis("profile", range(columns), offset=0)
            writer.field("counts", axes=("feature", "profile"), shape=(rows, columns), quantity="count", units="count")
            with table.open() as stream:
                for offset in range(0, rows, width):
                    batch = [next(stream).rstrip("\n").split("\t") for _ in range(width)]
                    writer.axis("feature", tuple(int(row[0]) for row in batch), offset=offset,
                        context=tuple({"symbol": row[1]} for row in batch),
                        origins=({"source": "axis", "selector": {"rows": (offset, offset + width)}},))
                    writer.dense("counts", values[offset:offset + width], np.ones((width, columns), bool),
                        offset=(offset, 0), origins=({"source": "counts", "selector": {"rows": (offset, offset + width)}},))
            writer.finish()
        del values
        paths.append(target)
    limits = Limits(catalog_bytes=64 * 1024)
    additional_data_bytes = 128 * 1024**2
    data_line = next(line for line in Path("/proc/self/status").read_text().splitlines() if line.startswith("VmData:"))
    data_bytes = int(data_line.split()[1]) * 1024
    old_limit = resource.getrlimit(resource.RLIMIT_DATA)
    proposed = data_bytes + additional_data_bytes
    if old_limit[0] != resource.RLIM_INFINITY:
        proposed = min(proposed, old_limit[0])
    resource.setrlimit(resource.RLIMIT_DATA, (proposed, old_limit[1]))
    with ExitStack() as stack:
        sources = [stack.enter_context(DataStore(path, limits=limits)) for path in paths]
        for source in sources:
            getter = source.store.get_version
            def checked(identifier, *args, getter=getter, **kwargs):
                if identifier.startswith("fields/"):
                    raise AssertionError("lookup decoded a numerical field")
                return getter(identifier, *args, **kwargs)
            stack.enter_context(patch.object(source.store, "get_version", checked))
        build_axis_index(((source, "feature") for source in sources), root / "index", context_roles=("symbol",), limits=limits)
        with AxisIndex(root / "index", limits=limits) as index:
            result = index.verify_index(datasets=sources)
            assert result == {"locators": 256, "original_positions": 65536}
            expected = {(source.identity, i) for source in sources for i in range(0, rows, 1024)}
            actual = {(entry["source"]["dataset"]["dataset"], entry["coordinate"])
                      for entry in index.find("TP53", datasets=sources)}
            assert actual == expected
            for number, source in enumerate(sources):
                with (root / (str(number) + ".tsv")).open() as stream:
                    for block in source.manifest["axes"]["feature"]:
                        labels, contexts, _ = source.axis("feature", range(block["offset"], block["offset"] + block["count"]))
                        for label, context in zip(labels, contexts, strict=True):
                            original_id, original_symbol = next(stream).rstrip("\n").split("\t")
                            assert label == int(original_id) and context == {"symbol": original_symbol}
                    assert not stream.read()
                assert source._cache_bytes <= limits.catalog_bytes
    report = {"status": "PASS", **result, "matches": len(actual), "unopened_numerical_coordinates": 2 * rows * columns,
              "metadata_cache_bytes": limits.catalog_bytes, "additional_data_bytes": additional_data_bytes,
              "seconds": time.monotonic() - started,
              "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (root / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    print(json.dumps(execute(parser.parse_args().root.resolve()), indent=2))
