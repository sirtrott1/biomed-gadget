"""Content identities for the application I/O implementation and providers."""
from importlib import metadata
from pathlib import Path
import hashlib


def runtime_identity():
    root = Path(__file__).resolve().parent
    files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(root.glob("*.py"))}
    adapter = root.parent / "gadget_rcdb.py"
    if adapter.is_file():
        files[adapter.name] = hashlib.sha256(adapter.read_bytes()).hexdigest()
    combined = "".join(f"{name}:{digest}\n" for name, digest in sorted(files.items()))
    providers = {}
    for name in ("rexgraph", "rexgraph-rcdb", "rexgraph-rcql", "rexgraph-cpp-bindings"):
        try:
            distribution = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue
        record = distribution.read_text("RECORD")
        providers[name] = {"version": distribution.version,
                          "record_sha256": None if record is None else hashlib.sha256(record.encode()).hexdigest()}
    return {"source_sha256": hashlib.sha256(combined.encode()).hexdigest(), "files": files, "providers": providers}
