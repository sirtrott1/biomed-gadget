"""Explicit export of native datasets and their linked provenance closure."""
from pathlib import Path
import json
import os
import shutil

from ._records import fingerprint
from .store import DataStore


def export_bundle(source, output, *, dependencies=(), resolvers=None, include_originals=False):
    """Copy a checked native closure without rewriting store or record identities.

    Dataset directories retain RCDB's own storage format and complete history.
    The JSON file is only a location index. It contains no scientific values.
    """
    from rexgraph.io.publication import publication_lock, staged_publication
    output = Path(output).absolute()
    if os.path.lexists(output):
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stores = {store.identity: store for store in (source, *dependencies)}
    resolvers = dict(resolvers or {})
    pending, closure = [source.identity], {}
    while pending:
        identity = pending.pop()
        if identity in closure:
            continue
        if identity not in stores:
            raise KeyError("missing linked dataset " + identity)
        store = stores[identity]
        closure[identity] = store
        for reference in store.manifest["parents"]:
            parent = stores.get(reference["dataset"])
            if parent is None or any(parent.reference[key] != reference[key]
                                     for key in ("store_id", "version", "state_digest")):
                raise ValueError("linked parent dataset identity or manifest version differs")
            pending.append(parent.identity)
    index = {"schema": "gadget.native-data.locations.v1", "root": source.identity, "datasets": {}, "originals": {}}
    with publication_lock(output.parent):
        if os.path.lexists(output):
            raise FileExistsError(output)
        with staged_publication(output, directory=True) as staged:
            (staged / ".gitignore").write_text("*\n")
            for identity, store in closure.items():
                store.verify()
                relative = Path("datasets") / identity
                target = staged / relative
                target.mkdir(parents=True)
                shutil.copytree(store.path / "data.rcdb", target / "data.rcdb")
                with DataStore(target, limits=store.limits, dataset=identity) as copied:
                    if copied.reference != store.reference:
                        raise ValueError("exported native dataset manifest changed")
                    copied.verify()
                index["datasets"][identity] = relative.as_posix()
                if include_originals:
                    resolver = resolvers.get(identity)
                    if resolver is None and store.manifest["sources"]:
                        raise KeyError("no original source resolver for " + identity)
                    index["originals"][identity] = {}
                    for key, reference in store.manifest["sources"].items():
                        path = resolver.resolve(reference)
                        relative_original = Path("originals") / reference["sha256"]
                        copied = staged / relative_original
                        copied.parent.mkdir(exist_ok=True)
                        if not copied.exists():
                            shutil.copyfile(path, copied)
                        if fingerprint(copied) != {name: reference[name] for name in ("sha256", "bytes")}:
                            raise ValueError("exported original source content differs")
                        index["originals"][identity][key] = relative_original.as_posix()
            (staged / "locations.json").write_text(json.dumps(index, indent=2) + "\n")
    return index
