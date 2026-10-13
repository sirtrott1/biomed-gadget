"""Native format setup preserves matching readers and refuses occupied bindings."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

from rexgraph import io


def load_adapter():
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    return adapter


def registry_state():
    return tuple((name, id(value), getattr(value, "name", None),
                  id(getattr(value, "load", None)), id(getattr(value, "save", None)),
                  getattr(value, "extensions", None), io._FORMATS.meta(name))
                 for name, value in io._FORMATS.items())


def refuses(action):
    before = registry_state(), io.format_extensions()
    try:
        action()
    except ValueError:
        pass
    else:
        raise AssertionError("Expected format registration refusal")
    assert (registry_state(), io.format_extensions()) == before


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gadget", type=Path)
    parser.add_argument("publication", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--publisher", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--format", choices=("gadget-native", "gadget-native-compressed"))
    args = parser.parse_args()
    binary, publication, root = (path.resolve() for path in
                                 (args.gadget, args.publication, args.root))
    raw = publication.read_bytes()
    adapter = load_adapter()

    def install(name, gadget=binary, encoded=len(raw), decoded=len(raw)):
        if name == "gadget-native":
            return adapter.install_format(gadget=gadget)
        return adapter.install_compressed_format(gadget=gadget,
                                                 max_encoded_bytes=encoded, max_output_bytes=decoded)

    if args.case:
        name = args.format
        extension = ".rfcxs" if name == "gadget-native" else ".rcblob"

        def foreign(*positional, **keywords):
            raise AssertionError("Occupied reader/saver must not execute during setup")

        if args.case == "foreign-reader":
            occupied = io.register_format(name, load=foreign, extensions=(extension,))
        elif args.case == "save-only":
            occupied = io.register_format(name, save=foreign, extensions=(extension,))
        elif args.case == "missing-reader":
            occupied = io.register_format(name, extensions=(extension,))
        elif args.case in ("extra-save", "different-extension"):
            original = install(name)
            occupied = io.register_format(name, load=original.load,
                                           save=foreign if args.case == "extra-save" else None,
                                           extensions=(extension,) if args.case == "extra-save" else (".other",))
        elif args.case in ("extension-conflict", "late-extension-conflict"):
            original = install(name) if args.case == "late-extension-conflict" else None
            occupied = io.register_format("foreign-format", load=foreign, extensions=(extension,))
            refuses(lambda: install(name))
            assert io._FORMATS.get("foreign-format") is occupied
            assert io._FORMATS.get(name) is original
            print(name, args.case, "PASS")
            return
        else:
            raise AssertionError("Unknown format registration control")
        refuses(lambda: install(name))
        assert io._FORMATS.get(name) is occupied
        print(name, args.case, "PASS")
        return

    root.mkdir(parents=True)
    for value in (0, -1, True, "100", None, 1.5):
        refuses(lambda: install("gadget-native-compressed", encoded=value))
        refuses(lambda: install("gadget-native-compressed", decoded=value))
    installed = {name: install(name) for name in ("gadget-native", "gadget-native-compressed")}
    for name, original in installed.items():
        assert install(name) is original
    adapter = load_adapter()
    alias = binary.parent/".."/binary.parent.name/binary.name
    for name, original in installed.items():
        assert install(name, gadget=alias) is original
        refuses(lambda: install(name, gadget=Path(sys.executable)))
    for encoded, decoded in ((len(raw)+1, len(raw)), (len(raw), len(raw)+1)):
        refuses(lambda: install("gadget-native-compressed", encoded=encoded, decoded=decoded))
    for value in (0, -1, True, "100", None, 1.5):
        refuses(lambda: install("gadget-native-compressed", encoded=value))
        refuses(lambda: install("gadget-native-compressed", decoded=value))
    assert io.load(publication) == raw
    for codec in ("zstd", "zlib"):
        output = root/(codec+".rcblob")
        adapter.compress_publication(publication, output, gadget=binary, codec=codec,
                                     publisher=args.publisher)
        assert output.stat().st_size <= len(raw)
        assert io.load(output) == raw
        assert io.load(output, format="gadget-native-compressed") == raw
    for case in ("foreign-reader", "save-only", "missing-reader", "extra-save", "different-extension",
                 "extension-conflict", "late-extension-conflict"):
        for name in installed:
            process = subprocess.run([sys.executable, "-B", __file__, str(binary), str(publication), str(root),
                                      "--case", case, "--format", name], capture_output=True, text=True, timeout=30)
            assert process.returncode == 0, (case, name, process.stdout, process.stderr)
            print(process.stdout.strip(), flush=True)
    for name, original in installed.items():
        assert io._FORMATS.get(name) is original
    io.unregister_format("gadget-native-compressed")
    replacement = install("gadget-native-compressed", encoded=2**80, decoded=len(raw))
    assert replacement is not installed["gadget-native-compressed"]
    assert install("gadget-native-compressed", encoded=2**80, decoded=len(raw)) is replacement
    assert io.load(root/"zstd.rcblob") == raw
    assert publication.read_bytes() == raw
    receipt = {"status": "PASS", "foreign_controls": 14, "invalid_budget_refusals": 24,
               "matching_reader_identity_after_reload": True, "resolved_path_alias_reused": True,
               "different_binary_and_budget_refusals": 4, "registry_unchanged_on_refusal": True,
               "explicit_unregister_before_reconfiguration": True,
               "complete_original_bytes_equal": True, "source_bytes": len(raw),
               "scope": "Serial application initialization through the existing canonical format registry"}
    (root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("Native format reuse, unchanged occupied providers, budgets and exact bytes PASS", flush=True)


if __name__ == "__main__":
    main()
