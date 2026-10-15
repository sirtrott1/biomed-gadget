"""Native codec setup reuses matching capabilities and preserves occupied ones."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys

from rcdb.codecs import available_record_codecs, record_codec, register_record_codec


def load_adapter():
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    return adapter


def main():
    binary, publication = (Path(p).resolve() for p in sys.argv[1:3])
    raw = publication.read_bytes()
    adapter = load_adapter()
    reference = adapter.install_codec(gadget=binary, max_bytes=len(raw))
    original = record_codec(reference)

    def refuses(action):
        before = available_record_codecs()
        installed = record_codec(reference)
        try:
            action()
        except (TypeError, ValueError):
            pass
        else:
            raise AssertionError("Expected codec registration refusal")
        assert available_record_codecs() == before
        assert record_codec(reference) is installed

    for provider in (adapter, load_adapter()):
        alias = binary.parent/".."/binary.parent.name/binary.name
        for path in (binary, alias):
            assert provider.install_codec(gadget=path, max_bytes=len(raw)) == reference
            assert record_codec(reference) is original
    assert original.open(raw) == raw
    refuses(lambda: adapter.install_codec(gadget=sys.executable, max_bytes=len(raw)))
    refuses(lambda: adapter.install_codec(gadget=binary, max_bytes=len(raw)+1))
    for maximum in (0, -1, True, "100", None, 1.5, 2**80):
        refuses(lambda: adapter.install_codec(gadget=binary, max_bytes=maximum))
    for replacement in (None, 0, 1, "yes"):
        refuses(lambda: adapter.install_codec(gadget=binary, max_bytes=len(raw),
                                              replace=replacement))

    def foreign(*args):
        raise AssertionError("Provider setup must not execute an occupied capability")

    for changes in ({"encode": foreign}, {"decode": foreign}, {"admission": foreign},
                    {"object_type": "OtherPublication"}, {"max_bytes": len(raw)+1}):
        occupied = replace(original, **changes)
        register_record_codec(occupied, replace=True)
        refuses(lambda: adapter.install_codec(gadget=binary, max_bytes=len(raw)))
        assert record_codec(reference) is occupied
    assert adapter.install_codec(gadget=binary, max_bytes=len(raw), replace=True) == reference
    restored = record_codec(reference)
    assert restored is not original and restored.open(raw) == raw
    assert adapter.install_codec(gadget=binary, max_bytes=len(raw)) == reference
    assert record_codec(reference) is restored
    print("Native codec reload/path reuse, 18 refusals, occupied provider preservation and explicit replacement PASS")


if __name__ == "__main__":
    main()
