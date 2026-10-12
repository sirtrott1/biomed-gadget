"""RCQL native field owners at exact versions without Python numerical conversion."""
from contextlib import closing
import importlib.util
from pathlib import Path
import subprocess
import sys

from rcdb import LocalStore, RecordCodec
from rcdb.header import DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result

binary, selector, verifier, publication, request, root = map(lambda x: Path(x).resolve(), sys.argv[1:])
root.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location(
    "gadget_rcdb", Path(__file__).resolve().parents[1] / "tools/rcdb/gadget_rcdb.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
control = root / "request.rfcxs"
subprocess.run([str(binary), "native-record", "--input", str(request), "--out", str(control)], check=True)
codec = adapter.install_codec(gadget=binary, max_bytes=RecordCodec.__dataclass_fields__["max_bytes"].default)
name = adapter.install_rcql_field_selector(selector=selector, gadget=binary)
assert name == adapter.install_rcql_field_selector(selector=selector, gadget=binary)
with closing(LocalStore(root / "fields.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
    adapter.put_publication(store, "source", publication, expected_version=0)
    adapter.put_publication(store, "source", publication, expected_version=1)
    adapter.put_publication(store, "control", control, expected_version=0)
    declaration = query(source("db"), call(name, "source", 1, "control", 1))
    restricted = BoundSource(store, SourcePolicy.allow("identity", "records", "files"))
    try:
        run(declaration, sources={"db": restricted})
    except (TypeError, ValueError, PermissionError):
        pass
    else:
        raise AssertionError("Field selection admitted a missing execution capability")
    admitted = BoundSource(store, SourcePolicy.allow("identity", "records", "files", "gadget_field_select"))
    result = run(declaration, sources={"db": admitted})
    assert type(result.values[0]) is bytes
    method = result.execution[0]["methods"][0]
    assert method["method"] == "gadget-native-field-selection"
    assert method["publication_version"] == method["request_version"] == 1
    restored = unpack_result(pack_result(result))
    assert restored.values == result.values and restored.execution == result.execution
    output = root / "selection.rfcxs"
    output.write_bytes(result.values[0])
    subprocess.run([str(verifier), str(publication), str(output), str(request)], check=True)
    adapter.put_publication(store, "selection", output, expected_version=0)
    assert adapter.read_publication(store, "selection", version=1).value == result.values[0]
print("RCQL exact selected versions/capabilities/native field owners and portable receipts PASS")
