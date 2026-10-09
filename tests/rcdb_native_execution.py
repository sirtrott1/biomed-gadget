"""Version-bound RCQL execution of a binary GADGET assay request."""
from contextlib import closing
import importlib.util
from pathlib import Path
import sys

from rcdb import LocalStore, RecordCodec
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result

module_path = Path(__file__).resolve().parents[1]/"tools"/"rcdb"/"gadget_rcdb.py"
spec = importlib.util.spec_from_file_location("gadget_rcdb", module_path)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
binary, workspace, request, root = map(lambda x: Path(x).resolve(), sys.argv[1:5])
compression = BlobCodecSpec(name=sys.argv[5] if len(sys.argv) > 5 else "none")
root.mkdir(parents=True, exist_ok=True)
codec = adapter.install_codec(gadget=binary, max_bytes=RecordCodec.__dataclass_fields__["max_bytes"].default)
adapter.install_rcql_executor(gadget=binary)
with closing(LocalStore(root/"execution.rcdb", compression=compression,
                        record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
    adapter.put_publication(store, "workspace", workspace, expected_version=0)
    adapter.put_publication(store, "request", request, expected_version=0)
    declaration = query(source("db"), call("GADGET_ASSAY_QUERY", "workspace", 1, "request", 1))
    restricted = BoundSource(store, SourcePolicy.allow("identity", "records", "files"))
    try:
        run(declaration, sources={"db": restricted})
    except (TypeError, PermissionError, ValueError):
        pass
    else:
        raise AssertionError("Executor admitted a query without its execution capability")
    admitted = BoundSource(store, SourcePolicy.allow("identity", "records", "files", "gadget_execute"))
    result = run(declaration, sources={"db": admitted})
    payload = result.values[0]
    assert adapter.validate_publication(payload, gadget=binary) == payload
    method = result.execution[0]["methods"][0]
    assert method["method"] == "gadget-native-assay-query"
    assert method["workspace_version"] == method["request_version"] == 1
    restored = unpack_result(pack_result(result))
    assert restored.values == (payload,) and restored.execution == result.execution
    path = root/"operation.rfcxs"
    path.write_bytes(payload)
    adapter.put_publication(store, "operation", path, expected_version=0)
    assert adapter.read_publication(store, "operation", version=1).value == payload
print("RCDB workspace/request selection, policy admission, C++ execution and binary result persistence PASS")
