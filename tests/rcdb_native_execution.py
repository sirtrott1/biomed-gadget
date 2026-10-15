"""Version-bound RCQL execution of a binary GADGET assay request."""
from contextlib import closing
import importlib.util
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

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
workspace_bytes, request_bytes = workspace.read_bytes(), request.read_bytes()


def rejects(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("Expected native execution refusal")


for invalid in (b"invalid", bytearray(request_bytes), request_bytes[:-1]):
    rejects(lambda: adapter.execute_publication(workspace_bytes, invalid, gadget=binary))
    rejects(lambda: adapter.execute_publication(invalid, request_bytes, gadget=binary))

with patch.object(subprocess, "run", wraps=subprocess.run) as dispatched:
    adapter.execute_publication(workspace_bytes, request_bytes, gadget=binary)
assert [Path(c.args[0][c.args[0].index("--input")+1]).name
        for c in dispatched.call_args_list if c.args[0][1] == "native-record"] == ["result.rfcxs"]

# A canonical file replacement must not change the already captured query answer.
original_run = subprocess.run
replacements = []


def replace_answer(command, **kwargs):
    if command[1] == "native-record":
        path = Path(command[command.index("--input") + 1])
        if path.name == "result.rfcxs":
            path.write_bytes(request_bytes)
            replacements.append(path)
    return original_run(command, **kwargs)


subprocess.run = replace_answer
try:
    rejects(lambda: adapter.execute_publication(workspace_bytes, request_bytes, gadget=binary))
finally:
    subprocess.run = original_run
assert len(replacements) == 1
assert workspace.read_bytes() == workspace_bytes and request.read_bytes() == request_bytes
codec = adapter.install_codec(gadget=binary, max_bytes=RecordCodec.__dataclass_fields__["max_bytes"].default)
name = adapter.install_rcql_executor(gadget=binary)
assert name == adapter.install_rcql_executor(gadget=binary)
rejects(lambda: adapter.install_rcql_executor(gadget=workspace))
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
