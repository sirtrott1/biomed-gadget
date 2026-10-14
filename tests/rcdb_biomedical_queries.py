"""Exact RCQL record versions and complete native biomedical field owners."""
from contextlib import closing
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

from rcdb import LocalStore, RecordCodec
from rcdb.header import DEFAULT_RECORD_CODECS
from rcql import BoundSource, SourcePolicy, call, query, run, source
from rcql.result_codec import pack_result, unpack_result


def main(binary, selector, factory, verifier, consumer, root):
    root.mkdir(parents=True, exist_ok=False)
    repository = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("gadget_rcdb", repository/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    checks = 0

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def invoke(*args, json_output=True):
        result = subprocess.run(list(map(str, args)), text=True, capture_output=True, timeout=300)
        if result.returncode:
            raise RuntimeError(result.stderr)
        if not result.stdout:
            return None
        return json.loads(result.stdout) if json_output else result.stdout

    def refuses(action):
        try:
            action()
        except (ValueError, TypeError, PermissionError):
            check(True, "declared source or capability refused")
        else:
            raise AssertionError("invalid source, cutoff, carrier or capability was admitted")

    fixture_parent = root/"fixtures"
    fixture_receipt = invoke(sys.executable, "-I", "-B", repository/"tests/native_biomedical_fields.py",
                            binary, consumer, fixture_parent)
    check(fixture_receipt["status"] == "PASS", "independent native producer controls")
    fixtures = list(fixture_parent.glob("biomedical-fields-*"))
    check(len(fixtures) == 1, "isolated native fixture owner")
    fixture = fixtures[0]
    requests = json.loads((fixture/"requests.json").read_text())
    codec = adapter.install_codec(gadget=binary, max_bytes=RecordCodec.__dataclass_fields__["max_bytes"].default)
    executor = adapter.install_rcql_executor(gadget=binary)
    field_selector = adapter.install_rcql_field_selector(selector=selector, gadget=binary)
    field_count = 0
    results = []
    with closing(LocalStore(root/"queries.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as store:
        allowed = BoundSource(store, SourcePolicy.allow("identity", "records", "files",
                                                       "gadget_execute", "gadget_field_select"))
        restricted = BoundSource(store, SourcePolicy.allow("identity", "records", "files"))

        def execute(name, *arguments, policy=allowed):
            return run(query(source("db"), call(name, *arguments)), sources={"db": policy})

        def portable(result):
            check(type(result.values[0]) is bytes, "native bytes remain distinct from Python mathematical state")
            restored = unpack_result(pack_result(result))
            check(restored.values == result.values and restored.execution == result.execution,
                  "portable typed execution retains complete bytes and source trace")
            return result.execution[0]["methods"][0]

        for name, request in requests.items():
            case = root/name
            case.mkdir()
            control_json = case/"request.json"
            control_json.write_text(json.dumps(request)+"\n")
            control = case/"request.rfcxs"
            invoke(binary, "native-record", "--input", control_json, "--out", control)
            workspace = fixture/"native"/(name+".workspace.rfcxs")
            workspace_id, request_id = name+"-workspace", name+"-request"
            adapter.put_publication(store, workspace_id, workspace, expected_version=0)
            adapter.put_publication(store, request_id, control, expected_version=0)
            # New latest versions deliberately have a different source contract.
            adapter.put_publication(store, workspace_id, control, expected_version=1)
            adapter.put_publication(store, request_id, workspace, expected_version=1)
            refuses(lambda: execute(executor, workspace_id, 1, request_id, 1, policy=restricted))
            refuses(lambda: execute(executor, workspace_id, 2, request_id, 1))
            refuses(lambda: execute(executor, workspace_id, 1, request_id, 2))
            result = execute(executor, workspace_id, 1, request_id, 1)
            trace = portable(result)
            check(trace["workspace_version"] == trace["request_version"] == 1 and
                  trace["workspace_state_digest"] == adapter.read_publication(store, workspace_id, version=1).state_digest and
                  trace["request_state_digest"] == adapter.read_publication(store, request_id, version=1).state_digest,
                  "query uses exact historical workspace and request owners")
            artifact = case/"operation.rfcxs"
            artifact.write_bytes(result.values[0])
            observed = invoke(binary, "record-view", "--input", artifact, "--member", "result")
            expected = json.loads((fixture/"native"/(name+".json")).read_text())
            check(observed == expected, "typed native query agrees with the independent complete biomedical result")
            check(invoke(consumer, "--replay", artifact)["status"] == "PASS",
                  "typed producer replays from embedded source after original intake relocation")
            operation_id = name+"-operation"
            adapter.put_publication(store, operation_id, artifact, expected_version=0)
            adapter.put_publication(store, operation_id, control, expected_version=1)
            request_root = case/"fields"
            generated = invoke(factory, "--all-fields", artifact, request_root)
            manifest = invoke(binary, "record-view", "--input", request_root/"manifest.rfcxs")
            keys = manifest["fields"]
            check(len(keys) == len(set(keys)) == generated["fields"] and bool(keys),
                  "every actual native tensor has one original-source selector")
            for key in keys:
                field_count += 1
                controls = request_root/str(key)
                for kind in ("full", "window"):
                    field_id = name+"-"+str(key)+"-"+kind
                    adapter.put_publication(store, field_id, controls/(kind+".rfcxs"), expected_version=0)
                    adapter.put_publication(store, field_id, controls/"future.rfcxs", expected_version=1)
                    refuses(lambda: execute(field_selector, operation_id, 1, field_id, 1, policy=restricted))
                    refuses(lambda: execute(field_selector, operation_id, 1, field_id, 2))
                    selected = execute(field_selector, operation_id, 1, field_id, 1)
                    trace = portable(selected)
                    check(trace["publication_version"] == trace["request_version"] == 1 and
                          trace["source_state_digest"] == adapter.read_publication(store, operation_id, version=1).state_digest and
                          trace["request_state_digest"] == adapter.read_publication(store, field_id, version=1).state_digest,
                          "field query binds exact original record versions and execution plan")
                    output = controls/(kind+"-selection.rfcxs")
                    output.write_bytes(selected.values[0])
                    invoke(verifier, artifact, output, controls/(kind+".rfcxs"), json_output=False)
                    check(True, "C++ compares whole source prefix, every tensor/carrier/support and window coordinate")
                for kind in ("foreign", "wrong-carrier"):
                    bad_id = name+"-"+str(key)+"-"+kind
                    adapter.put_publication(store, bad_id, controls/(kind+".rfcxs"), expected_version=0)
                    refuses(lambda: execute(field_selector, operation_id, 1, bad_id, 1))
            results.append(dict(name=name, fields=len(keys), sha256=hashlib.sha256(result.values[0]).hexdigest()))
    with closing(LocalStore(root/"queries.rcdb", record_codecs=(*DEFAULT_RECORD_CODECS, codec))) as reopened:
        for item in results:
            raw = adapter.read_publication(reopened, item["name"]+"-operation", version=1).value
            check(hashlib.sha256(raw).hexdigest() == item["sha256"],
                  "reopened history retains every exact selected producer publication")
    receipt = dict(status="PASS", checks=checks, actual_fields=field_count, producers=results,
                   fixture_scope=fixture_receipt,
                   scope="Complete version-bound native biomedical queries and all tensor windows; generated independent controls, not a biological population or socket admission qualification")
    (root/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main(*[Path(value).resolve() for value in sys.argv[1:]])
