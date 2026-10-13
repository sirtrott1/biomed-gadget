"""CLI admission and native durable study roundtrip on the build filesystem."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

binary, fixture, root, schema = map(Path, sys.argv[1:])
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)
shutil.copytree(fixture, root / "input")


def invoke(*arguments, code=0):
    completed = subprocess.run([str(binary), *map(str, arguments)], capture_output=True, text=True)
    assert completed.returncode == code, (arguments, completed.returncode, completed.stderr)
    return json.loads(completed.stdout) if code == 0 else completed.stderr


completion = invoke("snapshot", "--input", root / "input", "--measurements", root / "input/measurements.json", "--out", root / "study")
initial = invoke("store-init", "--snapshot", root / "study", "--out", root / "store")
assert initial["root"]["sequence"] == "0"
contract = json.loads(schema.read_text())
assert contract["properties"]["expected"]["properties"]["source_namespace"]["const"] == initial["root"]["source_namespace"]
request = dict(schema="gadget.store.publish.v1", transaction="18446744073709551615", expected=initial["root"])
(root / "request.json").write_text(json.dumps(request))
published = invoke("store-publish", "--store", root / "store", "--snapshot", root / "study", "--request", root / "request.json")
assert published["root"]["sequence"] == "1"
retry = invoke("store-publish", "--store", root / "store", "--snapshot", root / "study", "--request", root / "request.json")
assert retry["receipt"] == published["receipt"]
assert invoke("store-status", "--store", root / "store", "--transaction", request["transaction"])["status"] == "committed"
assert invoke("store-status", "--store", root / "store", "--transaction", "17")["status"] == "absent"
audit = invoke("audit", "--store", root / "store")
assert audit.pop("store")["root"] == published["root"]
assert audit == invoke("audit", "--snapshot", root / "study")
shutil.copytree(root / "study", root / "spaced")
metadata = subprocess.run(
    [str(binary), "record-view", "--input", str(root / "study/study.rfcxs"),
     "--member", "study", "--format", "compact"],
    check=True, capture_output=True,
).stdout
for name, view in (("snapshot.json", json.dumps(completion).encode()), ("metadata.json", metadata)):
    (root / "spaced" / name).write_bytes(b"\n" + view + b"\n")
spaced = invoke("store-init", "--snapshot", root / "spaced", "--out", root / "spaced-store")
assert invoke("audit", "--store", root / "spaced-store").pop("store")["root"] == spaced["root"]
for name in ("snapshot.json", "metadata.json", "study.rfcxs"):
    assert (root / "spaced-store/studies/0000000000000000" / name).read_bytes() == (root / "spaced" / name).read_bytes()
shutil.copytree(root / "spaced", root / "native-only")
for name in ("snapshot.json", "metadata.json"):
    (root / "native-only" / name).unlink()
invoke("store-init", "--snapshot", root / "native-only", "--out", root / "native-only-store")
assert invoke("audit", "--store", root / "native-only-store")["source_digest"] == audit["source_digest"]
for arguments in [
    ("audit", "--store", root / "store", "--input", root / "input"),
    ("audit", "--store", root / "store", "--measurements", root / "input/measurements.json"),
    ("snapshot", "--store", root / "store", "--out", root / "recaptured"),
    ("store-status", "--store", root / "store", "--transaction", "01"),
    ("store-status", "--store", root / "store", "--transaction", "0"),
    ("store-status", "--store", root / "store", "--out", root / "irrelevant"),
    ("store-status", "--store", root / "store", "--store", root / "store"),
    ("store-publish", "--store", root / "store", "--snapshot", root / "study"),
    ("store-init", "--snapshot", root / "study", "--out", root / "store"),
]:
    invoke(*arguments, code=2)
print("PASS: durable CLI roundtrip, optional view bytes, native-only studies, retained history, uint64 maximum, exact retry, status and strict flags")
