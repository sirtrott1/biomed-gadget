"""Strict importer, source identity, missingness and canonical CPU action contracts."""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    binary, root = sys.argv[1], Path(sys.argv[2])
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    raw = "10,0,1,0\n20,1,1,?\n30,?,2,3\n"
    dictionary = "Synthetic numeric age and cholesterol; sex and chest pain are categorical codes.\n"
    (root / "table.txt").write_text(raw)
    (root / "dictionary.txt").write_text(dictionary)
    columns = [dict(id=name, channel="covariate" if name in ("sex", "cp") else "observation",
                    quantity="synthetic_" + name, units="synthetic_unit",
                    value_type="categorical" if name in ("sex", "cp") else "quantitative",
                    dictionary_start=0, dictionary_end=len(dictionary))
               for name in ("age", "sex", "cp", "chol")]
    spec = dict(schema="gadget.table.import.v1", id="table-contract", source_namespace="45099",
                source=dict(path="table.txt", sha256="sha256:" + hashlib.sha256(raw.encode()).hexdigest(),
                            delimiter=",", header=False),
                dictionary=dict(path="dictionary.txt", sha256="sha256:" + hashlib.sha256(dictionary.encode()).hexdigest()),
                columns=columns, missing_tokens=["?"], relation_columns=["sex", "cp"],
                question="Synthetic importer contract", population="three synthetic records",
                selection="every original row", limitations=["synthetic nonclinical fixture"],
                source_url="synthetic://table-contract", attribution="GADGET authored fixture", license="synthetic")
    spec_path = root / "spec.json"
    spec_path.write_text(json.dumps(spec) + "\n")

    def command(*arguments, success=True):
        result = subprocess.run([binary, *map(str, arguments)], capture_output=True, text=True, timeout=60)
        assert (result.returncode == 0) == success, (arguments, result.stderr)
        return json.loads(result.stdout) if success else None

    study = root / "study"
    command("import-table", "--spec", spec_path, "--out", study)
    assert (study / "table.txt").read_text() == raw
    declaration = json.loads((study / "study.json").read_text())
    aliases = {field["id"]: field["key"] for field in declaration["fields"]}
    assert aliases["all.chol"] == aliases["sex.chol"] == aliases["cp.chol"]

    def query(value, success=True, packet=study):
        path = root / "query.json"
        path.write_text(json.dumps(value) + "\n")
        return command("query", "--study", packet / "study.json", "--query", path, success=success)

    coverage = query(dict(operation="coverage", field="all.chol"))["result"]
    assert coverage["reported"] == 2 and coverage["missing"] == 1
    assert coverage["observations"][0]["value"]["rational"] == "0"
    assert coverage["observations"][1]["value"] is None
    summary = query(dict(operation="field_summary", field="all.chol", support="reported_only"))["result"]
    assert summary["mean"]["rational"] == "3/2" and summary["sum_of_squares"]["rational"] == "9"
    assert summary["coordinates"] == [0, 2]
    query(dict(operation="typed_moment", context="all", field="all.age", other_field="all.chol", order=0), success=False)
    query(dict(operation="field_summary", field="all.sex", support="reported_only"), success=False)
    native = query(dict(operation="native_boundary", field_set="category.memberships", backend="cpu_scheduled", tile_rows=1))
    assert [[value["rational"] for value in row["values"]] for row in native["result"]] == [["1", "1"], ["1", "1"], ["0", "1"]]
    assert native["execution"]["cpu_resources"]["worker_capacity"] >= native["execution"]["selected_workers"]
    before = (study / "study.rfcxs").read_bytes()
    command("import-table", "--spec", spec_path, "--out", study, success=False)
    assert (study / "study.rfcxs").read_bytes() == before
    broken = dict(spec)
    broken["source"] = {**spec["source"], "sha256": "sha256:" + "0" * 64}
    spec_path.write_text(json.dumps(broken))
    command("import-table", "--spec", spec_path, "--out", root / "broken", success=False)
    assert not (root / "broken").exists()
    raw = "10,0,1,?\n20,1,1,?\n30,?,2,?\n"
    (root / "table.txt").write_text(raw)
    spec["source"]["sha256"] = "sha256:" + hashlib.sha256(raw.encode()).hexdigest()
    spec_path.write_text(json.dumps(spec))
    command("import-table", "--spec", spec_path, "--out", root / "absent")
    absent = query(dict(operation="field_summary", field="all.chol", support="reported_only"), packet=root / "absent")["result"]
    assert absent["reported"] == 0 and absent["mean"] is None and absent["missing"] == 3
    print("PASS: exact table admission, shared tensor aliases, complete rows, zero/absence, P5 membership action and refusal boundaries")


if __name__ == "__main__":
    main()
