"""Original source support, ordered provenance and native publication contracts."""
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import tempfile


def main(binary, consumer, root):
    root.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="source-fields-", dir=root))
    path = Path(__file__).resolve().parents[1]/"examples/bioinformatics/relational-ensembles/replay.py"
    spec = importlib.util.spec_from_file_location("source_fields_fixture", path)
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    fields = dict(A=[6, 2, 4, 8], zero=[0]*4, missing=[None]*4, partial=[None, 0, 8, None])
    fixture.authored(root/"inputs", field_values=fields)
    subprocess.run([str(binary), "import-assays", "--spec", str(root/"inputs/import.json"),
                    "--out", str(root/"workspace")], check=True, capture_output=True, text=True)
    study = next((root/"workspace").rglob("study.json"))
    declaration = json.loads(study.read_text())
    queries = {"inspect": dict(operation="inspect")}
    for name in fields:
        queries["coverage_"+name] = dict(operation="coverage", field="all."+name)
    checks = 0

    def execute(query):
        control = root/"control.json"
        control.write_text(json.dumps(query))
        process = subprocess.run([str(binary), "query", "--study", str(study), "--query", str(control)],
                                 check=True, capture_output=True, text=True)
        return json.loads(process.stdout)["result"]

    catalog = execute(queries["inspect"])
    assert catalog["study_id"] == declaration["id"]
    assert len(catalog["operations"]) == 34 and len(catalog["fields"]) == len(declaration["fields"])
    checks += 3
    for name, values in fields.items():
        actual = execute(queries["coverage_"+name])
        assert actual["coordinates"] == len(values) and actual["reported"] == sum(v is not None for v in values)
        for i, (row, expected) in enumerate(zip(actual["observations"], values, strict=True)):
            assert row["coordinate"] == i and row["reported"] == (expected is not None)
            assert row["value"] is None if expected is None else Fraction(row["value"]["rational"]) == expected
            checks += 3
    for index, cell in enumerate(declaration.get("cells", [])):
        query = dict(operation="provenance", entity=cell["id"])
        queries["provenance_"+str(index)] = query
        actual = execute(query)
        assert actual["entity"] == cell
        assert all(int(item["key"]) > 0 and Fraction(item["magnitude"]["rational"]) ==
                   abs(Fraction(item["coefficient"]["rational"])) for item in actual["occurrences"])
        checks += len(actual["occurrences"])+1
    control = root/"queries.json"
    control.write_text(json.dumps(queries))
    (root/"native").mkdir()
    for name, query in queries.items():
        request = root/("proposal-"+name+".json")
        request.write_text(json.dumps(dict(schema="gadget.scientific.agent.proposal.v1", action="query", query=query)))
        native = root/"native"/("cli-"+name+".rfcxs")
        process = subprocess.run([str(binary), "proposal", "--study", str(study), "--request", str(request), "--out", str(native)],
                                 check=True, capture_output=True, text=True)
        answer = json.loads(process.stdout)
        assert answer["status"] == "executed" and answer["certificate"]["result"] == execute(query) and native.is_file()
        checks += 3
    source = subprocess.run([str(binary), "source-index", "--study", str(study)],
                            check=True, capture_output=True, text=True)
    plan = dict(schema="gadget.scientific.plan.v1", id="proposal-source-plan", question="Retain every source result",
                source_digest=json.loads(source.stdout)["source_digest"],
                steps=[dict(id=name, query=query) for name, query in queries.items()])
    request = root/"proposal-plan.json"
    request.write_text(json.dumps(dict(schema="gadget.scientific.agent.proposal.v1", action="plan", plan=plan)))
    native = root/"native/cli-plan.rfcxs"
    process = subprocess.run([str(binary), "proposal", "--study", str(study), "--request", str(request), "--out", str(native)],
                             check=True, capture_output=True, text=True)
    answer = json.loads(process.stdout)
    assert answer["status"] == "executed" and native.is_file()
    assert [step["id"] for step in answer["certificate"]["steps"]] == list(queries)
    checks += 3
    process = subprocess.run([str(consumer), str(study), str(control), str(root/"native"), "study"],
                             capture_output=True, text=True)
    assert process.returncode == 0, (process.stdout, process.stderr)
    print(process.stdout.strip())
    print(json.dumps(dict(status="PASS", independent_checks=checks,
                          scope="Complete original authored support, source declarations and ordered incidence provenance")))


if __name__ == "__main__":
    main(*[Path(value).resolve() for value in sys.argv[1:]])
