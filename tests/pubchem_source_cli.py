"""Exercise the source adapter through its native query and result commands."""
import importlib.util
import json
import shutil
from pathlib import Path
import subprocess
import sys


def main(reader, parent, binary, witness):
    root = parent / "typed"
    def run(mode, *args):
        result = subprocess.run([str(reader), mode, *map(str, args)], capture_output=True, text=True, timeout=90)
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout)
    output = root / "query-result.rfcxs"; output.unlink(missing_ok=True)
    query = run("pubchem-query", root / "evidence.rfcxs", root / "query.json", output)
    assert query["selected"] == 6
    output = root / "context-result.rfcxs"; output.unlink(missing_ok=True)
    context = run("pubchem-result", root / "evidence.rfcxs", root / "dispatch.json", root / "analysis.rfcxs", output)
    assert context["features"] == json.loads((root / "expected.json").read_text())
    def gadget(*args):
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=90)
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout) if result.stdout.strip() else {}
    repository = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("linear_example", repository / "examples/bioinformatics/linear-models/replay.py")
    example = importlib.util.module_from_spec(spec); spec.loader.exec_module(example)
    for name in ("inputs", "workspace"):
        if (root / name).exists(): shutil.rmtree(root / name)
    declarations = example.stage(root / "inputs")
    gadget("import-assays", "--spec", root / "inputs/import.json", "--out", root / "workspace")
    defaults = dict(schema="gadget.pipeline.default.v1", profile="custom", question="Exact response fields and explicit typed assertions",
        goals=[dict(id="fit", method="linear_model_comparison", assay="measurements",
            parameters={key: value for key, value in declarations["categories"].items() if key != "operation"})])
    parameters = json.loads((root / "dispatch.json").read_text())
    parameters.pop("schema"); parameters.pop("operation")
    parameters["query"]["feature_path"] = ["result", "fields", "view", "strata", 0, "models", 1, "responses"]
    parameters["query"]["feature_id_path"] = ["field"]
    defaults["goals"].append(dict(id="evidence", method="feature_evidence_context", result_input="fit", parameters=parameters))
    plan = root / "pipeline.json"; plan.write_text(json.dumps(defaults) + "\n")
    directory = root / "pipeline"
    if directory.exists(): shutil.rmtree(directory)
    execution = run("pubchem-pipeline", root / "evidence.rfcxs", plan, root / "workspace/workspace.rfcxs", directory)
    fit, mapped = execution["results"]
    assert mapped["analysis_sha256"] == fit["analysis_record"]["sha256"]
    case = dict(directory=str(directory.resolve()), method="feature_evidence_context",
        result=dict(path=mapped["result"], sha256=mapped["result_sha256"]),
        request=dict(path=mapped["request"], sha256=mapped["request_sha256"]), analysis_sha256=mapped["analysis_sha256"])
    path = root / "witness.json"; path.write_text(json.dumps(case) + "\n")
    checked = subprocess.run([str(witness), "--witness", str(path)], capture_output=True, text=True, timeout=90)
    assert checked.returncode == 0, checked.stderr
    fields = json.loads(checked.stdout)["result"]["features"]
    assert [row["status"] for row in fields] == ["mapped", "unmapped", "mapped"]
    assert fields[0]["correspondences"][0]["evidence_identity"] == "pubchem:sid:100"
    print("PASS: native source query, feature correspondence and pipeline witnesses")


if __name__ == "__main__": main(*map(Path, sys.argv[1:]))
