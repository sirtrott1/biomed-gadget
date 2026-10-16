"""Native biomedical owners with independent original coordinate comparisons."""
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main(binary, consumer, parent):
    parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="biomedical-fields-", dir=parent))
    example = Path(__file__).resolve().parents[1] / "examples/bioinformatics/relational-gluing"

    def module(name):
        spec = importlib.util.spec_from_file_location(name, example / (name + ".py"))
        value = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(value)
        return value

    def save(name, value):
        path = root / name
        path.write_text(value if isinstance(value, str) else json.dumps(value) + "\n")
        return path

    def run(*args, fail=False):
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=300)
        assert bool(result.returncode) == fail, (args, result.stdout, result.stderr)
        return result.stderr if fail else json.loads(result.stdout)

    def native(*args):
        result = subprocess.run([str(consumer), *map(str, args)], capture_output=True, text=True, timeout=300)
        assert result.returncode == 0, (result.stdout, result.stderr)
        return json.loads(result.stdout)

    spec, gluing = module("make_example").make_example(root / "inputs", repetitions=2)
    raw = root / "inputs/classes.txt"
    raw.write_text("3\n3\nx\ny\nz\n3\t0\t1\t2\t0\n2\t1\t2\t4\n1\t0\t2\n")
    spec["resources"].append(dict(id="classes", format="salmon_eq", method="Original zero and full arity counts",
        source=dict(path=raw.name, sha256="sha256:" + hashlib.sha256(raw.read_bytes()).hexdigest())))
    (root / "inputs/import.json").write_text(json.dumps(spec))
    run("import-assays", "--spec", root / "inputs/import.json", "--out", root / "workspace")
    workspace = root / "workspace/workspace.rfcxs"
    reference = json.loads((root / "inputs/reference.json").read_text())
    left = reference["orders"]["left"]
    right = reference["orders"]["common"]
    mapping = dict(domain_quantity="declared_observation_signal", domain_units="1",
        codomain_quantity="declared_observation_signal", codomain_units="1", method="Original rectangular map",
        evidence="Every signed coefficient supplied", entries=[
            dict(left_id=left[0], right_id=right[2], coefficient="-1/2"),
            dict(left_id=left[1], right_id=right[0], coefficient="0"),
            dict(left_id=left[2], right_id=right[1], coefficient="2")])
    transfer = dict(schema="gadget.assays.query.v1", operation="field_transfer", metric="identity", map=mapping,
        left=dict(assay="left", partition=0, field="all.signal-0"),
        right=dict(assay="common", partition=0, field="all.signal-0"))
    allocation = dict(schema="gadget.assays.query.v1", operation="allocation_bounds", resource="classes",
        targets=[dict(id="x", candidates=["x"]), dict(id="yz", candidates=["y", "z"]),
                 dict(id="all", candidates=["x", "y", "z"])], kernel_certificate=True,
        candidate_search=dict(method="Original equivalence classes", completeness="Declared export only",
            reference_sha256="sha256:" + "0" * 64), scale="3/2", scale_provenance="Declared sample scale")
    requests = dict(gluing=gluing, transfer=transfer, allocation=allocation)
    requests["same_dimensions"] = copy.deepcopy(transfer)
    requests["same_dimensions"]["right"]["assay"] = "right"
    requests["same_dimensions"]["map"]["entries"] = [dict(left_id=l, right_id=r, coefficient="1")
        for l, r in zip(left, reference["orders"]["right"])]
    queries = save("requests.json", requests)
    receipt = native(workspace, queries, root / "native")
    checks = receipt["native_checks"]
    glued = json.loads((root / "native/gluing.json").read_text())["result"]
    module("verify").verify(glued, reference, workspace.parent / "workspace.json")
    mapped = json.loads((root / "native/transfer.json").read_text())["result"]
    expected = [Q(-10), Q(0), Q(20), Q(0)]
    l = [Q(0), Q(10), Q(20), Q(30)]
    assert list(map(Q, mapped["transported_right"])) == expected
    assert list(map(Q, mapped["residual"])) == [a - b for a, b in zip(l, expected)]
    assert Q(mapped["residual_quadrance"]) == sum((a - b)**2 for a, b in zip(l, expected))
    assert Q(mapped["pairing"]) == sum(a*b for a, b in zip(l, expected))
    assert mapped["map_rank"] == 2 and mapped["map_kernel_dimension"] == 1
    bounds = json.loads((root / "native/allocation.json").read_text())["result"]["samples"][0]
    assert Q(bounds["total"]) == 6
    assert [(Q(b["lower"]), Q(b["upper"])) for b in bounds["bounds"]] == [(Q(2), Q(2)), (Q(4), Q(4)), (Q(6), Q(6))]
    checks += 8
    for change in [lambda q: q["map"]["entries"][0].update(right_id=left[0]),
                   lambda q: q["map"].update(domain_units="wrong"),
                   lambda q: q["map"]["entries"].append(q["map"]["entries"][0])]:
        invalid = copy.deepcopy(transfer)
        change(invalid)
        run("assay-query", "--workspace", workspace, "--request", save("invalid.json", invalid), fail=True)
        checks += 1
    question = "Retain the declared glued field with its complete source maps"
    intent = dict(schema="gadget.pipeline.intent.v1", profile="multi_omics",
        recipes=[dict(id="glue", purpose="Original source pushout", goal=dict(method="relational_gluing",
            parameters={k: v for k, v in gluing.items() if k not in {"schema", "operation"}}))],
        answer_contract=dict(schema="gadget.pipeline.answer.contract.v1", question=question,
            requirements=[dict(id="source-bound-glue", observable="relational_gluing", bindings=["glue"])]))
    intent_path = save("intent.json", intent)
    catalog = run("recipe-catalog", "--workspace", workspace, "--intent", intent_path)
    proposal = dict(schema="gadget.pipeline.recipe.proposal.v1", action="recipes", goals=["glue"],
        source_digest=catalog["source_digest"], intent_sha256=catalog["intent_sha256"])
    planning = run("recipe-proposal", "--workspace", workspace, "--intent", intent_path,
        "--request", save("proposal.json", proposal), "--prompt", question, "--out", root / "pipeline")
    assert planning["status"] == "admitted"
    planning_path = root / "pipeline" / planning["execution"]["planning"]
    checks += native("--planning", planning_path)["native_checks"]
    # Archives contain all original sources; remove former intake locations from
    # view before a fresh process reopens and reexecutes every native result.
    (root / "inputs").rename(root / "moved-inputs")
    (root / "workspace").rename(root / "moved-workspace")
    checks += native("--planning", planning_path)["native_checks"]
    for name in requests:
        checks += native("--replay", root / "native" / (name + ".rfcxs"))["native_checks"]
    print(json.dumps(dict(status="PASS", checks=checks,
        scope="Original native pushout, rectangular and equal dimension independent carriers, zero counts, complete archives and relocated replay")))


if __name__ == "__main__":
    main(*[Path(value).resolve() for value in sys.argv[1:]])
