"""Original source comparisons and native data to exact model execution."""
import csv
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools/rcdb"))
spec = importlib.util.spec_from_file_location("source_intake", ROOT / "examples/source-formats/import_native.py")
intake = importlib.util.module_from_spec(spec); spec.loader.exec_module(intake)
from gadget_data import DataStore, Limits, SourceResolver
from gadget_data.bridge import export_gadget
from rexgraph.value import Absent


def main(reader, binary, importer, witness, root):
    root.mkdir(parents=True, exist_ok=True)
    def save(name, value):
        path = root / name
        path.write_text(value if isinstance(value, str) else json.dumps(value) + "\n")
        return path
    def sha(path): return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()
    def run(*args, fail=False):
        result = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True, timeout=90)
        assert bool(result.returncode) == fail, result.stderr
        return None if fail else json.loads(result.stdout) if result.stdout.strip() else {}
    for name in ("counts", "binding", "invalid", "late"):
        if (root / name).exists(): shutil.rmtree(root / name)
    samples = ("GTEX-AA-1", "GTEX-AA-2", "GTEX-BB-1", "GTEX-CC-1")
    matrix = save("counts.gct", "#1.2\n3\t4\nName\tDescription\t" + "\t".join(samples) + "\n" +
        "ENSG0001.1\tfirst\t1\t3\t5\t7\nENSG0002.1\tsecond\t0\t2\tNA\t6\n" +
        "ENSG0003.1\tlarge\t9007199254740993\t9007199254740995\t9007199254740997\t9007199254740999\n")
    attributes = save("samples.tsv", "SAMPID\tSMTSD\tannotation\n" +
        "".join(sample + "\tdeclared tissue\tNA\n" for sample in samples))
    subjects = save("donors.tsv", "SUBJID\tSEX\tAGE\nGTEX-AA\t1\t20-29\nGTEX-BB\t2\tNA\nGTEX-CC\t1\t60-69\n")
    request = dict(layout="gct12", feature="gene", measurement="integer_counts",
        sample_attributes=str(attributes), subject_phenotypes=str(subjects),
        limits=dict(source_bytes=1048576, rows=10000, columns=10000, record_bytes=1048576,
                    field_bytes=65536, cells=1000000))
    limits = Limits(record_bytes=1048576, decode_bytes=8388608)
    intake.import_measurements(matrix, root / "counts", reader=reader, profile="gtex", request=request,
                               name="authored counts", limits=limits, tile_rows=2, tile_columns=2)
    lines = list(csv.reader(matrix.open(), delimiter="\t"))
    truth = [[None if token == "NA" else Q(token) for token in row[2:]] for row in lines[3:]]
    with DataStore(root / "counts", limits=limits) as store:
        assert store.verify_sources(SourceResolver(root))
        selected = store.select("measurements", rows=(0, 1, 2), columns=(3, 0, 2, 1))
        for row in range(3):
            for column, original in enumerate((3, 0, 2, 1)):
                assert bool(selected.observed[row, column]) == (truth[row][original] is not None)
                if selected.observed[row, column]: assert selected.values[row, column] == truth[row][original]
        assert selected.row_labels == ("ENSG0001.1", "ENSG0002.1", "ENSG0003.1")
        assert [row["donor"] for row in selected.column_context] == ["GTEX-CC", "GTEX-AA", "GTEX-BB", "GTEX-AA"]
        assert selected.column_context[1]["donor/AGE"] == "20-29"
        data = root / "data.rfcxs"
        data.unlink(missing_ok=True)
        export_gadget(selected, data, importer=importer, namespace=806201, max_bytes=8388608)
    binding = dict(schema="gadget.data.study.v1", data_sha256=sha(data), id="authored count model",
        design=dict(question="Which explicitly declared field represents these counts?", population="four original observations",
                    selection="all four coordinates", limitations=["Repeated donors retained; algebraic fit has no sampling independence claim"]),
        responses=[dict(id="first", row=0), dict(id="partial", row=1), dict(id="large", row=2)],
        covariates=[dict(id=name, values=["3", "0", "2", "1"], quantity="authored design", units="dimensionless", value_type="quantitative")
                    for name in ("x", "x_copy")], max_bytes=8388608)
    study = root / "study.rfcxs"; study.unlink(missing_ok=True)
    request_path = save("binding.json", binding)
    run("bind-data-study", "--input", data, "--request", request_path, "--out-record", study)
    query = dict(operation="linear_model_comparison", responses=["first", "partial", "large"],
        terms=[dict(id=name, factors=[name]) for name in ("x", "x_copy")], controls=[], selection=[],
        support="joint_reported", metric="identity", intercept=True,
        coefficient_contrasts=[dict(id="combined", weights=[dict(column="x", value="1"), dict(column="x_copy", value="1")]),
                              dict(id="single", weights=[dict(column="x", value="1")])])
    model = root / "model.rfcxs"; model.unlink(missing_ok=True)
    answer = run("query", "--study", study, "--query", save("model.json", query), "--out", model)
    comparison = answer["result"]; full = comparison["strata"][0]["models"][1]
    assert comparison["strata"][0]["coordinates"] == [0, 1, 3]
    assert full["rank"] == 2 and full["coefficient_nullity"] == 1
    for response in full["responses"]:
        assert response["residual_quadrance"] == "0"
        assert response["coefficient_contrasts"][0]["value"] == "2"
        assert response["coefficient_contrasts"][1]["value"] is None
    assert list(map(Q, full["responses"][2]["input"])) == [Q(9007199254740999), Q(9007199254740993), Q(9007199254740995)]
    driver_spec = importlib.util.spec_from_file_location("native_count_analysis", ROOT / "examples/source-formats/analyze_counts.py")
    driver = importlib.util.module_from_spec(driver_spec); driver_spec.loader.exec_module(driver)
    driver_output = root / "driver"
    if driver_output.exists(): shutil.rmtree(driver_output)
    selected_binding = {name: binding[name] for name in ("id", "design", "responses", "covariates")}
    selected_binding.update(rows=[0, 1, 2], columns=[3, 0, 2, 1])
    retained = driver.analyze(root / "counts", driver_output, binding=selected_binding, query=query,
        binary=binary, importer=importer, namespace=806202, limits=limits)
    driven = run("record-view", "--input", retained, "--member", "result")
    assert driven["result"]["strata"][0]["models"][1]["responses"] == full["responses"]
    wrong = dict(binding, data_sha256="sha256:" + "0" * 64)
    run("bind-data-study", "--input", data, "--request", save("wrong.json", wrong), "--out-record", root / "wrong.rfcxs", fail=True)
    assert not (root / "wrong.rfcxs").exists()
    limited = dict(binding, max_bytes=1048576, covariates=[dict(binding["covariates"][0], values=["9" * 65536] * 4)])
    run("bind-data-study", "--input", data, "--request", save("limited.json", limited), "--out-record", root / "limited.rfcxs", fail=True)
    assert not (root / "limited.rfcxs").exists()
    # Native model response identities are selected through array coordinates.
    gmt = "one\tauthored family\tA\tB\ntwo\tother family\tC\n"
    mapping = "feature_id\tmember\tevidence\nfirst\tA\tdeclared\npartial\tB\tdeclared\nlarge\tC\tdeclared\n"
    sets = dict(schema="gadget.pathways.source.v1", id="authored", release="one", identifier_namespace="declared:member",
        format="gmt", data=gmt, sha256="sha256:" + hashlib.sha256(gmt.encode()).hexdigest(),
        mapping=dict(id="authored", release="one", origin="explicit_declared", feature_namespace="model:response",
                     identifier_namespace="declared:member", format="tsv", data=mapping,
                     sha256="sha256:" + hashlib.sha256(mapping.encode()).hexdigest()))
    source_record = root / "sets.rfcxs"; source_record.unlink(missing_ok=True)
    run("native-record", "--input", save("sets.json", sets), "--out", source_record)
    biological = dict(schema="gadget.pathways.overrepresentation.v1", operation="pathway_overrepresentation",
        gene_sets_sha256=sha(source_record), feature_namespace="model:response",
        features_path=["result", "fields", "view", "strata", 0, "models", 1, "responses"], feature_id_path=["field"],
        background_indices=[0, 1, 2], selected_indices=[0, 1],
        mapping_policy=dict(ambiguous="refuse", unmapped="refuse", many_to_one="refuse", duplicate_occurrences="retain_count_once"),
        set_ids=["one", "two"], size_filter=dict(minimum=1, maximum=3, basis="members_in_background"),
        membership_policy=dict(duplicates="retain_count_once", empty="retain_exclude", outside_background="retain_exclude"),
        correction="bonferroni", selection_method="Explicit fixture response selection",
        null_assumptions=["Uniform subset of the mapped background conditional on selection size"],
        limits=dict(source_bytes=1048576, mapping_rows=100, membership_cells=1000, sets=100, features=100,
                    field_bytes=65536, kernel_work=100000, archive_bytes=8388608, declared_buffer_bytes=8388608))
    operation = dict(schema="gadget.result-operation.v1", operation="pathway_overrepresentation", analysis_binding="actual_result",
        analysis_bytes=8388608, source_record=dict(path=str(source_record), sha256=sha(source_record), max_bytes=8388608), query=biological)
    output = root / "pathway.rfcxs"; output.unlink(missing_ok=True)
    pathway = run("analyze-record", "--analysis", model, "--request", save("pathway.json", operation), "--out-record", output)
    assert Q(pathway["tests"][0]["p_value"]) == Q(1, 3)
    assert Q(pathway["tests"][0]["adjusted_p_value"]) == Q(2, 3)
    # Replay uses retained native owners while foreign source files are offline.
    replayed = root / "replayed.rfcxs"; replayed.unlink(missing_ok=True)
    for original in (matrix, attributes, subjects): original.rename(original.with_suffix(original.suffix + ".offline"))
    try:
        run("replay-scientific", "--record", model, "--out-record", replayed)
        restored = run("record-view", "--input", replayed, "--member", "result")
        assert restored["result"]["strata"][0]["models"][1]["responses"] == full["responses"]
    finally:
        for original in (matrix, attributes, subjects): original.with_suffix(original.suffix + ".offline").rename(original)
    # Reuse the existing assay pipeline and select the actual native fit leaf.
    demo_spec = importlib.util.spec_from_file_location("linear_example", ROOT / "examples/bioinformatics/linear-models/replay.py")
    demo = importlib.util.module_from_spec(demo_spec); demo_spec.loader.exec_module(demo)
    for directory in (root / "pipeline-inputs", root / "workspace", root / "pipeline"):
        if directory.exists(): shutil.rmtree(directory)
    declarations = demo.stage(root / "pipeline-inputs")
    run("import-assays", "--spec", root / "pipeline-inputs/import.json", "--out", root / "workspace")
    population = ["all.negative", "all.zero", "all.positive"]
    pipeline_sets = dict(sets)
    pipeline_sets["mapping"] = dict(sets["mapping"], data="feature_id\tmember\tevidence\n" +
        "".join(feature + "\t" + member + "\tdeclared\n" for feature, member in zip(population, "ABC", strict=True)))
    pipeline_sets["mapping"]["sha256"] = "sha256:" + hashlib.sha256(pipeline_sets["mapping"]["data"].encode()).hexdigest()
    family = root / "pipeline-sets.rfcxs"; family.unlink(missing_ok=True)
    run("native-record", "--input", save("pipeline-sets.json", pipeline_sets), "--out", family)
    pathway_parameters = dict(operation)
    pathway_parameters.pop("schema"); pathway_parameters.pop("operation")
    pathway_parameters["source_record"] = dict(path=str(family.resolve()), sha256=sha(family), max_bytes=8388608)
    pathway_parameters["query"] = dict(biological, gene_sets_sha256=sha(family))
    defaults = dict(schema="gadget.pipeline.default.v1", profile="custom", question="Exact field fit and declared set family",
        goals=[dict(id="fit", method="linear_model_comparison", assay="measurements", parameters={key: value for key, value in declarations["categories"].items() if key != "operation"}),
               dict(id="sets", method="pathway_overrepresentation", result_input="fit", parameters=pathway_parameters)])
    compiled = run("pipeline-prepare", "--workspace", root / "workspace/workspace.rfcxs", "--request", save("defaults.json", defaults))
    executed = run("pipeline", "--workspace", root / "workspace/workspace.rfcxs", "--request", save("pipeline.json", compiled), "--out", root / "pipeline")
    fit_receipt, set_receipt = executed["results"]
    assert fit_receipt["analysis_scope"] == "one_native_result"
    assert fit_receipt["analysis_record"]["sha256"] == set_receipt["analysis_sha256"]
    carried = run("record-view", "--input", root / "pipeline/sets/result.rfcxs", "--member", "result")
    assert Q(carried["tests"][0]["p_value"]) == Q(1, 3)
    witness_case=dict(directory=str((root/"pipeline").resolve()), method="pathway_overrepresentation",
        result=dict(path=set_receipt["result"],sha256=set_receipt["result_sha256"]),
        request=dict(path=set_receipt["request"],sha256=set_receipt["request_sha256"]),
        analysis_sha256=set_receipt["analysis_sha256"])
    checked=subprocess.run([str(witness),"--witness",str(save("witness.json",witness_case))],capture_output=True,text=True,timeout=60)
    assert checked.returncode==0,checked.stderr
    # A late invalid original row must discard every staged native block.
    bad = save("invalid.gct", matrix.read_text().replace("9007199254740999", "-1"))
    try:
        intake.import_measurements(bad, root / "late", reader=reader, profile="gtex", request=request,
                                   name="late invalid source", limits=limits, tile_rows=1, tile_columns=1)
    except ValueError: pass
    else: raise AssertionError("late invalid source published")
    assert not (root / "late").exists()
    header = ["BindingDB Reactant_set_id", "BindingDB MonomerID", "Target Name",
              "Target Source Organism According to Curator or DataSource", "Ki (nM)", "Kd (nM)", "IC50 (nM)", "EC50 (nM)"]
    rows = [["7", "1", "reported construct", "organism", "0", "", ">=1e-40", "~2"],
            ["7", "1", "reported construct", "organism", "0", "", ">=1e-40", "~2"]]
    bdb = save("binding.tsv", "\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in rows))
    intake.import_measurements(bdb, root / "binding", reader=reader, profile="bindingdb",
        request=dict(limits=dict(source_bytes=1048576, rows=100, columns=100, field_bytes=65536)),
        name="authored thresholds", limits=limits, tile_rows=1, tile_columns=2)
    with DataStore(root / "binding", limits=limits) as store:
        values = store.select("measurements", rows=(0, 1), columns=(0, 1, 2, 3))
        comparisons = store.select("comparison", rows=(0, 1), columns=(0, 1, 2, 3))
        assert values.row_labels == (0, 1)
        assert list(values.observed[0]) == [True, False, True, True]
        assert list(values.values[0, [0, 2, 3]]) == [0, Q(1, 10**40), 2]
        assert list(comparisons.values[0, [0, 2, 3]]) == [1, 5, 6]
        assert values.row_context[0]["compound"] == "bindingdb:monomer:1"
        assert values.row_context[0]["reported_construct"] != values.row_context[1]["reported_construct"]
        assert len(tuple(store.find_axis("features", context={"compound": "bindingdb:monomer:1"}))) == 2
        assert store.verify()
    try: intake.coefficient("1e99999999999999999999999999999", 100)
    except MemoryError: pass
    else: raise AssertionError("unbounded exact decimal allocation admitted")
    assert intake.coefficient("0e99999999999999999999999999999", 100) == 0
    assert intake.coefficient("9" * 5000, 6000, integer=True) == 10**5000 - 1
    print(json.dumps(dict(status="PASS", scope="Original source coordinates, exact counts, ambiguity, support, pathway law, source pins and atomic failure")))


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
