"""Run existing exact GADGET analyses on selected native measurement fields."""
import argparse
import json
from pathlib import Path
import subprocess

from gadget_data import DataStore, Limits
from gadget_data._records import fingerprint
from gadget_data.bridge import export_gadget
from rexgraph.io.publication import staged_publication


def analyze(store_path, output, *, binding, query, binary, importer, namespace,
            field="measurements", limits=None):
    """The binding declares responses, observation order and covariate values."""
    limits = limits or Limits()
    if query.get("operation") not in ("linear_model_comparison", "observed_response_model",
            "constrained_response_model", "finite_reference_calibration", "relational_field_ensemble"):
        raise ValueError("select an existing exact scientific operation")
    if set(binding) != {"id", "design", "rows", "columns", "responses", "covariates"}:
        raise ValueError("explicit measurement binding members required")
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with staged_publication(output, directory=True) as stage:
        stage = Path(stage)
        with DataStore(store_path, limits=limits) as source:
            selected = source.select(field, rows=tuple(binding["rows"]), columns=tuple(binding["columns"]))
            export_gadget(selected, stage / "data.rfcxs", importer=importer, namespace=namespace,
                          max_bytes=limits.selection_bytes)
        request = {"schema": "gadget.data.study.v1", "data_sha256": "sha256:" + fingerprint(stage / "data.rfcxs")["sha256"],
            **{name: binding[name] for name in ("id", "design", "responses", "covariates")},
            "max_bytes": limits.selection_bytes}
        (stage / "binding.json").write_text(json.dumps(request) + "\n")
        (stage / "query.json").write_text(json.dumps(query) + "\n")
        def run(*arguments):
            with (stage / "display.json").open("w") as display:
                subprocess.run([str(binary), *map(str, arguments)], check=True, stdout=display)
        run("bind-data-study", "--input", stage / "data.rfcxs", "--request", stage / "binding.json",
            "--out-record", stage / "study.rfcxs")
        run("query", "--study", stage / "study.rfcxs", "--query", stage / "query.json", "--out", stage / "result.rfcxs")
    return output / "result.rfcxs"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--binding", type=Path, required=True); parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True); parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--namespace", type=int, required=True)
    args = parser.parse_args()
    print(analyze(args.store, args.output, binding=json.loads(args.binding.read_text()),
        query=json.loads(args.query.read_text()), binary=args.binary, importer=args.importer, namespace=args.namespace))


if __name__ == "__main__":
    main()
