"""Independent original row oracles for native legacy analysis."""
import csv
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main(binary, consumer, fixture, parent):
    parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="native-analysis-", dir=parent))
    shutil.copytree(fixture, root/"input")
    source = root/"input"
    checks = 0
    def run(*args):
        p = subprocess.run(list(map(str, args)), text=True, capture_output=True)
        assert p.returncode == 0, (p.stdout, p.stderr)
        return p.stdout
    run(binary, "analyze", "--input", source, "--measurements", source/"measurements.json", "--out", root/"analysis")
    print(run(consumer, root/"analysis/analysis.rfcxs", source, source/"measurements.json").strip())
    context_paths = {"AD": "pag data/Alzheimers", "HD": "pag data/Huntingtons/huntington_all_pager_files"}
    get_paths = {"AD": "gene lists/Alzheimer_Disease_gene_analysis 1.csv", "HD": "gene lists/Huntington_Disease_gene_analysis_new 1.csv"}
    original_genes = set()
    participants, observed = {}, {}
    for context in ["AD", "HD"]:
        observed[context] = list(csv.DictReader((source/get_paths[context]).open()))
        original_genes.update(row["Gene"] for row in observed[context])
        participants[context] = []
        for path in sorted((source/context_paths[context]).glob("*.csv")):
            grouped = {}
            for row in csv.DictReader(path.open()):
                gene, value = row["Gene symbol"], Q(row["RP_score"])
                assert gene not in grouped or grouped[gene] == value
                grouped[gene] = value
                original_genes.add(gene)
            participants[context].append(grouped)
    genes = sorted(original_genes)
    def dot(x, y): return sum((a*b for a, b in zip(x, y)), Q(0))
    def apply(a, x): return [dot(row, x) for row in a]
    def close(actual, expected):
        nonlocal checks
        assert math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=1e-10), (actual, expected)
        checks += 1
    fields = {c: [[Q(any(row["Gene"] == g and channel in row["GET"] for row in observed[c])) for g in genes] for channel in "GET"] for c in ["AD", "HD"]}
    operators = {}
    for realization in ["support", "rp"]:
        for context in ["AD", "HD"]:
            columns = [[Q(g in p) if realization == "support" else p.get(g, Q(0)) for g in genes] for p in participants[context]]
            metrics = [dot(col, col) or Q(1) for col in columns]
            operator = [[sum((col[i]*col[j]/metric for col, metric in zip(columns, metrics)), Q(0)) for j in range(len(genes))] for i in range(len(genes))]
            operators[realization, context] = operator
            powers = fields[context]
            prefix = ("alzheimers" if context == "AD" else "huntingtons") + "_" + realization
            for power in range(5):
                moment = [[dot(fields[context][i], powers[j]) for j in range(3)] for i in range(3)]
                for normalized in [False, True]:
                    filename = prefix + ("_typed_moment_cosine_p" if normalized else "_typed_moment_p") + str(power) + ".csv"
                    rows = list(csv.reader((root/"analysis/tables"/filename).open()))[1:]
                    for i, row in enumerate(rows):
                        for j, actual in enumerate(row[1:]):
                            denominator = math.sqrt(float(moment[i][i]*moment[j][j]))
                            close(actual, float(moment[i][j])/denominator if normalized and denominator else 0 if normalized else moment[i][j])
                powers = [apply(operator, f) for f in powers]
        ad, hd = operators[realization, "AD"], operators[realization, "HD"]
        rows = list(csv.DictReader((root/"analysis/tables"/("gene_role_differential_"+realization+".csv")).open()))
        assert [row["gene"] for row in rows] == genes
        for gene, row in enumerate(rows):
            a, h = [r[gene] for r in ad], [r[gene] for r in hd]
            for name, expected in {"mu1_AD": ad[gene][gene], "mu1_HD": hd[gene][gene], "mu2_AD": dot(a,a), "mu2_HD": dot(h,h), "cross_AD_HD": dot(a,h), "delta_L_squared": dot([x-y for x,y in zip(a,h)], [x-y for x,y in zip(a,h)])}.items():
                close(row[name], expected)
    for item in json.loads((root/"analysis/manifests/artifacts.json").read_text())["artifacts"]:
        path = root/"analysis"/item["path"]
        with path.open("rb") as stream: digest = hashlib.file_digest(stream, "sha256").hexdigest()
        assert item["bytes"] == path.stat().st_size and item["sha256"] == "sha256:"+digest
        checks += 2
    assert not (root/"analysis.partial").exists()
    print(json.dumps(dict(status="PASS", independent_checks=checks, scope="All original authored gene rows, repeated primary channels, five complete moment orders, normalized matrices, full differential sweeps and streamed artifact identities")))


if __name__ == "__main__":
    main(*[Path(value).resolve() for value in sys.argv[1:]])
