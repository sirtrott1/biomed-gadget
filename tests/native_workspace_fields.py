"""Independent original coordinate oracles for native workspace publications."""
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from research_cli import rref


def main(binary, consumer, parent):
    parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="workspace-fields-", dir=parent))
    checks = 0

    def save(name, value):
        path = root/name
        path.write_text(json.dumps(value) + "\n")
        return path

    def file(name, value):
        path = root/name
        path.write_text(value)
        return dict(path=name, sha256="sha256:"+hashlib.sha256(path.read_bytes()).hexdigest())

    def run(*args, fail=False):
        result = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True)
        if fail:
            assert result.returncode, result.stdout
            return
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout)

    participant = lambda cell, value: dict(cell=cell, exists=True, negative=Q(value)<0, magnitude=str(abs(Q(value))))
    edges = [dict(id=f"e{i}", meaning="Original oriented square edge",
                  participants=[participant(f"v{i}", -1), participant(f"v{(i+1)%4}", 1)]) for i in range(4)]
    relations = edges + [dict(id="face", grade=2, meaning="Twice the original square cycle",
                             participants=[participant(f"e{i}", 2) for i in range(4)])]
    design = dict(question="Do actual owners survive workspace selection and persistence?",
                  population="Authored exact source contracts", selection="Every original row and occurrence",
                  limitations=["Computational witness qualification"])
    assay = dict(id="torsion", kind="quantitative", format="dense_tsv",
                 source=file("matrix.tsv", "id\tpartial\tzero\n"+"".join(f"v{i}\t{['NA','0','8','NA'][i]}\t0\n" for i in range(4))),
                 features=file("features.tsv", "id\nv0\nv1\nv2\nv3\n"),
                 observations=file("observations.tsv", "id\npartial\nzero\n"),
                 relations=file("relations.json", json.dumps(dict(relations=relations))),
                 source_namespace="942000", carrier="features", quantity="authored_value", units="units",
                 channel="observation", method="Original exact control rows", arithmetic="exact_rational",
                 implicit_entries="unobserved", missing_tokens=["NA"])
    partitioned = copy.deepcopy(assay)
    partitioned.update(id="partitioned", source_namespace="943000", partition_size=2)
    partitioned.pop("relations")
    run("import-assays", "--spec", save("import.json", dict(schema="gadget.assays.import.v1", id="native-workspace",
        design=design, assays=[assay, partitioned], resources=[], links=[])), "--out", root/"workspace")
    workspace = root/"workspace/workspace.rfcxs"
    births = [dict(grade=0, index=i, level=0) for i in range(4)] + [dict(grade=1, index=i, level=1) for i in range(4)] + [dict(grade=2,index=0,level=2)]
    persistence = dict(schema="gadget.assays.query.v1", operation="relational_persistence", assay="torsion", partition=0,
                       levels=["0","1","2"], births=births, lattice_snapshots=True,
                       filtration_meaning="Original vertices then edges then twice-loop filling", evidence="All original cell births")
    native = lambda query: dict(schema="gadget.assays.query.v1", operation="native_query", assay="torsion", query=query)
    queries = dict(inspect=dict(schema="gadget.assays.query.v1", operation="inspect"),
                   coverage=native(dict(operation="coverage",field="all.partial")),
                   summary=native(dict(operation="field_summary",field="all.zero",support="reported_only")),
                   persistence=persistence)
    queries["partition_order"]=dict(schema="gadget.assays.query.v1", operation="native_query", assay="partitioned",
                                   partitions=[1,0], query=dict(operation="coverage",field="all.partial"))
    empty = copy.deepcopy(persistence)
    empty["levels"]=["-1","0","1","2","3"]
    for b in empty["births"]: b["level"]+=1
    queries["empty_stages"]=empty
    queries["batch"]=dict(schema="gadget.assays.query.v1",operation="query_batch",concurrent_jobs=2,
                          queries=[dict(id=name,query=queries[name]) for name in ["coverage","inspect","persistence"]])
    queries["partition_run"] = dict(queries["partition_order"], concurrent_jobs=2)
    control=save("queries.json",queries)
    process=subprocess.run([str(consumer),str(workspace),str(control),str(root/"native")],text=True,capture_output=True)
    assert process.returncode==0, (process.stdout,process.stderr)
    print(process.stdout.strip())
    for name in ["persistence","empty_stages"]:
        answer=json.loads((root/"native"/(name+".json")).read_text())["result"]
        request=queries[name]
        boundary={0:[],1:[[-1,0,0,1],[1,-1,0,0],[0,1,-1,0],[0,0,1,-1]],2:[[2],[2],[2],[2]]}
        for t,stage in enumerate(answer["stages"]):
            selected={k:[b["index"] for b in request["births"] if b["grade"]==k and b["level"]<=t] for k in range(3)}
            for k,grade in enumerate(stage["grades"]):
                expected=[[Q(boundary[k][r][c]) for c in selected[k]] for r in selected[k-1]] if k else []
                n=len(selected[k]); upper=[[Q(boundary[k+1][r][c]) for c in selected[k+1]] for r in selected[k]] if k<2 else [[] for _ in selected[k]]
                rank=lambda a:len(rref(a)[1]) if a and a[0] else 0
                assert grade["betti"]==n-rank(expected)-rank(upper)
                assert grade["ambient_coordinates"]==selected[k]
                entries=[[i,j,str(x)] for i,row in enumerate(expected) for j,x in enumerate(row) if x]
                assert grade["boundary"]==dict(rows=len(expected),columns=n,entries=entries)
                checks+=3
        assert answer["stages"][-1]["grades"][1]["integral_homology"]["torsion_invariants"]==["2"]
        assert all(b["multiplicity"]>0 for b in answer["barcode"])
        checks+=2
    coverage=json.loads((root/"native/coverage.json").read_text())["result"]["answers"][0]["certificate"]["result"]
    assert [r["reported"] for r in coverage["observations"]]==[False,True,True,False]
    assert coverage["observations"][1]["value"]["rational"]=="0"
    checks+=2
    for mutation in [lambda q:q["births"][-1].update(level=0),lambda q:q["levels"].reverse(),lambda q:q.update(partition=1)]:
        bad=copy.deepcopy(persistence);mutation(bad)
        run("assay-query","--workspace",workspace,"--request",save("invalid.json",bad),fail=True)
        checks+=1
    print(json.dumps(dict(status="PASS",independent_checks=checks,scope="Original square coordinates, missing versus zero, rational stages and integer torsion")))


if __name__ == "__main__":
    main(*[Path(value).resolve() for value in sys.argv[1:]])
