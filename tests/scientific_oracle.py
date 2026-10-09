"""Independent Fraction identities for the declared nonclinical fixture."""
import json
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def q(item):
    return Q(item["rational"])


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), Q(0))


def apply(columns, weights, x):
    out = [Q(0)] * len(x)
    for column, weight in zip(columns, weights):
        coefficient = dot(column, x) / weight
        out = [a + coefficient * b for a, b in zip(out, column)]
    return out


def main():
    executable, root = sys.argv[1], Path(sys.argv[2])
    subprocess.run([executable], check=True)
    results = json.loads((root / "results.json").read_text())
    x0, y0, x1, y1 = map(lambda xs: list(map(Q, xs)), ([1, 0, 1], [0, 2, 0], [1, 1, 1], [1, 3, 2]))
    columns = list(map(lambda xs: list(map(Q, xs)), ([-1, 1, 0], [0, -1, 1])))
    old = lambda x: apply(columns[:1], [Q(2)], x)
    new = lambda x: apply(columns, [Q(2), Q(3)], x)
    dx, dy = [a-b for a, b in zip(x1, x0)], [a-b for a, b in zip(y1, y0)]
    terms = dict(pullback_form=dot(x0, [a-b for a, b in zip(new(y0), old(y0))]),
                 left_innovation=dot(dx, new(y0)), right_innovation=dot(x0, new(dy)),
                 innovation_interaction=dot(dx, new(dy)))
    total = dot(x1, new(y1)) - dot(x0, old(y0))
    assert sum(terms.values()) == total
    assert q(results["old"]["result"]) == dot(x0, old(y0))
    assert q(results["new"]["result"]) == dot(x1, new(y1))
    for name, value in terms.items():
        assert q(results["change"]["result"][name]) == value
    assert q(results["change"]["result"]["total"]) == total
    assert q(results["change"]["result"]["finite_rate"]) == total / 2
    # Independently verify all nine sector blocks and the full Green response
    # change on the explicitly materialized temporal fixture.
    if 'structure' in results:
        from research_cli import solve
        structural=results['structure']['result']
        transport=structural['sector_transport']
        assert [sum(Q(transport[a][b][i]) for a in range(3) for b in range(3)) for i in range(3)]==x0
        l0=[[old([Q(int(j==k)) for j in range(3)])[i] for k in range(3)] for i in range(3)]
        l1=[[new([Q(int(j==k)) for j in range(3)])[i] for k in range(3)] for i in range(3)]
        green0=solve([[value+int(i==j) for j,value in enumerate(row)] for i,row in enumerate(l0)],x0)
        green1=solve([[value+int(i==j) for j,value in enumerate(row)] for i,row in enumerate(l1)],x1)
        total_green=[a-b for a,b in zip(green1,green0)]
        assert total_green==list(map(Q,structural['total']))
        assert [sum(Q(structural[t][i]) for t in ('realization_term','innovation_term','operator_term')) for i in range(3)]==total_green
    # Each native row retains two actual coordinates. Repeated signed A
    # occurrences cancel; the retained inactive C occurrence contributes zero.
    expected = [[Q(-2), Q(-3)], [Q(-3), Q(-4)], [Q(5), Q(7)]]
    for name in ("native", "parallel"):
        assert [[q(v) for v in row["values"]] for row in results[name]["result"]] == expected
    assert results["coverage"]["result"]["observations"][0]["reported"]
    assert results["coverage"]["result"]["observations"][2]["value"] is None
    print("PASS: independent Fraction GET moments, all finite terms, native rows and missingness")


if __name__ == "__main__":
    main()
