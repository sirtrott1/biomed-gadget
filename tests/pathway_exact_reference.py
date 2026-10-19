"""Enumerate original finite subsets independently of the native law construction."""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main(binary):
    def execute(*args):
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout)
    template = execute("--request")
    randomizer = random.Random(763)
    checks = 0
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "source.json"
        for N in range(1, 9):
            population = set(range(N))
            sets = [set(), population, set(randomizer.sample(range(N), randomizer.randrange(N + 1)))]
            ids = ["source:" + str(i) for i in range(N)]
            gmt = "".join("set" + str(j) + "\tdeclared\t" + "\t".join("member" + str(i) for i in sorted(members)) + "\n"
                          for j, members in enumerate(sets))
            mapping = "feature_id\tmember\tevidence\n" + "".join(identity + "\tmember" + str(i) + "\texplicit\n" for i, identity in enumerate(ids))
            for draws in range(N + 1):
                chosen = set(randomizer.sample(range(N), draws))
                subsets = tuple(map(set, itertools.combinations(range(N), draws)))
                tails = [Q(sum(len(sample & members) >= len(chosen & members) for sample in subsets), len(subsets)) for members in sets]
                for correction in ("none", "bonferroni", "benjamini_hochberg", "benjamini_yekutieli"):
                    query = json.loads(json.dumps(template))
                    query.update(background_indices=list(range(N)), selected_indices=sorted(chosen),
                                 set_ids=["set" + str(j) for j in range(len(sets))], correction=correction)
                    query["size_filter"].update(minimum=0, maximum=N)
                    metadata = dict(id="authored", release="one", identifier_namespace="declared:member",
                        sha256="sha256:" + hashlib.sha256(gmt.encode()).hexdigest(), mapping_id="authored", mapping_release="one",
                        mapping_origin="explicit_declared", feature_namespace=query["feature_namespace"],
                        mapping_sha256="sha256:" + hashlib.sha256(mapping.encode()).hexdigest())
                    path.write_text(json.dumps(dict(features=ids, gmt=gmt, mapping=mapping, source=metadata, request=query)))
                    result = execute("--evaluate", path)
                    assert result["tested_family_size"] == len(sets)
                    harmonic = sum((Q(1, j) for j in range(1, len(sets) + 1)), Q(0))
                    for test, members, tail in zip(result["tests"], sets, tails, strict=True):
                        assert Q(test["p_value"]) == tail
                        distribution = Counter(len(sample & members) for sample in subsets)
                        counts = list(map(int, test["intersection_counts"]))
                        masses = list(map(Q, test["reference_masses"]))
                        assert dict(zip(counts, masses, strict=True)) == {k: Q(count, len(subsets)) for k, count in distribution.items()}
                        if correction == "none": expected = tail
                        elif correction == "bonferroni": expected = min(1, len(sets) * tail)
                        else:
                            expected = min(Q(1), *(len(sets) * value / sum(p <= value for p in tails) *
                                (harmonic if correction == "benjamini_yekutieli" else 1) for value in tails if value >= tail))
                        assert Q(test["adjusted_p_value"]) == expected
                        checks += 3
    print(json.dumps(dict(status="PASS", exact_comparisons=checks, source="complete original subset enumeration")))


if __name__ == "__main__": main(Path(sys.argv[1]))
