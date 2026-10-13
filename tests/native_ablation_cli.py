"""Native ablation sources, complete field families and original owner replay."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    path = Path(__file__).resolve().parents[1] / 'examples/bioinformatics/structural-ablations/replay.py'
    spec = importlib.util.spec_from_file_location('native_ablation_demo', path)
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    case = root / 'source'
    case.mkdir()
    query = demo.authored(case, workers=2)
    demo.execute(binary, case, query)
    ranks = copy.deepcopy(query)
    ranks['homology_output'] = 'ranks'
    lattice = copy.deepcopy(query)
    lattice['lattice_snapshots'] = True
    compatible = demo.compatible_cover(query)
    empty = copy.deepcopy(query)
    empty.pop('sheaf')
    empty.pop('comparisons')
    original_study = next((case / 'workspace').rglob('study.json'))
    declaration = json.loads(original_study.read_text())
    empty['roles'] = [demo.cell_role('all', [cell['id'] for cell in declaration['cells']],
                                   'All explicitly identified source cells')]
    empty['ablations'] = [dict(id='zero-carrier', remove_roles=['all'])]
    demo.bind_evidence(empty, case / 'workspace')
    queries = {'complete': query, 'ranks': ranks, 'lattice': lattice,
               'compatible': compatible, 'empty': empty}
    weighted = root / 'weighted'
    shutil.copytree(original_study.parent, weighted)
    declaration['contexts'][0]['metrics'] = [dict(kind='diagonal', values=values)
        for values in (['1', '2', '3', '4'], ['2', '3', '5', '7'], ['11'])]
    (weighted / 'study.json').write_text(json.dumps(declaration))
    path = root / 'queries.json'
    path.write_text(json.dumps(queries))
    result = subprocess.run([str(consumer), str(case / 'workspace/workspace.rfcxs'),
                             str(path), str(root / 'results'), 'ablation'],
                            capture_output=True, text=True, timeout=600)
    if result.returncode:
        raise AssertionError(result.stderr)
    print(result.stdout.strip())
    result = subprocess.run([str(consumer), str(weighted / 'study.json'), str(path),
                             str(root / 'weighted-results'), 'ablation'],
                            capture_output=True, text=True, timeout=600)
    if result.returncode:
        raise AssertionError(result.stderr)
    print(result.stdout.strip())


if __name__ == '__main__':
    main(*[Path(value).resolve() for value in sys.argv[1:]])
