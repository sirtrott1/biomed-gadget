"""Native model and ensemble ownership, complete fields and portable replay."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    examples = Path(__file__).resolve().parents[1] / 'examples/bioinformatics'
    observed = load('native_observed', examples / 'observed-responses/replay.py')
    constrained = load('native_constrained', examples / 'constrained-responses/replay.py')
    ensemble = load('native_ensemble', examples / 'relational-ensembles/replay.py')
    receipts = []
    for name, mode in [('observed', 'response'), ('constrained', 'response'), ('ensemble', 'ensemble')]:
        case = root / name
        if name == 'ensemble':
            query = ensemble.authored(case / 'inputs')
            ensemble.execute(binary, case, query)
        else:
            query = observed.authored(case, 2)
            if name == 'constrained':
                query = constrained.constrained(query)
                constrained.execute(binary, case, query)
            else:
                observed.execute(binary, case, query)
        queries = {'complete': query}
        if name == 'observed':
            empty = copy.deepcopy(query)
            empty['model'] = dict(intercept=False, terms=[])
            empty.pop('dose_comparisons')
            queries['empty'] = empty
        path = root / (name + '-queries.json')
        path.write_text(json.dumps(queries))
        result = subprocess.run([str(consumer), str(case / 'workspace/workspace.rfcxs'),
                                 str(path), str(root / (name + '-results')), mode],
                                capture_output=True, text=True, timeout=300)
        if result.returncode:
            raise AssertionError(name + ': ' + result.stderr)
        receipts.append(dict(case=name, **json.loads(result.stdout)))
    for size in (131, 3):
        case = root / ('green-' + str(size))
        ids = ['coordinate-' + str(i) for i in range(size)]
        signal = list(range(size))
        signal[0] -= 1
        signal[-1] += 1
        relations = [dict(id='edge-' + str(i), grade=1, meaning='Declared line boundary',
                          participants=[ensemble.participant(ids[i], -1), ensemble.participant(ids[i + 1], 1)])
                     for i in range(size - 1)]
        ensemble.stage(case / 'inputs', ids, dict(signal=signal), relations, [])
        ensemble.run(binary, 'import-assays', '--spec', case / 'inputs/import.json', '--out', case / 'workspace')
        result = subprocess.run([str(consumer), str(case / 'workspace/workspace.rfcxs'),
                                 str(root / 'observed-queries.json'), str(root / ('green-results-' + str(size))), 'green'],
                                capture_output=True, text=True, timeout=300)
        if result.returncode:
            raise AssertionError('green: ' + result.stderr)
        receipts.append(dict(case='green-' + str(size), **json.loads(result.stdout)))
    print(json.dumps(dict(status='PASS', cases=receipts)))


if __name__ == '__main__':
    main(*[Path(v).resolve() for v in sys.argv[1:]])
