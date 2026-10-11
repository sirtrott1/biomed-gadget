"""Retained model fields, complete coefficient families and relocated replay."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path


def load(name):
    path = Path(__file__).resolve().parents[1] / 'examples/bioinformatics' / name / 'replay.py'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    def command(executable, *args):
        result = subprocess.run([str(executable), *map(str, args)], capture_output=True, text=True, timeout=240)
        if result.returncode:
            raise RuntimeError(result.stderr + result.stdout)
        return result.stdout
    linear = load('linear-models')
    queries = linear.stage(root/'models/inputs')
    queries['partial'] = copy.deepcopy(queries['continuous'])
    queries['partial']['terms'].append(dict(id='partial', factors=['all.partial_x']))
    linear.save(root/'models/queries.json', queries)
    command(binary, 'import-assays', '--spec', root/'models/inputs/import.json', '--out', root/'models/workspace')
    print(command(consumer, root/'models/workspace/workspace.rfcxs', root/'models/queries.json', root/'models/results', 'models'))
    finite = load('finite-calibration')
    observations = [dict(id=f'specimen-{i}', block=f'pair-{i//2}', treated=str(i%2)) for i in range(8)]
    pairs = [[row['id'] for row in observations[i:i+2]] for i in range(0, 8, 2)]
    fields = dict(A=['10','12','20','24','30','36','40','48'], B=['7','10','15','12','20','21','30','29'],
                  flat=['10','10','20','20','30','30','40','40'])
    query = finite.stage(root/'finite/inputs', observations, fields, pairs)
    command(binary, 'import-assays', '--spec', root/'finite/inputs/import.json', '--out', root/'finite/workspace')
    finite.bind_evidence(query, root/'finite/workspace')
    queries = {'undefined': query}
    complete = copy.deepcopy(query)
    complete['model']['responses'] = complete['model']['responses'][:2]
    complete['hypotheses'] = complete['hypotheses'][:2]
    for state in complete['law']['states']:
        state['responses'] = state['responses'][:2]
    queries['complete'] = complete
    finite.save(root/'finite/queries.json', queries)
    print(command(consumer, root/'finite/workspace/workspace.rfcxs', root/'finite/queries.json', root/'finite/results', 'finite'))
    source = root/'models/results/near_collinear.rfcxs'
    command(binary, 'replay-scientific', '--record', source, '--out-record', root/'cli-replay.rfcxs')
    view = json.loads(command(binary, 'record-view', '--input', root/'cli-replay.rfcxs'))
    assert view['result']['result']['strata'][0]['models'][1]['rank'] == 3
    refused = subprocess.run([str(binary), 'replay-scientific', '--record', str(root/'models/results/changed-field.rfcxs'),
                              '--out-record', str(root/'refused.rfcxs')], capture_output=True, text=True, timeout=240)
    assert refused.returncode and not (root/'refused.rfcxs').exists()
    print(json.dumps(dict(status='PASS', native_cli_replay=True)))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve())
