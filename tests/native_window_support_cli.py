"""Native windows and ordered eligibility on complete original support."""
from pathlib import Path
from fractions import Fraction
import json
import subprocess
import sys
import shutil
import hashlib

binary, window_consumer, support_consumer, root = [Path(v).resolve() for v in sys.argv[1:]]
shutil.rmtree(root, ignore_errors=True)
inputs = root / 'inputs'
inputs.mkdir(parents=True)
values = {'a': [2] + [0]*299, 'c': [0]*300, 'd': [0]*200 + [None] + [0]*99,
          'missing': [None]*300}

def raw(name, text):
    (inputs/name).write_text(text)
    return dict(path=name, sha256='sha256:'+hashlib.sha256(text.encode()).hexdigest())

assay = dict(id='x', kind='rna_seq', format='dense_tsv', source_namespace='976701',
    carrier='features', quantity='read_count', units='counts', channel='observation',
    method='authored counts', arithmetic='exact_rational', implicit_entries='measured_zero',
    missing_tokens=['NA'], partition_size=300)
assay['source'] = raw('matrix.tsv', 'id\ta\tc\td\tmissing\n'+''.join(
    '\t'.join(['g'+str(i)]+['NA' if values[f][i] is None else str(values[f][i]) for f in values])+'\n'
    for i in range(300)))
assay['features'] = raw('features.tsv', 'id\n'+''.join('g'+str(i)+'\n' for i in range(300)))
assay['observations'] = raw('observations.tsv', 'id\na\nc\nd\nmissing\n')
spec = dict(schema='gadget.assays.import.v1', id='native-eligibility',
    design=dict(question='Retain exact declared support and windows', population='authored observations',
                selection='all', limitations=['authored fixture']), assays=[assay], resources=[], links=[])
(inputs/'import.json').write_text(json.dumps(spec))
subprocess.run([str(binary), 'import-assays', '--spec', str(inputs/'import.json'),
                '--out', str(root/'workspace')], check=True, capture_output=True, text=True)
study = next((root/'workspace').rglob('study.json'))

def query(name, declaration, good=True, out=None):
    path = root/(name+'.json')
    path.write_text(json.dumps(declaration))
    command = [str(binary), 'query', '--study', str(study), '--query', str(path)]
    if out is not None:
        command += ['--out', str(out)]
    p = subprocess.run(command, capture_output=True, text=True)
    assert (p.returncode == 0) == good, p.stderr
    return json.loads(p.stdout) if good else None

source = root/'source.rfcxs'
query('source', dict(operation='field_summary', field='all.a', support='reported_only'), out=source)
windows = {}
checks = 0
for field in values:
    for lo, hi in [(0, 0), (0, 300), (199, 202), (200, 201)]:
        for preparation in ['reuse', 'one_shot']:
            name = field+'-'+str(lo)+'-'+str(hi)+'-'+preparation
            declaration = dict(operation='tensor_window', field='all.'+field, lower=[lo], upper=[hi],
                               available_sequence='0', preparation=preparation)
            windows[name] = declaration
            result = query(name, declaration)['result']
            retained = values[field][lo:hi]
            assert result['reported'] == sum(v is not None for v in retained)
            assert result['missing'] == sum(v is None for v in retained)
            assert [r['source_coordinate'] for r in result['observations']] == list(range(lo, hi))
            actual = [None if r['value'] is None else Fraction(r['value']['rational'])
                      for r in result['observations']]
            assert actual == retained
            total = sum(v for v in retained if v is not None)
            assert result['sum'] is None if result['reported'] == 0 else Fraction(result['sum']['rational']) == total
            assert result['preparation_retained'] == (preparation == 'reuse')
            checks += len(retained)+6
windows['reused-full'] = dict(windows['a-0-300-reuse'])
(root/'windows.json').write_text(json.dumps(windows))
subprocess.run([str(window_consumer), str(source), str(root/'windows.json'), str(root/'native-windows')],
                check=True, capture_output=True, text=True)
support = {}
for fields in [['d'], ['a', 'd'], ['d', 'a'], ['a', 'c', 'd'], ['d', 'missing'], ['missing']]:
    for mode in ['all_reported', 'any_reported', 'none_reported']:
        name = '-'.join(fields)+'-'+mode
        declaration = dict(operation='support_mask', fields=['all.'+f for f in fields], mode=mode,
                           backend='cpu_tiles', workers=1, tile_words=1)
        support[name] = declaration
        result = query(name, declaration)
        expected = []
        for i in range(300):
            reports = [values[f][i] is not None for f in fields]
            selected = all(reports) if mode == 'all_reported' else any(reports) if mode == 'any_reported' else not any(reports)
            if selected:
                expected.append(i)
        assert result['result']['selected'] == len(expected)
        assert [r['coordinate'] for r in result['result']['selection']] == expected
        assert result['execution']['exact_cpu_reference_passed']
        checks += 302
(root/'support.json').write_text(json.dumps(support))
subprocess.run([str(support_consumer), str(source), str(root/'support.json'), str(root/'native-support')],
                check=True, capture_output=True, text=True)
query('bad-bounds', dict(operation='tensor_window', field='all.a', lower=[0], upper=[301], available_sequence='0'), False)
query('reversed-bounds', dict(operation='tensor_window', field='all.a', lower=[3], upper=[2], available_sequence='0'), False)
query('duplicate-fields', dict(operation='support_mask', fields=['all.a', 'all.a'], mode='all_reported', backend='cpu_tiles'), False)
query('unknown-mode', dict(operation='support_mask', fields=['all.a'], mode='unknown', backend='cpu_tiles'), False)
print(json.dumps(dict(status='PASS', independent_checks=checks, window_queries=len(windows), support_queries=len(support))))
