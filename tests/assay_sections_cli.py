"""Declared inter-assay sections; independent rational projection and portable lineage."""
import copy
import hashlib
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)
example = Path(__file__).resolve().parent.parent / 'examples/bioinformatics/assay-sections'
shutil.copytree(example, root / 'inputs')
checks = 0


def check(condition):
    global checks
    assert condition
    checks += 1


def run(*args, fail=False):
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    assert bool(p.returncode) == fail, (args, p.stdout, p.stderr)
    check(True)
    if not fail:
        return json.loads(p.stdout)


def save(name, value):
    p = root / name
    p.write_text(json.dumps(value, indent=2) + '\n')
    return p


def rref(rows):
    rows = [list(map(Q, row)) for row in rows]
    k = 0
    for j in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        factor = rows[k][j]
        rows[k] = [x / factor for x in rows[k]]
        for i in range(len(rows)):
            if i != k:
                factor = rows[i][j]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[k])]
        k += 1
        if k == len(rows):
            break
    return rows, k


def solve(matrix, rhs):
    n = len(rhs)
    rows, rank = rref([list(row) + [b] for row, b in zip(matrix, rhs)])
    check(rank == n)
    return [row[-1] for row in rows]


def matrix(record):
    rows = [[Q(0)] * record['columns'] for _ in range(record['rows'])]
    for e in record['entries']:
        rows[e['row']][e['column']] = Q(e['value'])
    return rows


def payload(result):
    result = copy.deepcopy(result)
    result.pop('execution', None)
    return result


run('import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
manifest = root/'workspace/workspace.json'
query = json.loads((root/'inputs/query.json').read_text())
result = run('assay-query', '--workspace', manifest, '--request', root/'inputs/query.json')['result']
D = [[1, 0, -1, 0, 0, 0], [0, 0, 0, 1, Q(-1, 2), 0],
     [0, 1, 0, 1, 0, -1], [0]*6, [0]*6, [0]*6]
dependencies = [{0, 2}, {3, 4}, {1, 3, 5}, {0}, set(), set()]
weights = list(map(Q, range(1, 7)))
check(matrix(result['coboundary']) == D)
check(result['interface_grade'] == 1 and result['interface_chain_verified'])
check([c['observable'] for c in result['constraint_coordinates']] == ['signal'] * 6)
check(result['execution']['prepared_solves'] == 3)
check(result['execution']['completed_tiles'] == 5)
check(len(result['interfaces'][2]['participants']) == 3)
check(len(result['interfaces'][3]['participants']) == 2)
check(not result['interfaces'][4]['participants'][0]['exists'])
check(result['interfaces'][5]['participants'][0]['magnitude'] == '0')
for section in result['sections']:
    x = [None if v is None else Q(v) for f in section['fields'] for v in f['input']]
    observed = [i for i, v in enumerate(x) if v is not None]
    assessed = [i for i, dep in enumerate(dependencies) if all(x[j] is not None for j in dep)]
    rows = [[Q(D[i][j]) for j in observed] for i in assessed]
    basis, rank = rref(rows)
    basis = basis[:rank]
    gram = [[sum(a*b/weights[j] for a, b, j in zip(u, v, observed)) for v in basis] for u in basis]
    defect = [sum(a*x[j] for a, j in zip(u, observed)) for u in basis]
    coeff = solve(gram, defect) if rank else []
    residual = [sum(u[c]*a for u, a in zip(basis, coeff))/weights[j] for c, j in enumerate(observed)]
    compatible = [x[j]-r for j, r in zip(observed, residual)]
    actual_residual = [None if v is None else Q(v) for f in section['fields'] for v in f['residual']]
    actual_compatible = [None if v is None else Q(v) for f in section['fields'] for v in f['compatible']]
    for j, r, c in zip(observed, residual, compatible):
        check(actual_residual[j] == r and actual_compatible[j] == c)
    for j in set(range(6))-set(observed):
        check(actual_residual[j] is None and actual_compatible[j] is None)
    check(Q(section['residual_quadrance']) == sum(weights[j]*r*r for j, r in zip(observed, residual)))
    check(section['orthogonality'] == '0' and section['reconstruction_verified'])
    for i, c in enumerate(section['constraints']):
        check(c['assessed'] == (i in assessed))
        check(c['missing_coordinates'] == sorted(j for j in dependencies[i] if x[j] is None))
        check(c['input_defect'] == (str(sum(Q(D[i][j])*x[j] for j in observed)) if i in assessed else None))
    group = result['prepared_groups'][section['prepared_group']]
    check(group['coordinates'] == observed and group['assessed_rows'] == assessed)
    check(group['constraint_rank'] == rank and group['compatible_dimension'] == len(observed)-rank)
    F = matrix(group['residual_frame'])
    G = matrix(group['gram'])
    check(G == [[sum(weights[j]*F[k][a]*F[k][b] for k, j in enumerate(observed))
                  for b in range(rank)] for a in range(rank)])
# Equal labels across sources require the supplied maps, not an inferred identity.
modified = copy.deepcopy(query)
modified['interfaces'][0]['participants'][1]['map']['entries'][0]['coefficient'] = '2'
changed = run('assay-query', '--workspace', manifest, '--request', save('changed-map.json', modified))['result']
check(payload(changed) != payload(result))
# Joined serial and multi-worker actions agree on the complete mathematical output.
serial = copy.deepcopy(query); serial['workers'] = 1
check(payload(run('assay-query', '--workspace', manifest, '--request', save('serial.json', serial))['result']) == payload(result))
# Unequal interface stalk dimensions and multivariate restriction maps.
rectangular = copy.deepcopy(query)
rectangular['interfaces'][0]['axes'].append('total')
for participant in rectangular['interfaces'][0]['participants']:
    participant['map']['entries'].extend([
        dict(observable_id='total', source_id='f1', coefficient='2'),
        dict(observable_id='total', source_id='f2', coefficient='1')])
wide_interface = run('assay-query', '--workspace', manifest,
                     '--request', save('rectangular.json', rectangular))['result']
expected_D = [D[0], [2, 1, -2, -1, 0, 0]] + D[1:]
check(matrix(wide_interface['coboundary']) == expected_D)
first = wide_interface['sections'][0]
x = [Q(v) for field in first['fields'] for v in field['input']]
basis, rank = rref(expected_D); basis = basis[:rank]
gram = [[sum(a*b/w for a, b, w in zip(u, v, weights)) for v in basis] for u in basis]
coeff = solve(gram, [sum(a*b for a,b in zip(u,x)) for u in basis])
expected_r = [sum(u[j]*a for u,a in zip(basis,coeff))/weights[j] for j in range(6)]
check([Q(v) for f in first['fields'] for v in f['residual']] == expected_r)
check([Q(v) for f in first['fields'] for v in f['compatible']] == [a-b for a,b in zip(x,expected_r)])
check(Q(first['residual_quadrance']) == sum(w*r*r for w,r in zip(weights,expected_r)))

# All declaration admission, including inactive occurrence maps, stays explicit.
for name, edit in [
    ('metric-zero', lambda q: q['views'][0]['metric']['diagonal'].__setitem__(0, '0')),
    ('metric-units', lambda q: q['views'][0]['metric'].__setitem__('quadrance_units', 'other')),
    ('domain-channel', lambda q: q['interfaces'][0]['participants'][0]['map'].__setitem__('domain_channel', 'G')),
    ('inactive-domain', lambda q: q['interfaces'][4]['participants'][0]['map'].__setitem__('domain_units', 'other')),
    ('foreign-id', lambda q: q['interfaces'][0]['participants'][0]['map']['entries'][0].__setitem__('source_id', 'foreign')),
    ('duplicate-occurrence', lambda q: q['interfaces'][1]['participants'][0].__setitem__('key', '1')),
    ('missing-field', lambda q: q['views'][0]['fields'][0].__setitem__('field', 'all.foreign')),
    ('negative-weight', lambda q: q['interfaces'][0]['participants'][0].__setitem__('magnitude', '-1')),
]:
    bad = copy.deepcopy(query); edit(bad)
    run('assay-query', '--workspace', manifest, '--request', save(name+'.json', bad), fail=True)
# Default pipelines and bound local-model recipe selection use the same contract.
plan = run('pipeline-prepare', '--workspace', manifest, '--request', root/'inputs/pipeline.json')
run('pipeline', '--workspace', manifest, '--request', save('plan.json', plan), '--out', root/'pipeline')
derived = root/'pipeline/sections/workspace/workspace.rfcxs'
check(derived.exists())
run('replay-derivation', '--workspace', derived, '--assay', 'rna-compatible', '--out', root/'replay')
catalog = run('recipe-catalog', '--workspace', manifest, '--intent', root/'inputs/intent.json')
proposal = dict(schema='gadget.pipeline.recipe.proposal.v1', action='recipes',
                source_digest=catalog['source_digest'], intent_sha256=catalog['intent_sha256'], goals=['support'])
admitted = run('recipe-proposal', '--workspace', manifest, '--intent', root/'inputs/intent.json',
               '--request', save('proposal.json', proposal), '--prompt', 'Assess support after compatibility')
check(admitted['status'] == 'admitted')
# Previous local-section ceilings: 258 direct-sum coordinates and 33 sections.
study = json.loads((root/'workspace/assays/0/partitions/0/study.json').read_text())
local = dict(operation='local_sections', context='all', construction='cellular_sheaf',
             stalks=[dict(grade=0,key=c['key'],axes=[f'a{i}' for i in range(129)])
                     for c in study['cells'] if c['grade']==0], maps=[],
             sections=[dict(id=f's{i}',grade=0,values=[str(i)]*258,quantity='declared_test',units='test')
                       for i in range(33)], metric='identity')
local_query = dict(schema='gadget.assays.query.v1',operation='native_query',assay='rna',query=local)
wide = run('assay-query', '--workspace', manifest, '--request', save('wide-local.json',local_query))
local_result = wide['result']['answers'][0]['certificate']['result']
check(local_result['reconstruction_exact'] and len(local_result['sections']) == 33)
check(local_result['cohomology'][0]['dimension'] == 258)
for i, section in enumerate(local_result['sections']):
    check(section['glued'] == [str(i)]*258 and section['residual'] == ['0']*258)
# The archive alone suffices after the source workspace and original source files move.
shutil.copytree(derived.parent, root/'portable')
(root/'workspace').rename(root/'workspace-unavailable')
(root/'inputs').rename(root/'inputs-unavailable')
replayed = run('replay-derivation', '--workspace', root/'portable/workspace.rfcxs',
               '--assay', 'rna-compatible', '--out', root/'detached-replay')
check(replayed['status'] == 'PASS')
print(json.dumps(dict(status='PASS', independent_checks=checks, sections=5, source_coordinates=6,
                      constraint_coordinates=6, support_patterns=3, portable_replay=True)))
