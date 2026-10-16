"""Declared inter assay sections; independent rational projection and portable lineage."""
import copy
import hashlib
import importlib.util
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


def complete_families(result):
    for group in result['prepared_groups']:
        action = matrix(group['action'])
        kernel = matrix(group['kernel_frame'])
        dependencies = matrix(group['row_dependencies'])
        frame = matrix(group['residual_frame'])
        n, m = group['action']['columns'], group['action']['rows']
        rank = rref(action)[1]
        check(group['kernel_frame']['rows'] == n)
        check(group['kernel_frame']['columns'] == n-rank == group['compatible_dimension'])
        check(rref(kernel)[1] == n-rank)
        check(group['row_dependencies']['rows'] == m)
        check(group['row_dependencies']['columns'] == m-rank)
        check(rref(dependencies)[1] == m-rank)
        check(all(sum(action[i][j]*kernel[j][k] for j in range(n)) == 0
                  for i in range(m) for k in range(n-rank)))
        check(all(sum(action[i][j]*dependencies[i][k] for i in range(m)) == 0
                  for j in range(n) for k in range(m-rank)))
        check(rref([kernel[j]+frame[j] for j in range(n)])[1] == n)


def paired_capture_contracts():
    tool = example.parents[2]/'tools/prepare_paired_assay_sections.py'
    module_spec = importlib.util.spec_from_file_location('paired_sections', tool)
    preparation = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(preparation)
    folder = root/'paired-capture'
    folder.mkdir()

    def artifact(name, body):
        (folder/name).write_text(body)
        return dict(path=name, sha256='sha256:'+hashlib.sha256(body.encode()).hexdigest())

    barcodes = artifact('barcodes.tsv', 'cell-a\ncell-b\ncell-c\n')
    features = artifact('features.tsv', 'g1\tGene1\tGene Expression\ng2\tGene2\tGene Expression\na1\tAntibody1\tAntibody Capture\na2\tAntibody2\tAntibody Capture\n')
    original = artifact('matrix.mtx', '%%MatrixMarket matrix coordinate integer general\n4 3 8\n1 1 2\n2 1 3\n3 1 1\n4 1 2\n1 2 7\n3 2 7\n2 3 4\n4 3 4\n')
    observations = artifact('observations.tsv', 'id\ncell-a\ncell-b\ncell-c\n')
    assays = []
    for name, ids, entries, namespace, channel in [
        ('rna', ['g1','g2'], [(1,1,2),(2,1,3),(1,2,7),(2,3,4)], '978101', 'E'),
        ('antibody', ['a1','a2'], [(1,1,1),(2,1,2),(1,2,7),(2,3,4)], '978102', 'observation')]:
        dictionary = artifact(name+'-features.tsv', 'id\n'+'\n'.join(ids)+'\n')
        matrix_source = artifact(name+'.mtx', '%%MatrixMarket matrix coordinate integer general\n2 3 4\n'+''.join(f'{i} {j} {v}\n' for i,j,v in entries))
        assays.append(dict(id=name, kind='single_cell', format='matrix_market', source=matrix_source,
            features=dictionary, observations=observations, source_namespace=namespace, carrier='features',
            quantity='umi_count', units=name+'_assigned_umi', channel=channel, method='Authored captured counts',
            arithmetic='exact_decimal', implicit_entries='measured_zero', missing_tokens=[], provenance=[barcodes,features,original]))
    source = dict(schema='gadget.assays.import.v1', id='paired-source-controls', design=dict(
        question='Check original matched source columns', population='Three authored columns', selection='All',
        limitations=['Constructed source control']), assays=assays, resources=[], links=[])
    (folder/'paired-import.json').write_text(json.dumps(source)+'\n')
    contract = preparation.source_contract(folder)
    check([p['section'] for p in contract['correspondence']['pairs']] == ['cell-a','cell-b','cell-c'])
    destination = root/'paired-workflow'
    p = subprocess.run([sys.executable, str(tool), '--binary', str(binary), '--source', str(folder),
        '--out', str(destination), '--source-namespace', '978103', '--workers', '2'], capture_output=True, text=True)
    check(p.returncode == 0)
    paired = json.loads((destination/'analysis/sections/result.json').read_text())['result']
    complete_families(paired)
    check([len(paired['sections'][0]['fields'][i]['input']) for i in range(2)] == [1,2])
    for section, expected in zip(paired['sections'], ([5,1,2],[7,7,0],[4,0,4])):
        check([Q(v) for f in section['fields'] for v in f['input']] == expected)
        defect = Q(expected[0]-expected[1]-expected[2])
        residual = [defect/3, -defect/3, -defect/3]
        check([Q(v) for f in section['fields'] for v in f['residual']] == residual)
        check(Q(section['residual_quadrance']) == defect*defect/3)
    for filename, changed in [('barcodes.tsv','cell-a\ncell-a\ncell-c\n'),
                              ('observations.tsv','id\ncell-b\ncell-a\ncell-c\n')]:
        original_bytes = (folder/filename).read_bytes()
        (folder/filename).write_text(changed)
        bad = copy.deepcopy(source)
        for assay in bad['assays']:
            for record in [assay['observations'], *assay['provenance']]:
                if record['path'] == filename:
                    record['sha256'] = 'sha256:'+hashlib.sha256(changed.encode()).hexdigest()
        (folder/'paired-import.json').write_text(json.dumps(bad))
        try:
            preparation.source_contract(folder)
        except ValueError:
            check(True)
        else:
            check(False)
        (folder/filename).write_bytes(original_bytes)
    (folder/'paired-import.json').write_text(json.dumps(source))


run('import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
manifest = root/'workspace/workspace.json'
query = json.loads((root/'inputs/query.json').read_text())
query['complete_families'] = True
result = run('assay-query', '--workspace', manifest, '--request', save('complete-query.json', query),
             '--out', root/'section owner #1')['result']
complete_families(result)
published = root/'section owner #1/result.rfcxs'
published_hash = hashlib.sha256(published.read_bytes()).hexdigest()
run('assay-query', '--workspace', manifest, '--request', root/'complete-query.json',
    '--out', root/'section owner #1', fail=True)
check(hashlib.sha256(published.read_bytes()).hexdigest() == published_hash)
paired_capture_contracts()
legacy = run('assay-query', '--workspace', manifest, '--request', root/'inputs/query.json')['result']
explicit_legacy = copy.deepcopy(query)
explicit_legacy['complete_families'] = False
old_profile = run('assay-query', '--workspace', manifest,
                  '--request', save('legacy-query.json', explicit_legacy))['result']
check(payload(legacy) == payload(old_profile))
check(all('kernel_frame' not in g and 'row_dependencies' not in g for g in legacy['prepared_groups']))
for wrong in ['complete', 1, None, [], {}]:
    bad = copy.deepcopy(query)
    bad['complete_families'] = wrong
    run('assay-query', '--workspace', manifest,
        '--request', save('invalid-families.json', bad), fail=True)
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
# Joined serial and multi worker actions agree on the complete mathematical output.
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
complete_families(wide_interface)
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
    refused = root / ('refused-output-' + name)
    run('assay-query', '--workspace', manifest, '--request', root/(name+'.json'), '--out', refused, fail=True)
    check(not refused.exists())
# Default pipelines and bound local model recipe selection use the same contract.
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
# Previous local section ceilings: 258 direct sum coordinates and 33 sections.
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
renderer_path = example.parents[2] / 'tools/render_assay_sections.py'
spec = importlib.util.spec_from_file_location('section_reading', renderer_path)
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)
reading = root / 'reading with spaces'
operation = root/'native owner #1/assay sections.rfcxs'
operation.parent.mkdir()
shutil.copyfile(root/'section owner #1/result.rfcxs', operation)
receipt = renderer.render(binary, operation, reading)
check(receipt['sections'] == len(result['sections']))
check(json.loads((reading/'view.json').read_text())['result'] == result)
check(receipt['native_sha256'] == 'sha256:' + hashlib.sha256(operation.read_bytes()).hexdigest())
check('Missing coordinates' in (reading/'report.html').read_text())
check('native%20owner%20%231/assay%20sections.rfcxs' in (reading/'report.html').read_text())
try:
    renderer.render(binary, root/'complete-query.json', root/'invalid reading')
except ValueError:
    check(not (root/'invalid reading').exists())
else:
    check(False)
print(json.dumps(dict(status='PASS', independent_checks=checks, sections=5, source_coordinates=6,
                      constraint_coordinates=6, support_patterns=3, portable_replay=True)))
