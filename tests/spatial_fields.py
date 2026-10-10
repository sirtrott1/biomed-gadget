#!/usr/bin/env python3
"""Check spatial carriers, measured fields, native charts and local restrictions."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

gadget, verifier, sequence, section_verifier, root = map(Path, sys.argv[1:])
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)


def run(*args, success=True):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    if (result.returncode == 0) != success:
        raise AssertionError(f'{args}: {result.stdout}\n{result.stderr}')
    return result


def native(name, value):
    source, target = root / (name + '.json'), root / (name + '.rfcxs')
    source.write_text(json.dumps(value) + '\n')
    run(gadget, 'native-record', '--input', source, '--out', target)
    source.unlink()
    return target


def artifact(path):
    return dict(path=path.name, sha256='sha256:' + hashlib.sha256(path.read_bytes()).hexdigest())


design = dict(question='Qualify original spatial measurements and local sections',
              population='authored fixture', selection='every original row', limitations=['synthetic data'])


def prepare(name, positions, diagonal, radius, implicit='measured_zero', relations=None):
    folder = root / name
    folder.mkdir()
    columns = ['x' + str(i) for i in range(len(diagonal))]
    observations = folder / 'observations.tsv'
    observations.write_text('id\t' + '\t'.join(columns) + '\tcoordinate_frame\tcoordinate_units\n' +
                            ''.join(f'p{i}\t' + '\t'.join(map(str, p)) + '\tsection\tum\n'
                                    for i, p in enumerate(positions)))
    features = folder / 'features.tsv'
    features.write_text('id\tsymbol\ng1\tfirst\ng2\tsecond\n')
    matrix = folder / 'matrix.mtx'
    entries = [(1, i + 1, i * 3 + 2) for i in range(len(positions))]
    entries += [(2, 1, 0), (2, len(positions), 7)]
    matrix.write_text('%%MatrixMarket matrix coordinate integer general\n' +
                      f'2 {len(positions)} {len(entries)}\n' +
                      ''.join(f'{g} {p} {v}\n' for g, p, v in entries))
    assay = dict(id='measurements', kind='spatial_expression', format='matrix_market', source=artifact(matrix),
                 features=artifact(features), observations=artifact(observations), source_namespace='750001',
                 carrier='observations', quantity='umi_count', units='gene_assigned_umi', channel='E',
                 method='Original measured counts', arithmetic='exact_rational', implicit_entries=implicit,
                 missing_tokens=[])
    if relations:
        path = folder / 'relations.json'
        path.write_text(json.dumps(relations))
        assay['relations'] = artifact(path)
    spec = native(name + '-import', dict(schema='gadget.assays.import.v1', id=name, design=design,
                                       assays=[assay], resources=[], links=[]))
    # Resolve paths against the same directory as the authored specification.
    shutil.move(spec, folder / 'import.rfcxs')
    run(gadget, 'import-assays', '--spec', folder / 'import.rfcxs', '--out', folder / 'import')
    workspace = folder / 'import/workspace.rfcxs'
    value = dict(schema='gadget.assays.query.v1', operation='spatial_relations', assay='measurements',
                 coordinates=columns, frame='section', units='um', metric_diagonal=list(map(str, diagonal)),
                 radius_squared=str(radius), relation_type='physical.neighborhood', evidence='Declared coordinate protocol',
                 execution_certificate=True, charts='grade_one_stars')
    request = native(name + '-spatial', value)
    direct = folder / 'direct.rfcxs'
    run(gadget, 'assay-query', '--workspace', workspace, '--request', request, '--out-record', direct)
    if implicit == 'measured_zero':
        verified = run(verifier, direct, matrix, observations, features, workspace, request)
        assert f'genes\t2\n' in verified.stdout
        wrong = native(name + '-wrong-request', dict(value, units='different'))
        assert 'original bytes' in run(verifier, direct, matrix, observations, features, workspace, wrong,
                                      success=False).stderr
    materialize = native(name + '-materialize', dict(schema='gadget.assays.materialize.v1', id='geometry',
                         source_namespace='750003', design=design, recipe=value, channel='E'))
    run(gadget, 'derive-assays', '--workspace', workspace, '--request', materialize, '--out', folder / 'derived')
    derived = folder / 'derived/workspace.rfcxs'
    if implicit == 'measured_zero':
        run(verifier, '--assay', 'geometry', derived, matrix, observations, features, workspace, materialize)
    run(gadget, 'replay-derivation', '--workspace', derived, '--assay', 'geometry', '--out', folder / 'replayed')
    if implicit == 'measured_zero':
        run(verifier, '--assay', 'geometry', folder / 'replayed/workspace.rfcxs', matrix, observations,
            features, workspace, materialize)
    # All further queries read a standalone publication with no compatibility side files.
    portable = folder / 'portable.rfcxs'
    shutil.copyfile(derived, portable)
    shutil.rmtree(folder / 'derived')
    shutil.rmtree(folder / 'import')
    return folder, portable


def query(workspace, name, value, success=True):
    request = native(name, dict(schema='gadget.assays.query.v1', operation='native_query', assay='geometry', query=value))
    result = run(gadget, 'assay-query', '--workspace', workspace, '--request', request, success=success)
    if success:
        return json.loads(result.stdout)['result']['answers'][0]['certificate']['result']
    return result


folder, workspace = prepare('chain', [(0,), (1,), (2,), (7,)], [1], 1)
adjoint = query(workspace, 'chain-adjoint', dict(operation='adjoint', context='star:1', field='all.g1'))
assert [Q(row['value']['rational']) for row in adjoint] == [Q(3)]
sheaf = dict(operation='local_sections', context='star:1', construction='cellular_sheaf', metric='identity',
             stalks=[dict(grade=0, key='1', axes=['count']), dict(grade=0, key='2', axes=['count']),
                     dict(grade=1, key='1', axes=['difference'])],
             maps=[dict(grade=1, key='1', matrix=[['1']]), dict(grade=1, key='2', matrix=[['1']])],
             sections=[dict(id='g1', grade=0, fields=['all.g1'], axes=['count'], quantity='umi_count', units='gene_assigned_umi')])
local = query(workspace, 'chain-sections', sheaf)
assert [x['dimension'] for x in local['cohomology']] == [1, 0]
assert list(map(Q, local['sections'][0]['glued'])) == [Q(7, 2), Q(7, 2)]
assert list(map(Q, local['sections'][0]['coboundary_defect'])) == [Q(3)]
assert Q(local['sections'][0]['residual_quadrance']) == Q(9, 2)
green = query(workspace, 'chain-green', dict(operation='green_action', context='star:1', field='all.g1', **{'lambda': '1'}))
assert [Q(row['value']['rational']) for row in green] == [Q(3), Q(4)]
isolated = query(workspace, 'chain-isolated', dict(operation='hodge', context='star:4', field='all.g1'))
assert [Q(row['value']['rational']) for row in isolated] == [Q(0)]
program = native('chain-sequence', dict(schema='gadget.scientific.sequence.v1', assay='geometry',
    context_template=dict(fields=['all.g1', 'all.g2'], methods=['local_sections', 'green_action', 'hodge'],
                          restriction='scalar_identity', metric='identity', green_lambda='1')))
run(sequence, workspace, program, folder / 'sequence.rfcxs')
verified = run(section_verifier, folder / 'sequence.rfcxs', workspace, 'geometry', folder / 'summary.tsv')
assert 'sections\t8\n' in verified.stdout and 'Green_identities\t8\n' in verified.stdout
old = dict(relations=[dict(id=name, grade=1, type='authored.branch', meaning='Original multiway observation',
    participants=[dict(cell='p0', exists=True, negative=True, magnitude='1/2'),
                  dict(cell='p1', exists=True, negative=True, magnitude='1/2'),
                  dict(cell='p2', exists=True, negative=False, magnitude='1'),
                  dict(cell='p0', exists=False, negative=False, magnitude='2'),
                  dict(cell='p2', exists=True, negative=False, magnitude='0')]) for name in ['a', 'b']] +
    [dict(id='closure', grade=2, type='authored.filling', meaning='Original cancelling boundary',
          participants=[dict(cell='a', exists=True, negative=False, magnitude='1'),
                        dict(cell='b', exists=True, negative=True, magnitude='1')])],
    fields=[dict(id='flux', grade=1, channel='T', quantity='transport', units='authored_units',
                 method='Original relation field', values=[dict(relation='a', value='2'), dict(relation='b', value='0')])])
folder, inherited = prepare('inherited', [(0,), (1,), (2,), (7,)], [1], 1, relations=old)
summary = query(inherited, 'inherited-flux', dict(operation='field_summary', field='all.flux', support='reported_only'))
assert summary['reported'] == 2 and Q(summary['sum']['rational']) == 2
prepare('coincident', [(0,), (0,), (4,)], [1], 0)
prepare('anisotropic', [('1/2', 0, 0), ('3/2', 0, 0), (0, 1, 0), (0, 0, 1)], ['1/2', 2, 3], 2)
folder, partial = prepare('partial', [(0,), (1,), (2,), (7,)], [1], 1, implicit='unobserved')
summary = query(partial, 'partial-support', dict(operation='field_summary', field='all.g2', support='reported_only'))
assert summary['reported'] == 2 and Q(summary['sum']['rational']) == 7
query(partial, 'partial-refusal', dict(operation='adjoint', context='star:1', field='all.g2'), success=False)
print('PASS: original matrices, coordinates, occurrence order, measured zeros, absence, local sections, Green actions and complete replay')
