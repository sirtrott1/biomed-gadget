#!/usr/bin/env python3
"""Replay declared native observations without constructing an inverse."""
import argparse
import copy
import csv
import importlib.util
import json
import subprocess
from fractions import Fraction as F
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / 'declared-responses/replay.py'
_spec = importlib.util.spec_from_file_location('observation_sources', _path)
sources = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sources)
save, run, identity, matrix = sources.save, sources.run, sources.identity, sources.matrix


def specification(name, source, grade, n, accession, operator=None):
    fields = [sources.field('x', 100, 10, ['1'] + ['0'] * (n-1) if n else []),
              sources.field('x2', 101, 10, [F(i-2, i+1) for i in range(n)]),
              sources.field('zero', 102, 10, ['0'] * n),
              sources.field('partial', 103, 10, [None] * n)]
    category = sources.field('category', 104, 10, ['1'] * n)
    category.update(value_type='categorical', codebook={'1': 'declared_category'})
    fields.append(category)
    operation = sources.op('R', 1, 'resolving', 10, 10, operator,
                           **({} if operator is not None else dict(hodge_grade=grade, identity_shift='1')))
    return dict(id=name, namespace='9311', sequence='7', source=source,
                carriers=[dict(key='10', grade=grade, size=n),
                          dict(key='20', coordinate_system='200', size=accession['rows'])],
                operators=[operation, sources.op('S', 2, 'accession', 10, 20, accession)], fields=fields)


def execute(binary, producer, root, spec, observation=None, original=None):
    root.mkdir(parents=True, exist_ok=False)
    save(root/'spec.json', spec)
    command = [str(producer), str(root/'spec.json'), str(root/'study')]
    if original:
        command.append(str(original))
    subprocess.run(command, check=True, capture_output=True, text=True)
    query = dict(operation='linear_response_observation', publication='self',
                 source_sequence=spec['sequence'], operator='R', accession='S',
                 fields=['all.time', 'all.age'] if original else ['x', 'x2', 'zero', 'x'],
                 ordering=dict(kind='minimum_degree'), observation=observation or dict(kind='linear'),
                 workers=2)
    save(root/'query.json', query)
    save(root/'result.json', run(binary, 'query', '--study', root/'study/study.json',
                               '--query', root/'query.json'))


def examples(binary, producer, out, clinical=None, sparse_size=257):
    old, new, *_ = sources.square()
    accession = matrix([[1, 1, 0, 0], [0, 1, -1, 0], [1, 1, 0, 0],
                        [0, 0, 0, 0], [2, 2, 0, 0]])
    for name, source in [('unfilled-square', old['source']), ('filled-square', new['source'])]:
        execute(binary, producer, out/name, specification(name, source, 1, 4, accession))
    selected = [3, 1, 3, 0, 1, 3]
    selector = dict(rows=len(selected), columns=4,
                    entries=[dict(row=i, column=j, value='1') for i, j in enumerate(selected)])
    execute(binary, producer, out/'ordered-coordinates', specification('ordered-coordinates', old['source'], 1, 4, selector),
            dict(kind='coordinates', coordinates=selected))
    execute(binary, producer, out/'complete-reading', specification('complete-reading', old['source'], 1, 4, identity(4)))
    execute(binary, producer, out/'nonsymmetric', specification('nonsymmetric', sources.source(3), 0, 3,
        matrix([[2, -1, 0], [0, 0, 1]]), matrix([[2, 1, 0], [0, 3, 1], [1, 0, 4]])))
    case = out/'nonsymmetric'
    query = json.loads((case/'query.json').read_text())
    query['ordering'] = dict(kind='declared', rows=[2, 0, 1], columns=[1, 2, 0])
    save(case/'declared-query.json', query)
    save(case/'declared-result.json', run(binary, 'query', '--study', case/'study/study.json',
                                        '--query', case/'declared-query.json'))
    execute(binary, producer, out/'empty-output', specification('empty-output', sources.source(3), 0, 3,
                                                             dict(rows=0, columns=3, entries=[]), identity(3)))
    execute(binary, producer, out/'empty-source', specification('empty-source', sources.source(0), 0, 0, identity(0)))
    n = sparse_size
    assert n > 0
    tridiagonal = dict(rows=n, columns=n, entries=[dict(row=i, column=j, value=str(3 if i == j else -1))
                      for i in range(n) for j in range(max(0, i-1), min(n, i+2))])
    read = dict(rows=5, columns=n, entries=[dict(row=i, column=j, value=str(v))
                for i, row in enumerate([{0: 1}, {n-1: 1}, {0: 1}, {}, {0: 2}]) for j, v in row.items()])
    spec = specification('sparse-source', sources.source(n), 0, n, read, tridiagonal)
    # Independently planted complete source solutions qualify a larger sparse action.
    for field in spec['fields'][:2]:
        x = [F(i+1) if field['id'] == 'x' else F((-1)**i, i+1) for i in range(n)]
        field['values'] = [str(3*x[i] - (x[i-1] if i else 0) - (x[i+1] if i+1 < n else 0)) for i in range(n)]
    execute(binary, producer, out/'sparse-source', spec)
    if clinical:
        original = clinical/'p0/study.json'
        declaration = json.loads(original.read_text())
        with (original.parent/'table.txt').open(newline='') as f:
            rows = list(csv.DictReader(f))
        groups = [[i for i, row in enumerate(rows) if F(row['sex']) == sex] for sex in [1, 2]]
        entries = []
        for i, group in enumerate(groups):
            entries += [dict(row=i, column=j, value=str(F(1, len(group)))) for j in group]
            entries += [dict(row=2, column=j, value=str(F(1 if i == 0 else -1, len(group)))) for j in group]
        entries += [dict(row=3, column=j, value=str(F(1, len(groups[0])))) for j in groups[0]]
        spec = dict(id='recorded-followup-observations', namespace=declaration['native']['source_namespace'],
                    sequence=run(binary, 'query', '--study', original, '--query', save(out/'clinical-inspect.json', dict(operation='inspect')))['source_sequence'], carriers=[dict(key='20000', coordinate_system='20000', size=5)],
                    operators=[sources.op('R', 20001, 'resolving', 1, 1, hodge_grade=0, identity_shift='1',
                        method='Identity plus grade-zero Hodge action on the retained sex-membership source'),
                        sources.op('S', 20002, 'accession', 1, 20000, dict(rows=5, columns=len(rows), entries=entries),
                                   method='Sex-category mean readings, category difference, repeated first mean and zero row')], fields=[])
        execute(binary, producer, out/'recorded-followup', spec, original=original)
    from verify import verify
    records = [verify(p) for p in sorted(out.iterdir()) if (p/'result.json').is_file()]
    save(out/'verification.json', dict(status='PASS', cases=records,
                                      output_coordinates=sum(r['output_coordinates'] for r in records)))
    return records


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['binary', 'producer', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--clinical', type=Path)
    p.add_argument('--sparse-size', type=int, default=257)
    a = p.parse_args()
    records = examples(a.binary.resolve(), a.producer.resolve(), a.out.resolve(), a.clinical, a.sparse_size)
    print(json.dumps(dict(status='PASS', cases=len(records))))
