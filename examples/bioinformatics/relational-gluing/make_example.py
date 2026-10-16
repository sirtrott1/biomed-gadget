"""Author a declared gluing diagram; no source identity is inferred from labels."""
import argparse
import hashlib
import json
from pathlib import Path


def make_example(root, vertices=3, repetitions=1, modules=1):
    if modules < 1 or vertices < 2 * modules or repetitions < 1:
        raise ValueError('Each module needs two interface vertices; module and section counts must be positive')
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)

    def save(name, value):
        path = root / name
        path.write_text(value if isinstance(value, str) else json.dumps(value, indent=2) + '\n')
        return dict(path=name, sha256='sha256:' + hashlib.sha256(path.read_bytes()).hexdigest())

    def leg(cell, negative=False, magnitude='1', exists=True):
        return dict(cell=cell, exists=exists, negative=negative, magnitude=magnitude)

    def relation(id, grade, legs):
        return dict(id=id, grade=grade, meaning='Authored structural relation ' + id, participants=legs)

    roles = ['common', 'left', 'right']
    def name(role, part, module):
        base = 'interface' if role == 'common' and part == 'ab' else role + '-' + part
        return base + (f'-{module}' if module else '')
    shared = {role: [f'{role}-shared-{i}' for i in range(vertices)] for role in roles}
    orders = {role: shared[role][:] for role in roles}
    orders['left'] += [name('left', 'unique', m) for m in range(modules)]
    # Deliberately change local order, names, and native grade zero keys.
    orders['right'] = list(reversed(orders['right'])) + [name('right', 'unique', m) for m in range(modules)]
    relations, assays, source_fields = {}, [], {}
    for role in roles:
        records, faces, edge_values, face_values = [], [], {}, {}
        for m in range(modules):
            a, b = shared[role][2*m:2*m+2]
            ab = name(role, 'ab', m)
            legs = [leg(a, True), leg(b), leg(a, False, '1/3'), leg(a, True, '1/3'),
                    leg(b, False, '2', False), leg(b, False, '0')]
            records.append(relation(ab, 1, legs)); edge_values[ab] = '2'
            if role != 'common':
                c = name(role, 'unique', m)
                bc, ac, multi, face = [name(role, part, m) for part in ['bc', 'ac', 'multiway', 'face']]
                records += [relation(bc, 1, [leg(b, True), leg(c)]),
                            relation(ac, 1, [leg(a, True), leg(c)]),
                            relation(multi, 1, [leg(a, True), leg(b, True), leg(c, magnitude='2')])]
                faces.append(relation(face, 2, [leg(ab), leg(bc), leg(ac, True)]))
                edge_values.update({bc: '3', ac: '5', multi: '7'})
                face_values[face] = '7' if role == 'left' else '11'
        # The importer retains the original increasing grade declaration order.
        records += faces
        fields = []
        for repeat in range(repetitions):
            if role != 'common':
                fields += [dict(id=f'edges-{repeat}', grade=1, channel='observation',
                                quantity='declared_relation_signal', units='1',
                                method='Authored scalar section',
                                values=[dict(relation=id, value=v) for id, v in edge_values.items()]),
                           dict(id=f'faces-{repeat}', grade=2, channel='observation',
                                quantity='declared_relation_signal', units='1',
                                method='Authored scalar section',
                                values=[dict(relation=id, value=v) for id, v in face_values.items()])]
        relations[role] = records
        features = save(role + '-features.tsv', 'id\n' + '\n'.join(orders[role]) + '\n')
        observations = save(role + '-observations.tsv', 'id\n' + '\n'.join(f'signal-{i}' for i in range(repetitions)) + '\n')
        values = {id: str(i * 10) for i, id in enumerate(shared[role])}
        if role != 'common':
            values.update({name(role, 'unique', m): '30' if role == 'left' else '40' for m in range(modules)})
        source_fields[role] = {f'all.signal-{i}': dict(grade=0, values=values) for i in range(repetitions)}
        if role != 'common':
            source_fields[role].update({f'all.{kind}-{i}': dict(grade=g, values=v)
                                       for i in range(repetitions)
                                       for kind, g, v in [('edges', 1, edge_values), ('faces', 2, face_values)]})
        matrix = save(role + '-signal.tsv', 'id\t' + '\t'.join(f'signal-{i}' for i in range(repetitions)) + '\n' +
                      ''.join(id + '\t' + '\t'.join([values[id]] * repetitions) + '\n' for id in orders[role]))
        rel = save(role + '-relations.json', dict(relations=records, fields=fields))
        assays.append(dict(id=role, kind='declared_structural_assay', format='dense_tsv',
                           features=features, observations=observations, source=matrix,
                           source_namespace=str(10801000 + roles.index(role)), carrier='features',
                           quantity='declared_observation_signal', units='1', channel='observation',
                           method='Constructed gluing example', arithmetic='exact_rational',
                           implicit_entries='unobserved', missing_tokens=['NA'], relations=rel))
    design = dict(question='Combine two declared overlapping complexes and compatible scalar fields',
                  population='Constructed interface and two participant structures', selection='Every authored cell and occurrence',
                  limitations=['Constructed method example; these values are not biological measurements'])
    spec = dict(schema='gadget.assays.import.v1', id='declared-gluing', design=design, assays=assays, resources=[], links=[])
    save('import.json', spec)
    sources = {role: dict(key=str(i + 1), assay=role, partition=0, available_sequence='0') for i, role in enumerate(roles)}
    embeddings = {}
    for i, role in enumerate(['left', 'right']):
        cells = [dict(grade=0, **{'from': x, 'to': y}) for x, y in zip(shared['common'], shared[role])]
        cells += [dict(grade=1, **{'from': 'relation:' + name('common', 'ab', m),
                                  'to': 'relation:' + name(role, 'ab', m)}) for m in range(modules)]
        embeddings['common_' + role] = dict(key=str(11 + i), cells=cells,
            occurrences=[dict(grade=1, **{'from': str(6*m + k), 'to': str(13*m + k)})
                         for m in range(modules) for k in range(1, 7)])
    sections = []
    for i in range(repetitions):
        for grade, field, quantity in [(0, f'signal-{i}', 'declared_observation_signal'),
                                       (1, f'edges-{i}', 'declared_relation_signal'),
                                       (2, f'faces-{i}', 'declared_relation_signal')]:
            sections.append(dict(id=f'grade-{grade}-{i}', key=str(101 + len(sections)),
                                 left_field='all.' + field, right_field='all.' + field,
                                 quantity=quantity, units='1', channel='observation',
                                 method='Literal same-field gluing', evidence='The authored common-interface maps'))
    query = dict(schema='gadget.assays.query.v1', operation='relational_gluing',
                 declaration=dict(key='21', revision='1', pushout='22', method='Declared pushout',
                                  evidence='Complete interface embeddings in this request'),
                 sources=sources, embeddings=embeddings, sections=sections,
                 boundary_policy='to_grade_zero', workers=2)
    save('query.json', query)
    parameters = {k: v for k, v in query.items() if k not in ['schema', 'operation']}
    save('pipeline.json', dict(schema='gadget.pipeline.default.v1', profile='multi_omics',
                               question=design['question'], goals=[dict(id='glue', method='relational_gluing', parameters=parameters)]))
    save('reference.json', dict(orders=orders, relations=relations, source_fields=source_fields, vertices=vertices,
                                modules=modules, repetitions=repetitions, expected_cells=[vertices + 2*modules, 7*modules, 2*modules],
                                expected_occurrences=[0, 20*modules, 6*modules]))
    return spec, query


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--vertices', type=int)
    parser.add_argument('--repetitions', type=int, default=1)
    parser.add_argument('--modules', type=int, default=1)
    args = parser.parse_args()
    make_example(args.out, args.vertices if args.vertices is not None else max(3, 2*args.modules),
                 args.repetitions, args.modules)
