#!/usr/bin/env python3
"""Reproduce source declared structural and field span ablations with GADGET."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


methods = Path(__file__).resolve().parents[1]
source = module('ablation_source', methods/'relational-ensembles/replay.py')
exchange = module('ablation_expression', methods/'structural-calibration/replay.py')
save, digest, run = source.save, source.digest, source.run


def cell_role(name, cells, meaning):
    return dict(id=name, kind='source_cells', cells=cells, meaning=meaning, evidence=[])


def span_role(name, fields, meaning):
    return dict(id=name, kind='field_span', fields=fields, meaning=meaning, evidence=[])


def query(fields, roles, ablations, comparisons, workers):
    return dict(operation='source_structure_ablation', context='all', fields=fields,
                roles=roles, ablations=[dict(id=name, remove_roles=removed) for name, removed in ablations],
                comparisons=[dict(id=name, left=left, right=right) for name, left, right in comparisons],
                homology_output='representatives', **{'lambda': '1'}, workers=workers,
                interpretation=dict(population='The recorded fields and explicitly retained source cells',
                                    meaning='Compare declared source and span deletions on their certified common carriers', evidence=[]))


def authored(root, workers=1):
    source.authored(root/'inputs', field_values=dict(A=[6,2,4,8], B=[2,4,8,6],
                                                    probe=[1,-1,0,0], same_probe=[2,-2,0,0], zero=[0,0,0,0]))
    roles = [cell_role('higher-filling', ['relation:exchange'], 'The declared filling of the four repeated assertions'),
             cell_role('first-assertion', ['relation:assertion-0'], 'The first retained assertion of the multiway interface'),
             span_role('measurement-direction', ['all.probe','all.same_probe'], 'Remove the declared redundant contrast span from each original observed field')]
    return query(['all.A','all.B'], roles,
                 [('remove-filling',['higher-filling']), ('remove-assertion-and-filling',['first-assertion','higher-filling']),
                  ('remove-span',['measurement-direction']), ('remove-span-and-filling',['measurement-direction','higher-filling'])],
                 [('same-source-span','remove-filling','remove-span-and-filling'),
                  ('closed-source-overlap','remove-filling','remove-assertion-and-filling'),
                  ('cover','remove-filling','remove-span')], workers)


def expression(root, packet, panel, workers=1):
    exchange.expression(root, packet, panel, workers)
    selected = json.loads((root/'selected-features.json').read_text())
    relations = json.loads((root/'inputs/relations.json').read_text())['relations']
    pairs = ['relation:'+r['id'] for r in relations if r['id'].startswith('pair-')]
    conditions = ['relation:'+r['id'] for r in relations if r['id'].startswith('condition-')]
    roles = [cell_role('condition-membership', conditions, 'Remove the two recorded condition membership relations'),
             cell_role('donor-pairing', pairs, 'Remove the recorded within-donor comparison relations'),
             span_role('condition-contrast', ['all.treatment_contrast'], 'Remove the fixed dimensionless treated-minus-untreated direction on the original specimen carrier')]
    return query(['all.'+g['id'] for g in selected], roles,
                 [('remove-condition',['condition-membership']), ('remove-pairing',['donor-pairing']),
                  ('remove-all-relations',['condition-membership','donor-pairing']),
                  ('remove-contrast',['condition-contrast']), ('remove-contrast-and-condition',['condition-contrast','condition-membership'])],
                 [('design-overlap','remove-condition','remove-pairing'),
                  ('span-on-fixed-source','remove-condition','remove-contrast-and-condition'),
                  ('condition-cover','remove-condition','remove-contrast')], workers)


def sheaf_declaration(workspace, fields, authored_case):
    paths = list(workspace.rglob('study.json'))
    if len(paths) != 1:
        raise ValueError('This walkthrough requires one complete native source partition')
    study = json.loads(paths[0].read_text())
    grade_count = len(study['contexts'][0]['metrics'])
    if authored_case:
        bindings = [['all.object_A'], ['all.A'], ['all.top_B']]
    else:
        bindings = [fields] + [None]*(grade_count-1)
    stalks = []
    for cell in study['cells']:
        axes = bindings[cell['grade']] or fields
        stalks.append(dict(grade=cell['grade'], key=cell['key'], axes=axes))
    dimension = len(fields) if not authored_case else 1
    identity = [['1' if i == j else '0' for j in range(dimension)] for i in range(dimension)]
    maps = [dict(grade=o['grade'], key=o['key'], matrix=identity) for o in study['occurrences']]
    sections = []
    for grade in range(grade_count):
        section = dict(id='grade-'+str(grade), grade=grade, quantity='declared_measurement' if authored_case else 'read_count',
                       units='dimensionless' if authored_case else 'counts')
        if bindings[grade]:
            section['fields'] = bindings[grade]
        else:
            # An explicit zero cochain supplies no unmeasured relation amplitude.
            section['values'] = ['0']*sum(len(s['axes']) for s in stalks if s['grade'] == grade)
        sections.append(section)
    return dict(stalks=stalks, maps=maps, sections=sections, metric='identity')


def bind_evidence(q, workspace):
    source.bind_evidence(q, workspace)
    for role in q['roles']:
        role['evidence'] = q['interpretation']['evidence']


def compatible_cover(q):
    control=copy.deepcopy(q)
    control['sheaf']['sections'][0].pop('fields')
    control['sheaf']['sections'][0]['values']=['1']*sum(len(s['axes']) for s in control['sheaf']['stalks'] if s['grade']==0)
    control['sheaf']['sections'][1]['fields']=['all.probe']
    for stalk in control['sheaf']['stalks']:
        if stalk['grade']==1:stalk['axes']=['all.probe']
    return control


def execute(binary, root, q, authored_case=True):
    imported, seconds = run(binary, 'import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    save(root/'import-result.json', imported)
    bind_evidence(q, root/'workspace')
    q['sheaf'] = sheaf_declaration(root/'workspace', q['fields'], authored_case)
    save(root/'inputs/query.json', q)
    request = dict(schema='gadget.assays.query.v1', operation='native_query', assay='measurements', query=q)
    result, query_seconds = run(binary, 'assay-query', '--workspace', root/'workspace/workspace.json',
                                '--request', save(root/'request.json', request))
    save(root/'result.json', result)
    save(root/'certificate.json', result['result']['answers'][0]['certificate'])
    save(root/'timing.json', dict(import_seconds=seconds, query_seconds=query_seconds))
    if authored_case:
        control=compatible_cover(q)
        control_request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=control)
        controlled,_=run(binary,'assay-query','--workspace',root/'workspace/workspace.json',
                         '--request',save(root/'compatible-cover-request.json',control_request))
        save(root/'compatible-cover-certificate.json',controlled['result']['answers'][0]['certificate'])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--source', type=Path); p.add_argument('--panel', type=Path)
    p.add_argument('--workers', type=int, default=1)
    args = p.parse_args(); root = args.out.resolve(); root.mkdir(parents=True, exist_ok=False)
    if args.source:
        if not args.panel: p.error('--source requires --panel')
        q = expression(root, args.source.resolve(), args.panel.resolve(), args.workers)
    else:
        q = authored(root, args.workers)
    execute(args.binary.resolve(), root, q, not args.source)
    print(json.dumps(dict(status='completed', output=str(root))))


if __name__ == '__main__': main()
