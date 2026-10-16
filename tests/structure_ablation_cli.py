#!/usr/bin/env python3
"""Contract and source admissibility checks for structural role ablations."""
import copy
import importlib.util
import json
import subprocess
import sys
import shutil
from fractions import Fraction as Q
from pathlib import Path

path = Path(__file__).resolve().parents[1]/'examples/bioinformatics/structural-ablations/replay.py'
spec = importlib.util.spec_from_file_location('ablation_demo', path)
demo = importlib.util.module_from_spec(spec); spec.loader.exec_module(demo)


def main(binary, root):
    root.mkdir(parents=True, exist_ok=True)
    # These are generated fixture directories owned by this contract suite.
    for name in ('authored', 'weighted', 'wide'):
        fixture = root/name
        if fixture.is_dir():
            shutil.rmtree(fixture)
    checks = refusals = 0
    def check(value):
        nonlocal checks
        assert value
        checks += 1
    def command(*args):
        p = subprocess.run([str(binary), *map(str,args)], capture_output=True, text=True)
        if p.returncode: raise AssertionError(p.stderr)
        return json.loads(p.stdout)
    case = root/'authored'; case.mkdir(exist_ok=True)
    q = demo.authored(case)
    out = demo.execute(binary, case, q)['result']['answers'][0]['certificate']['result']
    def execute(query, refuse=False, workspace=None):
        nonlocal refusals
        request = dict(schema='gadget.assays.query.v1', operation='native_query', assay='measurements', query=query)
        p = subprocess.run([str(binary),'assay-query','--workspace',str(workspace or case/'workspace/workspace.json'),
                            '--request',str(demo.save(root/'request.json',request))], capture_output=True, text=True)
        if refuse:
            check(p.returncode != 0); refusals += 1; return
        if p.returncode: raise AssertionError(p.stderr)
        return json.loads(p.stdout)['result']['answers'][0]['certificate']['result']
    check([g['betti'] for g in out['original_grades']] == [3,2,0])
    check(out['execution']['source_preparations'] == 3)
    check(out['execution']['harmonic_basis_requested'])
    unfilled = out['ablations'][0]['grades'][1]
    check(unfilled['betti'] == 3 and unfilled['homology_image_rank'] == 2 and unfilled['killed_by_inclusion'] == 1)
    check(out['ablations'][0]['grades'][2]['cells'] == 0)
    check(unfilled['ambient_fillings']['columns'] == 1)
    check(out['ablations'][0]['grades'][2]['relative_betti'] == 1)
    check(out['ablations'][1]['grades'][1]['betti'] == 2 and out['ablations'][1]['grades'][1]['killed_by_inclusion'] == 0)
    check(out['ablations'][2]['removed_span_rank'] == 1)
    check(list(map(Q,out['ablations'][2]['fields'][0]['removed_span_field'])) == [2,-2,0,0])
    check(list(map(Q,out['ablations'][2]['fields'][1]['removed_span_field'])) == [-1,1,0,0])
    check(out['comparisons'][0]['sheaf_overlap_compatible'])
    check(not out['comparisons'][2]['sheaf_overlap_compatible'] and out['comparisons'][2]['union_covers_ambient'])
    check(any(not o['exists'] for o in out['source']['occurrences']))
    for a in out['ablations']:
        for f in a['fields']:
            x=list(map(Q,f['sectors']['full']))
            check(all(sum(Q(f['sectors'][s][i]) for s in ('gradient','curl','harmonic')) == v for i,v in enumerate(x)))
            check(sum(Q(f['quadrances'][s]) for s in ('gradient','curl','harmonic')) == Q(f['quadrances']['full']))
            check(all(Q(s)+Q(t)==Q(v) for s,t,v in zip(f['green']['structure_term'],f['green']['field_term'],f['green']['change'])))
    # Complete native section gluing, with a declared compatible control.
    control=copy.deepcopy(q)
    control['sheaf']['sections'][0].pop('fields');control['sheaf']['sections'][0]['values']=['1']*4
    control['sheaf']['sections'][1]['fields']=['all.probe']
    for s in control['sheaf']['stalks']:
        if s['grade']==1:s['axes']=['all.probe']
    glued=execute(control)
    check(glued['comparisons'][2]['sheaf_overlap_compatible'] and glued['comparisons'][2]['global_cochain_verified'])
    check(len(glued['comparisons'][2]['glued_sections']) == 9)
    # Ranks mode and native worker count change neither numerical fields nor images.
    ranks=copy.deepcopy(q);ranks['homology_output']='ranks';ranks['workers']=2
    simple=execute(ranks);check(not simple['execution']['harmonic_basis_requested'])
    for a,b in zip(simple['ablations'],out['ablations']):
        check(a['fields']==b['fields'])
        for g,h in zip(a['grades'],b['grades']):
            check(g['betti']==h['betti'] and g['homology_image_rank']==h['homology_image_rank'])
            check('harmonic_representatives' not in g)
    integral=copy.deepcopy(q);integral['lattice_snapshots']=True
    li=execute(integral)
    check(all(g['lattice']['free_rank']==g['betti'] for g in li['original_grades']))
    check(all(g['lattice']['free_rank']==g['betti'] for a in li['ablations'] for g in a['grades']))
    # The inherited positive diagonal form changes field projections while
    # preserving the source inclusion and rational homology ranks.
    original_study=next((case/'workspace').rglob('study.json'));weighted=root/'weighted'
    shutil.copytree(original_study.parent,weighted)
    declaration=json.loads((weighted/'study.json').read_text())
    declaration['contexts'][0]['metrics']=[dict(kind='diagonal',values=values) for values in (['1','2','3','4'],['2','3','5','7'],['11'])]
    demo.save(weighted/'study.json',declaration)
    answer=command('query','--study',weighted/'study.json','--query',demo.save(root/'weighted-query.json',q))
    weighted_result=answer['certificate']['result'] if 'certificate' in answer else answer['result']
    check(list(map(Q,weighted_result['ablations'][2]['fields'][0]['removed_span_field']))==[Q(6,5),Q(-6,5),0,0])
    check([g['betti'] for g in weighted_result['original_grades']]==[3,2,0])
    for a,b in zip(weighted_result['ablations'],out['ablations']):
        check([g['homology_image_rank'] for g in a['grades']]==[g['homology_image_rank'] for g in b['grades']])
    variants=[]
    def bad(mutator):
        v=copy.deepcopy(q);mutator(v);variants.append(v)
    bad(lambda v:v.update(fields=[])); bad(lambda v:v.update(fields=['all.A','all.A']))
    bad(lambda v:v.update(workers=0)); bad(lambda v:v.update(**{'lambda':'0'}))
    bad(lambda v:v.update(homology_output='embedded'));bad(lambda v:v.update(roles=[]))
    bad(lambda v:v['roles'][0].update(cells=['missing']))
    bad(lambda v:v['roles'][0].update(cells=['relation:exchange','relation:exchange']))
    bad(lambda v:v['roles'][0].update(kind='inferred_biology'))
    bad(lambda v:v['roles'][0].update(fields=['all.A']))
    bad(lambda v:v['roles'][2].update(fields=['all.top_B']))
    bad(lambda v:v['roles'][2].update(fields=['all.probe','all.probe']))
    bad(lambda v:v['roles'][2].update(fields=['all.missing']))
    bad(lambda v:v['ablations'][1].update(remove_roles=['first-assertion']))
    bad(lambda v:v['roles'][1].update(cells=['feature-2'])) # Inactive retained attachments still require closure.
    bad(lambda v:v['roles'][1].update(cells=['feature-1'])) # Cancelling occurrences still require closure.
    bad(lambda v:v['ablations'][0].update(remove_roles=['not-a-role']))
    bad(lambda v:v['ablations'][0].update(remove_roles=['higher-filling','higher-filling']))
    bad(lambda v:v['ablations'][1].update(id=v['ablations'][0]['id']))
    bad(lambda v:v['roles'].append(copy.deepcopy(v['roles'][0])))
    bad(lambda v:v['roles'].append(demo.span_role('unused',['all.zero'],'Unused declaration')))
    bad(lambda v:v['comparisons'][0].update(left='absent'))
    bad(lambda v:v['comparisons'][1].update(id=v['comparisons'][0]['id']))
    bad(lambda v:v['interpretation'].update(evidence=[dict(document='absent',start=0,end=1,text='x')]))
    bad(lambda v:v['sheaf']['maps'][0].update(matrix=[['2']])) # Original upward maps fail the complete cochain law.
    bad(lambda v:v['sheaf']['maps'].pop())
    bad(lambda v:v['sheaf']['sections'].pop())
    bad(lambda v:v['sheaf']['sections'][1].update(grade=0))
    bad(lambda v:v['sheaf'].update(metric='imputed'))
    bad(lambda v:v.update(lattice_snapshots='yes'))
    bad(lambda v:v.update(unknown='member'))
    for v in variants:execute(v,True)
    # All grade deletion admits the zero carrier, rather than inventing a cell.
    empty=copy.deepcopy(q);empty.pop('sheaf');empty.pop('comparisons')
    empty['roles']=[demo.cell_role('all',[c['id'] for grade in out['source']['cells'] for c in grade],'All explicitly identified source cells')]
    empty['ablations']=[dict(id='zero-carrier',remove_roles=['all'])]
    demo.bind_evidence(empty,case/'workspace')
    z=execute(empty);check(all(g['cells']==0 and g['betti']==0 for g in z['ablations'][0]['grades']))
    check(all(not f['restricted_input'] for f in z['ablations'][0]['fields']))
    # Sparse exact actions with more than 4096 coordinates; no harmonic frame.
    wide=root/'wide';demo.source.authored(wide/'inputs',arity=4101,filled=False)
    command('import-assays','--spec',wide/'inputs/import.json','--out',wide/'workspace')
    w=demo.query(['all.object_A'],[demo.cell_role('all-relations',[f'relation:assertion-{i}' for i in range(4)],'Remove the four declared primary assertions')],
                 [('remove-relations',['all-relations'])],[],2);w['homology_output']='ranks';demo.bind_evidence(w,wide/'workspace')
    v=execute(w,workspace=wide/'workspace/workspace.json');check(v['ablations'][0]['grades'][0]['betti']==4101)
    check(len(v['ablations'][0]['fields'][0]['sectors']['harmonic'])==4101)
    print(json.dumps(dict(status='PASS',checks=checks,refusals=refusals)))


if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
