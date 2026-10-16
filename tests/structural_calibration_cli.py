"""Independent exact assignment, structural observable, law and sensitivity checks."""
import copy
import importlib.util
import itertools
import json
import math
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def main(binary,root):
    shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    path=Path(__file__).resolve().parents[1]/'examples/bioinformatics/structural-calibration/replay.py'
    spec=importlib.util.spec_from_file_location('structural_demo',path);demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)
    q=demo.authored(root,workers=2);checks=refusals=0
    def check(c):
        nonlocal checks
        assert c
        checks+=1
    def command(*args,fail=False):
        nonlocal refusals
        p=subprocess.run([str(binary),*map(str,args)],text=True,capture_output=True,timeout=180)
        if bool(p.returncode)!=fail:raise AssertionError(p.stderr)
        if fail:refusals+=1;return None
        return json.loads(p.stdout)
    command('import-assays','--spec',root/'inputs/import.json','--out',root/'workspace');demo.bind_evidence(q,root/'workspace')
    def execute(query,fail=False,workspace=root/'workspace/workspace.json'):
        request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=query)
        demo.save(root/'request.json',request)
        value=command('assay-query','--workspace',workspace,'--request',root/'request.json',fail=fail)
        if value:return value['result']['answers'][0]['certificate']
    readings=dict(A=[6,2,4,8],B=[2,4,8,6],A_shift=[7,1,4,8],B_shift=[1,5,8,6])
    permutations=list(itertools.permutations(range(4)))
    def verify(cert,diagonal=(1,1,1,1),low=(1,1,1,1),top=1):
        out=cert['result'];d=list(map(Q,diagonal));z=[Q(1),Q(1),Q(-1),Q(-1)];u=[Q(1),Q(-1),Q(0),Q(0)]
        dot=lambda x,y:sum(a*b*c for a,b,c in zip(d,x,y))
        parts=lambda x:[list(x),[sum(x)/sum(1/c for c in d)/c for c in d],[c*dot(z,x)/dot(z,z) for c in z]]
        def hodge(x):
            p=parts(x);p.append([a-b-c for a,b,c in zip(*p)])
            return p
        expected=[]
        for scenario,axes in zip(out['scenarios'],(('A','B'),('A_shift','B_shift'))):
            scores=[]
            for state,perm in zip(scenario['statistics'],permutations):
                ps=[hodge([Q(readings[a][i]) for i in perm]) for a in axes]
                probe=hodge(u)[3]
                proj=lambda x:[v*dot(probe,x)/dot(probe,probe) for v in probe] if dot(probe,probe) else [Q(0)]*4
                A,B=ps;full=dot(A[0],A[0]);g,c,h=[dot(x,x) for x in A[1:]]
                b=[Q(-1),Q(1,3),Q(1,3),Q(1,3)]
                lg=sum(x*x*Q(m) for x,m in zip(b,low))*sum(1/v for v in d);lc=dot(z,z)/Q(top)
                green=[x/(1+lg)+y/(1+lc)+v for x,y,v in zip(A[1],A[2],A[3])]
                change=[x-y for x,y in zip(A[0],green)]
                row=[g/full,c/full,h/full,dot(proj(A[3]),proj(A[3]))/h,dot(change,change)/full,dot(proj(B[3]),proj(B[3]))/dot(B[3],B[3])]
                check(list(map(Q,state['statistics']))==row);scores.append(row)
                for record,fields in zip(state['readings']['fields'],ps):
                    for key,v in zip(('full','gradient','curl','harmonic'),fields):check(list(map(Q,record[key]))==v)
                    check(record['harmonic_kernel_exact'] and record['reconstruction_exact'])
                detail=state['readings']['hypotheses'][3];check(list(map(Q,detail['captured']))==proj(A[3]))
                check(list(map(Q,state['readings']['hypotheses'][4]['green']))==green)
            maxima=[max(row) for row in scores]
            for j,hyp in enumerate(scenario['hypotheses']):
                values=[row[j] for row in scores];obs=values[0]
                tail=Q(sum(x>=obs for x in values),24);ties=Q(values.count(obs),24);adjusted=Q(sum(x>=obs for x in maxima),24)
                check(Q(hyp['observed'])==obs and Q(hyp['reference_tail_mass'])==tail)
                check(Q(hyp['tie_mass'])==ties and Q(hyp['strict_tail_mass'])==tail-ties)
                check(Q(hyp['max_statistic_tail_mass'])==adjusted)
                for value in set(values):
                    p=Q(sum(x>=value for x in values),24)
                    check(Q(sum(Q(sum(y>=x for y in values),24)<=p for x in values),24)<=p)
            check(Q(scenario['joint_family']['global_tail_mass'])==Q(sum(x>=maxima[0] for x in maxima),24))
            expected.append(scores)
        for j,summary in enumerate(out['scenario_sensitivity']):
            values=[v[0][j] for v in expected];mean=sum(values)/2
            check(Q(summary['weighted_mean'])==mean and Q(summary['weighted_variance'])==sum((x-mean)**2 for x in values)/2)
            check(Q(summary['minimum'])==min(values) and Q(summary['maximum'])==max(values))
        check(out['execution']['completed_tiles']==48 and out['execution']['projector_preparations']==2)
        check([s['block_source_coordinates'][0] for s in out['states']]==[list(p) for p in permutations])
        check(all(Q(s['mass'])==Q(1,24) for s in out['states']))
        check(out['probe_frames'][0]['rank']==1 and len(out['coordinate_carrier'])==4)
    baseline=execute(q);verify(baseline)
    serial=copy.deepcopy(q);serial['workers']=1;one=execute(serial);verify(one)
    a,b=copy.deepcopy(baseline['result']),copy.deepcopy(one['result']);a.pop('execution');b.pop('execution');check(a==b)
    compact=copy.deepcopy(q);compact['reference_output']='statistics';v=execute(compact)['result']
    check(v['hypotheses']==baseline['result']['hypotheses'] and v['joint_family']==baseline['result']['joint_family'])
    check(all('readings' not in s for c in v['scenarios'] for s in c['statistics']))
    plan=copy.deepcopy(q);plan['execution']='plan';v=execute(plan)['result'];check(v['status']=='planned' and str(v['cost']['reference_states'])=='24' and 'states' not in v)
    undefined=copy.deepcopy(q);undefined.pop('uncertainty');undefined['fields']=['all.zero'];undefined['hypotheses']=[dict(id='zero',field='all.zero',statistic='sector_fraction',sector='harmonic')]
    v=execute(undefined)['result'];check(v['status']=='undefined_family_statistic' and v['hypotheses'][0]['observed'] is None)
    spec=json.loads((root/'inputs/import.json').read_text());spec['assays'][0]['metrics']=[dict(kind='diagonal',values=['2','3','4','5']),dict(kind='diagonal',values=['1','2','3','4']),dict(kind='diagonal',values=['5'])]
    command('import-assays','--spec',demo.save(root/'inputs/weighted.json',spec),'--out',root/'weighted-workspace')
    weighted=copy.deepcopy(q);demo.bind_evidence(weighted,root/'weighted-workspace');verify(execute(weighted,workspace=root/'weighted-workspace/workspace.json'),(1,2,3,4),(2,3,4,5),5)
    mutations=[
        lambda q:q.update(fields=['all.A','all.A']),lambda q:q.update(fields=['all.A','all.object_A']),
        lambda q:q.update(workers=0),lambda q:q.update(alpha='0'),lambda q:q.update(execution='sampled'),lambda q:q.update(reference_output='reduced'),
        lambda q:q['law'].update(kind='inferred_exchange'),lambda q:q['law'].update(evidence=[]),
        lambda q:q['law']['evidence'][0].update(text='invented'),lambda q:q['law'].update(assumptions=[]),
        lambda q:q['law']['blocks'][0].update(relation='relation:assertion-0'),lambda q:q['law']['blocks'][0].update(relation='missing'),
        lambda q:q['law']['blocks'][0].update(cells=['relation:assertion-0']),
        lambda q:q['law']['blocks'][0].update(cells=['relation:assertion-0']*4),
        lambda q:q['law']['blocks'].append(dict(id='overlap',relation='relation:exchange',cells=[f'relation:assertion-{i}' for i in range(4)])),
        lambda q:q['hypotheses'][0].update(sector='latent'),lambda q:q['hypotheses'][0].update(field='all.zero'),
        lambda q:q['hypotheses'][4].update(**{'lambda':'-1'}),lambda q:q['hypotheses'][3].update(probes=[]),
        lambda q:q['hypotheses'][3].update(probes=['all.object_A']),lambda q:q['hypotheses'][1].update(id='gradient-A'),
        lambda q:q['uncertainty']['scenarios'][0].update(fields=['all.A_shift','all.B']),
        lambda q:q['uncertainty']['scenarios'][1].update(weight='0'),lambda q:q['uncertainty'].update(observed_scenario='missing'),
        lambda q:q['uncertainty'].update(evidence=[]),lambda q:q['uncertainty']['scenarios'][1].update(fields=['all.A_shift'])]
    for mutate in mutations:
        bad=copy.deepcopy(q);mutate(bad);execute(bad,True)
    order=copy.deepcopy(q);order['uncertainty']['scenarios'].reverse();v=execute(order)['result']
    check(v['observed_scenario']=='original' and v['hypotheses']==baseline['result']['hypotheses'])
    check(v['scenario_sensitivity']==baseline['result']['scenario_sensitivity'])
    check(v['scenario_observable_ensemble']['mean']==baseline['result']['scenario_observable_ensemble']['mean'])
    # Declared block order is an assignment coordinate choice, not a source rewrite.
    reversed_block=copy.deepcopy(q);reversed_block['law']['blocks'][0]['cells'].reverse();v=execute(reversed_block)['result']
    check(v['hypotheses']==baseline['result']['hypotheses'] and v['states'][0]['block_source_coordinates'][0]==[3,2,1,0])
    # Retained occurrence participants differ from net boundary support.
    altered=json.loads((root/'inputs/relations.json').read_text())
    altered['relations'][-1]['participants']=[demo.source.participant('assertion-0',1),demo.source.participant('assertion-1',1),demo.source.participant('assertion-2',-2),demo.source.participant('assertion-3',1),demo.source.participant('assertion-3',-1),demo.source.participant('assertion-2',7,False)]
    path=demo.save(root/'inputs/cancellation.json',altered)
    spec=json.loads((root/'inputs/import.json').read_text());spec['assays'][0]['relations']=dict(path=path.name,sha256=demo.digest(path))
    command('import-assays','--spec',demo.save(root/'inputs/cancellation-import.json',spec),'--out',root/'cancel-workspace')
    cancellation=copy.deepcopy(q);cancellation.pop('uncertainty');cancellation['hypotheses']=[dict(id='full',field='all.A',statistic='sector_fraction',sector='full')];demo.bind_evidence(cancellation,root/'cancel-workspace')
    v=execute(cancellation,workspace=root/'cancel-workspace/workspace.json')['result']
    check(len(v['states'])==24 and len(v['exchange_blocks'][0]['source_occurrences'])==6)
    check(not v['exchange_blocks'][0]['source_occurrences'][-1]['exists'])
    missing_cancel=copy.deepcopy(cancellation);missing_cancel['law']['blocks'][0]['cells']=missing_cancel['law']['blocks'][0]['cells'][:3]
    execute(missing_cancel,True,root/'cancel-workspace/workspace.json')
    # A supported zero coefficient is still an existing participant.
    altered['relations'][-1]['participants'][-1]=demo.source.participant('assertion-3',0)
    path=demo.save(root/'inputs/supported-zero.json',altered);spec['assays'][0]['relations']=dict(path=path.name,sha256=demo.digest(path))
    command('import-assays','--spec',demo.save(root/'inputs/supported-zero-import.json',spec),'--out',root/'zero-occ-workspace')
    demo.bind_evidence(cancellation,root/'zero-occ-workspace');v=execute(cancellation,workspace=root/'zero-occ-workspace/workspace.json')['result']
    check(v['exchange_blocks'][0]['source_occurrences'][-1]['exists'] and v['exchange_blocks'][0]['source_occurrences'][-1]['coefficient']=='0')
    # Missing measured values are refused rather than filled by a law constructor.
    altered=json.loads((root/'inputs/relations.json').read_text());altered['fields'][0]['values'].pop()
    path=demo.save(root/'inputs/partial.json',altered);spec['assays'][0]['relations']=dict(path=path.name,sha256=demo.digest(path))
    command('import-assays','--spec',demo.save(root/'inputs/partial-import.json',spec),'--out',root/'partial-workspace')
    partial=copy.deepcopy(q);demo.bind_evidence(partial,root/'partial-workspace');execute(partial,True,root/'partial-workspace/workspace.json')
    # Zero probe span is defined as zero captured energy when the sector is nonzero.
    z=copy.deepcopy(q);z['hypotheses']=[dict(id='zero-probe',field='all.A',statistic='probe_fraction',sector='harmonic',probes=['all.zero'])]
    v=execute(z)['result'];check(v['hypotheses'][0]['observed']=='0' and v['hypotheses'][0]['reference_tail_mass']=='1' and v['probe_frames'][0]['rank']==0)
    # Repeated probe declarations in different hypotheses share their native solve.
    repeated=copy.deepcopy(q);repeated['hypotheses'][3]['probes']=['all.probe'];v=execute(repeated)['result']
    check(v['execution']['probe_projector_preparations']==1 and v['execution']['probe_field_decompositions']==1)
    # A denominator can vanish in only part of the complete law. Do not drop states.
    unstable=copy.deepcopy(q);unstable.pop('uncertainty');unstable['fields']=['all.probe'];unstable['hypotheses']=[dict(id='curl-probe',field='all.probe',statistic='probe_fraction',sector='curl',probes=['all.A'])]
    v=execute(unstable)['result'];check(v['status']=='undefined_family_statistic' and len(v['states'])==24 and v['hypotheses'][0]['reference_tail_mass'] is None)
    # Cardinality beyond 64 bit addressability can be inspected without enumeration.
    demo.source.authored(root/'wide/inputs',arity=21)
    command('import-assays','--spec',root/'wide/inputs/import.json','--out',root/'wide/workspace')
    huge=demo.query(['all.object_A'],[dict(id='all',relation='relation:assertion-0',cells=[f'feature-{i}' for i in range(21)])],[dict(id='full',field='all.object_A',statistic='sector_fraction',sector='full')])
    huge['execution']='plan';demo.bind_evidence(huge,root/'wide/workspace')
    v=execute(huge,workspace=root/'wide/workspace/workspace.json')['result'];check(v['cost']['reference_states']==str(math.factorial(21)))
    huge['execution']='exhaustive';execute(huge,True,root/'wide/workspace/workspace.json')
    print(json.dumps(dict(status='PASS',checks=checks,refusals=refusals)))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
