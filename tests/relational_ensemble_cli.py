"""Independent rational source, Hodge, moment, group and refusal contracts."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def main(binary, root):
    shutil.rmtree(root,ignore_errors=True)
    path=Path(__file__).resolve().parents[1]/'examples/bioinformatics/relational-ensembles/replay.py'
    spec=importlib.util.spec_from_file_location('ensemble_demo',path)
    demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)
    readings=dict(A=[6,2,4,8],B=[2,4,8,6],C=[3,7,5,1],zero=[0,0,0,0],same_A=[6,2,4,8],partial=[1,None,2,3])
    q=demo.authored(root/'inputs',field_values=readings)
    checks=refused=0
    def check(condition):
        nonlocal checks
        assert condition
        checks+=1
    def command(*args,fail=False):
        nonlocal refused
        p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True,timeout=180)
        check(bool(p.returncode)==fail)
        if fail:refused+=1;return None
        return json.loads(p.stdout)
    command('import-assays','--spec',root/'inputs/import.json','--out',root/'workspace')
    demo.bind_evidence(q,root/'workspace')
    workspace=root/'workspace/workspace.json'
    def execute(value,fail=False,workspace_=workspace):
        request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=value)
        demo.save(root/'request.json',request)
        answer=command('assay-query','--workspace',workspace_,'--request',root/'request.json',fail=fail)
        if answer:return answer['result']['answers'][0]['certificate']
    def vec(xs):return list(map(Q,xs))
    def dot(x,y):return sum((a*b for a,b in zip(x,y)),Q(0))
    def verify(cert,filled=True):
        result=cert['result'];request=cert['query'];weights=vec(result['weights'])
        channels={};z=[Q(1),Q(1),Q(-1),Q(-1)]
        metric_record=result['metric'];diagonal=[Q(0)]*4
        for e in metric_record['entries']:
            check(e['row']==e['column']);diagonal[e['row']]=Q(e['value'])
        def inner(x,y):return dot(x,[a*b for a,b in zip(diagonal,y)])
        for item in result['samples']:
            x=list(map(Q,readings[item['field'].removeprefix('all.')]))
            grad=[sum(x)/sum(1/d for d in diagonal)/d for d in diagonal]
            curl=[v*inner(z,x)/inner(z,z) for v in z] if filled else [Q(0)]*4
            harmonic=[v-g-c for v,g,c in zip(x,grad,curl)]
            parts=[x,grad,curl,harmonic];channels[item['field']]=parts
            for name,v in zip(('full','gradient','curl','harmonic'),parts):check(vec(item[name])==v)
            check(all(inner(a,b)==0 for a,b in ((grad,curl),(grad,harmonic),(curl,harmonic))))
            arity=len(result['coordinate_carriers'][0]);b=[Q(-1)]+[Q(1,arity-1)]*(arity-1)
            check(vec(item['boundary'])==[v*sum(x) for v in b])
            mg0=cert['context']['metrics'][0];md0=[Q(1)]*arity if mg0['kind']=='identity' else list(map(Q,mg0['values']))
            mg2=cert['context']['metrics'][2] if filled else None
            md2=Q(1) if not mg2 or mg2['kind']=='identity' else Q(mg2['values'][0])
            down=[dot(b,[v*d for v,d in zip(b,md0)])*sum(x)/d for d in diagonal]
            up=[v*inner(z,x)/md2 for v in z] if filled else [Q(0)]*4
            check(vec(item['hodge_down'])==down and vec(item['hodge_up'])==up)
            for response in item['green']:
                l=Q(response['lambda']);lg=dot(b,[v*d for v,d in zip(b,md0)])*sum(1/d for d in diagonal)
                lc=inner(z,z)/md2 if filled else Q(0)
                expected=[g/(1+l*lg)+c/(1+l*lc)+h for g,c,h in zip(grad,curl,harmonic)]
                check(vec(response['field'])==expected and response['residual_exact'])
            check(item['dirac_square_exact'] and item['reconstruction_exact'])
        for group in result['groups']:
            w=vec(group['conditional_weights']);ids=group['fields'];var=[];dep=[]
            mass=sum(weights[request['fields'].index(name)] for name in ids)
            check(Q(group['original_mass'])==mass and w==[weights[request['fields'].index(name)]/mass for name in ids])
            for ch,name in enumerate(('full','gradient','curl','harmonic')):
                fields=[channels[id][ch] for id in ids]
                mean=[sum(weight*v[i] for weight,v in zip(w,fields)) for i in range(4)]
                departures=[[x-m for x,m in zip(v,mean)] for v in fields];dep.append(departures)
                V=sum(weight*inner(d,d) for weight,d in zip(w,departures))
                S2=sum(w[a]*w[b]*inner(x,y)**2 for a,x in enumerate(departures) for b,y in enumerate(departures))
                summary=group['ensembles'][name];var.append(V)
                check(vec(summary['mean'])==mean and Q(summary['variance'])==V and Q(summary['second_contraction'])==S2)
                check(summary['contraction_space']==('samples' if len(ids)<4 else 'source_coordinates'))
                check((Q(summary['effective_dimension_2'])==V**2/S2) if V else summary['effective_dimension_2'] is None)
            check(var[0]==sum(var[1:]))
            for pair in group['sector_coupling']:
                a=('full','gradient','curl','harmonic').index(pair['left']);b=('full','gradient','curl','harmonic').index(pair['right'])
                if var[a] and var[b]:
                    numerator=sum(w[i]*w[j]*inner(x,y)*inner(dep[b][i],dep[b][j]) for i,x in enumerate(dep[a]) for j,y in enumerate(dep[a]))
                    check(Q(pair['value'])==numerator/(var[a]*var[b]))
                else:check(pair['value'] is None and pair['status']=='zero_sector_variance')
        check(result['execution']['projector_preparations']==2 and result['execution']['completed_tiles']==len(request['fields']))
    baseline=execute(q);verify(baseline)
    serial=copy.deepcopy(q);serial['workers']=1;v=execute(serial);verify(v)
    a=copy.deepcopy(baseline['result']);b=copy.deepcopy(v['result']);a.pop('execution');b.pop('execution');check(a==b)
    check(baseline['result']['comparisons'][0]['shared_fields']==['all.B'])
    many=copy.deepcopy(q);many['fields']=['all.A','all.B','all.C','all.zero','all.same_A'];many['weights']=['1/5']*5
    verify(execute(many))
    uneven=copy.deepcopy(q);uneven['weights']=['1/2','1/3','1/6'];uneven['interpretation']['kind']='declared_scenarios';verify(execute(uneven))
    singleton=demo.query(['A'],workers=1);singleton['interpretation']['evidence']=q['interpretation']['evidence'];verify(execute(singleton))
    constant=demo.query(['A','same_A']);constant['interpretation']['evidence']=q['interpretation']['evidence'];verify(execute(constant))
    zero=demo.query(['zero']);zero['interpretation']['evidence']=q['interpretation']['evidence'];verify(execute(zero))
    for mutate in [
        lambda v:v.update(fields=['all.A','all.partial','all.C']),
        lambda v:v.update(fields=['all.A','all.A','all.C']),
        lambda v:v.update(fields=['all.A','all.missing','all.C']),
        lambda v:v.update(fields=['all.A','all.object_A','all.C']),
        lambda v:v.update(context='unknown'),lambda v:v.update(workers=0),
        lambda v:v.update(weights=['1/2','1/2','1/2']),lambda v:v.update(weights=['1','0','0']),
        lambda v:v.update(weights=['-1','1','1']),lambda v:v.update(weights=['1/0','1/3','1/3']),
        lambda v:v.update(weights=['1/2','1/2']),lambda v:v.update(green_parameters=['0']),
        lambda v:v.update(green_parameters=['-1']),lambda v:v.update(green_parameters=['1','2/2']),
        lambda v:v['interpretation'].update(kind='population_null'),lambda v:v['interpretation'].update(evidence=[]),
        lambda v:v['interpretation']['evidence'][0].update(text='invented'),
        lambda v:v.update(groups=[dict(id='all',fields=['all.A'])]),
        lambda v:v.update(groups=[dict(id='empty',fields=[])]),
        lambda v:v.update(groups=[dict(id='bad',fields=['all.zero'])]),
        lambda v:v.update(comparisons=[dict(id='bad',left='all',right='unknown')]),
        lambda v:v.update(alpha='1/20')]:
        bad=copy.deepcopy(q);mutate(bad);execute(bad,True)
    # Same topology, positive nonidentity metrics and original tensor carrier.
    modified=json.loads((root/'inputs/import.json').read_text())
    modified['assays'][0]['metrics']=[dict(kind='diagonal',values=['2','3','4','5']),dict(kind='diagonal',values=['1','2','3','4']),dict(kind='diagonal',values=['5'])]
    command('import-assays','--spec',demo.save(root/'inputs/weighted.json',modified),'--out',root/'weighted-workspace')
    weighted=copy.deepcopy(q);demo.bind_evidence(weighted,root/'weighted-workspace')
    verify(execute(weighted,workspace_=root/'weighted-workspace/workspace.json'))
    top=demo.query(['top_A','top_B','top_C']);top['interpretation']['evidence']=q['interpretation']['evidence']
    v=execute(top)['result'];check(v['grade']==2 and Q(v['groups'][0]['ensembles']['full']['variance'])==Q(2,3))
    check([v['frames'][s]['rank'] for s in ('gradient','curl','harmonic')]==[1,0,0])
    objects=demo.query(['object_A','object_B']);objects['interpretation']['evidence']=q['interpretation']['evidence']
    v=execute(objects)['result'];check(v['grade']==0 and [v['frames'][s]['rank'] for s in ('gradient','curl','harmonic')]==[0,1,3])
    # An open complex retains the same primary assertions and has a larger harmonic space.
    open_q=demo.authored(root/'open/inputs',filled=False,arity=7,field_values=readings)
    command('import-assays','--spec',root/'open/inputs/import.json','--out',root/'open/workspace')
    demo.bind_evidence(open_q,root/'open/workspace');verify(execute(open_q,workspace_=root/'open/workspace/workspace.json'),False)
    # A primary assertion source supports arbitrary selected field axes; no small
    # fixed application field count is used by this producer.
    values={f'field_{i}':[i,i+1,i+2,i+3] for i in range(137)}
    large_q=demo.authored(root/'many/inputs',field_values=values)
    large_q=demo.query(list(values),workers=4)
    command('import-assays','--spec',root/'many/inputs/import.json','--out',root/'many/workspace')
    demo.bind_evidence(large_q,root/'many/workspace')
    v=execute(large_q,workspace_=root/'many/workspace/workspace.json')['result'];check(len(v['samples'])==137 and v['execution']['projector_preparations']==2)
    wide_q=demo.authored(root/'wide/inputs',arity=4101)
    wide_q=demo.query(['object_A','object_B'],workers=2)
    command('import-assays','--spec',root/'wide/inputs/import.json','--out',root/'wide/workspace')
    demo.bind_evidence(wide_q,root/'wide/workspace')
    v=execute(wide_q,workspace_=root/'wide/workspace/workspace.json')['result']
    check(v['frames']['harmonic']['rank']==4100 and not v['frames']['harmonic']['basis_materialized'])
    check(all(len(f['harmonic'])==4101 and f['harmonic_kernel_exact'] for f in v['samples']))
    check(Q(v['groups'][0]['ensembles']['full']['variance'])==Q(4101,4))
    print(json.dumps(dict(status='PASS',checks=checks,refusals=refused)))


if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
