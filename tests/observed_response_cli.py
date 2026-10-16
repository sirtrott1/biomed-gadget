#!/usr/bin/env python3
"""Admission, exact model and source contract checks for observed responses."""
import copy,csv,importlib.util,json,shutil,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'examples/bioinformatics/observed-responses/replay.py'
spec=importlib.util.spec_from_file_location('observed_demo',p);demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)


def main(binary,root):
    root.mkdir(parents=True,exist_ok=True)
    for name in ('authored','weighted','recoded','invalid-exposure','invalid-trials','invalid-count','missing-support','grade-one'):
        path=root/name
        if path.exists():shutil.rmtree(path)
    case=root/'authored';q=demo.authored(case,1)
    original=demo.execute(binary,case,q)['result']['answers'][0]['certificate']['result'];checks=refusals=0
    def check(v):
        nonlocal checks
        assert v;checks+=1
    def execute(query,refuse=False,study=None):
        nonlocal refusals
        request=demo.save(root/'query.json',query)
        p=subprocess.run([str(binary),'query','--study',str(study or next((case/'workspace').rglob('study.json'))),'--query',str(request)],capture_output=True,text=True)
        if refuse:check(p.returncode!=0);refusals+=1;return
        if p.returncode:raise AssertionError(p.stderr)
        d=json.loads(p.stdout);return d.get('certificate',d)['result']
    check(original['execution']['model_carriers']==3);check(original['execution']['model_preparations']==9)
    check(original['execution']['completed_tiles']==30)
    check(original['hodge_dimensions']==dict(gradient=0,curl=5,harmonic=3))
    check(original['preparations'][0]['models'][1]['rank']==3)
    check(all(p['models'][1]['coefficient_nullity']==1 for p in original['preparations']))
    for state in original['states']:
        for f in state['fits']:
            check(all(Q(x)==0 for x in f['normal_equations']))
            for group in ('observed','represented','residual'):
                s=f[group]['sectors'];qs=f[group]['quadrances'];n=len(s['full'])
                check(all(sum(Q(s[k][i]) for k in ('gradient','curl','harmonic'))==Q(s['full'][i]) for i in range(n)))
                check(sum(Q(qs[k]) for k in ('gradient','curl','harmonic'))==Q(qs['full']))
        for c in state['nested_comparisons']:check(Q(c['gain'])>=0)
    changed=copy.deepcopy(q);changed['workers']=2;parallel=execute(changed)
    check(parallel['states']==original['states']);check(parallel['preparations']==original['preparations'])
    changed=copy.deepcopy(q);changed.pop('reference_law');changed['responses'][0]['denominator']='all.exposure_minutes';changed['responses'][0]['rate_units']='events per minute'
    scaled=execute(changed)
    original_events=[f for f in original['states'][0]['fits'] if f['response']=='events']
    for a,b in zip(original_events,[f for f in scaled['states'][0]['fits'] if f['response']=='events']):
        check(a['objective_represented']==b['objective_represented']);check(a['residual_quadrance']==b['residual_quadrance'])
        check(all(Q(v)/60==Q(w) for v,w in zip(a['represented']['sectors']['full'],b['represented']['sectors']['full'])))
    # No intercept and empty column models retain valid native zero spans.
    zero=copy.deepcopy(q);zero['model']=dict(intercept=False,terms=[]);zero.pop('dose_comparisons');empty=execute(zero)
    check(all(p['models'][1]['rank']==0 for p in empty['preparations']))
    check(all(Q(f['residual_quadrance'])>=0 for f in empty['states'][0]['fits']))
    # Native source metric declaration; no response rescaling or inferred weights.
    study=next((case/'workspace').rglob('study.json'));weighted=root/'weighted';shutil.copytree(study.parent,weighted)
    d=json.loads((weighted/'study.json').read_text());d['contexts'][0]['metrics'][0]=dict(kind='diagonal',values=list(map(str,range(1,9))))
    demo.save(weighted/'study.json',d);out=execute(q,study=weighted/'study.json')
    check(out['hodge_dimensions']==original['hodge_dimensions']);check(out['states'][0]['fits'][1]['objective_represented']!=original['states'][0]['fits'][1]['objective_represented'])
    variants=[]
    def bad(change):
        v=copy.deepcopy(q);change(v);variants.append(v)
    bad(lambda v:v.update(workers=0));bad(lambda v:v.update(extra=1));bad(lambda v:v.update(context='missing'))
    bad(lambda v:v.update(responses=[]));bad(lambda v:v['responses'].append(copy.deepcopy(v['responses'][0])))
    bad(lambda v:v['responses'][0].update(kind='poisson'));bad(lambda v:v['responses'][0].pop('denominator'))
    bad(lambda v:v['responses'][0].update(denominator='all.events'));bad(lambda v:v['responses'][0].update(denominator='all.group0'))
    bad(lambda v:v['responses'][1].update(rate_units='counts'));bad(lambda v:v['responses'][2].update(field='all.dose'))
    bad(lambda v:v['responses'][2].update(levels=v['responses'][2]['levels'][:1]))
    bad(lambda v:v['responses'][2]['levels'][1].update(code='40'));bad(lambda v:v['responses'][2]['levels'][1].update(label='low'))
    bad(lambda v:v['responses'][2]['levels'][2].update(code='999'));bad(lambda v:v['responses'][3].update(levels=[]))
    bad(lambda v:v['model']['terms'][0].update(factors=['all.ordinal']));bad(lambda v:v['model']['terms'][0].update(factors=['all.events']))
    bad(lambda v:v['model']['terms'][0].update(factors=[]));bad(lambda v:v['model']['terms'][1].update(id='dose'))
    bad(lambda v:v['model']['ablations'][0].update(remove=['missing']));bad(lambda v:v['model'].update(intercept='true'))
    bad(lambda v:v['dose_comparisons'][0]['groups'][0].update(indicator='all.dose'))
    bad(lambda v:v['dose_comparisons'][0]['groups'][0].update(indicator='all.group1'))
    bad(lambda v:v['dose_comparisons'][0]['groups'].pop());bad(lambda v:v['dose_comparisons'][0]['groups'].reverse())
    bad(lambda v:v['dose_comparisons'][0].update(direction='predicted'));bad(lambda v:v['reference_law'].update(observed_state='missing'))
    bad(lambda v:v['reference_law']['states'][1].update(mass='0'));bad(lambda v:v['reference_law']['states'][1].update(mass='2/3'))
    bad(lambda v:v['reference_law']['states'][1].update(id='recorded'));bad(lambda v:v['reference_law']['states'][0]['response_fields'].pop('state'))
    bad(lambda v:v['reference_law']['states'][0]['response_fields'].update(events='all.events_ref'))
    bad(lambda v:v['reference_law']['states'][1]['response_fields'].update(events='all.dose'))
    bad(lambda v:v['reference_law'].update(alpha='0'));bad(lambda v:v.update(green_parameter='0'))
    bad(lambda v:v['interpretation'].update(evidence=[dict(document='missing',start=0,end=1,text='x')]))
    for v in variants:execute(v,True)
    # Categorical relabeling is an input declaration, not a numeric distance change.
    recoded=root/'recoded';qq=demo.authored(recoded,1);qq.pop('reference_law')
    path=recoded/'inputs/observations.tsv'
    rows=list(csv.reader(path.open(),delimiter='\t'));conversion={'40':'-7','10':'9000','70':'2'}
    for row in rows[1:]:
        for field in ('ordinal','ordinal_ref'):row[rows[0].index(field)]=conversion[row[rows[0].index(field)]]
    path.write_text(''.join('\t'.join(row)+'\n' for row in rows))
    spec=json.loads((recoded/'inputs/import.json').read_text());spec['assays'][0]['observations']['sha256']=demo.digest(path);demo.save(recoded/'inputs/import.json',spec)
    for l,code in zip(qq['responses'][2]['levels'],('-7','9000','2')):l['code']=code
    rr=demo.execute(binary,recoded,qq)['result']['answers'][0]['certificate']['result']
    for a,b in zip([f for f in original['states'][0]['fits'] if f['response']=='state'],[f for f in rr['states'][0]['fits'] if f['response']=='state']):
        check(a['objective_input']==b['objective_input']);check(a['objective_represented']==b['objective_represented'])
    # Actual observed data admission, rather than only malformed query syntax.
    for name,table,column,value in [('invalid-exposure','observations.tsv','exposure','0'),
                                    ('invalid-trials','measurements.tsv','trials','3/2'),
                                    ('invalid-count','measurements.tsv','events','-1'),
                                    ('missing-support','observations.tsv','exposure','NA')]:
        target=root/name;v=demo.authored(target,1);path=target/'inputs'/table
        rows=list(csv.reader(path.open(),delimiter='\t'))
        if table=='observations.tsv':rows[1][rows[0].index(column)]=value
        else:next(row for row in rows[1:] if row[0]==column)[1]=value
        path.write_text(''.join('\t'.join(row)+'\n' for row in rows))
        spec=json.loads((target/'inputs/import.json').read_text());spec['assays'][0]['observations' if table=='observations.tsv' else 'source']['sha256']=demo.digest(path);demo.save(target/'inputs/import.json',spec)
        demo.run(binary,'import-assays','--spec',target/'inputs/import.json','--out',target/'workspace')
        demo.source.bind_evidence(v,target/'workspace')
        for d in v['responses']:d['evidence']=v['interpretation']['evidence']
        v.pop('dose_comparisons');v.pop('reference_law')
        execute(v,True,next((target/'workspace').rglob('study.json')))
    # Count fields on grade one primary relations retain gradient, curl and
    # harmonic information, including the original higher filling.
    target=root/'grade-one';demo.authored(target,1)
    path=target/'inputs/relations.json';data=json.loads(path.read_text());ids=[d['id'] for d in data['relations'] if d['grade']==1]
    data['fields']=[dict(id=name,grade=1,channel='observation',quantity=quantity,units=units,
                         method='Authored primary-relation measurements',values=[dict(relation=rel,value=str(value)) for rel,value in zip(ids,values)])
                    for name,values,quantity,units in [('relationcounts',range(1,9),'event_count','counts'),('relationexposure',[1]*8,'exposure','days'),('relationfactor',list(range(4))*2,'declared_factor','dimensionless')]]
    demo.save(path,data);spec=json.loads((target/'inputs/import.json').read_text());spec['assays'][0]['relations']['sha256']=demo.digest(path);demo.save(target/'inputs/import.json',spec)
    qq=demo.query([demo.declaration('relationcounts','count','relationcounts',denominator='all.relationexposure',rate_units='events per day')],
                  [dict(id='factor',factors=['all.relationfactor'])],[],2)
    output=demo.execute(binary,target,qq)['result']['answers'][0]['certificate']['result']
    check(output['grade']==1);check(output['hodge_dimensions']==dict(gradient=5,curl=1,harmonic=2))
    check(list(map(Q,output['states'][0]['fits'][1]['observed']['sectors']['curl']))==[0,0,0,0,-1,-1,1,1])
    result=dict(status='PASS',checks=checks,refusals=refusals);demo.save(root/'contracts.json',result);print(json.dumps(result))


if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
