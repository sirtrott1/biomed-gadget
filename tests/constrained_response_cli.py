#!/usr/bin/env python3
"""Independent physical optimum and source admission contracts."""
import copy,csv,importlib.util,json,shutil,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path

examples=Path(__file__).resolve().parents[1]/'examples/bioinformatics/constrained-responses'
def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
demo=load('physical_demo',examples/'replay.py');verify=load('physical_verify',examples/'verify.py')

def main(binary,root):
    shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    checks=refusals=0
    def check(value):
        nonlocal checks
        assert value;checks+=1
    case=root/'authored';q=demo.constrained(demo.source.authored(case,1));out=demo.execute(binary,case,q)['result']['answers'][0]['certificate']['result']
    checks+=verify.verify(case)['checks']
    study=next((case/'workspace').rglob('study.json'))
    def execute(query,fail=False,owner=study):
        nonlocal refusals
        p=subprocess.run([str(binary),'query','--study',str(owner),'--query',str(demo.save(root/'query.json',query))],capture_output=True,text=True)
        if fail:check(p.returncode!=0);refusals+=1;return
        if p.returncode:raise AssertionError(p.stderr)
        certificate=json.loads(p.stdout);certificate=certificate.get('certificate',certificate);demo.save(root/'last-direct-certificate.json',certificate);return certificate['result']
    # Source ordinal violations are fixed jointly, without coordinate clipping.
    check(any(not c['admissible'] for c in out['states'][0]['ordinal_level_readings']))
    check(all(Q(x)>=0 for f in out['physical_shape']['results'] if 'ordinal_level_fields' in f for row in f['ordinal_level_fields'] for x in row))
    changed=copy.deepcopy(q);changed['workers']=2;parallel=execute(changed)
    check(parallel['physical_shape']['results']==out['physical_shape']['results'])
    check(parallel['physical_shape']['nested_comparisons']==out['physical_shape']['nested_comparisons'])
    weighted=root/'weighted';shutil.copytree(case,weighted)
    weighted_study=next((weighted/'workspace').rglob('study.json'));metadata=json.loads(weighted_study.read_text())
    metadata['contexts'][0]['metrics'][0]=dict(kind='diagonal',values=[str(i) for i in range(1,9)])
    demo.save(weighted_study,metadata);weighted_result=execute(q,owner=weighted_study)
    certificate=json.loads((root/'last-direct-certificate.json').read_text());demo.save(weighted/'certificate.json',certificate)
    checks+=verify.verify(weighted)['checks']
    check(weighted_result['physical_shape']['results']!=out['physical_shape']['results'])
    changed=copy.deepcopy(q);changed['physical_shape']['dose_constraints']=[dict(response='events',comparison='measured-dose')]
    dose=execute(changed);check(dose['physical_shape']['completed_tiles']==len(dose['physical_shape']['results']))
    for change in [lambda v:v.pop('physical_shape'),lambda v:v['physical_shape'].update(extra=True),
                   lambda v:v['physical_shape']['dose_constraints'].append(dict(response='missing',comparison='measured-dose')),
                   lambda v:v['physical_shape']['dose_constraints'].append(dict(response='events',comparison='missing')),
                   lambda v:v['physical_shape']['dose_constraints'].extend([dict(response='events',comparison='measured-dose')]*2),
                   lambda v:v['dose_comparisons'][0].update(direction='unrestricted'),
                   lambda v:v['physical_shape'].update(evidence=[dict(document='missing',start=0,end=1,text='x')])]:
        bad=copy.deepcopy(changed);change(bad);execute(bad,True)
    # Analytic count example: Q=sum_i(y_i beta0 beta1*d_i)^2, beta0>=0.
    analytic=root/'analytic';qq=demo.source.authored(analytic,1);qq.pop('reference_law');qq.pop('dose_comparisons');qq.pop('green_parameter')
    qq['responses']=qq['responses'][:1];qq['model']=dict(intercept=True,terms=[dict(id='dose',factors=['all.dose']),dict(id='alias',factors=['all.dose_alias'])])
    table=analytic/'inputs/measurements.tsv';rows=list(csv.reader(table.open(),delimiter='\t'))
    next(row for row in rows[1:] if row[0]=='events')[1:]=['0']*6+['100']*2;table.write_text(''.join('\t'.join(row)+'\n' for row in rows))
    metadata=analytic/'inputs/observations.tsv';rows=list(csv.reader(metadata.open(),delimiter='\t'))
    for i,row in enumerate(rows[1:]):
        row[rows[0].index('exposure')]='1';row[rows[0].index('exposure_minutes')]='60'
        row[rows[0].index('dose')]=str(i//2);row[rows[0].index('dose_alias')]=str(i//2)
    metadata.write_text(''.join('\t'.join(row)+'\n' for row in rows))
    spec=json.loads((analytic/'inputs/import.json').read_text());spec['assays'][0]['source']['sha256']=demo.digest(table);spec['assays'][0]['observations']['sha256']=demo.digest(metadata);demo.save(analytic/'inputs/import.json',spec)
    result=demo.execute(binary,analytic,demo.constrained(qq))['result']['answers'][0]['certificate']['result'];checks+=verify.verify(analytic)['checks']
    fitted=next(f for f in result['physical_shape']['results'] if f['model']=='full')
    check(list(map(Q,fitted['fits'][0]['objective_represented']))==[Q(150,7)*d for d in range(4) for _ in range(2)])
    check(Q(fitted['joint_residual_quadrance'])==Q(50000,7));check(Q(fitted['constraint_cost'])==Q(8000,7))
    check(fitted['fits'][0]['particular_coefficients']==['0','150/7','0'])
    scaled=copy.deepcopy(qq);scaled['responses'][0]['denominator']='all.exposure_minutes';scaled['responses'][0]['rate_units']='events per minute'
    scaled=demo.constrained(scaled)
    scaled['physical_shape']['evidence']=qq['interpretation']['evidence']
    scaled=execute(scaled,owner=next((analytic/'workspace').rglob('study.json')))
    other=next(f for f in scaled['physical_shape']['results'] if f['model']=='full')
    check(other['fits'][0]['objective_represented']==fitted['fits'][0]['objective_represented'])
    check(other['joint_residual_quadrance']==fitted['joint_residual_quadrance'])
    # Empty spans remain feasible; no positive rank special case assumptions.
    zero=copy.deepcopy(q);zero['model']=dict(intercept=False,terms=[]);zero.pop('reference_law');zero.pop('dose_comparisons')
    empty=execute(zero);check(all(all(Q(x)==0 for part in f['fits'] for x in part['objective_represented']) for f in empty['physical_shape']['results']))
    grade_case=root/'grade-one';demo.source.authored(grade_case,1)
    relations=grade_case/'inputs/relations.json';data=json.loads(relations.read_text());ids=[d['id'] for d in data['relations'] if d['grade']==1]
    data['fields']=[dict(id=name,grade=1,channel='observation',quantity=quantity,units=units,
        method='Authored primary relation measurement',values=[dict(relation=rel,value=str(value)) for rel,value in zip(ids,values)])
        for name,values,quantity,units in [('relationcounts',range(1,9),'event_count','counts'),('relationexposure',[1]*8,'exposure','days'),('relationfactor',list(range(4))*2,'declared_factor','dimensionless')]]
    demo.save(relations,data);spec=json.loads((grade_case/'inputs/import.json').read_text());spec['assays'][0]['relations']['sha256']=demo.digest(relations);demo.save(grade_case/'inputs/import.json',spec)
    grade_query=demo.source.query([demo.source.declaration('relationcounts','count','relationcounts',denominator='all.relationexposure',rate_units='events per day')],
                                 [dict(id='factor',factors=['all.relationfactor'])],[],1)
    grade_result=demo.execute(binary,grade_case,demo.constrained(grade_query))['result']['answers'][0]['certificate']['result']
    check(grade_result['grade']==1);check(grade_result['hodge_dimensions']==dict(gradient=5,curl=1,harmonic=2))
    for fit in grade_result['physical_shape']['results']:
        check(all(Q(x)>=0 for x in fit['fits'][0]['objective_represented']))
        sectors=fit['fits'][0]['represented']['sectors'];check(all(sum(Q(sectors[name][i]) for name in ('gradient','curl','harmonic'))==Q(sectors['full'][i]) for i in range(8)))
    result=dict(status='PASS',checks=checks,refusals=refusals);demo.save(root/'contracts.json',result);print(json.dumps(result))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
