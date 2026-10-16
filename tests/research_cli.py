"""Independent exact contracts for reusable endpoint and trajectory queries."""
import copy
import hashlib
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def rref(a):
    a=[list(map(Q,r)) for r in a]; pivots=[];row=0
    for col in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row];v=a[row][col];a[row]=[x/v for x in a[row]]
        for i in range(len(a)):
            if i!=row:
                v=a[i][col];a[i]=[x-v*y for x,y in zip(a[i],a[row])]
        pivots.append(col);row+=1
        if row==len(a):break
    return a,pivots


def solve(a,b):
    reduced,pivots=rref([list(row)+[value] for row,value in zip(a,b)])
    n=len(a);assert pivots==list(range(n))
    return [row[-1] for row in reduced]


def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def mv(a,x):return [dot(r,x) for r in a]
def columns(a):return [list(x) for x in zip(*a)]
def projection(cols,x):
    if not cols:return [Q(0)]*len(x)
    frame=columns(cols);_,pivots=rref(frame);basis=[cols[i] for i in pivots]
    coeff=solve([[dot(a,b) for b in basis] for a in basis],[dot(a,x) for a in basis])
    return [sum((a[i]*c for a,c in zip(basis,coeff)),Q(0)) for i in range(len(x))]


def matrix(record):return [list(map(Q,row)) for row in record['values']]


def check_trajectory(result):
    for interval in result['intervals']:
        b1=matrix(interval['B1']);b2=matrix(interval['B2']);n=len(interval['state_indices'])
        assert interval['H1_rank']==n-len(rref(b1)[1])-len(rref(b2)[1])
        l=[[sum((b1[k][i]*b1[k][j] for k in range(len(b1))),Q(0))+
            sum((b2[i][k]*b2[j][k] for k in range(interval['B2']['columns'])),Q(0))
            for j in range(n)] for i in range(n)]
        for field in interval['fields']:
            if not field.get('complete_support',True):
                assert 'harmonic' not in field;continue
            x=list(map(Q,field['input']));g=list(map(Q,field['gradient']));c=list(map(Q,field['curl']));h=list(map(Q,field['harmonic']));green=list(map(Q,field['green']))
            assert [a+b+d for a,b,d in zip(g,c,h)]==x
            assert dot(g,c)==dot(g,h)==dot(c,h)==0
            expected_g=projection(b1,x);expected_c=projection(columns(b2),x)
            assert g==expected_g and c==expected_c
            assert all(v==0 for v in mv(l,h))
            assert [a+b for a,b in zip(green,mv(l,green))]==x
            if 'controlled' in field:
                adjusted=field['controlled'];patient_y=list(map(Q,adjusted['patient_input']));residual_y=list(map(Q,adjusted['patient_residual']))
                assert residual_y==[a-b for a,b in zip(patient_y,projection(columns(matrix(result['control_frame'])),patient_y))]
                local={state:i for i,state in enumerate(interval['state_indices'])}
                assignments=[event['state'] for event in result['events'] if event['interval']==interval['interval']]
                state_y=[Q(0)]*n
                for state,value in zip(assignments,residual_y):state_y[local[state]]+=value
                assert state_y==list(map(Q,adjusted['input']))
                hy=[a-b-c for a,b,c in zip(state_y,projection(b1,state_y),projection(columns(b2),state_y))]
                assert hy==list(map(Q,adjusted['harmonic']))
                assert [hy[local[state]] for state in assignments]==list(map(Q,adjusted['patient_harmonic_evaluation']))
                stable=columns(matrix(interval['repeated_support_harmonic_frame']))
                repeated_y=projection(stable,hy)
                assert repeated_y==list(map(Q,adjusted['repeated_support_harmonic']))
                for feature in adjusted['fields']:
                    rx=list(map(Q,feature['patient_residual']));assert all(dot(col,rx)==0 for col in columns(matrix(result['control_frame'])))
                    state_x=[Q(0)]*n
                    for state,value in zip(assignments,rx):state_x[local[state]]+=value
                    hx=[a-b-c for a,b,c in zip(state_x,projection(b1,state_x),projection(columns(b2),state_x))]
                    assert hx==list(map(Q,feature['harmonic'])) and dot(hx,hy)==Q(feature['pairing'])
                    assert dot(projection(stable,hx),repeated_y)==Q(feature['repeated_support_pairing'])
    for comparison in result['temporal_comparisons']:
        ha=matrix(comparison['left_harmonic_images']);hb=matrix(comparison['right_harmonic_images'])
        a=len(rref(ha)[1]);b=len(rref(hb)[1]);joint=len(rref(ha+hb)[1])
        assert comparison['common_quotient_dimension']==a+b-joint
        for field in comparison['fields']:
            terms=[list(map(Q,field[k])) for k in ('realization_term','innovation_term','operator_term')]
            assert [sum(v) for v in zip(*terms)]==list(map(Q,field['total']))


def run(binary,root):
    if root.exists():shutil.rmtree(root)
    root.mkdir(parents=True);cases=[]
    first=[(0,0,2),(0,2,0),(0,2,2),(2,0,0),(2,0,2),(2,2,0)]
    first=[p for p,n in zip(first,(3,3,3,3,2,2)) for _ in range(n)]
    second=[(a,b,c) for a in (0,2) for b in (0,2) for c in (0,2) for _ in range(2)]
    header=['patient','group','size','size_copy','y','partial_y','x','z','missing_x']+[f'{f}{t}' for f in 'ABC' for t in range(4)]
    rows=[]
    for i in range(16):
        size=i+1;y=i%2;x=17+3*size+2*y;z=17+3*size
        values=[f'person-{i}','batch-X',str(size),str(size),str(y),'NA' if i==0 else str(y),str(x),str(z),'NA' if i==1 else str(x)]
        for f in range(3):
            a=100+i;b=a+(-1 if first[i][f]==0 else 1);c=b+(-1 if second[i][f]==0 else 1);d=c+(-1 if first[i][f]==0 else 1);values+=list(map(str,(a,b,c,d)))
        rows.append(values)
    raw='\n'.join(['\t'.join(header)]+['\t'.join(r) for r in rows])+'\n';(root/'table.txt').write_text(raw)
    dictionary='Authored nonclinical exact fixture: measured units and categorical identities.\n';(root/'dictionary.txt').write_text(dictionary)
    spec=dict(schema='gadget.table.import.v1',id='general-research-contract',source_namespace='871003',source=dict(path='table.txt',sha256='sha256:'+hashlib.sha256(raw.encode()).hexdigest(),delimiter='\t',header=True),dictionary=dict(path='dictionary.txt',sha256='sha256:'+hashlib.sha256(dictionary.encode()).hexdigest()),columns=[dict(id=name,channel='observation',quantity=name,units='synthetic_measurement' if name not in ('y','partial_y','patient','group') else 'dimensionless',value_type='categorical' if name in ('patient','group') else 'indicator' if name in ('y','partial_y') else 'quantitative',dictionary_start=0,dictionary_end=len(dictionary)) for name in header],identity_column='patient',contrasts=[dict(id='A_change',left='A1',right='A0')],missing_tokens=['NA'],relation_columns=['group'],question='General exact research contracts',population='sixteen authored nonclinical records',selection='all source records',limitations=['synthetic fixture'],source_url='synthetic://general-research',attribution='GADGET fixture',license='synthetic')
    (root/'spec.json').write_text(json.dumps(spec))
    def command(*args,success=True):
        p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True,timeout=90)
        assert (p.returncode==0)==success,(args,p.stderr)
        return json.loads(p.stdout) if success else None
    command('import-table','--spec',root/'spec.json','--out',root/'study')
    def query(value,success=True):
        path=root/'query.json';path.write_text(json.dumps(value));return command('query','--study',root/'study/study.json','--query',path,success=success)
    controlled=dict(operation='controlled_pairings',fields=['all.x','all.z'],endpoint='all.y',controls=['all.size','all.size_copy'],selection=[],support='joint_reported',metric='identity',intercept=True)
    receipt=query(controlled);r=receipt['result'];assert r['n']==16 and r['control_rank']==2
    assert r['entities']==[f'person-{i}' for i in range(16)]
    assert Q(r['fields'][0]['endpoint_coefficient'])==2 and Q(r['fields'][1]['pairing'])==0
    ry=list(map(Q,r['endpoint_residual']));raw_y=[Q(i%2) for i in range(16)]
    expected_y=projection([[Q(1)]*16,list(map(Q,range(1,17)))],raw_y)
    assert ry==[a-b for a,b in zip(raw_y,expected_y)]
    for f in r['fields']:
        x=list(map(Q,f['input']));p=projection([[Q(1)]*16,list(map(Q,range(1,17)))],x)
        residual=[a-b for a,b in zip(x,p)];assert residual==list(map(Q,f['residual'])) and dot(residual,ry)==Q(f['pairing'])
    parallel=query({**controlled,'workers':2})['result']
    assert {k:v for k,v in parallel.items() if k!='execution'}=={k:v for k,v in r.items() if k!='execution'}
    assert parallel['execution']['owner']=='rfc::numeric::PreparedExactSolve'
    assert parallel['execution']['gram_preparations']==1 and parallel['execution']['applications']==3
    assert parallel['execution']['worker_slots']==2 and parallel['execution']['completed_tiles']==2
    reverse=query({**controlled,'workers':2,'gram_fields':['all.z','all.x']})['result']
    assert reverse['fields']==r['fields'] and reverse['gram_fields']==['all.z','all.x']
    assert matrix(reverse['residual_gram'])==[list(reversed(row)) for row in reversed(matrix(r['residual_gram']))]
    empty=query({**controlled,'gram_fields':[]})['result']
    assert empty['fields']==r['fields'] and empty['residual_gram']==dict(rows=0,columns=0,values=[])
    partial={**controlled,'fields':['all.missing_x'],'endpoint':'all.partial_y'};assert query(partial)['result']['n']==14
    constant={**controlled,'selection':['all.y']};c=query(constant)['result'];assert not c['endpoint_identified'] and c['fields'][0]['endpoint_coefficient'] is None
    for edit in ({'controls':['all.y']},{'fields':['all.group']},{'selection':['all.size']},{'metric':'invented'},{'support':'impute'},{'extra':True},{'workers':0},{'gram_fields':['all.y']},{'gram_fields':['all.x','all.x']}):query({**controlled,**edit},False);cases.append(edit)
    trajectory=dict(operation='trajectory_structure',series=[dict(id=f,fields=['all.'+f+str(t) for t in range(4)]) for f in 'ABC'],selection=[],controls=['all.size','all.size_copy'],support='joint_reported',time=dict(axis='ordered_assay_visits',unit='ordinal_visit',positions=['0','1','2','3']),faces='coordinate_exchange_identities',semantics='observed_feature_change_relations',**{'lambda':'1'},minimum_witnesses=2,endpoint_fields=['all.y','all.partial_y'])
    tr=query(trajectory);check_trajectory(tr['result']);ints=tr['result']['intervals'];assert ints[0]['H1_rank']==2 and ints[1]['H1_rank']==0
    comp=tr['result']['temporal_comparisons'][0];assert comp['left_image_rank']==2 and comp['right_image_rank']==0 and comp['common_quotient_dimension']==0
    assert comp['boundary_compatible'] and not tr['result']['temporal_comparisons'][1]['boundary_compatible'];assert ints[0]['repeated_support_H1_image_rank']==2
    bad_series=copy.deepcopy(trajectory['series']);bad_series[0]['fields'][1]='all.y'
    for edit in ({'lambda':'0'},{'minimum_witnesses':1},{'faces':'infer_biology'},{'series':bad_series},{'time':dict(axis='visit',unit='ordinal_visit',positions=['0','2','1','3'])}):query({**trajectory,**edit},False);cases.append(edit)
    replay=query(trajectory);assert replay==tr
    if len(sys.argv)>3:
        requests=root/'trajectory-requests.json';requests.write_text(json.dumps({'complete':trajectory,'controlled':{**trajectory,'controls':['all.size','all.size_copy']}}))
        subprocess.run([sys.argv[3],str(root/'study/study.json'),str(requests),str(root/'native-trajectory')],check=True)
    coverage=query(dict(operation='coverage',field='all.A_change'))['result'];assert all(Q(v['value']['rational']) in (Q(-1),Q(1)) for v in coverage['observations'])
    # Batch extension must reuse immutable certificates and refuse changed
    # sources or already completed queries without damaging prior receipts.
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'tools'))
    from scientific_batch import execute
    batch_protocol=dict(schema='gadget.scientific.batch.v1',id='fixture-batch',question='Exact nonclinical field checks',jobs=[dict(id='controlled',query=controlled,description='Measured fixture')])
    protocol_file=root/'protocol.json';protocol_file.write_text(json.dumps(batch_protocol));batch=root/'batch'
    execute(binary,root/'study/study.json',protocol_file,batch)
    batch_protocol['jobs'].append(dict(id='coverage',query=dict(operation='coverage',field='all.A_change'),description='Joint contrast support'))
    protocol_file.write_text(json.dumps(batch_protocol));assert len(execute(binary,root/'study/study.json',protocol_file,batch,resume=True)['jobs'])==2
    # Display files do not own completion or the original protocol. Resume
    # directly from the checkpoint and source bound native plan publications.
    for name in ('manifest.json','protocol.json','catalog.json','job-0.json','job-1.json'):
        (batch/name).unlink()
    batch_protocol['jobs'].append(dict(id='inspection',query=dict(operation='inspect'),description='Native source binding'))
    protocol_file.write_text(json.dumps(batch_protocol))
    assert len(execute(binary,root/'study/study.json',protocol_file,batch,resume=True)['jobs'])==3
    batch_many={**batch_protocol,'jobs':[dict(id='long_identity_'+'x'*129 if i==0 else 'job_'+str(i),query=dict(operation='inspect'),description='Whole source') for i in range(17)]}
    many_file=root/'many-protocol.json';many_file.write_text(json.dumps(batch_many))
    assert len(execute(binary,root/'study/study.json',many_file,root/'many-batch',plan_steps=17)['jobs'])==17
    previous=(batch/'manifest.json').read_bytes()
    changed_protocol=copy.deepcopy(batch_protocol);changed_protocol['jobs'][0]['query']['controls']=[];protocol_file.write_text(json.dumps(changed_protocol))
    try:execute(binary,root/'study/study.json',protocol_file,batch,resume=True)
    except ValueError:pass
    else:raise AssertionError('Modified completed query resumed')
    assert (batch/'manifest.json').read_bytes()==previous
    protocol_file.write_text(json.dumps(batch_protocol));study_path=root/'study/study.json';original_study=study_path.read_bytes();changed_study=json.loads(original_study);changed_study['design']['question']='Different source declaration';study_path.write_text(json.dumps(changed_study))
    try:execute(binary,study_path,protocol_file,batch,resume=True)
    except ValueError:pass
    else:raise AssertionError('Modified admitted study resumed')
    assert (batch/'manifest.json').read_bytes()==previous;study_path.write_bytes(original_study)
    (root/'controlled.json').write_text(json.dumps(receipt,indent=2)+'\n');(root/'trajectory.json').write_text(json.dumps(tr,indent=2)+'\n')
    from research_scale import run as qualify_scale
    scale=qualify_scale(binary,root/"large-carrier")
    print(json.dumps({'status':'PASS','refusals':len(cases),'conditional_contrast':'2','spurious_size_pairing_removed':True,'persistent_H1_filled':True,'defective_transport_detected':True,'exact_replay':True,'large_carrier_comparisons':scale['comparisons']}))


if __name__=='__main__':run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
