"""Independent paired-count, gluing, solve and dependent workflow contracts."""
import copy, hashlib, json, shutil, subprocess, sys
from fractions import Fraction as F
from pathlib import Path


def main():
    binary,root=sys.argv[1],Path(sys.argv[2]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    rows=[['gA','gene','4','12','3','7','-1','0.01'],['gB','gene','0','0','3','3','0','NA'],['gC','gene','5','1','1','9','1','0.1'],['gD','gene','2','2','1','NA','0','NA']]
    header=['gene','kind','a0','a1','b0','b1','effect','padj'];raw='\n'.join(['\t'.join(header)]+['\t'.join(r) for r in rows])+'\n'
    dictionary='Authored count fixture with two explicitly matched units and reported zero.\n'
    (root/'table.tsv').write_text(raw);(root/'dictionary.txt').write_text(dictionary)
    spec=dict(schema='gadget.table.import.v1',id='paired-contract',source_namespace='890101',source=dict(path='table.tsv',sha256='sha256:'+hashlib.sha256(raw.encode()).hexdigest(),delimiter='\t',header=True),dictionary=dict(path='dictionary.txt',sha256='sha256:'+hashlib.sha256(dictionary.encode()).hexdigest()),columns=[dict(id=s,channel='E' if i>1 else 'observation',quantity='read_count' if i>1 else s,units='counts' if i>1 else 'dimensionless',value_type='quantitative' if i>1 else 'categorical',dictionary_start=0,dictionary_end=len(dictionary)) for i,s in enumerate(header)],identity_column='gene',missing_tokens=['NA'],relation_columns=['kind'],question='Exact matched field contracts',population='four authored features',selection='all records',limitations=['authored nonclinical data'],source_url='synthetic://matched-fields',attribution='GADGET contracts',license='synthetic')
    for i,quantity in [(6,'model_log2_fold_change'),(7,'adjusted_p_value')]:
        spec['columns'][i].update(channel='observation',quantity=quantity,units='dimensionless')
    def command(*args,success=True):
        p=subprocess.run([binary,*map(str,args)],capture_output=True,text=True,timeout=60)
        assert (p.returncode==0)==success,p.stderr
        return json.loads(p.stdout) if success else None
    (root/'spec.json').write_text(json.dumps(spec));command('import-table','--spec',root/'spec.json','--out',root/'study')
    declaration=json.loads((root/'study/study.json').read_text());study=root/'study/study.json'
    def query(q,success=True):
        (root/'query.json').write_text(json.dumps(q));return command('query','--study',study,'--query',root/'query.json',success=success)
    evidence=[dict(document='dictionary',start=0,end=len(dictionary),text=dictionary)]
    paired=dict(operation='paired_field_contrasts',pairs=[dict(id='A',left='all.a1',right='all.a0',left_scale='1',right_scale='1'),dict(id='B',left='all.b1',right='all.b0',left_scale='1',right_scale='1')],selection=[],support='joint_reported',normalization=dict(method='identity',evidence=evidence),output_units='counts',workers=1)
    serial=query(paired);paired['workers']=4;parallel=query(paired)
    assert serial['result']['features']==parallel['result']['features']
    assert parallel['result']['reported']==3 and parallel['result']['missing']==1
    for row,record in zip(rows,parallel['result']['features']):
        d=[F(row[3])-F(row[2]),F(row[5])-F(row[4])];mean=sum(d)/2;res=[v-mean for v in d]
        assert list(map(F,record['delta']))==d and F(record['mean_difference'])==mean
        assert list(map(F,record['residual']))==res and F(record['heterogeneity_quadrance'])==sum(v*v for v in res)
    assert parallel['result']['features'][1]['zero_pairs']==2
    compared=copy.deepcopy(paired);compared['benchmark']=dict(effect='all.effect',adjusted_p='all.padj',alpha='0.05')
    comparison=query(compared)['result'];assert comparison['reported']==3
    assert comparison['features'][0]['benchmark']['significant'] and not comparison['features'][0]['benchmark']['direction_agrees']
    assert comparison['features'][1]['benchmark']['adjusted_p'] is None
    bad=copy.deepcopy(compared);bad['benchmark']['adjusted_p']='all.effect';query(bad,False)
    broken=copy.deepcopy(paired);broken['pairs'][0]['left_scale']='1/0';query(broken,False)
    broken=copy.deepcopy(paired);broken['pairs'][1]['right']='all.a0';query(broken,False)
    broken=copy.deepcopy(paired);broken['normalization']['evidence'][0]['text']='invented';query(broken,False)
    cells=declaration['cells'];occurrences=declaration['occurrences'];lower=[c for c in cells if c['grade']==0];upper=[c for c in cells if c['grade']==1]
    # Exact differences on a four-way relation. Every flag map is explicit.
    maps=[];coeffs=[[4,0,0],[-4,4,0],[0,-4,4],[0,0,-4]]
    for o,c in zip(occurrences,coeffs):maps.append(dict(grade=1,key=o['key'],matrix=[[str(v)] for v in c]))
    local=dict(operation='local_sections',context='all',construction='cellular_sheaf',stalks=[dict(grade=0,key=c['key'],axes=['all.a0']) for c in lower]+[dict(grade=1,key=upper[0]['key'],axes=['difference.0','difference.1','difference.2'])],maps=maps,sections=[dict(id='counts',grade=0,fields=['all.a0'],quantity='read_count',units='counts')],metric='identity')
    section=query(local)['result'];assert section['cohomology'][0]['dimension']==1
    v=[F(r[2]) for r in rows];mean=sum(v)/4
    assert list(map(F,section['sections'][0]['glued']))==[mean]*4
    assert section['incompatible_sections']==1 and section['sections'][0]['input_origin']=='admitted source fields'
    bad=copy.deepcopy(local);bad['maps'].pop();query(bad,False)
    bad=copy.deepcopy(local);bad['stalks'][0]['axes']=['unrelated_axis'];query(bad,False)
    bad=copy.deepcopy(local);bad['sections'][0]['fields']=['all.b1'];query(bad,False)
    fibers=copy.deepcopy(local);fibers['construction']='boundary_fibers';fibers['maps']=[dict(grade=1,key=o['key'],matrix=[['1','0','0']]) for o in occurrences];fibers['sections']=[dict(id='source',grade=1,values=['2','3','4'],quantity='declared_stalk_amplitude',units='counts')]
    f=query(fibers)['result'];assert len(f['sections'][0]['relations'][0]['legs'])==4
    assert all(leg['value']==['2'] for leg in f['sections'][0]['relations'][0]['legs'])
    green=dict(operation='green_response',context='all',field='all.a0',**{'lambda':'1'},interior=[0,2],coordinates=dict(id='890999',cells=[c['key'] for c in lower]))
    y=[F(r['value']['rational']) for r in query(green)['result']['field']]
    # B=(1,...,1), hence (I+BB*)y=x.
    assert [z+sum(y) for z in y]==v
    plan=dict(schema='gadget.scientific.plan.v1',id='reuse',question='Prepared source reuse',source_digest=serial['source_digest'],steps=[dict(id=str(i),query=green) for i in range(2)])
    (root/'plan.json').write_text(json.dumps(plan));reuse=command('plan','--study',study,'--request',root/'plan.json')
    assert reuse['steps'][1]['certificate']['result']['factorization_reused']
    window=dict(operation='tensor_window',field='all.b1',lower=[1],upper=[4],available_sequence='0')
    w=query(window)['result'];assert w['reported']==2 and w['missing']==1 and w['sum']['rational']=='12'
    window['lower']=[3];window['upper']=[3];assert query(window)['result']['sum'] is None
    condition=dict(path=['incompatible_sections'],value_type='count',comparison='greater',value=0)
    inspect=dict(id='inspect',type='query',query=dict(operation='inspect'))
    branch=dict(id='branch',type='branch',**{'from':'sections'},condition=condition,on_true=inspect,on_false=dict(id='wrong',type='query',query=dict(operation='coverage',field='all.a0')))
    claims=dict(id='checked',type='check_claims',**{'from':'sections'},claims=[condition,dict(path=['sections',0,'compatible'],value_type='boolean',comparison='equal',value=True),dict(path=['sections',0,'residual_quadrance'],value_type='rational',comparison='greater',value='0')])
    entry=dict(id='graph',type='graph',workers=2,steps=[dict(id='sections',after=[],node=dict(id='sections',type='query',query=local)),dict(id='branch',after=['sections'],node=branch),dict(id='checked',after=['sections'],node=claims)])
    workflow=dict(schema='gadget.scientific.workflow.v1',id='dependent-contract',source_digest=serial['source_digest'],limits=dict(depth=12,nodes=32,model_calls=0,parallel_workers=2),definitions=[],entry=entry)
    (root/'workflow.json').write_text(json.dumps(workflow));run=command('workflow','--study',study,'--request',root/'workflow.json')
    assert run['status']=='complete';receipt=[e['result'] for e in run['trace'] if e['node_id']=='checked'][0]
    assert receipt['matched']==2 and receipt['contradicted']==1
    repair=dict(id='repair',type='branch',**{'from':'checked'},condition=dict(path=['contradicted'],value_type='count',comparison='greater',value=0),on_true=dict(id='repair_inspect',type='query',query=dict(operation='inspect')),on_false=dict(id='clean_inspect',type='query',query=dict(operation='inspect')))
    workflow['entry']['steps'].append(dict(id='repair',after=['checked'],node=repair));(root/'workflow.json').write_text(json.dumps(workflow))
    repaired=command('workflow','--study',study,'--request',root/'workflow.json');assert any(e['node_id']=='repair_inspect' for e in repaired['trace'])
    # A native dependency cycle rejects before a leaf can run.
    workflow['entry']['steps'][0]['after']=['branch'];(root/'workflow.json').write_text(json.dumps(workflow));command('workflow','--study',study,'--request',root/'workflow.json',success=False)
    # A declared cohort exclusion is distinct from an unreported count.
    selected_rows=copy.deepcopy(rows)
    for i,row in enumerate(selected_rows):row.append('1' if i==0 else '0')
    selected_raw='\n'.join(['\t'.join(header+['selected'])]+['\t'.join(r) for r in selected_rows])+'\n'
    (root/'table.tsv').write_text(selected_raw);selected_spec=copy.deepcopy(spec)
    selected_spec['source']['sha256']='sha256:'+hashlib.sha256(selected_raw.encode()).hexdigest()
    selected_spec['columns'].append(dict(id='selected',channel='observation',quantity='declared_selection',units='dimensionless',value_type='indicator',dictionary_start=0,dictionary_end=len(dictionary)))
    (root/'selected-spec.json').write_text(json.dumps(selected_spec));command('import-table','--spec',root/'selected-spec.json','--out',root/'selected-study')
    selected_query=copy.deepcopy(paired);selected_query['selection']=['all.selected']
    (root/'query.json').write_text(json.dumps(selected_query));selected=command('query','--study',root/'selected-study/study.json','--query',root/'query.json')['result']
    assert (selected['selected_rows'],selected['reported'],selected['missing'],selected['excluded_by_selection'])==(1,1,0,3)
    print('PASS: independent exact paired fields, support, gluing, explicit fibers, coupled solve, tensor windows, dependent jobs, branches and contradictory claims')

if __name__=='__main__':main()
