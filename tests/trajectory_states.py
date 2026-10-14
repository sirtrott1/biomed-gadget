"""All observed trajectory states survive native publication beyond128states."""
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys


def run(binary, consumer, root):
    if root.exists():shutil.rmtree(root)
    root.mkdir(parents=True)
    patterns=[p for p in itertools.product(range(3),repeat=6) if sum(p)==6][:129]
    assert len(patterns)==129
    header=['unit','group']+[f'c{c}v{v}' for c in range(6) for v in range(2)]
    rows=[]
    for i,p in enumerate(patterns):
        rows.append([f'unit-{i}','all-rows']+[str(x) for c in range(6) for x in [10,10+p[c]-1]])
    raw='\n'.join(['\t'.join(header)]+['\t'.join(row) for row in rows])+'\n'
    dictionary='Complete authored ternary change states at fixed coordinate total.\n'
    (root/'source.tsv').write_text(raw);(root/'dictionary.txt').write_text(dictionary)
    spec=dict(schema='gadget.table.import.v1',id='trajectory-state-contract',source_namespace='871005',
        source=dict(path='source.tsv',sha256='sha256:'+hashlib.sha256(raw.encode()).hexdigest(),delimiter='\t',header=True),
        dictionary=dict(path='dictionary.txt',sha256='sha256:'+hashlib.sha256(dictionary.encode()).hexdigest()),
        columns=[dict(id=n,channel='observation',quantity=n,units='synthetic_measurement',value_type='categorical' if n in ['unit','group'] else 'quantitative',dictionary_start=0,dictionary_end=len(dictionary)) for n in header],
        identity_column='unit',contrasts=[],missing_tokens=['NA'],relation_columns=['group'],
        question='Complete observed change states',population='all129authored rows',selection='all rows',limitations=['synthetic source'],source_url='synthetic://trajectory-states',attribution='GADGET fixture',license='synthetic')
    (root/'spec.json').write_text(json.dumps(spec));subprocess.run([str(binary),'import-table','--spec',str(root/'spec.json'),'--out',str(root/'study')],check=True,capture_output=True)
    query=dict(operation='trajectory_structure',series=[dict(id=f'channel-{c}',fields=[f'all.c{c}v0',f'all.c{c}v1']) for c in range(6)],selection=[],support='joint_reported',time=dict(axis='ordinal',unit='authored_visit',positions=['0','1']),faces='coordinate_exchange_identities',semantics='observed_feature_change_relations',**{'lambda':'1'},minimum_witnesses=2,endpoint_fields=[])
    (root/'query.json').write_text(json.dumps(query));p=subprocess.run([str(binary),'query','--study',str(root/'study/study.json'),'--query',str(root/'query.json')],check=True,capture_output=True,text=True);result=json.loads(p.stdout)['result']
    assert result['states']==[list(p) for p in patterns] and result['faces']==[]
    interval=result['intervals'][0];assert interval['state_indices']==list(range(129)) and interval['fields'][0]['input']==['1']*129
    for c,field in enumerate(interval['fields'][1:]):assert field['input']==[str(p[c]-1) for p in patterns]
    expected=[[str(int(p[c]==sign)) for p in patterns] for c in range(6) for sign in range(3)]
    assert interval['B1']['values']==expected and interval['B2']['columns']==0
    for event,p in zip(result['events'],patterns):assert event['changes']==[str(x-1) for x in p]
    (root/'requests.json').write_text(json.dumps({'complete':query}));subprocess.run([str(consumer),str(root/'study/study.json'),str(root/'requests.json'),str(root/'native')],check=True)
    assert (root/'source.tsv').read_bytes()==raw.encode()
    print(json.dumps(dict(status='PASS',states=129,channels=6,complete_original_rows=True,native_replay=True)))


if __name__=='__main__':run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),Path(sys.argv[3]).resolve())
