"""Trajectory declarations exceed former population, channel and visit ceilings."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import shutil


def run(binary, consumer, root):
    if root.exists():shutil.rmtree(root)
    root.mkdir(parents=True)
    population, channels, visits = 4097, 7, 9
    header = ['unit','group']+[f'c{c}v{v}' for c in range(channels) for v in range(visits)]
    raw = '\n'.join(['\t'.join(header)]+['\t'.join([f'unit-{i}','all-rows']+
        [str((i % 3)+c+v) for c in range(channels) for v in range(visits)]) for i in range(population)])+'\n'
    dictionary = 'Authored complete scalar fields with declared ordered measurements.\n'
    (root/'source.tsv').write_text(raw)
    (root/'dictionary.txt').write_text(dictionary)
    spec = dict(schema='gadget.table.import.v1', id='trajectory-scope-contract', source_namespace='871004',
        source=dict(path='source.tsv',sha256='sha256:'+hashlib.sha256(raw.encode()).hexdigest(),delimiter='\t',header=True),
        dictionary=dict(path='dictionary.txt',sha256='sha256:'+hashlib.sha256(dictionary.encode()).hexdigest()),
        columns=[dict(id=name,channel='observation',quantity=name,units='synthetic_measurement',
            value_type='categorical' if name in ['unit','group'] else 'quantitative',dictionary_start=0,dictionary_end=len(dictionary)) for name in header],
        identity_column='unit', contrasts=[], missing_tokens=['NA'], relation_columns=['group'],
        question='Complete ordered native trajectory fixture',population='all authored rows',selection='all source rows',
        limitations=['synthetic source'],source_url='synthetic://trajectory-scope',attribution='GADGET fixture',license='synthetic')
    (root/'spec.json').write_text(json.dumps(spec))
    subprocess.run([str(binary),'import-table','--spec',str(root/'spec.json'),'--out',str(root/'study')],check=True,capture_output=True)
    query=dict(operation='trajectory_structure',series=[dict(id=f'channel-{c}',fields=[f'all.c{c}v{v}' for v in range(visits)]) for c in range(channels)],
        selection=[],support='joint_reported',time=dict(axis='ordinal',unit='authored_visit',positions=list(map(str,range(visits)))),
        faces='coordinate_exchange_identities',semantics='observed_feature_change_relations',**{'lambda':'1'},
        minimum_witnesses=5000,endpoint_fields=[])
    (root/'query.json').write_text(json.dumps(query))
    p=subprocess.run([str(binary),'query','--study',str(root/'study/study.json'),'--query',str(root/'query.json')],check=True,capture_output=True,text=True)
    result=json.loads(p.stdout)['result']
    assert result['n']==population and len(result['intervals'])==visits-1
    assert result['states']==[[2]*channels] and result['faces']==[]
    assert len(result['events'])==population*(visits-1)
    for interval in result['intervals']:
        assert interval['fields'][0]['input']==[str(population)]
        assert all(field['input']==[str(population)] for field in interval['fields'][1:])
        assert interval['H1_rank']==0 and interval['repeated_support_H1_image_rank']==0
    (root/'requests.json').write_text(json.dumps({'complete':query}))
    subprocess.run([str(consumer),str(root/'study/study.json'),str(root/'requests.json'),str(root/'native')],check=True)
    assert (root/'source.tsv').read_bytes()==raw.encode()
    print(json.dumps(dict(status='PASS',population=population,channels=channels,visits=visits,
                         original_rows_unchanged=True,native_replay=True)))


if __name__=='__main__':run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),Path(sys.argv[3]).resolve())
