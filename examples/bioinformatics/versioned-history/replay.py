#!/usr/bin/env python3
"""Stage publication indexed histories on original coordinates."""
import argparse,csv,hashlib,json,shutil,subprocess,time
from pathlib import Path


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n');return path


def sha(path):return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary,*args):
    start=time.perf_counter();p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
    if p.returncode:raise RuntimeError(p.stderr)
    return json.loads(p.stdout),time.perf_counter()-start


def inspect(binary,root):
    return run(binary,'query','--study',root/'study.json','--query',save(root/'inspect-query.json',dict(operation='inspect')))[0]


def connect(binary,root,ids,field,word):
    # The coordinator admits exact existing publications. It changes neither
    # their source bytes nor their original measurement evidence.
    metadata=json.loads((root/ids[0]/'study.json').read_text())
    metadata['native']['path']=ids[0]+'/'+metadata['native']['path']
    for d in metadata['documents']:d['path']=ids[0]+'/'+d['path']
    ev=metadata['fields'][0]['evidence'];metadata['publications']=[]
    for id in ids:
        pinned=inspect(binary,root/id)
        metadata['publications'].append(dict(id=id,study=dict(path=id+'/study.json',sha256=sha(root/id/'study.json')),source_digest=pinned['source_digest'],meaning='Exact retained publication '+id,evidence=ev))
    save(root/'study.json',metadata)
    grades=len(metadata['contexts'][0]['metrics'])
    maps=[[dict(grade=g,kind='retained_identity',meaning='Identity only on retained source keys in the declared shared namespace',evidence=ev) for g in range(grades)] for _ in ids[1:]]
    lineages=[dict(kind='shared_native_namespace',population_relation='declared_population_change',meaning='Original row and occurrence identities share one explicitly declared native namespace; retained populations differ',evidence=ev) for _ in ids[1:]]
    query=dict(operation='temporal_field_history',endpoints=[dict(publication=id,context='all',field=field) for id in ids],correspondences=maps,lineages=lineages,word=word,**{'lambda':'1'})
    save(root/'history-query.json',query)
    return query


def evaluate(binary,root,query,name='certificate'):
    path=save(root/(name+'-query.json'),query)
    output,elapsed=run(binary,'query','--study',root/'study.json','--query',path)
    save(root/(name+'.json'),output);save(root/(name+'-timing.json'),dict(query_seconds=elapsed));return output


def authored(binary,fixture,root,native_series=False):
    subprocess.run([str(fixture),str(root)]+(['--series'] if native_series else []),check=True,capture_output=True,text=True)
    q=connect(binary,root,['p0','p1','p2'],'x',[dict(kind='hodge'),dict(kind='green',**{'lambda':'1'})])
    output=evaluate(binary,root,q)
    relation=json.loads(json.dumps(q));relation['endpoints']=[dict(publication=id,context='all',field='r') for id in ['p0','p1','p2']];relation['word']=[dict(kind='hodge'),dict(kind='boundary')];evaluate(binary,root,relation,'relations')
    incomplete=json.loads(json.dumps(q));incomplete['endpoints'][1]['field']='partial';incomplete['endpoints'][1]['extension']=dict(meaning='Authored conditional completion, not an observed value',evidence=json.loads((root/'p1/study.json').read_text())['fields'][0]['evidence'],values=[dict(entity='row:30',value='4')]);evaluate(binary,root,incomplete,'completion')
    return q,output


def clinical(binary,source,root):
    root.mkdir(parents=True,exist_ok=False);shutil.copytree(source,root/'original')
    provenance=json.loads((root/'original/provenance.json').read_text())
    assert sha(root/'original/lung.csv').removeprefix('sha256:')==provenance['decoder']['sha256']
    rows=list(csv.DictReader((root/'original/lung.csv').open()));assert len(rows)==228
    # These are constructed data access publications, not longitudinal visits.
    # Each selected row keeps every original variable, missing token and status.
    cuts=[80,160,len(rows)]
    for t,size in enumerate(cuts):
        inputs=root/('input'+str(t));inputs.mkdir();table=inputs/'observations.csv'
        with table.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows[:size])
        dictionary=(root/'original/lung.Rd').read_text()+'\nConstructed cumulative data-access publication '+str(t)+'. Original variables and missing values are unchanged. time is observed follow-up duration, not an imputed death date. source_row is the original source-row identity. This publication axis is an ordering of data availability, not patient visit time.\n'
        (inputs/'dictionary.txt').write_text(dictionary)
        descriptions=[dict(id=c,channel='observation',quantity='followup_duration' if c=='time' else c,units='day' if c=='time' else 'dimensionless',value_type='categorical' if c in ['source_row','sex','status','inst','ph.ecog'] else 'quantitative',dictionary_start=0,dictionary_end=len(dictionary.encode())) for c in rows[0]]
        spec=dict(schema='gadget.table.import.v1',id='clinical-data-publication-'+str(t),source_namespace='996001',source=dict(path='observations.csv',sha256=sha(table),delimiter=',',header=True),dictionary=dict(path='dictionary.txt',sha256=sha(inputs/'dictionary.txt')),columns=descriptions,missing_tokens=['NA'],relation_columns=['sex'],identity_column='source_row',question='How do admitted source fields change as original records become available?',population=str(size)+' original NCCTG records',selection='Constructed prefix of original source rows',limitations=['Constructed publication ordering; unchanged clinical measurements'],source_url='https://github.com/therneau/survival',attribution='Terry Therneau; NCCTG lung cancer data distributed with survival',license='GPL-2, survival package distribution')
        run(binary,'import-table','--spec',save(inputs/'import.json',spec),'--out',root/('p'+str(t)))
        d=json.loads((root/('p'+str(t))/'study.json').read_text());d['contexts'][0]['time']=dict(axis='data_access_publication',position=str(t),unit='publication',population=d['design']['population']);save(root/('p'+str(t))/'study.json',d)
    q=connect(binary,root,['p0','p1','p2'],'all.time',[dict(kind='hodge'),dict(kind='green',**{'lambda':'1'})]);output=evaluate(binary,root,q)
    return q,output


def workflow(binary,root):
    cert=json.loads((root/'certificate.json').read_text());q=json.loads((root/'history-query.json').read_text())
    steps=[dict(id='history',type='query',query=q),dict(id='publication-catalog',type='query',query=dict(operation='inspect'))]
    request=dict(schema='gadget.scientific.workflow.v1',id='publication-indexed-fields',source_digest=cert['source_digest'],limits=dict(depth=2,nodes=1+len(steps),model_calls=0,parallel_workers=len(steps)),definitions=[],entry=dict(id='joined',type='parallel',workers=len(steps),steps=steps))
    out,elapsed=run(binary,'workflow','--study',root/'study.json','--request',save(root/'workflow-request.json',request));assert out['status']=='complete',out.get('error');save(root/'workflow.json',out);save(root/'workflow-timing.json',dict(seconds=elapsed))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixture',type=Path);p.add_argument('--source',type=Path);p.add_argument('--native-series',action='store_true');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.fixture:authored(a.binary.resolve(),a.fixture.resolve(),a.out.resolve(),a.native_series)
    elif a.source:clinical(a.binary.resolve(),a.source.resolve(),a.out.resolve())
    else:raise ValueError('Supply a fixture or the sealed clinical source packet')
    workflow(a.binary.resolve(),a.out.resolve())
