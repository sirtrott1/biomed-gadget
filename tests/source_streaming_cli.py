import copy,json,shutil,subprocess,sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'examples/bioinformatics/source-streaming'))
import replay
binary,producer,contracts,root=map(lambda p:Path(p).resolve(),sys.argv[1:])
shutil.rmtree(root,ignore_errors=True)
result=replay.examples(binary,producer,root)
study=root/'authored/study/study.json'
subprocess.run([str(contracts),str(study)],check=True)
index=json.loads((root/'authored-pages/index.json').read_text())
request=dict(schema='gadget.scientific.source.read.v1',source_digest=index['source_digest'],publication='self',source_sequence=index['source_sequence'],stream=dict(kind='field',id='x'),offset=0,limit=2)
refusals=[]
for mutate in [lambda q:q.update(source_digest='sha256:'+'0'*64),lambda q:q.update(publication='other'),lambda q:q.update(source_sequence='07'),lambda q:q.update(source_sequence='8'),lambda q:q.update(offset=-1),lambda q:q.update(offset=5),lambda q:q.update(limit=0),lambda q:q.update(limit=-1),lambda q:q.update(stream=dict(kind='field',id='missing')),lambda q:q.update(stream=dict(kind='catalog',section='missing')),lambda q:q.update(extra=True)]:
    q=copy.deepcopy(request);mutate(q)
    p=subprocess.run([str(binary),'source-read','--study',str(study),'--request',str(replay.save(root/'invalid.json',q))],capture_output=True,text=True)
    assert p.returncode,p.stdout
    refusals.append(p.stderr.strip())
doc=(study.parent/'declarations.json').read_bytes();position=doc.index('β'.encode())
for offset,limit in [(position+1,4),(position,1)]:
    q={**request,'stream':dict(kind='document',id='declarations'),'offset':offset,'limit':limit}
    p=subprocess.run([str(binary),'source-read','--study',str(study),'--request',str(replay.save(root/'invalid.json',q))],capture_output=True,text=True)
    assert p.returncode,p.stdout
    refusals.append(p.stderr.strip())
q={**request,'offset':4,'limit':2}
page=replay.run(binary,'source-read','--study',study,'--request',replay.save(root/'terminal.json',q))
assert page['complete'] and not page['records'] and page['next_request'] is None
proposal=dict(schema='gadget.scientific.agent.proposal.v1',action='retrieve',read=request)
r=replay.run(binary,'proposal','--study',study,'--request',replay.save(root/'proposal.json',proposal))
assert r['status']=='retrieved'
replay.save(root/'contracts.json',dict(**result,refusals=refusals,current_lease_and_immutable_owner=True,terminal_page=True,retrieval_admitted=True))

# Independent publications retain separate cutoffs despite equal local cell keys.
second=root/'filled';old,new,*_=replay.observations.sources.square()
replay.observations.execute(binary,producer,second,
    {**replay.observations.specification('filled',new['source'],1,4,replay.observations.identity(4)),'sequence':'8'})
coordinator=json.loads(study.read_text());coordinator['native']['path']='authored/study/'+coordinator['native']['path']
for d in coordinator['documents']:d['path']='authored/study/'+d['path']
coordinator['publications']=[]
for name,path in [('unfilled',study),('filled',second/'study/study.json')]:
    idx=replay.run(binary,'source-index','--study',path)
    coordinator['publications'].append(dict(id=name,study=dict(path=str(path.relative_to(root)),sha256=replay.sha(path.read_bytes())),source_digest=idx['source_digest'],meaning='Independent original publication '+name,evidence=coordinator['operators'][0]['evidence']))
coord=replay.save(root/'coordinator.json',coordinator)
idx=replay.run(binary,'source-index','--study',coord)
foreign={**request,'source_digest':idx['source_digest'],'publication':'filled','source_sequence':'8','stream':dict(kind='catalog',section='cells'),'limit':100}
page=replay.run(binary,'source-read','--study',coord,'--request',replay.save(root/'foreign.json',foreign))
assert len(page['records'])==9 and page['binding']['publication_digest']!=page['binding']['source_digest']
q={**foreign,'source_sequence':'7'}
p=subprocess.run([str(binary),'source-read','--study',str(coord),'--request',str(replay.save(root/'invalid.json',q))],capture_output=True,text=True)
assert p.returncode
# A portable original count assay selects an actual native partition, not a detached table.
pipeline=root/'assay';pipeline.mkdir()
features=pipeline/'features.tsv';features.write_text('id\ng1\ng2\ng3\ng4\n')
obs=pipeline/'observations.tsv';obs.write_text('id\nA\nB\n')
counts=pipeline/'counts.tsv';counts.write_text('id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\ng4\t1\t2\n')
artifact=lambda p:dict(path=p.name,sha256=replay.sha(p.read_bytes()))
assay=dict(id='counts',kind='expression',format='dense_tsv',source=artifact(counts),features=artifact(features),observations=artifact(obs),source_namespace='931101',carrier='features',quantity='read_count',units='counts',channel='E',method='Authored original counts',arithmetic='exact_decimal',implicit_entries='measured_zero',missing_tokens=[])
spec=dict(schema='gadget.assays.import.v1',id='source-reader-assay',design=dict(question='Read original counts',population='All authored features',selection='Complete source',limitations=[]),assays=[assay],resources=[],links=[])
workspace=pipeline/'workspace'
replay.run(binary,'import-assays','--spec',replay.save(pipeline/'spec.json',spec),'--out',workspace)
manifest=workspace/'workspace.json';w=json.loads(manifest.read_text())
wrapper=dict(schema='gadget.assays.source.read.v1',source_digest=w['manifest_sha256'],assay='counts',partition=0,read=dict(schema='gadget.scientific.source.index.v1'))
ix=replay.run(binary,'assay-source-read','--workspace',manifest,'--request',replay.save(pipeline/'index-request.json',wrapper))['index']
wrapper['read']={**request,'source_digest':ix['source_digest'],'source_sequence':ix['source_sequence'],'stream':dict(kind='catalog',section='fields'),'limit':1}
read=replay.run(binary,'assay-source-read','--workspace',manifest,'--request',replay.save(pipeline/'request.json',wrapper))
assert read['page']['total']==2 and len(read['page']['records'])==1
continued=replay.run(binary,'assay-source-read','--workspace',manifest,'--request',replay.save(pipeline/'next.json',read['next_request']))
assert continued['page']['complete']
for mutation in [lambda x:x.update(source_digest='sha256:'+'0'*64),lambda x:x.update(partition=1),lambda x:x['read'].update(source_digest=request['source_digest'])]:
    q=copy.deepcopy(wrapper);mutation(q)
    p=subprocess.run([str(binary),'assay-source-read','--workspace',str(manifest),'--request',str(replay.save(root/'invalid.json',q))],capture_output=True,text=True)
    assert p.returncode
    refusals.append(p.stderr.strip())
replay.save(root/'contracts.json',dict(**result,refusals=refusals,foreign_publication_and_cutoff=True,assay_partition_continuation=True,current_lease_and_immutable_owner=True,terminal_page=True,retrieval_admitted=True))
