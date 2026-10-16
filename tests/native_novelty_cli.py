import copy
import ctypes
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'examples/bioinformatics/metric-novelty'))
import replay
import method_answers
binary,producer,session,root=[Path(v).resolve() for v in sys.argv[1:]]
shutil.rmtree(root,ignore_errors=True)
records=replay.examples(binary,producer,root)
base=root/'algebra';q=json.loads((base/'query.json').read_text());refusals=[]
changes=[
 lambda x:x['old'].update(source_sequence='8'),lambda x:x['new'].update(source_sequence='7'),
 lambda x:x['old'].update(source_sequence='07'),lambda x:x['new'].update(publication='unknown'),
 lambda x:x['new'].update(field='partial'),lambda x:x['new'].update(field='category'),
 lambda x:x['new'].update(metric='unknown'),lambda x:x.update(causes=[]),
 lambda x:x.update(causes=['part']),lambda x:x.update(causes=['partial']),lambda x:x.update(causes=['category']),
 lambda x:x['comparison'].update(old_carrier='11'),lambda x:x['comparison']['map'].update(columns=3),
 lambda x:x['comparison']['map']['entries'].append(copy.deepcopy(x['comparison']['map']['entries'][0])),
 lambda x:x.update(preparation='silently-trim'),lambda x:x['comparison']['evidence'][0].update(text='foreign bytes')]
for i,change in enumerate(changes):
 changed=copy.deepcopy(q);change(changed)
 p=subprocess.run([str(binary),'query','--study',str(base/'study.json'),'--query',str(replay.save(root/f'refused-{i}.json',changed))],capture_output=True,text=True)
 assert p.returncode,(changed,p.stdout);refusals.append(dict(case=i,stderr=p.stderr))
# Recreate complete actual sources for hypothesis refusals.
for name,J,M in [('rank-loss',replay.matrix([[1,1],[0,0],[1,1]]),replay.identity(3)),
                 ('nonpositive-complement',replay.matrix([[1,0],[0,1],[0,0]]),replay.matrix([[1,0,0],[0,1,0],[0,0,-1]])),
                 ('asymmetric',replay.matrix([[1,0],[0,1],[1,1]]),replay.matrix([[1,1,0],[0,1,0],[0,0,1]])),
                 ('singular-metric',replay.matrix([[1,0],[0,1],[1,1]]),replay.matrix([[1,0,0],[0,1,0],[0,0,0]]))]:
 case=root/name;case.mkdir()
 admitted=True
 for p,spec in [('p0',replay.declaration('old',7,[2,3],replay.identity(2))),('p1',replay.declaration('new',8,[3,2,8],M))]:
  produced=subprocess.run([str(producer),str(replay.save(case/(p+'-spec.json'),spec)),str(case/p)],capture_output=True,text=True)
  if produced.returncode:
   assert name!='rank-loss';refusals.append(dict(case=name,stage='native source admission',stderr=produced.stderr));admitted=False;break
 if not admitted:continue
 try:replay.connect(binary,case,J)
 except RuntimeError as error:refusals.append(dict(case=name,stderr=str(error)))
 else:raise AssertionError(name+' was admitted')
p=subprocess.run([str(session),str(base)],check=True,capture_output=True,text=True)
session_receipt=json.loads(p.stdout)
# General assay native recipes over the original imported source, not a substituted fixture.
pipeline=root/'pipeline';pipeline.mkdir()
for name,text in [('features.tsv','id\ng1\ng2\ng3\n'),('observations.tsv','id\nA\nB\n'),('counts.tsv','id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\n')]:
 (pipeline/name).write_text(text)
record=lambda n:dict(path=n,sha256=replay.sources.digest(pipeline/n))
assay=dict(id='counts',kind='expression',format='dense_tsv',source=record('counts.tsv'),features=record('features.tsv'),observations=record('observations.tsv'),source_namespace='940100',carrier='features',quantity='read_count',units='counts',channel='E',method='Authored complete paired count columns',arithmetic='exact_decimal',implicit_entries='measured_zero',missing_tokens=[])
spec=dict(schema='gadget.assays.import.v1',id='novelty-pipeline',design=dict(question='Separate the complete change between declared count columns',population='authored three-feature source',selection='all',limitations=['Algebra qualification']),assays=[assay],resources=[],links=[])
original=pipeline/'imported';replay.run(binary,'import-assays','--spec',replay.save(pipeline/'import.json',spec),'--out',original)
source=original/'assays/0/partitions/0/study.json'
addition=dict(id='novelty-counts',namespace='940100',operators=[replay.sources.op('M',940001,'metric',1,1,replay.identity(3))],fields=[])
subprocess.run([str(producer),str(replay.save(pipeline/'addition.json',addition)),str(original/'extended'),str(source)],check=True,capture_output=True,text=True)
envelope=json.loads((original/'workspace.json').read_text());packet=envelope['assays'][0]['partitions'][0]
packet.update(path='extended/study.json',sha256=replay.sources.digest(original/'extended/study.json'),source_digest=replay.run(binary,'query','--study',original/'extended/study.json','--query',replay.save(pipeline/'inspect.json',dict(operation='inspect')))['source_digest'])
envelope.pop('manifest_sha256')
codec=ctypes.CDLL('libjson-c.so.5');codec.json_tokener_parse.argtypes=[ctypes.c_char_p];codec.json_tokener_parse.restype=ctypes.c_void_p
codec.json_object_to_json_string_ext.argtypes=[ctypes.c_void_p,ctypes.c_int];codec.json_object_to_json_string_ext.restype=ctypes.c_char_p;codec.json_object_put.argtypes=[ctypes.c_void_p]
v=codec.json_tokener_parse(json.dumps(envelope).encode());envelope['manifest_sha256']='sha256:'+hashlib.sha256(codec.json_object_to_json_string_ext(v,2|16)).hexdigest();codec.json_object_put(v)
replay.save(original/'workspace.json',envelope)
d=json.loads((original/'extended/study.json').read_text());e=d['operators'][-1]['evidence']
query=dict(operation='metric_relative_novelty',old=dict(publication='self',source_sequence='0',field='all.A'),new=dict(publication='self',source_sequence='0',field='all.B',metric='M'),comparison=dict(old_carrier='1',new_carrier='1',map=replay.identity(3),meaning='Identity of the original three feature coordinates',evidence=e))
wrapped=dict(schema='gadget.assays.query.v1',operation='native_query',assay='counts',query=query)
result=replay.run(binary,'assay-query','--workspace',original/'workspace.json','--request',replay.save(pipeline/'query.json',wrapped))
assert result['result']['answers'][0]['certificate']['result']['split']['innovation']==['3','3','4']
assert result['result']['answers'][0]['certificate']['result']['split']['quadrance']=='34'
method_answers.MEANINGS['metric_relative_novelty']='metric_relative_novelty'
answer_checks=method_answers.qualify(binary,pipeline/'answers',original/'workspace.json',{'metric-change':(wrapped,result)})
# A bound local proposal uses exactly the same runtime independent admission.
proposal=dict(schema='gadget.scientific.agent.proposal.v1',action='query',query=q)
replay.save(base/'proposal.json',proposal)
p=subprocess.run([str(session),str(base),'workspace',str(original/'workspace.json')],check=True,capture_output=True,text=True)
replay.save(root/'contracts.json',dict(status='PASS',independent=records,refusals=refusals,session=session_receipt,workspace=json.loads(p.stdout),answer_checks=answer_checks))
print(json.dumps(dict(status='PASS',cases=len(records),refusals=len(refusals),answer_checks=answer_checks)))

# Native publications replay after both source workspaces are removed.
native_checks = []
for case in sorted(root.iterdir()):
    if not case.is_dir() or not (case/'study.json').exists() or not (case/'p0/study.json').exists():
        continue
    if not (case/'query.json').exists():
        continue
    if not (case/'result.json').exists():
        continue
    checked = subprocess.run([str(session), str(case), 'native'], capture_output=True, text=True)
    assert checked.returncode == 0, (case.name, checked.stdout, checked.stderr)
    native_checks.append(dict(case=case.name, **json.loads(checked.stdout)))
    for result in sorted((case/'native').glob('*/result.rfcxs')):
        replayed = result.with_name('replayed.rfcxs')
        subprocess.run([str(binary), 'replay-scientific', '--record', str(result), '--out-record', str(replayed)], check=True, capture_output=True, text=True)
        forged = result.with_name('changed-field.rfcxs')
        if forged.exists():
            target = result.with_name('refused.rfcxs')
            refusal = subprocess.run([str(binary), 'replay-scientific', '--record', str(forged), '--out-record', str(target)], capture_output=True, text=True)
            assert refusal.returncode and not target.exists(), refusal.stdout
replay.save(root/'native-contracts.json', dict(status='PASS', cases=native_checks))
print(json.dumps(dict(status='PASS', native_cases=len(native_checks), checks=sum(v['checks'] for v in native_checks))))
