import copy, importlib.util, json, shutil, subprocess, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'examples/bioinformatics/hodge-frames'))
import replay, verify
import method_answers
binary,producer,session,root=[Path(x).resolve() for x in sys.argv[1:]]
shutil.rmtree(root,ignore_errors=True)
replay.examples(binary,producer,root)
records=[verify.verify(x) for x in sorted(root.iterdir()) if (x/'result.json').exists()]
# A coordinator selects each actual foreign publication, never its own equal-size source.
coordinator=json.loads((root/'filled-square/study/study.json').read_text())
coordinator['native']['path']='filled-square/study/'+coordinator['native']['path']
for doc in coordinator['documents']:doc['path']='filled-square/study/'+doc['path']
coordinator['publications']=[]
for name in ['unfilled-square','filled-square']:
    path=root/name/'study/study.json';cert=json.loads((root/name/'result.json').read_text())
    coordinator['publications'].append(dict(id=name,study=dict(path=name+'/study/study.json',sha256=replay.sources.digest(path)),
        source_digest=cert['source_digest'],meaning='Selected actual native source publication',evidence=coordinator['fields'][0]['evidence']))
replay.save(root/'coordinator.json',coordinator)
for name in ['unfilled-square','filled-square']:
    query=json.loads((root/name/'query.json').read_text());query['publication']=name
    result=replay.run(binary,'query','--study',root/'coordinator.json','--query',replay.save(root/'foreign-query.json',query))
    expected=json.loads((root/name/'result.json').read_text())['result']
    assert result['result']==expected
# Exercise the production recipe and checked answer path on an admitted assay.
pipeline=root/'pipeline';pipeline.mkdir()
def source(name,contents):
    path=pipeline/name;path.write_text(contents)
    return dict(path=name,sha256=replay.sources.digest(path))
features=source('features.tsv','id\ng1\ng2\ng3\n')
observations=source('observations.tsv','id\nA\nB\n')
counts=source('counts.tsv','id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\n')
assay=dict(id='counts',kind='expression',format='dense_tsv',source=counts,features=features,
           observations=observations,source_namespace='921100',carrier='features',quantity='read_count',
           units='counts',channel='E',method='authored counts',arithmetic='exact_decimal',
           implicit_entries='measured_zero',missing_tokens=[])
spec=dict(schema='gadget.assays.import.v1',id='native-Hodge-pipeline',
          design=dict(question='Retain complete fields on the supplied source',population='three authored features',selection='all',limitations=['nonclinical fixture']),
          assays=[assay],resources=[],links=[])
replay.run(binary,'import-assays','--spec',replay.save(pipeline/'import.json',spec),'--out',pipeline/'workspace')
q_assay=dict(operation='hodge_frame_decomposition',publication='self',source_sequence='0',grade=0,
             grade_carriers=['1'],fields=['all.A','all.B'],metrics=dict(kind='unit'),workers=2)
wrapped=dict(schema='gadget.assays.query.v1',operation='native_query',assay='counts',query=q_assay)
result=replay.run(binary,'assay-query','--workspace',pipeline/'workspace/workspace.json',
                  '--request',replay.save(pipeline/'query.json',wrapped))
method_answers.MEANINGS['hodge_frame_decomposition']='hodge_frame_decomposition'
answer_checks=method_answers.qualify(binary,pipeline/'answers',pipeline/'workspace/workspace.json',{'fields':(wrapped,result)})
base=root/'positive-metric';q=json.loads((base/'query.json').read_text());refusals=0
for change in [lambda x:x.update(source_sequence='8'),lambda x:x.update(source_sequence='07'),
               lambda x:x.update(grade=0),lambda x:x.update(fields=['partial']),
               lambda x:x.update(fields=['category']),lambda x:x.update(fields=[]),
               lambda x:x.update(grade_carriers=['10','11']),
               lambda x:x.update(grade_carriers=['11','10','12']),
               lambda x:x['metrics'].update(tensors=['M0','M1']),
               lambda x:x['metrics'].update(tensors=['M1','M0','M2']),
               lambda x:x['metrics'].update(kind='diagonal'),
               lambda x:x.update(publication='unknown'),lambda x:x.update(workers=0)]:
    changed=copy.deepcopy(q);change(changed);replay.save(root/'refused-query.json',changed)
    p=subprocess.run([str(binary),'query','--study',str(base/'study/study.json'),
                      '--query',str(root/'refused-query.json')],capture_output=True,text=True)
    assert p.returncode!=0,changed;refusals+=1
p=subprocess.run([str(session),str(base)],capture_output=True,text=True)
assert p.returncode==0,p.stderr
replay.save(root/'session.json',json.loads(p.stdout))
# Full-family rejection is required even when the defective interface or
# metric is remote from the selected grade.
family_refusals=[]
for name,template in [('remote-boundary','grade-three'),('unselected-metric','positive-metric')]:
    spec=json.loads((root/template/'spec.json').read_text())
    if name=='remote-boundary':
        spec['source']['higher'][-1]['occurrences'][-1]['magnitude']='2'
        spec['fields'].append(replay.sources.field('v',105,10,['1']*len(spec['source']['vertices'])))
    else:spec['operators'][0]['map']['entries'][0]['value']='-2'
    path=replay.save(root/(name+'.json'),spec)
    made=subprocess.run([str(producer),str(path),str(root/name)],capture_output=True,text=True)
    if made.returncode:
        family_refusals.append(dict(case=name,stage='source admission'));continue
    query=json.loads((root/template/'query.json').read_text());query['fields']=['x']
    if name=='remote-boundary':query['grade']=0;query['fields']=['v']
    failed=subprocess.run([str(binary),'query','--study',str(root/name/'study.json'),
        '--query',str(replay.save(root/(name+'-query.json'),query))],capture_output=True,text=True)
    assert failed.returncode!=0,name
    family_refusals.append(dict(case=name,stage='native full-family preparation',reason=failed.stderr.strip()))
replay.save(root/'contracts.json',dict(status='PASS',refusals=refusals,independent=records,answer_checks=answer_checks,family_refusals=family_refusals))
print(json.dumps(dict(status='PASS',refusals=refusals,cases=len(records))))

# Reopen a result bound to a referenced publication after admission through a coordinator.
foreign_query = json.loads((root/'unfilled-square/query.json').read_text())
foreign_query['publication'] = 'unfilled-square'
replay.save(root/'foreign-native-query.json', foreign_query)
foreign_output = root/'foreign-native.rfcxs'
cli = subprocess.run([str(binary), 'query', '--study', str(root/'coordinator.json'),
                      '--query', str(root/'foreign-native-query.json'), '--out', str(foreign_output)],
                     capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr
cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(foreign_output),
                      '--out-record', str(root/'foreign-replay.rfcxs')], capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr

native_records = []
for case in sorted(root.iterdir()):
    if not (case/'query.json').is_file() or not (case/'study/study.json').is_file():
        continue
    checked = subprocess.run([str(session), str(case), 'native'], capture_output=True, text=True)
    assert checked.returncode == 0, (case.name, checked.stderr)
    native_records.append(dict(case=case.name, **json.loads(checked.stdout)))
    replayed = case/'native/cli-replay.rfcxs'
    cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(case/'native/result.rfcxs'),
                          '--out-record', str(replayed)], capture_output=True, text=True)
    assert cli.returncode == 0, cli.stderr
    changed = case/'native/changed-field.rfcxs'
    if changed.exists():
        refused = case/'native/refused.rfcxs'
        cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(changed),
                              '--out-record', str(refused)], capture_output=True, text=True)
        assert cli.returncode != 0 and not refused.exists(), case.name
replay.save(root/'native-contracts.json', native_records)
print(json.dumps(dict(status='PASS', native_cases=len(native_records), checks=sum(x['checks'] for x in native_records))))
