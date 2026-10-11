"""Native declared map contracts, independent rational terms and joined workflows."""
import copy, importlib.util, json, shutil, subprocess, sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1];sys.path.insert(0,str(repo/'examples/bioinformatics/declared-responses'))
import replay,verify
binary,producer,session,root=map(lambda s:Path(s).resolve(),sys.argv[1:]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
replay.authored(binary,producer,root/'square',replay.square());checks=[verify.verify(root/'square')]
rectangular=(replay.declaration('two-records',17,replay.matrix([[2,1],[0,3]]),replay.matrix([[1,0,2],[0,1,-1]]),replay.matrix([[1,2]]),[3,5,7]),replay.declaration('three-records',21,replay.matrix([[1,0,2],[0,2,1],[0,0,3]]),replay.matrix([[1,0],[2,1],[0,1]]),replay.matrix([[1,0,1],[0,1,2]]),[11,13]),replay.matrix([[1,0],[0,1],[1,-1]]),replay.matrix([[1,0,0],[0,1,0]]),replay.matrix([[1],[2]]))
replay.authored(binary,producer,root/'rectangular',rectangular);checks.append(verify.verify(root/'rectangular'))
# The native API's default dimension cap must not leak into GADGET.
n=65
large=(replay.declaration('many-original-coordinates',3,replay.identity(n),replay.matrix([[1] for _ in range(n)]),dict(rows=1,columns=n,entries=[dict(row=0,column=i,value='1') for i in range(n)]),[2]),replay.declaration('many-original-coordinates-new',4,replay.identity(n),replay.matrix([[2] for _ in range(n)]),dict(rows=1,columns=n,entries=[dict(row=0,column=i,value='1') for i in range(n)]),[3]),replay.identity(n),replay.identity(1),replay.identity(1))
replay.authored(binary,producer,root/'beyond-native-default',large);checks.append(verify.verify(root/'beyond-native-default'))
# Complete zero and empty source-coordinate carriers are distinct from missing support.
empty=dict(rows=0,columns=0,entries=[])
case=(replay.declaration('zero-carrier',0,empty,empty,empty,[],replay.source(0)),replay.declaration('zero-carrier-new',1,empty,empty,empty,[],replay.source(0)),empty,empty,empty)
replay.authored(binary,producer,root/'empty',case)
out=json.loads((root/'empty/observed.json').read_text())['result'];assert out['aggregate']==[] and out['pairing']['aggregate_quadrance']=='0'
q=json.loads((root/'square/observed-query.json').read_text());refusals=0;variants=0
def run(query,fail=False):
    p=subprocess.run([str(binary),'query','--study',str(root/'square/study.json'),'--query',str(replay.save(root/'query.json',query))],capture_output=True,text=True)
    assert bool(p.returncode)==fail,p.stderr
    return json.loads(p.stdout) if not fail else p.stderr
for change in [lambda x:x['old'].update(source_sequence='8'),lambda x:x['new'].update(source_sequence='7'),lambda x:x['old'].update(field='partial'),lambda x:x['old'].update(resolving='I'),lambda x:x['new'].update(realization='A'),lambda x:x['observation'].update(old='I'),lambda x:x['pairing'].update(metric='L'),lambda x:x['comparisons']['ambient'].update(old_carrier='20'),lambda x:x['comparisons']['source']['map'].update(rows=2),lambda x:x['comparisons']['ambient']['map']['entries'].append(dict(row=0,column=0,value='1')),lambda x:x['comparisons']['source']['evidence'][0].update(text='unanchored')]:
    x=copy.deepcopy(q);change(x);run(x,True);refusals+=1
for lam in ['0','-1/7','2/3']:
    x=copy.deepcopy(q);x['lambda']=lam;run(x);variants+=1
x=copy.deepcopy(q);x['lambda']='-1/2';run(x,True);refusals+=1
x=copy.deepcopy(q);x['old']['field']='zero';x['new']['field']='zero';out=run(x)['result'];assert all(v=='0' for v in out['aggregate']) and out['pairing']['aggregate_quadrance']=='0';variants+=1
p=subprocess.run([str(session),str(root/'square')],capture_output=True,text=True);assert p.returncode==0,p.stderr
replay.save(root/'session.json',json.loads(p.stdout));replay.save(root/'contracts.json',dict(status='PASS',independent=checks,refusals=refusals,parameter_cases=variants))
print(json.dumps(dict(status='PASS',refusals=refusals,parameter_cases=variants)))

# Native publications replay after both source workspaces are removed.
native_checks = []
for case in sorted(root.iterdir()):
    if not case.is_dir() or not (case/'study.json').exists() or not (case/'p0/study.json').exists():
        continue
    if not (case/'response-query.json').exists():
        continue
    if not (case/'response.json').exists():
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
