"""Native publication cutoffs, coupled law admission and independent oracles."""
import copy,importlib.util,json,shutil,subprocess,sys
from pathlib import Path
root_repo=Path(__file__).resolve().parents[1]
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
r=load(root_repo/'examples/bioinformatics/temporal-execution/replay.py','execution_replay');v=load(root_repo/'examples/bioinformatics/temporal-execution/verify.py','execution_verify')
binary,fixture,root=map(lambda s:Path(s).resolve(),sys.argv[1:]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True);checks=0;refusals=0

def check(value):
    global checks
    assert value;checks+=1

def execute(study,q,fail=False):
    global refusals
    p=subprocess.run([str(binary),'query','--study',str(study),'--query',str(r.history.save(root/'query.json',q))],capture_output=True,text=True)
    check(bool(p.returncode)==fail)
    if fail:refusals+=1;return p.stderr
    return json.loads(p.stdout)
q,out=r.history.authored(binary,fixture,root/'native',True);checks+=v.load_history().authored(root/'native')['checks'];r.history.workflow(binary,root/'native')
for t,size in enumerate([0,1,2]):
    inspect=r.history.inspect(binary,root/'native'/('p'+str(t)));binding=inspect['result']['native_history'];check(binding['cutoff']==str(t));check(binding['retained_segments']==size);check(len(out['result']['endpoints'][t]['field'])==[3,3,2][t]);check(out['result']['endpoints'][t]['declaration']['publication_binding']['source_sequence']==str(t))
# Same source and field meanings through independent archives and actual frames.
_,independent=r.history.authored(binary,fixture,root/'independent')
for a,b in zip(out['result']['steps'],independent['result']['steps']):check(a==b)
for a,b in zip(out['result']['endpoints'],independent['result']['endpoints']):
    for key in ['field','word_field','coordinates','grade']:check(a[key]==b[key])
path=root/'native/p1/study.json';original=json.loads(path.read_text())
for value in [None,'3','-1','01','18446744073709551616',1,'0','2']:
    d=copy.deepcopy(original)
    if value is None:d['native'].pop('cutoff')
    else:d['native']['cutoff']=value
    r.history.save(path,d);execute(path,dict(operation='inspect'),True)
r.history.save(path,original)
q,joint=r.joint(binary,root/'coupled');checks+=v.coupled(root/'coupled')['checks'];study=root/'coupled/study/study.json'
for mutate in [lambda x:x['states'][0].update(weight='0'),lambda x:x['states'][0].update(weight='1'),lambda x:x['states'][1].update(id=x['states'][0]['id']),lambda x:x['states'][0]['values'].pop(),lambda x:x['states'][0]['values'][0].update(entity='unknown'),lambda x:x['states'][0]['values'][1].update(entity=x['states'][0]['values'][0]['entity']),lambda x:x['states'][0]['values'][1].update(event_time='4'),lambda x:x['states'][0]['values'][0].update(event_time=None),lambda x:x.update(meaning=''),lambda x:x.update(evidence=[]),lambda x:x['constraints'][0].update(evidence=[]),lambda x:x['constraints'][0].update(allowed=[]),lambda x:x['constraints'][0].update(allowed=[['5']]),lambda x:x['constraints'][0].update(allowed=[['5','1']]),lambda x:x['constraints'][0].update(entities=['unknown']),lambda x:x['constraints'][0].update(entities=[x['constraints'][0]['entities'][0]]*2)]:
    bad=copy.deepcopy(q);mutate(bad['joint_family']);execute(study,bad,True)
check(execute(study,q)==joint)
# Entity assignment order has no effect on the numerical joint family.
permuted=copy.deepcopy(q)
for state in permuted['joint_family']['states']:state['values'].reverse()
changed=execute(study,permuted)['result']['joint_family'];check(changed['states']==joint['result']['joint_family']['states']);check(changed['event_moments']==joint['result']['joint_family']['event_moments'])
receipt=dict(status='PASS',checks=checks,refusals=refusals);r.history.save(root/'contracts.json',receipt);print(json.dumps(receipt))
