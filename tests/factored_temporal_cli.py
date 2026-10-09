"""Full factored map equivalence, source binding, weighting and noncommuting words."""
import copy,importlib.util,json,shutil,subprocess,sys
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(repo/'examples/bioinformatics/factored-temporal'))
from maps import Expansion,expand_result
from verify import verify_pair
spec=importlib.util.spec_from_file_location('history_replay',repo/'examples/bioinformatics/versioned-history/replay.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
binary,fixture,root=map(lambda s:Path(s).resolve(),sys.argv[1:]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
r.authored(binary,fixture,root/'native',True)
cases=0;coordinates=0;bindings=0;refusals=0

def run(query,fail=False):
    p=subprocess.run([str(binary),'query','--study',str(root/'native/study.json'),'--query',str(r.save(root/'query.json',query))],capture_output=True,text=True)
    assert bool(p.returncode)==fail,p.stderr
    return json.loads(p.stdout) if not fail else p.stderr


def pair(query):
    global cases,coordinates,bindings
    full=run(query);q=copy.deepcopy(query);q['map_output']='factored';factored=run(q)
    check=verify_pair(full,factored);cases+=1;coordinates+=check['complete_map_coordinates'];bindings+=check['source_binding_coordinates']
    explicit=copy.deepcopy(query);explicit['map_output']='enumerated';e=run(explicit)
    assert e['result']==full['result']
    r.save(root/('full-'+str(cases)+'.json'),full);r.save(root/('factored-'+str(cases)+'.json'),factored)
    return full,factored

for name in ['history-query','relations-query','completion-query']:
    pair(json.loads((root/'native'/(name+'.json')).read_text()))
base=json.loads((root/'native/history-query.json').read_text())
words=[ [{'kind':'adjoint'},{'kind':'hodge'},{'kind':'boundary'},{'kind':'green','lambda':'2/3'}],
        [{'kind':'diagonal_field','fields':['x','x','x']},{'kind':'hodge'},{'kind':'green','lambda':'7/5'}],
        [{'kind':'hodge'},{'kind':'diagonal_field','fields':['x','x','x']},{'kind':'hodge'},{'kind':'green','lambda':'1/1'}]]
for word in words:
    q=copy.deepcopy(base);q['word']=word;pair(q)
# Unequal grade metrics at each endpoint expose an incorrect transpose-as-adjoint.
for t in range(3):
    path=root/('native/p'+str(t)+'/study.json');d=json.loads(path.read_text())
    for k,metric in enumerate(d['contexts'][0]['metrics']):
        size=sum(c['grade']==k for c in d['cells']);metric.clear();metric.update(kind='diagonal',values=[str(2+t+k+i) for i in range(size)])
    r.save(path,d)
# Re-admit the changed declarations under their actual new hashes and digests.
r.connect(binary,root/'native',['p0','p1','p2'],'x',base['word'])
for word in words+[base['word']]:
    q=copy.deepcopy(base);q['word']=word;pair(q)
# Unsupported output contracts and missing observations remain rejected.
for mode in ['approximate','selected',1,None,{},[]]:
    bad=copy.deepcopy(base);bad['map_output']=mode;run(bad,True);refusals+=1
bad=json.loads((root/'native/completion-query.json').read_text());bad['map_output']='factored';bad['endpoints'][1].pop('extension');run(bad,True);refusals+=1
# The independent reader rejects altered factor dimensions, cycles and singular solves.
_,cert=pair(base)
archive=cert['result']['factored_maps']
for change in ['cycle','shape','singular']:
    a=copy.deepcopy(archive)
    if change=='cycle':a['nodes'][0].update(kind='transpose',input=0)
    elif change=='shape':a['nodes'][-1]['rows']+=1
    else:
        a['nodes']=[dict(kind='matrix',rows=1,columns=1,entries=[]),dict(kind='solve',rows=1,columns=1,system=0)]
    try:Expansion(a).all()
    except (ValueError,ZeroDivisionError):refusals+=1
    else:raise AssertionError('altered map accepted')
receipt=dict(status='PASS',cases=cases,complete_map_coordinates=coordinates,source_binding_coordinates=bindings,refusals=refusals)
r.save(root/'contracts.json',receipt);print(json.dumps(receipt))
