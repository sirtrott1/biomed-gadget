"""Temporal source changes, interval evidence and mathematical refusal contracts."""
import copy,importlib.util,json,shutil,subprocess,sys
from pathlib import Path
from fractions import Fraction as Q

examples=Path(__file__).resolve().parents[1]/'examples/bioinformatics/longitudinal-endpoints'
def load(name):
    spec=importlib.util.spec_from_file_location(name,examples/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
replay,verify=load('replay'),load('verify')
binary,fixture,root=map(lambda s:Path(s).resolve(),sys.argv[1:]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True);checks=refusals=0

def check(value):
    global checks
    assert value;checks+=1

def execute(study,query,fail=False):
    global refusals
    p=subprocess.run([str(binary),'query','--study',str(study/'study.json'),'--query',str(replay.save(root/'query.json',query))],capture_output=True,text=True)
    check(bool(p.returncode)==fail)
    if fail:refusals+=1;return
    return json.loads(p.stdout)

q,out=replay.history(binary,fixture,root/'history');checks+=verify.verify_history(root/'history')['checks'];study=root/'history/study'
# The same admitted endpoints support the original general cross moment owner.
moment=dict(operation='temporal_moment_change',context='baseline',other_context='followup',field='baseline.x',other_field='baseline.x',new_field='followup.x',new_other_field='followup.x',order=1,correspondence=q['correspondences'][0][1]);r=execute(study,moment)['result']
check(Q(r['total']['rational'])==Q(48));check(Q(r['total']['rational'])==sum(Q(r[k]['rational']) for k in ('pullback_form','left_innovation','right_innovation','innovation_interaction')))
without=copy.deepcopy(moment);without.pop('correspondence');execute(study,without,True)
for mutate in [lambda v:v['correspondences'][0][1]['entries'].append(copy.deepcopy(v['correspondences'][0][1]['entries'][0])),lambda v:v['correspondences'][0][1]['entries'][0].update(row='gene:A'),lambda v:v['correspondences'][0][1].update(grade=0),lambda v:v['correspondences'][0][1].update(evidence=[]),lambda v:v['endpoints'][1].update(context='baseline'),lambda v:v.update(word=[dict(kind='diagonal_field',fields=['baseline.x','followup.G','final.x'])]),lambda v:v.update(word=[dict(kind='boundary'),dict(kind='boundary')]),lambda v:v.update(**{'lambda':'0'})]:
    bad=copy.deepcopy(q);mutate(bad);execute(study,bad,True)
# Common terminal accumulation is independent of the selected word.
check(out['result']['terminal_channels']==[dict(label='initial',field=['4','-1']),dict(label='innovation_0',field=['5','-4']),dict(label='innovation_1',field=['4','-2'])])
# Opposing structural and amplitude changes cancel completely at the output.
cancel=copy.deepcopy(q);cancel['endpoints']=cancel['endpoints'][:2];cancel['correspondences']=cancel['correspondences'][:1];cancel['word']=[dict(kind='boundary')];j=cancel['correspondences'][0][0];j['kind']='explicit';j['entries']=[dict(row='gene:'+name,column='gene:A',value=str(value)) for name,value in zip('ABC',(3,-2,-1))]
c=execute(study,cancel)['result']['steps'][0];check(c['field_delta']==['0','0','0']);check(Q(c['field_delta_quadrance'])==0);check(Q(c['endpoint_quadrance_change']['total'])==12);check(any(Q(e['value'])!=0 for e in c['field_delta_channel_gram']['entries']))
# Metric changes require no structural births. All source maps may be incomplete.
metadata=json.loads((study/'study.json').read_text());weighted=copy.deepcopy(metadata);weighted['contexts'][1]['metrics'][1]=dict(kind='diagonal',values=['2','3','5']);replay.save(study/'study.json',weighted)
w=execute(study,q)['result'];check(w['steps'][0]['structure']['metric_pullback_defect']!=out['result']['steps'][0]['structure']['metric_pullback_defect']);replay.save(study/'study.json',metadata)
partial=copy.deepcopy(q);partial['correspondences'][0][1]['entries']=[];p=execute(study,partial)['result'];check(p['steps'][0]['structure']['field_innovation']==['3','1','6']);check(p['identity_verified'])
# Pure metric evolution leaves the native source and the observed field intact.
metric_only=copy.deepcopy(metadata);final=next(c for c in metric_only['contexts'] if c['id']=='final');final['block']='2';final['metrics'][1]=dict(kind='diagonal',values=['2','3','5']);same=copy.deepcopy(next(f for f in metadata['fields'] if f['id']=='followup.x'));same.update(id='final.same',context='final');metric_only['fields'].append(same);replay.save(study/'study.json',metric_only)
pure=copy.deepcopy(q);pure['endpoints']=[dict(context='followup',field='followup.x'),dict(context='final',field='final.same')];pure['word']=[dict(kind='hodge')];pure['correspondences']=[[dict(grade=k,kind='retained_identity',meaning='Retained identical source with authored metric change',evidence=q['correspondences'][0][k]['evidence']) for k in range(3)]]
m=execute(study,pure)['result']['steps'][0];check(all(g['structure_unchanged'] for g in m['source_changes']));check(m['structure']['field_innovation']==['0','0','0']);check(m['structure']['metric_pullback_defect']['entries']);replay.save(study/'study.json',metadata)
# Order 17 was formerly refused; this is a finite exact request, not a derivative.
high=copy.deepcopy(moment);high['order']=17;check(execute(study,high)['result']['identity_verified'])
q,out=replay.endpoints(binary,root/'intervals');checks+=verify.verify_endpoints(root/'intervals')['checks'];study=root/'intervals/study'
check(out['result']['unidentified_event_count']==5)
for mutate in [lambda v:v.update(times=['4','2']),lambda v:v.update(times=['-1']),lambda v:v.update(population='invented population'),lambda v:v.update(time_unit='week'),lambda v:v.update(evidence=[]),lambda v:v.update(upper_finite='all.lower'),lambda v:v.update(lower='all.upper')]:
    bad=copy.deepcopy(q);mutate(bad);execute(study,bad,True)
# Two references to the same group share every uncertain original event.
shared=copy.deepcopy(q);shared['contrasts']=[dict(id='identical',left='all.groupA',right='all.groupA',meaning='The same original readout cancels exactly')];z=execute(study,shared)['result'];check(not z['contrast_readout']['entries']);check(all(v['event_fraction_contrast_lower']==v['event_fraction_contrast_upper']==['0'] for v in z['readings']))
unknown=copy.deepcopy(q);unknown['contrasts'][0]['left']='invented';execute(study,unknown,True)
# Bounds for empty strata remain explicitly undefined, with no division by zero.
empty=copy.deepcopy(q);empty['selection']=['all.groupA'];empty['strata']=['all.groupB'];empty.pop('contrasts',None);e=execute(study,empty)['result'];check(not e['groups'][1]['defined'] and e['groups'][1]['size']=='0')
# Query certificates are immutable rerunnable artifacts.
check(execute(study,q)==out)
receipt=dict(status='PASS',checks=checks,refusals=refusals);replay.save(root/'contracts.json',receipt);print(json.dumps(receipt))
