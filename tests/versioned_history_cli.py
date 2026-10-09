"""Publication identity, support and changing-population contracts."""
import copy,importlib.util,json,shutil,subprocess,sys
from pathlib import Path
examples=Path(__file__).resolve().parents[1]/'examples/bioinformatics/versioned-history'

def load(name):
    spec=importlib.util.spec_from_file_location(name,examples/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
r,v=load('replay'),load('verify')
binary,fixture,root=map(lambda x:Path(x).resolve(),sys.argv[1:]);shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True);checks=0;refusals=0

def check(value):
    global checks
    assert value;checks+=1


def execute(query,fail=False):
    global refusals
    p=subprocess.run([str(binary),'query','--study',str(root/'authored/study.json'),'--query',str(r.save(root/'query.json',query))],capture_output=True,text=True)
    check(bool(p.returncode)==fail)
    if fail:refusals+=1;return p.stderr
    return json.loads(p.stdout)
q,out=r.authored(binary,fixture,root/'authored');checks+=v.authored(root/'authored')['checks'];r.workflow(binary,root/'authored')
check(execute(q)==out)
for mutate in [lambda x:x.pop('lineages'),lambda x:x['lineages'][0].update(kind='unknown'),lambda x:x['lineages'][0].update(population_relation='same_study_population'),lambda x:x['lineages'][0].update(evidence=[]),lambda x:x['endpoints'][1].update(publication='missing'),lambda x:x['endpoints'][1].update(field='missing'),lambda x:x['endpoints'][1].update(context='missing'),lambda x:x['correspondences'][0][0].update(grade=1),lambda x:x['correspondences'][0].pop(),lambda x:x['endpoints'][1].update(field='partial'),lambda x:x['endpoints'][1].update(field='partial',extension=dict(meaning='Missing entity',evidence=q['lineages'][0]['evidence'],values=[]))]:
    bad=copy.deepcopy(q);mutate(bad);execute(bad,True)
conditional=json.loads((root/'authored/completion-query.json').read_text())
for mutate in [lambda x:x['endpoints'][1]['extension']['values'].append(dict(entity='row:20',value='0')),lambda x:x['endpoints'][1]['extension']['values'].append(dict(entity='row:999',value='0')),lambda x:x['endpoints'][1]['extension']['values'].append(dict(entity='row:30',value='5')),lambda x:x['endpoints'][1]['extension'].update(evidence=[])]:
    bad=copy.deepcopy(conditional);mutate(bad);execute(bad,True)
# Equal keys in deliberately independent sources supply no implicit identity.
independent=copy.deepcopy(q);independent['lineages'][0]['kind']='independent_sources';execute(independent,True)
for g,m in enumerate(independent['correspondences'][0]):
    m['kind']='explicit';m['entries']=[]
    a=out['result']['endpoints'][0]['declaration']['source_grades'][g];b=out['result']['endpoints'][1]['declaration']['source_grades'][g]
    for left in a:
        for right in b:
            if left['key']==right['key']:m['entries'].append(dict(row=right['id'],column=left['id'],value='1'))
i=execute(independent)['result'];check(not i['steps'][0]['source_identity_comparable']);check(i['steps'][0]['source_changes']==[]);check('retained_observation_alignment' not in i['steps'][0]);check(i['endpoints']==out['result']['endpoints'] or i['endpoints'][1]['field']==out['result']['endpoints'][1]['field']);check(i['terminal_quadrance']==out['result']['terminal_quadrance'])
r.save(root/'authored/independent-query.json',independent);r.save(root/'authored/independent.json',execute(independent))
# A complete map is allowed to have an image smaller than the next field carrier.
empty=copy.deepcopy(q);empty['correspondences'][0][0].update(kind='explicit',entries=[]);e=execute(empty)['result'];check(e['steps'][0]['structure']['field_innovation']==out['result']['endpoints'][1]['field'])
metadata=json.loads((root/'authored/study.json').read_text());original=copy.deepcopy(metadata)
for mutate in [lambda x:x['publications'][0]['study'].update(sha256='sha256:'+'0'*64),lambda x:x['publications'][0].update(source_digest='sha256:'+'0'*64),lambda x:x['publications'].append(copy.deepcopy(x['publications'][0])),lambda x:x['publications'][0].update(id='self'),lambda x:x['publications'][0]['study'].update(path='../study.json'),lambda x:x['publications'][0].update(evidence=[])]:
    bad=copy.deepcopy(original);mutate(bad);r.save(root/'authored/study.json',bad);execute(q,True)
r.save(root/'authored/study.json',original)
# Altering an external measurement meaning never becomes a coordinate transform.
child=root/'authored/p1/study.json';child_bytes=child.read_text();d=json.loads(child_bytes);changed=copy.deepcopy(d);changed['fields'][0]['units']='other-unit';r.save(child,changed)
modified=copy.deepcopy(original);modified['publications'][1]['study']['sha256']=r.sha(child);modified['publications'][1]['source_digest']=r.inspect(binary,root/'authored/p1')['source_digest'];r.save(root/'authored/study.json',modified);execute(q,True);child.write_text(child_bytes);r.save(root/'authored/study.json',original)
# Publication-local labels can coincide while their exact source binding differs.
check(out['publication_bindings'][0]['field_semantics']['x']['id']==out['publication_bindings'][1]['field_semantics']['x']['id']);check(out['publication_bindings'][0]['source_digest']!=out['publication_bindings'][1]['source_digest'])
# Original occurrence changes remain visible under a zero aggregate boundary.
changes=out['result']['steps'][0]['source_changes'][1]['occurrence_changes'];check(any(c['key']=='1007' for c in changes));check(any(c['orientation_term']!='0' for c in changes));check(any(c['attachment_changed'] for c in changes));check(out['result']['steps'][1]['source_changes'][1]['removed_cell_keys']==['102'])
receipt=dict(status='PASS',checks=checks,refusals=refusals);r.save(root/'contracts.json',receipt);print(json.dumps(receipt))
