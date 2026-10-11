"""Dictionary-defined exact maps, native branch bindings and detached lineage."""
import copy
import hashlib
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path
binary=Path(sys.argv[1]).resolve(); root=Path(sys.argv[2]).resolve()
if root.exists():shutil.rmtree(root)
root.mkdir(parents=True)
public=Path(__file__).resolve().parents[1]/'examples/bioinformatics/staged-assays'
shutil.copytree(public,root/'inputs')
checks=0
def check(condition,message='contract'):
 global checks
 assert condition,message
 checks+=1
def save(name,value):
 p=root/name;p.write_text(json.dumps(value,indent=2)+'\n');return p
def run(*args,fail=False):
 p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
 if fail:
  check(p.returncode!=0,p.stdout);return p.stderr
 check(p.returncode==0,p.stderr)
 return json.loads(p.stdout) if p.stdout.strip() else None
def payload(v):
 if isinstance(v,list):return [payload(x) for x in v]
 if isinstance(v,dict):return {k:payload(x) for k,x in v.items() if k not in ['source_digest','certificate_sha256','gadget_sha256','gadget_sources_sha256','rexgraph_sha256','rexgraph_headers_sha256','native_sha256'] and not(k=='execution' and isinstance(x,dict) and x.get('owner','').startswith('rfc::'))}
 return v
run('import-assays','--spec',root/'inputs/import.json','--out',root/'source')
manifest=root/'source/workspace.json'; workspace=json.loads(manifest.read_text())
for name in ['rna','atac']:
 query=json.loads((root/'inputs'/f'{name}-design.json').read_text())
 result=run('assay-query','--workspace',manifest,'--request',root/'inputs'/f'{name}-design.json')['result']
 original=(public.parent/'assay-sections'/f'{name}.tsv').read_text().splitlines()
 check(result['rows']==[r.split('\t')[1:] for r in original[1:]])
 check(result['quantity']==query['count_policy'].replace('disjoint_count_sum','umi_count' if name=='rna' else 'read_count'))
 proof=result['design']['observations'];dictionary=(root/'inputs'/f'{name}-observations.tsv').read_bytes()
 check(proof['dictionary_sha256']=='sha256:'+hashlib.sha256(dictionary).hexdigest())
 check(proof['indexed_columns']==3)
 for row in proof['rows']:
  span=row['evidence'][0];check(dictionary[span['start']:span['end']].decode()==span['text'])
 check(proof['rows'][-1]['excluded'] and proof['rows'][-1]['missing_columns']==['state'])
 check(proof['rows'][-2]['excluded'] and not proof['rows'][-2]['missing_columns'])
 compiled=result['design']['compiled_request'];reference=run('assay-query','--workspace',manifest,'--request',save('compiled-'+name+'.json',compiled))['result']
 measured=copy.deepcopy(result);measured.pop('design');check(payload(measured)==payload(reference))
 serial=copy.deepcopy(query);serial['workers']=1
 a=payload(run('assay-query','--workspace',manifest,'--request',save('serial-'+name+'.json',serial))['result']);b=payload(result)
 a['design']['compiled_request'].pop('workers');b['design']['compiled_request'].pop('workers');check(a==b)
# Arbitrary authored overlapping, signed, multicolumn feature/observation maps.
q=json.loads((root/'inputs/rna-design.json').read_text());q.pop('count_policy')
q['feature_design']=dict(mode='groups',groups=[dict(id='response',where=[dict(column='pathway',values=['response'])],coefficient='1/2'),dict(id='gene-A',where=[dict(column='symbol',values=['gene-A'])],coefficient='-2')],missing_tokens=[],unassigned_policy='error',overlap_policy='allow')
q['observation_design']['groups']=[dict(id='two-states',where=[dict(column='state',values=['before','after']),dict(column='cell_type',values=['malignant'])],coefficient='3/2')]
r=run('assay-query','--workspace',manifest,'--request',save('signed-groups.json',q))['result']
check([[Q(v) for v in row] for row in r['rows']]==[[Q(51,4)],[Q(-30)]])
check(r['design']['features']['rows'][0]['groups']==['response','gene-A'])
q['feature_design']['groups'][0]['coefficient']='0'
zero=run('assay-query','--workspace',manifest,'--request',save('zero-map.json',q))['result'];check(zero['rows'][0]==['0'])
# Admission policies apply to dictionary membership independently of measured support.
base=json.loads((root/'inputs/rna-design.json').read_text())
for name,edit in [
 ('unassigned',lambda q:q['observation_design'].__setitem__('unassigned_policy','error')),
 ('unknown-column',lambda q:q['observation_design']['groups'][0]['where'][0].__setitem__('column','foreign')),
 ('empty-group',lambda q:q['observation_design']['groups'][0]['where'][0].__setitem__('values',['foreign'])),
 ('missing-value',lambda q:q['observation_design']['groups'][0]['where'][1].__setitem__('values',['NA'])),
 ('duplicate-clause',lambda q:q['observation_design']['groups'][0]['where'].append(copy.deepcopy(q['observation_design']['groups'][0]['where'][0]))),
 ('zero-count-map',lambda q:q['observation_design']['groups'][0].__setitem__('coefficient','0')),
 ('identity-policies',lambda q:q['feature_design'].__setitem__('overlap_policy','allow')),
 ('unknown-member',lambda q:q['observation_design'].__setitem__('inferred_groups',True)),
]:
 bad=copy.deepcopy(base);edit(bad);run('assay-query','--workspace',manifest,'--request',save(name+'.json',bad),fail=True)
for overlap_policy in ['error','allow']:
 bad=copy.deepcopy(base);other=copy.deepcopy(bad['observation_design']['groups'][0]);other['id']='duplicate';bad['observation_design']['groups'].append(other);bad['observation_design']['overlap_policy']=overlap_policy
 run('assay-query','--workspace',manifest,'--request',save('overlap-'+overlap_policy+'.json',bad),fail=True)
# The native work graph joins both generating branches and the selected source.
plan=run('pipeline-prepare','--workspace',manifest,'--request',root/'inputs/pipeline.json');save('plan.json',plan)
index=run('pipeline','--workspace',manifest,'--request',root/'plan.json','--out',root/'pipeline')
check([x['id'] for x in index['results']]==['rna-stage','atac-stage','combined','sections']+['support-'+str(i) for i in range(5)])
check(index['dispatcher']=='rfc::runtime::WorkGraph' and index['status']=='complete')
bound=root/'pipeline/combined/workspace';combined=run('record-view','--input',bound/'workspace.rfcxs','--member','manifest')
check({a['declaration']['id'] for a in combined['assays']}=={'rna','atac','protein'})
check(len(combined['links'])==5 and len(combined['assays'])==3)
for a in combined['assays']:
 alias=a['declaration']['id']; grouped=alias!='protein'
 parent=root/'pipeline'/(alias+'-stage')/'workspace' if grouped else root/'source'
 prior=(run('record-view','--input',parent/'workspace.rfcxs','--member','manifest') if grouped else json.loads((parent/'workspace.json').read_text()))['assays'][0 if grouped else 2]
 for p,old in zip(a['partitions'],prior['partitions']):
  check(p['source_digest']==old['source_digest'] and p['source_namespace']==old['source_namespace'])
  check(Path(p['path']).suffix=='.rfcxs','Native study binding')
# Independent six-coordinate before/after projection witnesses from the declaration.
section=json.loads((root/'pipeline/sections/result.json').read_text())['assays'][0]
# Full numerical certificate is retained in the published assay provenance.
derived_root=root/'pipeline/sections/workspace';derived=json.loads((derived_root/'workspace.json').read_text())
proof=json.loads((derived_root/derived['assays'][0]['originals']['provenance-0']['path']).read_text())['result']
check(Q(proof['sections'][0]['residual_quadrance'])==Q(1420,51))
check([Q(v) for f in proof['sections'][0]['fields'] for v in f['compatible']]==[Q(3),Q(84,17),Q(3),Q(191,51),Q(382,51),Q(443,51)])
check(proof['sections'][1]['residual_quadrance']=='0' and proof['sections'][1]['complete'])
check(not proof['sections'][2]['complete'] and proof['sections'][2]['fields'][0]['input'][0] is None)
# Model-bound selection closes through every explicitly selected input branch.
intent=json.loads((root/'inputs/intent.json').read_text());catalog=run('recipe-catalog','--workspace',manifest,'--intent',root/'inputs/intent.json')
check(next(x for x in catalog['recipes'] if x['id']=='combined')['dependencies']==['atac-stage','rna-stage'])
proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',action='recipes',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],goals=['support'])
admitted=run('recipe-proposal','--workspace',manifest,'--intent',root/'inputs/intent.json','--request',save('proposal.json',proposal),'--prompt',intent['answer_contract']['question'])
check(admitted['included_goals']==['rna-stage','atac-stage','combined','sections','support'])
for selected in ['later','combined']:
 bad=copy.deepcopy(json.loads((root/'inputs/pipeline.json').read_text()));bad['goals'][2]['parameters']['assays'][0]['input']=selected
 run('pipeline-prepare','--workspace',manifest,'--request',save('bad-input-'+selected+'.json',bad),fail=True)
for name,edit in [
 ('unknown-assay',lambda q:q['goals'][2]['parameters']['assays'][0].__setitem__('assay','foreign')),
 ('duplicate-alias',lambda q:q['goals'][2]['parameters']['assays'][1].__setitem__('id','rna')),
 ('bad-link',lambda q:q['goals'][2]['parameters']['links'][0].__setitem__('left_observation','foreign')),
 ('namespace-collision',lambda q:q['goals'][2]['parameters']['assays'][1].update(input='rna-stage',assay='rna-pseudobulk')),
]:
 bad=copy.deepcopy(json.loads((root/'inputs/pipeline.json').read_text()));edit(bad)
 prepared=run('pipeline-prepare','--workspace',manifest,'--request',save(name+'-default.json',bad));save(name+'-plan.json',prepared)
 run('pipeline','--workspace',manifest,'--request',root/(name+'-plan.json'),'--out',root/('failed-'+name),fail=True)
 check(not(root/('failed-'+name)).exists())
# A bound collection is portable and its replay declaration contains only local records.
shutil.copytree(bound,root/'portable-bound');shutil.copytree(derived_root,root/'portable-derived')
shutil.move(root/'pipeline',root/'removed-pipeline');shutil.move(root/'source',root/'removed-source')
run('native-workspace','--workspace',root/'portable-bound/workspace.rfcxs','--out',root/'bound-reopened.rfcxs')
check((root/'bound-reopened.rfcxs').read_bytes()==(root/'portable-bound/workspace.rfcxs').read_bytes(),'Complete bound archive reencoding')
run('replay-derivation','--workspace',root/'portable-derived/workspace.json','--assay','rna-compatible','--out',root/'replayed')
# Each grouped branch remains reproducible from its nested original source files.
for name in ['rna','atac']:
 run('replay-derivation','--workspace',root/'portable-bound/workspace.rfcxs','--assay',name,'--out',root/('replayed-'+name))
# Binding also retains a nontrivial three-grade source and selected physical bytes.
def source_record(name,value):
 p=root/name;p.write_text(value);return dict(path=name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())
spec=copy.deepcopy(json.loads((root/'inputs/import.json').read_text()));spec['assays']=spec['assays'][:1]
a=spec['assays'][0];a['id']='topological';a['carrier']='observations';a['source_namespace']='976100'
for key in ['source','features','observations']:a[key]['path']='inputs/'+a[key]['path']
vertices=['rna-before-1','rna-before-2','rna-after-1']
cycle=[dict(id=f'e{i}',grade=1,meaning='Explicit oriented comparison',participants=[dict(cell=vertices[i],exists=True,negative=True,magnitude='1'),dict(cell=vertices[(i+1)%3],exists=True,negative=False,magnitude='1')]) for i in range(3)]
cycle[0]['participants'] += [dict(cell=vertices[0],exists=False,negative=False,magnitude='7'),dict(cell=vertices[1],exists=True,negative=True,magnitude='0')]
face=dict(id='face',grade=2,meaning='Authored cycle filling',participants=[dict(cell=e['id'],exists=True,negative=False,magnitude='1') for e in cycle])
cell3=dict(id='retained-third-grade',grade=3,meaning='Authored zero-boundary primary identity with separate cancelling occurrences',participants=[dict(cell='face',exists=True,negative=n,magnitude='1') for n in [False,True]])
a['relations']=source_record('topology.json',json.dumps(dict(relations=cycle+[face,cell3])))
sequence=source_record('sequence.fa','>example\nACGTACGT\n');spec['resources']=[dict(id='sequence',format='fasta',source=sequence,method='Authored literal source sequence')]
run('import-assays','--spec',save('topology-import.json',spec),'--out',root/'topology-source')
binding=copy.deepcopy(json.loads((root/'inputs/bindings.json').read_text()));binding['assays']=[dict(input='source',assay='topological',id='retained')];binding['resources']=[dict(input='source',resource='sequence',id='bound-sequence')];binding['links']=[]
settings=dict(schema='gadget.pipeline.default.v1',profile='custom',question='Retain a complete authored chain family and selected physical source',goals=[dict(id='bind',method='bind_assays',parameters={k:v for k,v in binding.items() if k!='schema'}),dict(id='chain',method='chain_certificate',input='bind',assay='retained',parameters=dict(context='all',grade=3,coefficient_domain='rational',values=['1']))])
plan=run('pipeline-prepare','--workspace',root/'topology-source/workspace.json','--request',save('topology-default.json',settings))
run('pipeline','--workspace',root/'topology-source/workspace.json','--request',save('topology-plan.json',plan),'--out',root/'topology-pipeline')
original=json.loads((root/'topology-source/workspace.json').read_text());bound_path=root/'topology-pipeline/bind/workspace';w=json.loads((bound_path/'workspace.json').read_text())
old=original['assays'][0]['partitions'][0];new=w['assays'][0]['partitions'][0]
check((bound_path/new['path']).read_bytes()==(root/'topology-source'/old['path']).with_name('study.rfcxs').read_bytes())
check((bound_path/w['resources'][0]['original']['path']).read_bytes()==(root/'sequence.fa').read_bytes())
chain=json.loads((root/'topology-pipeline/chain/results/partition-0.json').read_text())['result']
check(chain['facts']['cycle'] and chain['facts']['cancelled_lower_coordinates']==1)
# Reimport retains original occurrences and complete chain geometry.
run('import-assays','--spec',bound_path/'reimport.json','--out',root/'topology-reimport')
request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='retained',query=dict(operation='chain_certificate',context='all',grade=3,coefficient_domain='rational',values=['1']))
r=run('assay-query','--workspace',root/'topology-reimport/workspace.json','--request',save('topology-query.json',request))
check(payload(r['result']['answers'][0]['certificate']['result'])==payload(chain))
print(json.dumps(dict(status='PASS',checks=checks)))
