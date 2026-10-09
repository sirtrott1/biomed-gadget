"""Researcher-bound answer coverage, native conflict cores and retained outputs."""
import copy
import hashlib
import itertools
import json
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)
checks = 0

def check(value):
    global checks
    assert value
    checks += 1

def run(*args, fail=False):
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    assert bool(p.returncode) == fail, (args,p.stdout,p.stderr)
    check(True)
    return p.stderr if fail else json.loads(p.stdout)

def save(name, value):
    path = root/name
    path.write_text(value if isinstance(value, str) else json.dumps(value))
    return path

def source(name, value):
    path = save(name, value)
    return dict(path=name, sha256='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest())

features=source('features.tsv','id\nf1\nf2\nf3\n')
observations=source('observations.tsv','id\nA\nB\n')
counts=source('counts.tsv','id\tA\tB\nf1\t2\t5\nf2\t0\tNA\nf3\t7\t11\n')
design=dict(question='Authored measured counts',population='three declared features',selection='all',limitations=['Constructed method fixture'])
assay=dict(id='counts',kind='expression',format='dense_tsv',source=counts,features=features,observations=observations,source_namespace='986100',carrier='features',quantity='read_count',units='counts',channel='E',method='Authored exact counts',arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA'],partition_size=2)
whole={**assay,'id':'whole','source_namespace':'987100','partition_size':3}
spec=dict(schema='gadget.assays.import.v1',id='scientific-answer-contracts',design=design,assays=[assay,whole],resources=[],links=[])
run('import-assays','--spec',save('import.json',spec),'--out',root/'workspace')
manifest=root/'workspace/workspace.json'
question='Inspect reported support in A and B and the declared Green field of A.'
recipes=[dict(id=name,purpose='Authored '+name,goal=dict(method=method,parameters=params)) for name,method,params in [
    ('support-A','coverage',dict(fields=['all.A'])),
    ('alternative-A','coverage',dict(fields=['all.A'])),
    ('support-B','coverage',dict(fields=['all.B'])),
    ('green-A','green_action',dict(context='all',field='all.A',**{'lambda':'3/2'}))]]
recipes[-1]['goal']['assay']='whole'
requirements=[dict(id='A-observed',observable='field_support',bindings=['support-A','alternative-A']),dict(id='B-observed',observable='field_support',bindings=['support-B']),dict(id='A-green',observable='green_action',bindings=['green-A'])]
intent=dict(schema='gadget.pipeline.intent.v1',profile='rna_seq',assay='counts',recipes=recipes,answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,requirements=requirements))

def prepare(value, workspace=manifest):
    path=save('active-intent.json',value)
    catalog=run('recipe-catalog','--workspace',workspace,'--intent',path)
    check(catalog['scientific_admission']['status'] in {'satisfied','unsatisfied','untyped'})
    return path,catalog

def proposal(catalog, names=None, reason=None):
    out=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],action='clarify' if reason else 'recipes')
    out.update(dict(reason=reason) if reason else dict(goals=names))
    return out

def admit(path,catalog,names=None,reason=None,out=None,prompt=question,workspace=manifest):
    args=['recipe-proposal','--workspace',workspace,'--intent',path,'--request',save('active-proposal.json',proposal(catalog,names,reason)),'--prompt',prompt]
    if out: args+=['--out',root/out]
    return run(*args)

path,catalog=prepare(intent)
check(catalog['scientific_admission']['status']=='satisfied')
check(catalog['scientific_admission']['native']['constraint_complete'])
check(catalog['scientific_admission']['native']['assignment']==[1,1,1])
check(next(x for x in catalog['observables'] if x['id']=='population_p_value')['methods']==[])
# Every numerical recipe has one reviewed producing observable. Workspace binding
# is orchestration, and population calibration has no numerical producer.
definitions=run('pipelines')['recipes']
expected_methods={r['method'] for r in definitions if not r['workspace_binding']}
reviewed=[m for o in catalog['observables'] if o['available'] for m in o['methods']]
check(set(reviewed)==expected_methods and len(reviewed)==len(set(reviewed)))
check('constrained_response_model' in expected_methods)
contract_schema=json.loads((Path(__file__).resolve().parents[1]/'schemas/pipeline-answer-contract.schema.json').read_text())
declared_observables=contract_schema['properties']['requirements']['items']['properties']['observable']['enum']
check(set(declared_observables)=={o['id'] for o in catalog['observables']})

# Exhaustively compare native coverage with an independent set-based oracle.
for n in range(1,len(recipes)+1):
    for names in itertools.combinations([r['id'] for r in recipes],n):
        answer=admit(path,catalog,list(names))
        covered=all(set(names)&set(r['bindings']) for r in requirements)
        check((answer['status']=='admitted')==covered)
        check(answer['scientific_admission']['native']['global_section_exists']==covered)
        if not covered:
            check('pipeline' not in answer)
            core=answer['scientific_admission']['native']['unsat_core']
            check(core['subset_minimal'] and len(core['charts'])==2 and len(core['cells'])==1)
            cell=core['cells'][0]
            check(set(core['charts'])=={'require/'+requirements[cell]['id'],'supply/'+requirements[cell]['id']})
            charts={r['name']:r['allowed'] for r in answer['scientific_admission']['native']['charts']}
            predicates=[charts[name] for name in core['charts']]
            check(not any(all(v==allowed for allowed in predicates) for v in (0,1)))
            for removed in range(len(predicates)):
                check(any(all(v==allowed for i,allowed in enumerate(predicates) if i!=removed) for v in (0,1)))
false_clarify=admit(path,catalog,reason='No recipe is available',out='false-clarify')
check(false_clarify['status']=='selection_required' and not (root/'false-clarify').exists())
check(false_clarify['available_answers']['requirements'][0]['eligible_bindings']==['support-A','alternative-A'])
answer=admit(path,catalog,['support-A','support-B','green-A'],out='complete')
check(answer['status']=='admitted')

def verify_receipts(directory, answer):
    index=answer['execution']
    receipt=run('record-view','--input',directory/index['answers'],'--member','result')
    check(index['answers_sha256']=='sha256:'+hashlib.sha256((directory/index['answers']).read_bytes()).hexdigest())
    check(receipt['question_sha256']=='sha256:'+hashlib.sha256(answer['question'].encode()).hexdigest())
    check(receipt['admission_sha256']==answer['scientific_admission']['admission_sha256'])
    unsigned={k:v for k,v in receipt.items() if k!='answers_sha256'}
    # json-c pretty encoding differs from Python's; artifact byte seals are checked above.
    check(unsigned['schema']=='gadget.pipeline.answers.v1')
    leaves=[]
    for req in receipt['answers']:
        check(req['artifacts'])
        for artifact in req['artifacts']:
            check(artifact['binding'] in {r['binding'] for r in req['scope_witnesses']})
            for ref in artifact['results']:
                data=(directory/ref['path']).read_bytes()
                check(ref['sha256']=='sha256:'+hashlib.sha256(data).hexdigest())
                certificate=run('record-view','--input',directory/ref['path'])
                if certificate.get('schema')=='gadget.operation.native.v1':certificate=certificate['result']
                check(certificate['source_digest']==ref['source_digest'])
                check(certificate['query']==ref['query'] and ref['pointer']=='/result')
                check(certificate['certificate_sha256']==ref['certificate_sha256'])
                leaves.append((req['id'],certificate))
    return receipt,leaves
receipt,leaves=verify_receipts(root/'complete',answer)
check(len(leaves)==5)  # Both support partitions and the complete Green carrier.
A=[c['result'] for name,c in leaves if name=='A-observed']
B=[c['result'] for name,c in leaves if name=='B-observed']
check(sum(c['reported'] for c in A)==3 and sum(c['reported'] for c in B)==2)
rows=[r for c in A for r in c['observations']]
check(any(r['reported'] and r['value']['rational']=='0' for r in rows))
check(any(not r['reported'] for c in B for r in c['observations']))
check([Fraction(r['value']['rational']) for name,c in leaves if name=='A-green' for r in c['result']]==[2,0,7])
# Literal bytes, rather than interpreted keyword similarity, bind the question.
for changed in [question+' ',question.lower(),'Compute population-level p-values for a treatment effect from these data.']:
    error=run('recipe-proposal','--workspace',manifest,'--intent',path,'--request',save('stale-question.json',proposal(catalog,['support-A','support-B','green-A'])),'--prompt',changed,'--out',root/'stale-question',fail=True)
    check('question differs' in error and not (root/'stale-question').exists())
# A tempting purpose cannot relabel the method's actual output meaning.
unsupported=copy.deepcopy(intent)
unsupported['answer_contract']['question']='Compute population-level p-values for a treatment effect from these data.'
unsupported['answer_contract']['requirements']=[dict(id='treatment-p-value',observable='population_p_value',bindings=['support-A'])]
unsupported['recipes'][0]['purpose']='Compute treatment p-values with population significance'
p,c=prepare(unsupported)
check(c['scientific_admission']['status']=='unsatisfied')
for names,reason in [(['support-A'],None),(None,'P-values require sampling and null declarations')]:
    result=admit(p,c,names,reason,out='unavailable',prompt=unsupported['answer_contract']['question'])
    check(result['status']=='unsupported' and not (root/'unavailable').exists())
    check(result['scientific_admission']['native']['unsat_core']['subset_minimal'])
# Equal dimensions/labels or support in another declared field is not the requested scope.
wrong=copy.deepcopy(intent);wrong['answer_contract']['requirements']=[dict(id='Green-of-A',observable='green_action',bindings=['support-A'])]
p,c=prepare(wrong)
check(c['scientific_admission']['status']=='unsatisfied')
check(admit(p,c,['support-A'])['status']=='unsupported')
p,c=prepare(intent)
check(admit(p,c,['alternative-A','green-A'],out='wrong-field')['status']=='incompatible')
check(not (root/'wrong-field').exists())
for label,change in [
    ('unknown-observable',lambda x:x['answer_contract']['requirements'][0].update(observable='invented_p_value')),
    ('unknown-binding',lambda x:x['answer_contract']['requirements'][0].update(bindings=['foreign'])),
    ('duplicate-binding',lambda x:x['answer_contract']['requirements'][0].update(bindings=['support-A','support-A'])),
    ('duplicate-requirement',lambda x:x['answer_contract']['requirements'].append(copy.deepcopy(x['answer_contract']['requirements'][0]))),
    ('empty-requirements',lambda x:x['answer_contract'].update(requirements=[])),
    ('injected-contract',lambda x:x['answer_contract'].update(parameters={'null':'invented'}))]:
    invalid=copy.deepcopy(intent);change(invalid)
    run('recipe-catalog','--workspace',manifest,'--intent',save(label+'.json',invalid),fail=True)
# Intent seals bind scopes and requested output meanings, not just recipe names.
stale=copy.deepcopy(intent);stale['answer_contract']['requirements'][0]['bindings']=['support-A']
run('recipe-proposal','--workspace',manifest,'--intent',save('stale-contract.json',stale),'--request',save('old-contract-proposal.json',proposal(c,['support-A','support-B','green-A'])),'--prompt',question,'--out',root/'stale-contract',fail=True)
check(not (root/'stale-contract').exists())
# Retain explicit manual compatibility, but do not call untyped provenance a semantic proof.
legacy=copy.deepcopy(intent);legacy.pop('answer_contract')
p,c=prepare(legacy)
check(c['scientific_admission']['status']=='untyped')
check(admit(p,c,['support-A'],out='manual-legacy')['scientific_admission']['status']=='untyped')
check(not (root/'manual-legacy/answers.json').exists())
# Long questions have no inherited 4096-byte policy; native string representation applies.
long=copy.deepcopy(intent);long['answer_contract']['question']='Observe declared support. '+('alpha beta '*500)
p,c=prepare(long)
check(admit(p,c,['support-A','support-B','green-A'],prompt=long['answer_contract']['question'])['status']=='admitted')
# UTF-8 question hashes bind original bytes, with no keyword translator.
unicode=copy.deepcopy(intent);unicode['answer_contract']['question']='Quels champs mesurés satisfont β ?'
p,c=prepare(unicode)
check(admit(p,c,['support-A','support-B','green-A'],prompt=unicode['answer_contract']['question'])['scientific_admission']['question_bytes']==len(unicode['answer_contract']['question'].encode()))
# Static answer coverage does not bypass runtime scientific parameter admission.
invalid=copy.deepcopy(intent);invalid['answer_contract']['requirements']=[dict(id='unavailable-field',observable='field_support',bindings=['support-A'])]
invalid['recipes'][0]['goal']['parameters']['fields']=['all.absent']
p,c=prepare(invalid)
run('recipe-proposal','--workspace',manifest,'--intent',p,'--request',save('invalid-field-proposal.json',proposal(c,['support-A'])),'--prompt',question,'--out',root/'invalid-field',fail=True)
check(not (root/'invalid-field').exists())

# Dynamic normalization evidence is retained in the effective native request.
paired=copy.deepcopy(intent)
paired['recipes']=[dict(id='paired',purpose='Declared matched difference',goal=dict(method='expression_contrast',normalization_evidence_field='all.A',parameters=dict(pairs=[dict(id='pair',left='all.A',right='all.B',left_scale='1',right_scale='1')],selection=[],support='joint_reported',normalization=dict(method='Authored identity scales'),output_units='counts',workers=1)))]
paired['answer_contract']['requirements']=[dict(id='difference',observable='paired_field_contrast',bindings=['paired'])]
p,c=prepare(paired)
result=admit(p,c,['paired'],out='paired')
_,leaves=verify_receipts(root/'paired',result)
check([Fraction(r['mean_difference']) for _,cert in leaves for r in cert['result']['features']]==[-3,-4])
check(all(cert['query']['normalization']['evidence'] for _,cert in leaves))
# Existing staged protocol: two preparation branches, explicit binding and native sheaf.
inputs=root/'inputs';shutil.copytree(Path(__file__).resolve().parents[1]/'examples/bioinformatics/staged-assays',inputs)
run('import-assays','--spec',inputs/'import.json','--out',root/'staged-source')
staged_manifest=root/'staged-source/workspace.json'
staged=json.loads((inputs/'intent.json').read_text())
staged['answer_contract']['requirements']=[dict(id='rna-groups',observable='linear_aggregation',bindings=['rna-stage']),dict(id='atac-groups',observable='linear_aggregation',bindings=['atac-stage']),dict(id='interfaces',observable='assay_sections',bindings=['sections']),dict(id='published-support',observable='field_support',bindings=['support'])]
p,c=prepare(staged,staged_manifest)
result=admit(p,c,['support'],out='staged',workspace=staged_manifest,prompt=staged['answer_contract']['question'])
check(result['included_goals']==['rna-stage','atac-stage','combined','sections','support'])
check(len(result['execution']['results'])==9)
_,leaves=verify_receipts(root/'staged',result)
check(len(leaves)==8)
section=next(cert for name,cert in leaves if name=='interfaces')['result']['sections'][0]
check(Fraction(section['residual_quadrance'])==Fraction(1420,51))
check([Fraction(v) for field in section['fields'] for v in field['compatible']]==[3,Fraction(84,17),3,Fraction(191,51),Fraction(382,51),Fraction(443,51)])
# Answer pointers remain valid after the generating source moves and the run is copied.
shutil.move(root/'staged-source',root/'moved-source')
shutil.copytree(root/'staged',root/'detached')
verify_receipts(root/'detached',result)
# No former 16/32 requirement cutoff; all obligations and source-scoped witnesses retained.
many=copy.deepcopy(intent)
many['answer_contract']['requirements']=[dict(id='observation-'+str(i),observable='field_support',bindings=['support-A']) for i in range(600)]
p,c=prepare(many)
check(c['scientific_admission']['native']['cells']==600)
check(admit(p,c,['support-A'])['scientific_admission']['native']['assignment']==[1]*600)
check(not list(root.glob('.gadget-pipeline-*')))
(root/'verification.json').write_text(json.dumps(dict(status='PASS',checks=checks,subset_selections=15,native_partitions=2,staged_steps=9,large_requirements=600),indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=checks)))
