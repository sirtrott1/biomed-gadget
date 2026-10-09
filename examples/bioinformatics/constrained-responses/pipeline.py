#!/usr/bin/env python3
"""Bind the physical-response recipe and eight numerical fact owners."""
import argparse,json
from pathlib import Path
from replay import save,digest,run


def execute(binary,root):
    q=json.loads((root/'inputs/query.json').read_text());workspace=root/'workspace/workspace.json'
    question='Which represented response fields satisfy their physical bounds and declared measured-dose directions, and what objective change does that require?'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='custom',assay='measurements',
                recipes=[dict(id='response-model',purpose=question,goal=dict(method='constrained_response_model',parameters={k:v for k,v in q.items() if k!='operation'}))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,
                requirements=[dict(id='response-model',observable='constrained_response_model',bindings=['response-model'])]))
    save(root/'intent.json',intent);catalog,_=run(binary,'recipe-catalog','--workspace',workspace,'--intent',root/'intent.json');save(root/'recipe-catalog.json',catalog)
    if catalog['scientific_admission']['status']!='satisfied':raise ValueError('Response question not admitted')
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],action='recipes',goals=['response-model'])
    receipt,_=run(binary,'recipe-proposal','--workspace',workspace,'--intent',root/'intent.json','--request',save(root/'proposal.json',proposal),'--prompt',question,'--out',root/'pipeline');save(root/'pipeline-receipt.json',receipt)
    answers=root/'pipeline/answers.json';artifact=json.loads(answers.read_text())['answers'][0]['artifacts'][0]
    observed=json.loads((root/'certificate.json').read_text())['result']
    observed_state=observed['states'][observed['observed_state_index']]['id']
    full=next(i for i,f in enumerate(observed['physical_shape']['results']) if f['state']==observed_state and f['model']=='full' and f['response']==q['responses'][0]['id'])
    selections=[('full-objective',['physical_shape','results',full,'joint_residual_quadrance'],'rational','Joint original-coordinate objective of the first full constrained response'),
                ('constraint-cost',['physical_shape','results',full,'constraint_cost'],'rational','Additional observed quadrance required by its declared physical constraints'),
                ('first-readout',['physical_shape','results',full,'fits',0,'represented','sectors','full',0],'rational','First original specimen readout of the full constrained response'),
                ('coefficient-nullity',['physical_shape','results',full,'fits',0,'coefficient_kernel','columns'],'count','Complete free coefficient dimension of that model'),
                ('active-rank',['physical_shape','results',full,'active_rank'],'count','Independent active constraints in the global optimum witness'),
                ('harmonic-dimension',['hodge_dimensions','harmonic'],'count','Harmonic dimension of the retained source grade'),
                ('worker-slots',['physical_shape','requested_workers'],'count','Declared joined native worker reservation'),
                ('completed-fits',['physical_shape','completed_tiles'],'count','Joined constrained model applications')]
    facts=dict(schema='gadget.scientific.facts.v1',source_digest=catalog['source_digest'],question=question,
               answers=dict(path='pipeline/answers.json',sha256=digest(answers)),facts=[dict(id=name,path=path,value_type=kind,meaning=meaning,
               **{'from':dict(kind='answer',answer='response-model',binding=artifact['binding'],step=artifact['step'],result=0)}) for name,path,kind,meaning in selections])
    save(root/'facts.json',facts)
    full,_=run(binary,'fact-catalog','--workspace',workspace,'--request',root/'facts.json','--artifacts',root)
    model,_=run(binary,'fact-catalog','--workspace',workspace,'--request',root/'facts.json','--artifacts',root,'--view','model')
    save(root/'facts-catalog.json',full);save(root/'facts-model-view.json',model)
    claims=dict(schema='gadget.scientific.claims.proposal.v1',action='claims',source_digest=full['source_digest'],facts_sha256=full['facts_sha256'],question_sha256=full['question_sha256'],
                claims=[dict(id=f['id'],fact=f['id'],comparison='equal',value=f['value']) for f in full['facts']])
    receipt,_=run(binary,'claim-proposal','--workspace',workspace,'--facts',root/'facts.json','--request',save(root/'claims.json',claims),'--artifacts',root);save(root/'claims-receipt.json',receipt)
    if receipt['status']!='matched' or receipt['certificate']['missing']:raise ValueError('Fact ownership check failed')
    return dict(status='PASS',facts=8)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--run',type=Path,required=True);a=p.parse_args();print(json.dumps(execute(a.binary.resolve(),a.run.resolve())))
