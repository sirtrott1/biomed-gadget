#!/usr/bin/env python3
"""Bind the observed-response recipe and eight numerical fact owners."""
import argparse,json
from pathlib import Path
from replay import save,digest,run


def execute(binary,root):
    q=json.loads((root/'inputs/query.json').read_text());workspace=root/'workspace/workspace.json'
    question='How do declared response models change observed fields, and which measured dose or condition comparisons are supported by these records?'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='custom',assay='measurements',
                recipes=[dict(id='response-model',purpose=question,goal=dict(method='observed_response_model',parameters={k:v for k,v in q.items() if k!='operation'}))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,
                requirements=[dict(id='response-model',observable='observed_response_model',bindings=['response-model'])]))
    save(root/'intent.json',intent);catalog,_=run(binary,'recipe-catalog','--workspace',workspace,'--intent',root/'intent.json');save(root/'recipe-catalog.json',catalog)
    if catalog['scientific_admission']['status']!='satisfied':raise ValueError('Response question not admitted')
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],action='recipes',goals=['response-model'])
    receipt,_=run(binary,'recipe-proposal','--workspace',workspace,'--intent',root/'intent.json','--request',save(root/'proposal.json',proposal),'--prompt',question,'--out',root/'pipeline');save(root/'pipeline-receipt.json',receipt)
    answers=root/'pipeline/answers.json';artifact=json.loads(answers.read_text())['answers'][0]['artifacts'][0]
    selections=[('model-rank',['preparations',0,'models',1,'rank'],'count','Rank of the first full declared model'),
                ('coefficient-nullity',['preparations',0,'models',1,'coefficient_nullity'],'count','Free coefficient dimension of that full model'),
                ('residual-quadrance',['states',0,'fits',1,'residual_quadrance'],'rational','Count-space residual quadrance of the first full fit'),
                ('curl-quadrance',['states',0,'fits',1,'residual','quadrances','curl'],'rational','Curl quadrance of its complete residual readout'),
                ('harmonic-quadrance',['states',0,'fits',1,'residual','quadrances','harmonic'],'rational','Harmonic quadrance of its complete residual readout'),
                ('harmonic-dimension',['hodge_dimensions','harmonic'],'count','Harmonic field dimension on the declared source grade'),
                ('source-preparations',['execution','model_carriers'],'count','Shared model carrier preparations for the response family'),
                ('completed-fits',['execution','completed_tiles'],'count','Joined native response-model applications')]
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
