#!/usr/bin/env python3
"""Bind a relational field question to a native pipeline and selected result facts."""
import argparse
import json
from pathlib import Path

from replay import digest, run, save


def execute(binary, root):
    query=json.loads((root/'inputs/query.json').read_text())
    workspace=root/'workspace/workspace.json'
    question='How are the complete declared field variations distributed across the relational complex?'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='custom',assay='measurements',
                recipes=[dict(id='relational-ensemble',purpose=question,
                    goal=dict(method='relational_field_ensemble',parameters={k:v for k,v in query.items() if k!='operation'}))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,
                    requirements=[dict(id='relational-fields',observable='relational_field_ensemble',bindings=['relational-ensemble'])]))
    save(root/'intent.json',intent)
    catalog,_=run(binary,'recipe-catalog','--workspace',workspace,'--intent',root/'intent.json')
    if catalog['scientific_admission']['status']!='satisfied':
        raise ValueError('The native relational-fields recipe does not cover the authored question')
    save(root/'recipe-catalog.json',catalog)
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],
                  intent_sha256=catalog['intent_sha256'],action='recipes',goals=['relational-ensemble'])
    receipt,_=run(binary,'recipe-proposal','--workspace',workspace,'--intent',root/'intent.json',
                 '--request',save(root/'proposal.json',proposal),'--prompt',question,'--out',root/'pipeline')
    save(root/'pipeline-receipt.json',receipt)
    answers=root/'pipeline/answers.json'
    artifact=json.loads(answers.read_text())['answers'][0]['artifacts'][0]
    selected=[
        ('full-spread',['groups',0,'ensembles','full','variance'],'rational','Weighted spread of the complete fields in the declared source metric'),
        ('full-dimension',['groups',0,'ensembles','full','effective_dimension_2'],'rational','Order-two effective dimension; no field coordinate is discarded'),
        ('gradient-spread',['groups',0,'ensembles','gradient','variance'],'rational','Weighted spread of the complete gradient components'),
        ('curl-spread',['groups',0,'ensembles','curl','variance'],'rational','Weighted spread of the complete curl components'),
        ('harmonic-spread',['groups',0,'ensembles','harmonic','variance'],'rational','Weighted spread of the complete harmonic components'),
        ('gradient-rank',['frames','gradient','rank'],'count','Dimension of the exact gradient frame for this source and metric'),
        ('curl-rank',['frames','curl','rank'],'count','Dimension of the declared upper-boundary image'),
        ('harmonic-rank',['frames','harmonic','rank'],'count','Dimension of the simultaneous boundary and upper-adjoint kernel')]
    facts=dict(schema='gadget.scientific.facts.v1',source_digest=catalog['source_digest'],question=question,
               answers=dict(path='pipeline/answers.json',sha256=digest(answers)),facts=[
                   dict(id=name,path=path,value_type=kind,meaning=meaning,
                        **{'from':dict(kind='answer',answer='relational-fields',binding=artifact['binding'],step=artifact['step'],result=0)})
                   for name,path,kind,meaning in selected])
    save(root/'facts.json',facts)
    full,_=run(binary,'fact-catalog','--workspace',workspace,'--request',root/'facts.json','--artifacts',root)
    view,_=run(binary,'fact-catalog','--workspace',workspace,'--request',root/'facts.json','--artifacts',root,'--view','model')
    save(root/'facts-catalog.json',full);save(root/'facts-model-view.json',view)
    claims=dict(schema='gadget.scientific.claims.proposal.v1',action='claims',source_digest=full['source_digest'],
                facts_sha256=full['facts_sha256'],question_sha256=full['question_sha256'],
                claims=[dict(id=f['id'],fact=f['id'],comparison='equal',value=f['value']) for f in full['facts']])
    checked,_=run(binary,'claim-proposal','--workspace',workspace,'--facts',root/'facts.json',
                  '--request',save(root/'claims.json',claims),'--artifacts',root)
    if checked['status']!='matched' or checked['certificate']['missing']:
        raise ValueError('The authored literal report did not cover the selected result facts')
    save(root/'claims-receipt.json',checked)
    return dict(status='PASS',selected_facts=len(selected),scientific_admission='satisfied')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    a=p.parse_args();print(json.dumps(execute(a.binary.resolve(),a.run.resolve())))
