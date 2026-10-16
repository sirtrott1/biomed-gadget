#!/usr/bin/env python3
"""Run the reviewed source ablation recipe and bind selected certificate facts."""
import argparse
import json
from pathlib import Path
from replay import save, digest, run


def execute(binary, root):
    q=json.loads((root/'inputs/query.json').read_text());workspace=root/'workspace/workspace.json'
    question='How do the declared source relations and observed field spans determine the surviving field and homology components?'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='custom',assay='measurements',
                recipes=[dict(id='source-ablation',purpose=question,
                              goal=dict(method='source_structure_ablation',parameters={k:v for k,v in q.items() if k!='operation'}))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,
                                     requirements=[dict(id='source-ablation',observable='source_structure_ablation',bindings=['source-ablation'])]))
    save(root/'intent.json',intent)
    catalog,_=run(binary,'recipe-catalog','--workspace',workspace,'--intent',root/'intent.json')
    if catalog['scientific_admission']['status']!='satisfied':raise ValueError('Source-ablation question was not admitted')
    save(root/'recipe-catalog.json',catalog)
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],
                  intent_sha256=catalog['intent_sha256'],action='recipes',goals=['source-ablation'])
    receipt,_=run(binary,'recipe-proposal','--workspace',workspace,'--intent',root/'intent.json',
                  '--request',save(root/'proposal.json',proposal),'--prompt',question,'--out',root/'pipeline')
    save(root/'pipeline-receipt.json',receipt)
    answers=root/'pipeline/answers.json';artifact=json.loads(answers.read_text())['answers'][0]['artifacts'][0]
    result=json.loads((root/'certificate.json').read_text())['result'];grade=result['grade']
    span=next(i for i,a in enumerate(result['ablations']) if a['removed_span_rank'])
    selected=[('original-homology',['original_grades',grade,'betti'],'count','Original rational homology dimension on the declared field grade'),
              ('retained-homology',['ablations',0,'grades',grade,'betti'],'count','Homology dimension after the first declared source deletion'),
              ('inclusion-image',['ablations',0,'grades',grade,'homology_image_rank'],'count','Dimension of retained classes represented in the original ambient homology'),
              ('filled-in-ambient',['ablations',0,'grades',grade,'killed_by_inclusion'],'count','Retained classes whose representatives have an original ambient filling'),
              ('original-harmonic-energy',['original_fields',0,'quadrances','harmonic'],'rational','Original metric quadrance of the first measured field harmonic component'),
              ('retained-harmonic-energy',['ablations',0,'fields',0,'quadrances','harmonic'],'rational','Metric quadrance of that field harmonic component after the first source deletion'),
              ('removed-span-rank',['ablations',span,'removed_span_rank'],'count','Independent dimension of the declared observed span removed from the fields'),
              ('prepared-source-owners',['execution','source_preparations'],'count','Distinct source Hodge owners prepared for the complete declared family')]
    facts=dict(schema='gadget.scientific.facts.v1',source_digest=catalog['source_digest'],question=question,
               answers=dict(path='pipeline/answers.json',sha256=digest(answers)),facts=[
                   dict(id=name,path=path,value_type=kind,meaning=meaning,
                        **{'from':dict(kind='answer',answer='source-ablation',binding=artifact['binding'],step=artifact['step'],result=0)})
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
    if checked['status']!='matched' or checked['certificate']['missing']:raise ValueError('Selected report facts were not covered')
    save(root/'claims-receipt.json',checked)
    return dict(status='PASS',selected_facts=len(selected),scientific_admission='satisfied')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    a=p.parse_args();print(json.dumps(execute(a.binary.resolve(),a.run.resolve())))
