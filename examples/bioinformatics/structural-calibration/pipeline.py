#!/usr/bin/env python3
"""Bind a relational field question to a native pipeline and selected result facts."""
import argparse
import json
from pathlib import Path

from replay import digest, run, save


def execute(binary, root):
    query=json.loads((root/'inputs/query.json').read_text())
    workspace=root/'workspace/workspace.json'
    question='Which declared structural observables change under complete joint exchange assignments?'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='custom',assay='measurements',
                recipes=[dict(id='structural-calibration',purpose=question,
                    goal=dict(method='structural_reference_calibration',parameters={k:v for k,v in query.items() if k!='operation'}))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,
                    requirements=[dict(id='structural-calibration',observable='structural_reference_calibration',bindings=['structural-calibration'])]))
    alignment=None
    if (root/'alignment-request.json').exists():
        alignment=json.loads((root/'alignment-request.json').read_text())['query']
        intent['recipes'].append(dict(id='condition-alignment',purpose='Compare the declared field-to-condition alignments jointly',
                                      goal=dict(method='structural_reference_calibration',parameters={k:v for k,v in alignment.items() if k!='operation'})))
        intent['answer_contract']['requirements'].append(dict(id='condition-alignment',observable='structural_reference_calibration',bindings=['condition-alignment']))
    save(root/'intent.json',intent)
    catalog,_=run(binary,'recipe-catalog','--workspace',workspace,'--intent',root/'intent.json')
    if catalog['scientific_admission']['status']!='satisfied':
        raise ValueError('The native structural-calibration recipe does not cover the authored question')
    save(root/'recipe-catalog.json',catalog)
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],
                  intent_sha256=catalog['intent_sha256'],action='recipes',goals=[r['id'] for r in intent['recipes']])
    receipt,_=run(binary,'recipe-proposal','--workspace',workspace,'--intent',root/'intent.json',
                 '--request',save(root/'proposal.json',proposal),'--prompt',question,'--out',root/'pipeline')
    save(root/'pipeline-receipt.json',receipt)
    answers=root/'pipeline/answers.json'
    artifact=json.loads(answers.read_text())['answers'][0]['artifacts'][0]
    probe=next(i for i,h in enumerate(query['hypotheses']) if h['statistic']=='probe_fraction')
    selected=[
        ('observed-sector',['hypotheses',0,'observed'],'rational','First declared Hodge-sector quadrance fraction'),
        ('observed-alignment',['hypotheses',probe,'observed'],'rational','Sector quadrance captured by the fixed declared probe span'),
        ('alignment-tail',['hypotheses',probe,'reference_tail_mass'],'rational','Inclusive conditional tail under the complete declared joint exchange law'),
        ('alignment-tie',['hypotheses',probe,'tie_mass'],'rational','Joint reference mass exactly equal to the observed alignment'),
        ('alignment-family-tail',['hypotheses',probe,'max_statistic_tail_mass'],'rational','Inclusive maximum-statistic tail for the declared family'),
        ('family-global-tail',['joint_family','global_tail_mass'],'rational','Inclusive global tail of the declared family maximum'),
        ('reference-states',['cost','reference_states'],'count','Complete number of coordinate assignments in the declared reference law'),
        ('harmonic-rank',['frames','harmonic_rank'],'count','Dimension of the retained harmonic field sector')]
    owners={name:artifact for name,_,_,_ in selected}
    if alignment is not None:
        alignment_artifact=next(a for a in json.loads(answers.read_text())['answers'] if a['id']=='condition-alignment')['artifacts'][0]
        target=next(i for i,h in enumerate(alignment['hypotheses']) if h['statistic']=='probe_fraction')
        selected=[(name,['hypotheses',target,path[-1]] if name in ('alignment-tail','alignment-tie','alignment-family-tail') else path,kind,meaning) for name,path,kind,meaning in selected]
        for name in ('alignment-tail','alignment-tie','alignment-family-tail','family-global-tail'):owners[name]=alignment_artifact
    facts=dict(schema='gadget.scientific.facts.v1',source_digest=catalog['source_digest'],question=question,
               answers=dict(path='pipeline/answers.json',sha256=digest(answers)),facts=[
                   dict(id=name,path=path,value_type=kind,meaning=meaning,
                        **{'from':dict(kind='answer',answer='condition-alignment' if alignment is not None and name in ('alignment-tail','alignment-tie','alignment-family-tail','family-global-tail') else 'structural-calibration',binding=owners[name]['binding'],step=owners[name]['step'],result=0)})
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
