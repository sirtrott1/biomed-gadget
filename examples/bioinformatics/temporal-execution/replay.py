#!/usr/bin/env python3
"""Reproduce owning historical cutoffs and a declared coupled observation law."""
import argparse,importlib.util,json
from pathlib import Path

def load(directory,name):
    path=Path(__file__).resolve().parents[1]/directory/(name+'.py');spec=importlib.util.spec_from_file_location(directory+'_'+name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
history=load('versioned-history','replay');intervals=load('longitudinal-endpoints','replay')

def joint(binary,root):
    q,out=intervals.endpoints(binary,root)
    values=[['2','5','1','4','7','10','21/2','12'],['2','20','5','8','7',None,'43/4','12'],['2',None,'1','4','7','10','21/2','12']]
    coordinates=out['result']['coordinate_carrier'];ids=[c['id'] for c in coordinates]
    q['joint_family']=dict(meaning='Three authored complete event-time assignments sharing one finite parameter across original records, evaluation times and sectors',evidence=q['evidence'],states=[dict(id='state-'+str(s),weight=w,values=[dict(entity=id,event_time=v) for id,v in zip(ids,row)]) for s,(row,w) in enumerate(zip(values,['1/4','1/4','1/2']))],constraints=[dict(id='shared-event-choice',entities=[ids[1],ids[2]],allowed=[[row[1],row[2]] for row in values],meaning='Only the three jointly authored choices are permitted on these original records',evidence=q['evidence'])])
    query=intervals.save(root/'joint-query.json',q);result,elapsed=intervals.run(binary,'query','--study',root/'study/study.json','--query',query)
    intervals.save(root/'joint.json',result);intervals.save(root/'joint-timing.json',dict(seconds=elapsed))
    steps=[dict(id='joint-observations',type='query',query=q),dict(id='source-catalog',type='query',query=dict(operation='inspect'))]
    request=dict(schema='gadget.scientific.workflow.v1',id='coupled-endpoint-family',source_digest=result['source_digest'],limits=dict(depth=2,nodes=1+len(steps),model_calls=0,parallel_workers=len(steps)),definitions=[],entry=dict(id='joined',type='parallel',workers=len(steps),steps=steps))
    workflow,elapsed=intervals.run(binary,'workflow','--study',root/'study/study.json','--request',intervals.save(root/'workflow-request.json',request));assert workflow['status']=='complete',workflow.get('error');intervals.save(root/'workflow.json',workflow);intervals.save(root/'workflow-timing.json',dict(seconds=elapsed))
    return q,result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    history.authored(a.binary.resolve(),a.fixture.resolve(),a.out.resolve()/'native',True);history.workflow(a.binary.resolve(),a.out.resolve()/'native');joint(a.binary.resolve(),a.out.resolve()/'coupled')
