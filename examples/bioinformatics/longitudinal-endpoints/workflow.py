#!/usr/bin/env python3
"""Execute two native queries through the existing joined C++ work graph."""
import argparse,json
from pathlib import Path
from replay import save,run
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
querypath=root/'study/history-query.json' if (root/'study/history-query.json').exists() else root/'inputs/query.json';q=json.loads(querypath.read_text());cert=json.loads((root/'certificate.json').read_text());field=q['endpoints'][0]['field'] if 'endpoints' in q else q['lower']
queries=[dict(id='analysis',type='query',query=q),dict(id='original-support',type='query',query=dict(operation='coverage',field=field))]
request=dict(schema='gadget.scientific.workflow.v1',id='reported-longitudinal-question',source_digest=cert['source_digest'],limits=dict(depth=2,nodes=1+len(queries),model_calls=0,parallel_workers=len(queries)),definitions=[],entry=dict(id='joined',type='parallel',workers=len(queries),steps=queries))
out,elapsed=run(a.binary.resolve(),'workflow','--study',root/'study/study.json','--request',save(root/'workflow-request.json',request));assert out['status']=='complete',out.get('error');save(root/'workflow.json',out);save(root/'workflow-timing.json',dict(seconds=elapsed));print(json.dumps(dict(status='PASS',seconds=elapsed)))
