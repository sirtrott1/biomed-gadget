#!/usr/bin/env python3
"""Run both complete-map contracts for an admitted temporal history."""
import argparse,json,subprocess
from pathlib import Path
from verify import verify_pair


def replay(binary, case, query_name, stem):
    case=case.resolve();binary=binary.resolve()
    query=json.loads((case/query_name).read_text());query.pop('map_output',None)
    results=[]
    for mode in ('enumerated','factored'):
        q=dict(query,map_output=mode)
        path=case/(stem+'-'+mode+'-query.json');path.write_text(json.dumps(q,indent=2)+'\n')
        result=subprocess.run([str(binary),'query','--study',str(case/'study.json'),'--query',str(path)],capture_output=True,text=True,check=True)
        cert=json.loads(result.stdout);results.append(cert)
        (case/(stem+'-'+mode+'.json')).write_text(result.stdout)
    receipt=verify_pair(*results)
    (case/(stem+'-expansion.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--case',type=Path,required=True);p.add_argument('--query',default='history-query.json');p.add_argument('--stem',default='maps');a=p.parse_args()
    print(json.dumps(replay(a.binary,a.case,a.query,a.stem)))
