#!/usr/bin/env python3
"""Measure complete native actions and two lossless cause word encodings."""
import argparse
import copy
import hashlib
import importlib.util
import json
import statistics
import subprocess
import sys
from pathlib import Path

_path=Path(__file__).resolve().parents[1]/'source-admission/measure.py'
_spec=importlib.util.spec_from_file_location('shared_process_measurement',_path)
recorder=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(recorder)


def measure(session,examples,out,repeat=3,iterations=3):
    if repeat<=0 or iterations<=0:raise ValueError('positive measurement counts required')
    out.mkdir(parents=True,exist_ok=False);actions=[]
    for name in ['authored/algebra','authored/sparse-257','clinical']:
        for mode in ['reuse','ephemeral']:
            samples=[]
            for i in range(repeat):
                stem=name.replace('/','-')+'-'+mode+'-'+str(i)
                target=out/(stem+'.json')
                with (out/(stem+'.stderr')).open('w') as error:
                    worker=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--record-child',
                        '--session',str(session),'--examples',str(examples/name),'--out',str(target),
                        '--iterations',str(iterations),'--mode',mode],stdout=subprocess.PIPE,stderr=error,text=True,check=True)
                observed=json.loads(worker.stdout);native=json.loads(target.read_text())
                samples.append(dict(**observed,native=native))
                assert int(native['preparations'])==1 if mode=='reuse' else int(native['preparations'])==iterations+1
                assert int(native['retained_plans'])==(1 if mode=='reuse' else 0)
            actions.append(dict(case=name,mode=mode,runs=repeat,first_ms=statistics.median(s['native']['first_call_microseconds']/1000 for s in samples),
                subsequent_ms=statistics.median(s['native']['reused_calls_microseconds']/iterations/1000 for s in samples),
                peak_rss_kib=statistics.median(s['peak_rss_kib'] for s in samples),preparations=samples[0]['native']['preparations'],
                retained_plans=samples[0]['native']['retained_plans'],numeric_entries=samples[0]['native']['numeric_entries'],samples=samples))
    c=json.loads((examples/'authored/sparse-257/result.json').read_text())['result']['causes']
    expanded={k:copy.deepcopy(v) for k,v in c.items() if k not in ['fields','word']}
    expanded['word']=[dict(field=v['field'],**copy.deepcopy(c['fields'][v['dictionary_index']])) for v in c['word']]
    restored=[dict(field=v['field'],**c['fields'][v['dictionary_index']]) for v in c['word']]
    assert expanded['word']==restored
    encode=lambda x:json.dumps(x,ensure_ascii=False,separators=(',',':')).encode()
    report=dict(status='PASS',scope='Sequential whole-process wait4 RSS; action timings include full certificate construction and serialization. No kernel-only or universal scaling claim.',
        session_sha256='sha256:'+hashlib.sha256(session.read_bytes()).hexdigest(),actions=actions,
        serialization=dict(dictionary_bytes=len(encode(c)),expanded_bytes=len(encode(expanded)),complete_equivalent_value_words=True,
            original_coordinates=len(c['fields'][0]['values']),ordered_slots=len(c['word']),full_Gram_extents=[c['pairings']['rows'],c['pairings']['columns']]))
    (out/'performance.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',required=True,type=Path);p.add_argument('--examples',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--repeat',type=int,default=3);p.add_argument('--iterations',type=int,default=3);p.add_argument('--mode',choices=['reuse','ephemeral'],default='reuse');p.add_argument('--record-child',action='store_true',help=argparse.SUPPRESS);a=p.parse_args()
    if a.record_child:
        # A small independent recorder must not inherit the orchestration heap.
        with a.out.open('w') as target:
            print(json.dumps(recorder.measure_command([a.session,a.examples,str(a.iterations),a.mode],stdout=target)))
    else:
        print(json.dumps(measure(a.session.resolve(),a.examples.resolve(),a.out.resolve(),a.repeat,a.iterations),indent=2))
