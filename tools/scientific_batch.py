#!/usr/bin/env python3
"""Replay a declared batch through GADGET's admitted scientific plans.

No numerical calculation or model call occurs in this runner. Each native plan
admits every leaf before executing it. Native plan publications own each result;
a native administrative checkpoint retains batch completion and original protocol bytes.
"""
import argparse
import hashlib
import json
import re
import subprocess
import tempfile
import time
from pathlib import Path


def sha(path):
    return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def execute(binary, study, protocol, destination, plan_steps=1, resume=False):
    if plan_steps<1:raise ValueError('Positive caller-selected plan size required')
    protocol_bytes=protocol.read_bytes().decode('utf-8')
    declaration=json.loads(protocol_bytes)
    if set(declaration)!= {'schema','id','question','jobs'} or declaration['schema']!='gadget.scientific.batch.v1':
        raise ValueError('Unsupported batch declaration')
    jobs=declaration['jobs']; names=[j['id'] for j in jobs]
    if not jobs or len(set(names))!=len(names) or any(not isinstance(s,str) or not re.fullmatch(r'[A-Za-z0-9_-]+',s) for s in names):
        raise ValueError('Unique nonempty batch identities required')
    if any(set(j)!= {'id','query','description'} for j in jobs):
        raise ValueError('Undeclared job member')
    if destination.exists() and not resume:
        raise ValueError('New batch destination required')
    def record(path, member=None):
        command=[str(binary),'record-view','--input',str(path)]
        if member is not None:command+=['--member',member]
        return json.loads(subprocess.run(command,capture_output=True,text=True,check=True).stdout)
    def publication_step(job):
        path=destination/job['publication']
        if sha(path)!=job['publication_sha256']:raise ValueError('Changed native plan bytes')
        return record(path,'result')['steps'][job['step']]['certificate']
    previous=None
    if resume:
        checkpoint=record(destination/'manifest.rfcxs')
        previous=checkpoint['manifest'];original=checkpoint['protocol']
        if json.loads(checkpoint['protocol_bytes'])!=original:
            raise ValueError('Checkpoint differs from retained original protocol bytes')
        if previous['status']!='complete' or len(previous['jobs'])!=len(original['jobs']):
            raise ValueError('Resume requires a completed, verified prefix')
        if any(original[k]!=declaration[k] for k in ('schema','id','question')) or original['jobs']!=jobs[:len(original['jobs'])]:
            raise ValueError('A resumed protocol can only append unchanged declared jobs')
        for job,declared in zip(previous['jobs'],original['jobs'],strict=True):
            if job['id']!=declared['id'] or publication_step(job)['query']!=declared['query']:
                raise ValueError('Changed completed query or receipt')
    else:
        destination.mkdir(parents=True)
    manifest=dict(schema='gadget.scientific.batch.receipt.v1',status='running',
                  protocol_sha256=sha(protocol),jobs=[] if previous is None else previous['jobs'].copy(),
                  elapsed_seconds=None,plan_steps=plan_steps)
    started=time.monotonic();started_run=previous is None
    def checkpoint():
        # These are batch-control declarations on an administrative witness.
        # Mathematical source/results stay owned by the referenced native plans.
        with tempfile.TemporaryDirectory(prefix='.checkpoint-',dir=destination) as stage:
            control=Path(stage)/'control.json';staged=Path(stage)/'control.rfcxs'
            write(control,dict(manifest=manifest,protocol=declaration,protocol_bytes=protocol_bytes))
            subprocess.run([str(binary),'native-record','--input',str(control),
                            '--out',str(staged)],
                           capture_output=True,text=True,check=True)
            staged.replace(destination/'manifest.rfcxs')
        write(destination/'manifest.json',manifest)  # Optional compatibility view.
    try:
        audit=subprocess.run([str(binary),'audit','--study',str(study)],capture_output=True,text=True,check=True)
        catalog=json.loads(audit.stdout)
        if previous is not None:
            if previous['source_digest']!=catalog['source_digest']:raise ValueError('Changed admitted study')
            inspection=destination.parent/(destination.name+'-inspect.json');write(inspection,dict(operation='inspect'))
            current=json.loads(subprocess.run([str(binary),'query','--study',str(study),'--query',str(inspection)],capture_output=True,text=True,check=True).stdout)
            inspection.unlink()
            prior=publication_step(previous['jobs'][0])
            if any(current[k]!=prior[k] for k in ('gadget_sources_sha256','rexgraph_headers_sha256','native_sha256')):
                raise ValueError('Changed compiled mathematical/application source')
            prefix=str(len(previous['jobs']))
            (destination/('prior-manifest-'+prefix+'.rfcxs')).write_bytes((destination/'manifest.rfcxs').read_bytes())
            write(destination/('prior-manifest-'+prefix+'.json'),previous)
            write(destination/('prior-protocol-'+prefix+'.json'),original)
            manifest['prior_manifest_sha256']=sha(destination/('prior-manifest-'+prefix+'.rfcxs'))
        started_run=True
        write(destination/'catalog.json',catalog)
        write(destination/'protocol.json',declaration)
        manifest['source_digest']=catalog['source_digest']
        for offset in range(len(manifest['jobs']),len(jobs),plan_steps):
            subset=jobs[offset:offset+plan_steps]
            plan=dict(schema='gadget.scientific.plan.v1',id=declaration['id']+'-'+str(offset),
                      question=declaration['question'],source_digest=catalog['source_digest'],
                      steps=[dict(id=j['id'],query=j['query']) for j in subset])
            request=destination/('plan-'+str(offset)+'.json');write(request,plan)
            begin=time.monotonic()
            native_result=destination/('plan-'+str(offset)+'-receipt.rfcxs')
            process=subprocess.run([str(binary),'plan','--study',str(study),'--request',str(request),'--out',str(native_result)],
                                   capture_output=True,text=True)
            if process.returncode:
                raise RuntimeError(process.stderr.strip() or process.stdout[:4000])
            receipt=record(native_result,'result')
            if receipt!=json.loads(process.stdout):raise RuntimeError('Changed native plan display')
            if receipt['status']!='complete':raise RuntimeError('Incomplete native plan')
            plan_result=destination/('plan-'+str(offset)+'-receipt.json');write(plan_result,receipt)
            for index,(job,step) in enumerate(zip(subset,receipt['steps'],strict=True)):
                if job['id']!=step['id']:raise RuntimeError('Native step identity mismatch')
                path=destination/('job-'+str(offset+index)+'.json');write(path,step['certificate'])
                manifest['jobs'].append(dict(id=job['id'],path=path.name,sha256=sha(path),publication=native_result.name,publication_sha256=sha(native_result),step=index))
            checkpoint()
            print(f"Completed {len(manifest['jobs'])}/{len(jobs)} queries in {time.monotonic()-begin:.2f} seconds",flush=True)
        manifest['status']='complete'
    except Exception as error:
        manifest['status']='failed';manifest['error']=str(error)
        raise
    finally:
        manifest['elapsed_seconds']=time.monotonic()-started+(previous['elapsed_seconds'] if previous else 0)
        if started_run:checkpoint()
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',required=True,type=Path)
    parser.add_argument('--study',required=True,type=Path)
    parser.add_argument('--protocol',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--plan-steps',type=int,default=1,help='Caller-selected leaves per native plan; default 1')
    parser.add_argument('--resume',action='store_true',help='Append to a completed batch after verifying immutable inputs and receipts')
    args=parser.parse_args()
    execute(args.binary.resolve(),args.study.resolve(),args.protocol.resolve(),args.out.resolve(),args.plan_steps,args.resume)
