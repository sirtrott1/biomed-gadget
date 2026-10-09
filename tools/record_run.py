#!/usr/bin/env python3
"""Record a local qualification process without mixing compilation and runtime."""
import argparse
import json
import resource
import subprocess
import time
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stem',type=Path,required=True)
    p.add_argument('--input',type=Path)
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args(); command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:p.error('a command is required')
    a.stem.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    with (a.stem.parent/(a.stem.name+'-stdout.json')).open('wb') as out, (a.stem.parent/(a.stem.name+'-stderr.log')).open('wb') as err:
        run=subprocess.run(command,input=a.input.read_bytes() if a.input else None,stdout=out,stderr=err)
    receipt=dict(command=command,returncode=run.returncode,elapsed_seconds=time.monotonic()-start,
                 peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (a.stem.parent/(a.stem.name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt));raise SystemExit(run.returncode)

if __name__=='__main__':main()
