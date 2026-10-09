#!/usr/bin/env python3
"""Run an authored comparison protocol on a portable copy of an admitted study."""
import argparse
import json
import shutil
from pathlib import Path

from replay import benchmark, digest, run, save


def execute(binary, study, protocol, root, with_benchmark):
    if root.exists():
        raise ValueError('Choose a new output directory')
    declaration = json.loads(study.read_text())
    request = json.loads(protocol.read_text())
    if request['schema'] != 'gadget.scientific.batch.v1':
        raise ValueError('Expected an authored scientific batch protocol')
    # Preserve complete source artifacts, including relative declared paths.
    records = declaration['documents'] + [declaration['native']]
    source_root = study.parent.resolve()
    destination = root/'source'
    destination.mkdir(parents=True)
    for record in records:
        relative = Path(record['path'])
        path = (source_root/relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(source_root) or not path.is_file():
            raise ValueError('Declared study artifact escapes its source root')
        if digest(path) != record['sha256']:
            raise ValueError('Declared source artifact hash changed')
        target = destination/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
    shutil.copyfile(study,destination/'study.json')
    save(root/'protocol.json',request)
    receipts, times, identities = {}, {}, set()
    for job in request['jobs']:
        name = job['id']
        if not isinstance(name,str) or not name or name in identities or '/' in name or '\\' in name or name in ('.','..'):
            raise ValueError('Invalid or duplicate protocol job identity')
        identities.add(name)
        if job['query']['operation'] != 'linear_model_comparison':
            raise ValueError('This comparison protocol expects the declared linear model method')
        query_path = save(root/f'{name}.query.json',job['query'])
        receipts[name], times[name] = run(binary,'query','--study',destination/'study.json','--query',query_path)
        save(root/f'{name}.certificate.json',receipts[name])
    if with_benchmark:
        save(root/'benchmark.json',benchmark(receipts))
    save(root/'timings.json',times)
    manifest = dict(schema='gadget.linear.comparison.reproduction.v1',
                    source_study_sha256=digest(destination/'study.json'),protocol_sha256=digest(root/'protocol.json'),
                    results=[dict(id=name,path=f'{name}.certificate.json',sha256=digest(root/f'{name}.certificate.json')) for name in receipts],
                    input_scope='Complete admitted source packet and authored comparisons; upstream experimental processing remains in the original source provenance')
    save(root/'manifest.json',manifest)
    return receipts


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--study',type=Path,required=True)
    p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--benchmark',action='store_true')
    a=p.parse_args()
    results=execute(a.binary.resolve(),a.study.resolve(),a.protocol.resolve(),a.out.resolve(),a.benchmark)
    print(json.dumps(dict(status='completed',cases=list(results),out=str(a.out.resolve()))))
