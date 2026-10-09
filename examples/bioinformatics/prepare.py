#!/usr/bin/env python3
"""Package caller-declared axes and matrix bytes into a portable import specification."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True);p.add_argument('--features',type=Path,required=True);p.add_argument('--observations',type=Path,required=True)
    p.add_argument('--format',choices=['dense_tsv','matrix_market','h5ad'],required=True)
    p.add_argument('--id',required=True);p.add_argument('--namespace',required=True)
    p.add_argument('--kind',required=True);p.add_argument('--quantity',required=True);p.add_argument('--units',required=True)
    p.add_argument('--method',required=True);p.add_argument('--question',required=True);p.add_argument('--population',required=True)
    p.add_argument('--carrier',choices=['features','observations'],default='features')
    p.add_argument('--channel',choices=['G','E','T','observation'],default='observation')
    p.add_argument('--implicit-entries',choices=['measured_zero','unobserved'],required=True)
    p.add_argument('--arithmetic',choices=['exact_decimal','exact_rational','ieee754_exact'],required=True)
    p.add_argument('--missing-token',action='append',default=[]);p.add_argument('--nan-policy',choices=['reject','unobserved'])
    p.add_argument('--relations',type=Path);p.add_argument('--limitation',action='append',default=[])
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():p.error('output already exists')
    a.out.mkdir(parents=True)
    try:
        def retain(path,name):
            target=a.out/name;shutil.copyfile(path,target)
            with target.open('rb') as stream:hash=hashlib.file_digest(stream,'sha256').hexdigest()
            return dict(path=name,sha256='sha256:'+hash)
        assay=dict(id=a.id,kind=a.kind,format=a.format,source=retain(a.matrix,'matrix.source'),features=retain(a.features,'features.tsv'),observations=retain(a.observations,'observations.tsv'),source_namespace=a.namespace,carrier=a.carrier,quantity=a.quantity,units=a.units,method=a.method,channel=a.channel,arithmetic=a.arithmetic,implicit_entries=a.implicit_entries,missing_tokens=a.missing_token)
        if a.relations:assay['relations']=retain(a.relations,'relations.json')
        if a.nan_policy:assay['nan_policy']=a.nan_policy
        spec=dict(schema='gadget.assays.import.v1',id=a.id,design=dict(question=a.question,population=a.population,selection='all rows and observations in the supplied axes',limitations=a.limitation),assays=[assay],resources=[],links=[])
        (a.out/'import.json').write_text(json.dumps(spec,indent=2)+'\n')
        print(a.out/'import.json')
    except BaseException:
        shutil.rmtree(a.out);raise

if __name__=='__main__':main()
