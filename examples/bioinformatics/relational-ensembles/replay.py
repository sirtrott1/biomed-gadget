#!/usr/bin/env python3
"""Stage retained relational fields and run the general C++ ensemble method."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import time
from fractions import Fraction as Q
from pathlib import Path


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value if isinstance(value, str) else json.dumps(value, indent=2)+'\n')
    return path


def digest(path):
    return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary, *args):
    start=time.perf_counter()
    result=subprocess.run([str(binary),*map(str,args)],text=True,capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr)
    return json.loads(result.stdout),time.perf_counter()-start


def stage(root, ids, measurements, relations, relation_fields, *, units='dimensionless',
          quantity='declared_measurement', namespace='993101', question=None):
    if root.exists():
        raise ValueError('Choose a new input directory')
    root.mkdir(parents=True)
    samples=list(measurements)
    paths={
        'values':save(root/'measurements.tsv','id\t'+'\t'.join(samples)+'\n'+
                      ''.join(name+'\t'+'\t'.join(str(measurements[s][i]) for s in samples)+'\n'
                              for i,name in enumerate(ids))),
        'features':save(root/'features.tsv','id\n'+''.join(name+'\n' for name in ids)),
        'observations':save(root/'observations.tsv','id\n'+''.join(name+'\n' for name in samples)),
        'relations':save(root/'relations.json',dict(relations=relations,fields=relation_fields))}
    def ref(path):return dict(path=path.name,sha256=digest(path))
    question=question or 'How do the complete declared fields inhabit the retained relational complex?'
    spec=dict(schema='gadget.assays.import.v1',id='relational-field-ensemble',
              design=dict(question=question,population='Declared identified measurements',selection='Original fields and retained relation declarations',
                          limitations=['Structural readings refer to the declared measurement relations']),
              assays=[dict(id='measurements',kind='expression',format='dense_tsv',source=ref(paths['values']),
                           features=ref(paths['features']),observations=ref(paths['observations']),relations=ref(paths['relations']),
                           source_namespace=namespace,carrier='features',quantity=quantity,units=units,channel='E',
                           method='Original declared scalar measurement fields',arithmetic='exact_rational',
                           implicit_entries='unobserved',missing_tokens=['NA'],partition_size=len(ids))],resources=[],links=[])
    save(root/'import.json',spec)
    return spec


def participant(cell, coefficient, exists=True):
    value=Q(coefficient)
    return dict(cell=cell,exists=exists,negative=value<0,magnitude=str(abs(value)))


def query(fields, *, weights=None, groups=None, comparisons=None, workers=2):
    result=dict(operation='relational_field_ensemble',context='all',fields=['all.'+s for s in fields],
                weights=weights or [str(Q(1,len(fields)))]*len(fields),workers=workers,
                interpretation=dict(kind='observed_fields',population='The explicitly selected retained measurement fields',
                                    weighting='Declared positive field weights; conditional weights within each selected group',evidence=[]),
                green_parameters=['1','3/2'])
    if groups is not None:result['groups']=groups
    if comparisons is not None:result['comparisons']=comparisons
    return result


def authored(root, *, filled=True, arity=4, field_values=None):
    ids=[f'feature-{i}' for i in range(arity)]
    branch=[participant(ids[0],'-1')]+[participant(name,str(Q(1,arity-1))) for name in ids[1:]]
    # These extra source flags cancel or are inactive; they remain occurrences.
    branch += [participant(ids[1],'1'),participant(ids[1],'-1'),participant(ids[2],'7',False)]
    relations=[dict(id=f'assertion-{i}',grade=1,meaning='Separate repeated assertion of the declared multi-participant assay interface',participants=branch) for i in range(4)]
    if filled:
        relations.append(dict(id='exchange',grade=2,meaning='Authored zero-boundary comparison of four repeated interface assertions',
                              participants=[participant(f'assertion-{i}',1 if i<2 else -1) for i in range(4)]))
    readings=field_values or dict(A=[6,2,4,8],B=[2,4,8,6],C=[3,7,5,1],zero=[0,0,0,0],same_A=[6,2,4,8])
    relation_fields=[dict(id=name,grade=1,channel='E',quantity='declared_measurement',units='dimensionless',
                         method='Authored scalar fields on retained primary assertions',
                         values=[dict(relation=f'assertion-{i}',value=str(v)) for i,v in enumerate(values) if v is not None]) for name,values in readings.items()]
    if filled:
        relation_fields += [dict(id='top_'+name,grade=2,channel='E',quantity='declared_measurement',units='dimensionless',
                                method='Authored scalar fields on retained primary assertions',values=[dict(relation='exchange',value=str(i))]) for i,name in enumerate(('A','B','C'))]
    stage(root,ids,dict(object_A=list(range(arity)),object_B=list(range(1,arity+1))),relations,relation_fields)
    q=query(['A','B','C'],groups=[dict(id='first',fields=['all.A','all.B']),dict(id='second',fields=['all.B','all.C'])],
            comparisons=[dict(id='overlapping-groups',left='first',right='second')])
    save(root/'query.json',q)
    return q


def expression(root, source, panel):
    root.mkdir(parents=True,exist_ok=False)
    original=root/'source';original.mkdir()
    receipt=json.loads((source/'source.json').read_text())
    required=list(dict.fromkeys(['source.json','counts.tsv','samples.tsv']+list(receipt['artifacts'])))
    for name in required:
        path=(source/name).resolve()
        if not path.is_relative_to(source.resolve()) or not path.is_file():raise ValueError('Invalid retained source artifact')
        if name in receipt['artifacts'] and digest(path).removeprefix('sha256:')!=receipt['artifacts'][name]:raise ValueError('Original source hash changed')
        shutil.copyfile(path,original/name)
    selected=json.loads(panel.read_text());names=[x['id'] for x in selected]
    if not names or len(set(names))!=len(names):raise ValueError('Declare unique selected features')
    with (original/'counts.tsv').open() as handle:
        rows=csv.reader(handle,delimiter='\t');header=next(rows)
        values={row[0]:row[1:] for row in rows if row[0] in names}
    if set(values)!=set(names):raise ValueError('Selected gene missing from original source')
    observations=list(csv.DictReader((original/'samples.tsv').open(),delimiter='\t'))
    by_id={row['sample']:row for row in observations}
    ids=header[1:]
    relations=[]
    for block in dict.fromkeys(by_id[name]['cell'] for name in ids):
        members=[name for name in ids if by_id[name]['cell']==block]
        untreated=[name for name in members if by_id[name]['dex']=='untrt']
        treated=[name for name in members if by_id[name]['dex']=='trt']
        if len(untreated)!=1 or len(treated)!=1:raise ValueError('Expected documented specimen pairs')
        relations.append(dict(id='pair-'+block,grade=1,meaning='Recorded donor pairing, oriented untreated to treated for comparison',participants=[participant(untreated[0],-1),participant(treated[0],1)]))
    for condition in ('untrt','trt'):
        members=[name for name in ids if by_id[name]['dex']==condition]
        relations.append(dict(id='condition-'+condition,grade=1,meaning='Recorded condition membership with declared equal positive shares',participants=[participant(name,Q(1,len(members))) for name in members]))
    stage(root/'inputs',ids,{name:values[name] for name in names},relations,[],units='counts',quantity='read_count',namespace='993201',
          question='How do the selected original gene-count fields vary in the declared donor-pair and condition-membership complex?')
    q=query(names)
    q['interpretation']['population']='The selected raw gene-count fields on eight recorded airway specimens'
    save(root/'inputs/query.json',q)
    save(root/'selected-features.json',selected)
    save(root/'source-witness.json',dict(files={name:digest(original/name) for name in required},
                                       selection=selected,normalization='none',
                                       relation_semantics='Specimen pairing and condition membership, not a regulatory-network assertion'))
    return q


def bind_evidence(q, workspace):
    paths=list(workspace.rglob('assay.json'))
    if len(paths)!=1:raise ValueError('This example requires one complete source carrier')
    raw=paths[0].read_text()
    q['interpretation']['evidence']=[dict(document='assay_declaration',start=0,end=len(raw.encode()),text=raw)]


def execute(binary, root, q):
    imported,import_seconds=run(binary,'import-assays','--spec',root/'inputs/import.json','--out',root/'workspace')
    save(root/'import-result.json',imported)
    bind_evidence(q,root/'workspace');save(root/'inputs/query.json',q)
    request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=q)
    result,seconds=run(binary,'assay-query','--workspace',root/'workspace/workspace.json','--request',save(root/'request.json',request))
    save(root/'result.json',result);save(root/'timing.json',dict(import_seconds=import_seconds,query_seconds=seconds))
    save(root/'certificate.json',result['result']['answers'][0]['certificate'])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--panel',type=Path)
    p.add_argument('--workers',type=int,default=1,help='Explicit native CPU worker budget')
    p.add_argument('--unfilled',action='store_true',help='Run the authored example without its declared grade-two filling')
    a=p.parse_args();binary=a.binary.resolve();root=a.out.resolve()
    if a.source:
        if not a.panel:p.error('--source requires a declared --panel')
        q=expression(root,a.source.resolve(),a.panel.resolve())
    else:
        root.mkdir(parents=True,exist_ok=False);q=authored(root/'inputs',filled=not a.unfilled)
    q['workers']=a.workers
    execute(binary,root,q)
    print(json.dumps(dict(status='completed',output=str(root))))


if __name__=='__main__':main()
