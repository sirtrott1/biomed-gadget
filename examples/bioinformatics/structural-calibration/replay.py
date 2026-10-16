#!/usr/bin/env python3
"""Stage source declared exchange laws and calibrate native structural readings."""
import argparse
import csv
import copy
import importlib.util
import json
import shutil
from fractions import Fraction as Q
from pathlib import Path

# Reuse the existing portable source staging protocol; C++ owns all analyses.
path=Path(__file__).resolve().parents[1]/'relational-ensembles/replay.py'
spec=importlib.util.spec_from_file_location('retained_source',path)
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
save,digest,run=source.save,source.digest,source.run


def query(fields,blocks,hypotheses,workers=1,profile='statistics'):
    return dict(operation='structural_reference_calibration',context='all',fields=fields,
                law=dict(kind='blocked_permutations',blocks=blocks,
                         null_hypothesis='The complete field axes are jointly exchangeable within the declared blocks, conditional on the retained design',
                         population='The identified source coordinates of the declared experiment',
                         assumptions=['Each within-block coordinate assignment has equal conditional mass',
                                      'One assignment acts jointly on every selected field axis',
                                      'The probe declarations and source complex are fixed before viewing reference scores'],
                         derivation='Uniform counting measure on the product of the declared block permutation groups',evidence=[]),
                hypotheses=hypotheses,workers=workers,execution='exhaustive',reference_output=profile,alpha='1/20')


def authored(root,workers=1):
    readings=dict(A=[6,2,4,8],B=[2,4,8,6],A_shift=[7,1,4,8],B_shift=[1,5,8,6],
                  probe=[1,-1,0,0],same_probe=[2,-2,0,0],zero=[0,0,0,0])
    source.authored(root/'inputs',field_values=readings)
    blocks=[dict(id='interface-exchange',relation='relation:exchange',cells=[f'relation:assertion-{i}' for i in range(4)])]
    hypotheses=[dict(id='gradient-A',field='all.A',statistic='sector_fraction',sector='gradient'),
                dict(id='curl-A',field='all.A',statistic='sector_fraction',sector='curl'),
                dict(id='harmonic-A',field='all.A',statistic='sector_fraction',sector='harmonic'),
                dict(id='harmonic-probe-A',field='all.A',statistic='probe_fraction',sector='harmonic',probes=['all.probe','all.same_probe']),
                dict(id='green-change-A',field='all.A',statistic='green_change_fraction',**{'lambda':'1'}),
                dict(id='harmonic-probe-B',field='all.B',statistic='probe_fraction',sector='harmonic',probes=['all.probe'])]
    q=query(['all.A','all.B'],blocks,hypotheses,workers,'full_fields')
    q['uncertainty']=dict(kind='declared_joint_scenarios',observed_scenario='original',
                         population='Two explicitly authored paired measurement scenarios',
                         weighting='Equal sensitivity weights; these are not an experimental null law or estimated measurement probabilities',
                         derivation='Original fields and a declared joint perturbation retaining both source axes',
                         assumptions=['The two authored scenarios are the entire selected sensitivity family'],evidence=[],
                         scenarios=[dict(id='original',weight='1/2',fields=['all.A','all.B']),
                                    dict(id='perturbed',weight='1/2',fields=['all.A_shift','all.B_shift'])])
    return q


def expression(root,packet,panel,workers=1):
    original=root/'source';original.mkdir(parents=True)
    receipt=json.loads((packet/'source.json').read_text())
    required=list(dict.fromkeys(['source.json','counts.tsv','samples.tsv']+list(receipt['artifacts'])))
    for name in required:
        path=(packet/name).resolve()
        if not path.is_relative_to(packet.resolve()) or not path.is_file():raise ValueError('Invalid retained source artifact')
        if name in receipt['artifacts'] and digest(path).removeprefix('sha256:')!=receipt['artifacts'][name]:raise ValueError('Original source hash changed')
        shutil.copyfile(path,original/name)
    selected=json.loads(panel.read_text());genes=[g['id'] for g in selected]
    with (original/'counts.tsv').open() as handle:
        rows=csv.reader(handle,delimiter='\t');samples=next(rows)[1:]
        counts={r[0]:r[1:] for r in rows if r[0] in genes}
    if set(counts)!=set(genes) or len(set(genes))!=len(genes):raise ValueError('Declare unique observed genes')
    metadata={r['sample']:r for r in csv.DictReader((original/'samples.tsv').open(),delimiter='\t')}
    relations=[];blocks=[]
    for donor in dict.fromkeys(metadata[s]['cell'] for s in samples):
        untreated=[s for s in samples if metadata[s]['cell']==donor and metadata[s]['dex']=='untrt']
        treated=[s for s in samples if metadata[s]['cell']==donor and metadata[s]['dex']=='trt']
        if len(untreated)!=1 or len(treated)!=1:raise ValueError('This experiment needs one recorded pair per donor')
        members=[untreated[0],treated[0]]
        relations.append(dict(id='pair-'+donor,grade=1,meaning='Recorded within-donor untreated-to-treated comparison',participants=[source.participant(members[0],-1),source.participant(members[1],1)]))
        blocks.append(dict(id='donor-'+donor,relation='relation:pair-'+donor,cells=members))
    for condition in ('untrt','trt'):
        members=[s for s in samples if metadata[s]['dex']==condition]
        relations.append(dict(id='condition-'+condition,grade=1,meaning='Recorded condition membership with equal shares',participants=[source.participant(s,Q(1,len(members))) for s in members]))
    spec=source.stage(root/'inputs',genes,{s:[counts[g][i] for g in genes] for i,s in enumerate(samples)},relations,[],
                      units='counts',quantity='read_count',namespace='993401',question='How do raw-count Hodge components align with the declared condition contrast under joint within-donor exchanges?')
    observations=save(root/'inputs/observations.tsv','id\tdex\tcell\ttreatment_contrast\n'+''.join(f'{s}\t{metadata[s]["dex"]}\t{metadata[s]["cell"]}\t{1 if metadata[s]["dex"]=="trt" else -1}\n' for s in samples))
    spec['assays'][0]['observations']['sha256']=digest(observations)
    spec['assays'][0]['carrier']='observations'
    spec['assays'][0]['partition_size']=len(samples)
    spec['assays'][0]['observation_fields']=[dict(id='treatment_contrast',column='treatment_contrast',kind='exact_value',
                                               quantity='declared_condition_contrast',units='dimensionless',
                                               method='Declared oriented contrast: treated +1, untreated -1; recorded labels retained',missing_tokens=['NA'])]
    save(root/'inputs/import.json',spec)
    hypotheses=[]
    for g in genes:
        for sector in ('curl','harmonic'):
            hypotheses.append(dict(id=g+'-'+sector,field='all.'+g,statistic='sector_fraction',sector=sector))
        hypotheses.append(dict(id=g+'-condition',field='all.'+g,statistic='probe_fraction',sector='curl',probes=['all.treatment_contrast']))
    q=query(['all.'+g for g in genes],blocks,hypotheses,workers,'full_fields')
    q['law']['population']='The eight recorded airway specimens; six preselected raw count fields share each donor exchange'
    save(root/'selected-features.json',selected)
    save(root/'source-witness.json',dict(files={name:digest(original/name) for name in required},selection=selected,
                                       normalization='none',relation_semantics='Recorded donor and treatment-condition membership, not a regulatory mechanism'))
    return q


def bind_evidence(q,workspace):
    paths=list(workspace.rglob('assay.json'))
    if len(paths)!=1:raise ValueError('This protocol needs one source assay')
    raw=paths[0].read_text();evidence=[dict(document='assay_declaration',start=0,end=len(raw.encode()),text=raw)]
    q['law']['evidence']=evidence
    if 'uncertainty' in q:q['uncertainty']['evidence']=evidence


def execute(binary,root,q):
    imported,seconds=run(binary,'import-assays','--spec',root/'inputs/import.json','--out',root/'workspace')
    save(root/'import-result.json',imported)
    bind_evidence(q,root/'workspace');save(root/'inputs/query.json',q)
    request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=q)
    result,elapsed=run(binary,'assay-query','--workspace',root/'workspace/workspace.json','--request',save(root/'request.json',request))
    save(root/'result.json',result);save(root/'certificate.json',result['result']['answers'][0]['certificate'])
    if (root/'selected-features.json').exists():
        alignment=copy.deepcopy(q)
        alignment['hypotheses']=[h for h in alignment['hypotheses'] if h['statistic']=='probe_fraction']
        request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=alignment)
        answer,alignment_seconds=run(binary,'assay-query','--workspace',root/'workspace/workspace.json','--request',save(root/'alignment-request.json',request))
        save(root/'alignment-certificate.json',answer['result']['answers'][0]['certificate'])
        save(root/'alignment-timing.json',dict(seconds=alignment_seconds))
    save(root/'timing.json',dict(import_seconds=seconds,query_seconds=elapsed))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--panel',type=Path);p.add_argument('--workers',type=int,default=1)
    a=p.parse_args();root=a.out.resolve();root.mkdir(parents=True,exist_ok=False)
    if a.source:
        if not a.panel:p.error('--source requires --panel')
        q=expression(root,a.source.resolve(),a.panel.resolve(),a.workers)
    else:q=authored(root,a.workers)
    execute(a.binary.resolve(),root,q);print(json.dumps(dict(status='completed',output=str(root))))
