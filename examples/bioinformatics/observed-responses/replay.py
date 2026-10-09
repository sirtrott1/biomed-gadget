#!/usr/bin/env python3
"""Reproduce observed count/exposure, ordered-category and measured dose models."""
import argparse, csv, importlib.util, itertools, json, shutil
from fractions import Fraction as Q
from pathlib import Path

p=Path(__file__).resolve().parents[1]/'relational-ensembles/replay.py'
spec=importlib.util.spec_from_file_location('response_source',p)
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
save,digest,run=source.save,source.digest,source.run


def covariate(name, quantity, units='dimensionless', kind='exact_value', level=None):
    d=dict(id=name,column=name,kind=kind,quantity=quantity,units=units,
           method='Original recorded observation metadata',missing_tokens=['NA'])
    if level is not None:d['level']=level
    return d


def stage(root, patients, values, metadata, covariates, relations, *, namespace):
    features=list(values)
    spec=source.stage(root/'inputs',features,{p:[values[f][i] for f in features] for i,p in enumerate(patients)},
                      relations,[],units='counts',quantity='event_count',namespace=namespace,
                      question='How do original count, ordered response and measured dose fields interact with the declared assay relations?')
    columns=list(metadata)
    table=save(root/'inputs/observations.tsv','id\t'+'\t'.join(columns)+'\n'+''.join(
        p+'\t'+'\t'.join(str(metadata[c][i]) for c in columns)+'\n' for i,p in enumerate(patients)))
    assay=spec['assays'][0];assay['observations']['sha256']=digest(table)
    assay['carrier']='observations';assay['partition_size']=len(patients)
    assay['observation_fields']=covariates
    save(root/'inputs/import.json',spec)


def declaration(name, kind, field, **kwargs):
    return dict(id=name,kind=kind,field='all.'+field,meaning='Declared observed '+name+' endpoint',evidence=[],**kwargs)


def query(responses, terms, ablations, workers):
    return dict(operation='observed_response_model',context='all',responses=responses,
                model=dict(intercept=True,terms=terms,ablations=ablations),workers=workers,green_parameter='1',
                interpretation=dict(population='All identified observations on the declared source carrier',
                                    meaning='Compare observed response fields, declared model spans and retained source relations',evidence=[]))


def authored(root,workers=2):
    root.mkdir(parents=True,exist_ok=False)
    ids=[f'specimen-{i+1}' for i in range(8)]
    swap=lambda x:[x[i^1] for i in range(8)]
    values=dict(events=[0,1,3,5,7,9,12,14],successes=[0,1,3,2,4,4,5,6],trials=[4,4,4,4,5,5,6,6])
    values.update(events_ref=swap(values['events']),successes_ref=swap(values['successes']))
    metadata=dict(exposure=[1,2]*4,exposure_minutes=[60,120]*4,
                  dose=[0,0,1,1,3,3,6,6],dose_alias=[0,0,1,1,3,3,6,6],
                  response_area=[100,92,70,80,30,45,15,18],ordinal=[40,40,10,40,10,70,70,70])
    metadata['ordinal_ref']=swap(metadata['ordinal']);metadata['response_area_ref']=swap(metadata['response_area'])
    cov=[covariate('exposure','observation_exposure','days'),covariate('exposure_minutes','observation_exposure','minutes'),
         covariate('dose','measured_dose','micromolar'),covariate('dose_alias','measured_dose','micromolar'),
         covariate('response_area','measured_area','square_millimeters'),covariate('response_area_ref','measured_area','square_millimeters'),
         covariate('ordinal','response_state',kind='categorical_code'),covariate('ordinal_ref','response_state',kind='categorical_code')]
    for g in range(4):
        name='group'+str(g);metadata[name]=[str(i//2) for i in range(8)]
        cov.append(covariate(name,'dose_group',kind='level_indicator',level=str(g)))
    relations=[dict(id='pair-'+str(g),grade=1,meaning='Retained within-dose replicate comparison',
                    participants=[source.participant(ids[2*g],-1),source.participant(ids[2*g+1],1)]) for g in range(4)]
    branch=[source.participant(ids[0],-1)]+[source.participant(i,Q(1,7)) for i in ids[1:]]
    branch += [source.participant(ids[1],1),source.participant(ids[1],-1),source.participant(ids[2],7,False)]
    relations += [dict(id='branch-'+str(g),grade=1,meaning='Distinct retained multiway assay comparison',participants=branch) for g in range(4)]
    relations += [dict(id='filling',grade=2,meaning='Declared zero-boundary comparison between repeated multiway assertions',
                       participants=[source.participant('branch-'+str(g),1 if g<2 else -1) for g in range(4)])]
    stage(root,ids,values,metadata,cov,relations,namespace='993701')
    responses=[declaration('events','count','events',denominator='all.exposure',rate_units='events per day'),
               declaration('successes','binomial_count','successes',denominator='all.trials',rate_units='dimensionless'),
               declaration('state','ordinal','ordinal',levels=[dict(label=l,code=str(c)) for l,c in [('low',40),('intermediate',10),('high',70)]]),
               declaration('area','continuous','response_area')]
    q=query(responses,[dict(id='dose',factors=['all.dose']),dict(id='dose_squared',factors=['all.dose','all.dose']),
                       dict(id='duplicate_dose',factors=['all.dose_alias'])],
            [dict(id='without_curvature',remove=['dose_squared'])],workers)
    q['dose_comparisons']=[dict(id='measured-dose',field='all.dose',groups=[dict(id='dose-'+str(g),indicator='all.group'+str(g)) for g in range(4)],
                               direction='nondecreasing',meaning='Measured replicate means and divided differences on the declared dose levels',evidence=[])]
    q['reference_law']=dict(states=[dict(id='recorded',mass='1/2',response_fields=dict(events='all.events',successes='all.successes',state='all.ordinal',area='all.response_area')),
                                     dict(id='authored-alternative',mass='1/2',response_fields=dict(events='all.events_ref',successes='all.successes_ref',state='all.ordinal_ref',area='all.response_area_ref'))],
                            observed_state='recorded',meaning='Explicit two-state authored reference law with fixed exposures, predictors and source; no biological sampling assertion',evidence=[],alpha='1/2')
    return q


def expression(root,packet,panel,workers=2):
    root.mkdir(parents=True,exist_ok=False);original=root/'source';original.mkdir()
    receipt=json.loads((packet/'source.json').read_text())
    required=list(dict.fromkeys(['source.json','counts.tsv','samples.tsv']+list(receipt['artifacts'])))
    for name in required:
        path=(packet/name).resolve()
        if not path.is_relative_to(packet.resolve()):raise ValueError('Source path escapes its packet')
        if name in receipt['artifacts'] and digest(path).removeprefix('sha256:')!=receipt['artifacts'][name]:raise ValueError('Original source digest changed')
        shutil.copyfile(path,original/name)
    selected=json.loads(panel.read_text());genes=[g['id'] for g in selected];values={}
    with (original/'counts.tsv').open() as f:
        rows=csv.reader(f,delimiter='\t');patients=next(rows)[1:]
        totals=[int(receipt['library_totals'][s]) for s in patients]
        for row in rows:
            if row[0] in genes:values[row[0]]=list(map(int,row[1:]))
    if set(values)!=set(genes):raise ValueError('Missing declared gene')
    samples={r['sample']:r for r in csv.DictReader((original/'samples.tsv').open(),delimiter='\t')}
    donors=list(dict.fromkeys(samples[s]['cell'] for s in patients))
    metadata=dict(retained_total=totals,treated=[1 if samples[s]['dex']=='trt' else 0 for s in patients])
    cov=[covariate('retained_total','retained_read_count','counts'),covariate('treated','recorded_treatment_condition')]
    terms=[];relations=[]
    for d,donor in enumerate(donors):
        name='donor'+str(d);metadata[name]=[samples[s]['cell'] for s in patients]
        cov.append(covariate(name,'recorded_donor',kind='level_indicator',level=donor))
        if d:terms.append(dict(id=name,factors=['all.'+name]))
        members=[s for s in patients if samples[s]['cell']==donor]
        members.sort(key=lambda s:samples[s]['dex']=='trt')
        relations.append(dict(id='pair-'+donor,grade=1,meaning='Recorded within-donor untreated-to-treated comparison',
                              participants=[source.participant(members[0],-1),source.participant(members[1],1)]))
    for condition in ('untrt','trt'):
        members=[s for s in patients if samples[s]['dex']==condition]
        relations.append(dict(id='condition-'+condition,grade=1,meaning='Recorded condition membership with equal declared shares',
                              participants=[source.participant(s,Q(1,len(members))) for s in members]))
    stage(root,patients,values,metadata,cov,relations,namespace='993702')
    spec=json.loads((root/'inputs/import.json').read_text());spec['assays'][0]['quantity']='read_count';save(root/'inputs/import.json',spec)
    terms.append(dict(id='condition',factors=['all.treated']))
    responses=[declaration(g,'count',g,denominator='all.retained_total',rate_units='dimensionless') for g in genes]
    q=query(responses,terms,[dict(id='without_condition',remove=['condition'])],workers)
    q['interpretation']['population']='The six predeclared gene fields on all eight recorded airway specimens'
    q['interpretation']['meaning']='Count-space objectives with complete retained-count exposure and authored donor/condition factors; report full raw and rate residuals'
    save(root/'selected-features.json',selected);save(root/'source-witness.json',dict(files={n:digest(original/n) for n in required},
        exposure='Exact column sum of all 63,677 retained original genes',totals=totals,selection=selected,
        relation_meaning='Recorded pairing and condition membership; no regulatory-network assertion',reference_law='none'))
    return q


def execute(binary,root,q):
    imported,seconds=run(binary,'import-assays','--spec',root/'inputs/import.json','--out',root/'workspace');save(root/'import-result.json',imported)
    source.bind_evidence(q,root/'workspace');e=q['interpretation']['evidence']
    for response in q['responses']:response['evidence']=e
    for dose in q.get('dose_comparisons',[]):dose['evidence']=e
    if 'reference_law' in q:q['reference_law']['evidence']=e
    save(root/'inputs/query.json',q)
    request=dict(schema='gadget.assays.query.v1',operation='native_query',assay='measurements',query=q)
    out,elapsed=run(binary,'assay-query','--workspace',root/'workspace/workspace.json','--request',save(root/'request.json',request))
    save(root/'result.json',out);save(root/'certificate.json',out['result']['answers'][0]['certificate'])
    save(root/'timing.json',dict(import_seconds=seconds,query_seconds=elapsed));return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--panel',type=Path);p.add_argument('--workers',type=int,default=2)
    a=p.parse_args();root=a.out.resolve()
    if a.source:
        if not a.panel:p.error('A declared panel is required')
        q=expression(root,a.source.resolve(),a.panel.resolve(),a.workers)
    else:q=authored(root,a.workers)
    execute(a.binary.resolve(),root,q);print(json.dumps(dict(status='PASS',output=str(root))))
