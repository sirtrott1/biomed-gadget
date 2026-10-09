#!/usr/bin/env python3
"""Run native temporal histories and interval-observed clinical endpoints."""
import argparse,csv,hashlib,json,shutil,subprocess,time
from pathlib import Path
from fractions import Fraction as Q

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(value if isinstance(value,str) else json.dumps(value,indent=2)+'\n');return path

def digest(path):return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()
def run(binary,*args):
    start=time.perf_counter();p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
    if p.returncode:raise RuntimeError(p.stderr)
    return json.loads(p.stdout),time.perf_counter()-start

def history(binary,fixture,root):
    subprocess.run([str(fixture),str(root/'study')],check=True,capture_output=True,text=True)
    query=json.loads((root/'study/history-query.json').read_text())
    out,elapsed=run(binary,'query','--study',root/'study/study.json','--query',root/'study/history-query.json')
    save(root/'certificate.json',out);save(root/'timing.json',dict(query_seconds=elapsed))
    return query,out

def clinical_rows(source):
    rows=list(csv.DictReader((source/'lung.csv').open()));out=[]
    for row in rows:
        status=int(Q(row['status']));assert status in (1,2)
        t=str(Q(row['time']));assert Q(t)>=0
        # Record original status coding, clinical variables and all original columns.
        out.append(dict(source_row=row['source_row'],lower=t,upper=t if status==2 else 'NA',
                        lower_closed='1' if status==2 else '0',upper_closed='1' if status==2 else '0',upper_finite='1' if status==2 else '0',
                        male=str(int(Q(row['sex'])==1)),female=str(int(Q(row['sex'])==2)),
                        **{k:row[k] for k in row if k!='source_row'}))
    return out

def authored_rows():
    intervals=[(2,2,1,1),(4,None,0,0),(0,5,1,1),(3,8,0,1),(7,7,1,1),(9,None,0,0),(10,11,0,0),(12,12,1,1)]
    return [dict(source_row=str(i),lower=str(l),upper=str(u) if u is not None else 'NA',lower_closed=str(lc),upper_closed=str(uc),upper_finite=str(int(u is not None)),groupA=str(int(i%2==0)),groupB=str(int(i%2!=0)),group='A' if i%2==0 else 'B') for i,(l,u,lc,uc) in enumerate(intervals)]

def endpoints(binary,root,source=None):
    root.mkdir(parents=True,exist_ok=False);inputs=root/'inputs';inputs.mkdir()
    if source:
        shutil.copytree(source,root/'source');receipt=json.loads((root/'source/provenance.json').read_text())
        if digest(root/'source/lung.csv').removeprefix('sha256:')!=receipt['decoder']['sha256']:raise ValueError('Clinical source digest changed')
        rows=clinical_rows(root/'source');strata=['male','female'];times=['90','180','365','730'];population='228 original NCCTG lung cancer records'
        dictionary=(root/'source/lung.Rd').read_text()+'\nEach original record remains in source/lung.csv. Derived fields: death status 2 gives [time,time]; censoring status 1 gives (time,infinity). Sex indicators preserve original sex codes 1 and 2. All original variables and missing values are retained. Row numbers are source identities, not supplied patient identifiers.\n'
        source_url='https://github.com/therneau/survival';attribution='Terry Therneau; NCCTG lung cancer data distributed with survival';license_name='GPL-2, survival package distribution'
    else:
        rows=authored_rows();strata=['groupA','groupB'];times=['0','2','4','7','10','14'];population='Eight authored interval observations'
        dictionary='Authored independent event-time intervals. lower/upper in days; upper NA and upper_finite=0 means unbounded above. The closed flags state endpoint inclusion. groupA/groupB are explicit alternating authored strata.\n';source_url='authored://event-intervals';attribution='GADGET method fixture';license_name='Authored nonclinical fixture'
    columns=list(rows[0]);raw=','.join(columns)+'\n'+''.join(','.join(row[c] for c in columns)+'\n' for row in rows)
    save(inputs/'observations.csv',raw);save(inputs/'dictionary.txt',dictionary)
    descriptions=[]
    for c in columns:
        is_flag=c in ('lower_closed','upper_closed','upper_finite',*strata)
        categorical=c in ('status','sex','inst','ph.ecog','source_row','group')
        descriptions.append(dict(id=c,channel='observation',quantity='event_time' if c in ('lower','upper') else c,units='day' if c in ('lower','upper','time') else 'dimensionless',value_type='indicator' if is_flag else ('categorical' if categorical else 'quantitative'),dictionary_start=0,dictionary_end=len(dictionary.encode())))
    # A category membership relation is an explicit comparison population, not
    # an assertion of a biological mechanism. No event time is supplied for NA.
    spec=dict(schema='gadget.table.import.v1',id='interval-observed-endpoints',source_namespace='995801',source=dict(path='observations.csv',sha256=digest(inputs/'observations.csv'),delimiter=',',header=True),dictionary=dict(path='dictionary.txt',sha256=digest(inputs/'dictionary.txt')),columns=descriptions,missing_tokens=['NA'],relation_columns=['sex' if source else 'group'],identity_column='source_row',question='What event and restricted-time bounds are determined by the reported observations?',population=population,selection='Every original row',limitations=['Bounds condition on recorded event intervals and independent admissible completions'],source_url=source_url,attribution=attribution,license=license_name)
    run(binary,'import-table','--spec',save(inputs/'import.json',spec),'--out',root/'study')
    study=json.loads((root/'study/study.json').read_text());evidence=next(f['evidence'] for f in study['fields'] if f['id']=='all.lower')
    query=dict(operation='censored_endpoint_bounds',context='all',lower='all.lower',upper='all.upper',lower_closed='all.lower_closed',upper_closed='all.upper_closed',upper_finite='all.upper_finite',selection=[],strata=['all.'+s for s in strata],times=times,origin='0',time_unit='day',population=population,event='Death' if source else 'Authored event',evidence=evidence,contrasts=[dict(id='group-difference',left='all.'+strata[0],right='all.'+strata[1],meaning='Difference between the two declared finite-population event fractions and restricted durations')])
    out,elapsed=run(binary,'query','--study',root/'study/study.json','--query',save(inputs/'query.json',query));save(root/'certificate.json',out);save(root/'timing.json',dict(query_seconds=elapsed));return query,out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixture',type=Path);p.add_argument('--source',type=Path);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    if args.fixture:history(args.binary.resolve(),args.fixture.resolve(),args.out.resolve())
    else:endpoints(args.binary.resolve(),args.out.resolve(),args.source.resolve() if args.source else None)
