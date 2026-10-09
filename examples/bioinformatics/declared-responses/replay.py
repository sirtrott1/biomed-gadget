#!/usr/bin/env python3
"""Produce source-bound response transitions and independently check every term."""
import argparse, copy, hashlib, json, subprocess
from pathlib import Path

def save(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2)+'\n'); return p
def digest(p): return 'sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()
def run(binary, *args):
    p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
    if p.returncode: raise RuntimeError(p.stderr)
    return json.loads(p.stdout)
def matrix(rows):
    return dict(rows=len(rows),columns=len(rows[0]) if rows else 0,
        entries=[dict(row=i,column=j,value=str(v)) for i,row in enumerate(rows) for j,v in enumerate(row) if v!=0])
def identity(n): return dict(rows=n,columns=n,entries=[dict(row=i,column=i,value='1') for i in range(n)])
def op(name,key,role,source,target,m=None,**extra):
    d=dict(id=name,key=str(key),role=role,input=str(source),output=str(target),quantity='declared_bilinear_form' if role=='metric' else 'declared_linear_action',units='dimensionless',method='Explicit authored positive bilinear form' if role=='metric' else 'Explicit authored linear map')
    d.update(extra)
    if m is not None:d['map']=m
    return d
def field(name,key,carrier,values):return dict(id=name,key=str(key),carrier=str(carrier),values=list(map(lambda x:None if x is None else str(x),values)),quantity='authored_input',units='dimensionless',value_type='quantitative')
def source(n):return dict(vertices=[str(i+1) for i in range(n)],relations=[])
def declaration(name,seq,L,I,A,u,src=None,grade=0,metric=None):
    carriers=[dict(key='10',grade=grade,size=L['rows']),dict(key='20',coordinate_system='200',size=I['columns']),dict(key='30',coordinate_system='300',size=A['rows'])]
    operators=[op('L',1,'resolving',10,10,L),op('I',2,'realization',20,10,I),op('A',3,'accession',10,30,A),op('M',4,'metric',30,30,metric or identity(A['rows']))]
    return dict(id=name,namespace='9001',sequence=str(seq),source=src or source(L['rows']),carriers=carriers,operators=operators,fields=[field('u',100,20,u),field('u2',101,20,[str(i+2) for i in range(I['columns'])]),field('zero',102,20,['0']*I['columns']),field('partial',103,20,[None]*I['columns'])])
def square():
    relations=[]
    for i in range(4):
        occ=[]
        for j,(cell,exists,negative,mag) in enumerate([(i+1,True,True,1),((i+1)%4+1,True,False,1),(i+1,True,True,2),(i+1,True,False,2),((i+2)%4+1,False,False,7)]):
            occ.append(dict(key=str(1001+5*i+j),cell=str(cell),exists=exists,negative=negative,magnitude=str(mag)))
        relations.append(dict(key=str(101+i),occurrences=occ))
    s=dict(vertices=['1','2','3','4'],relations=relations)
    ns=copy.deepcopy(s);ns['higher']=[dict(grade=2,key='201',occurrences=[dict(key=str(2001+i),cell=str(101+i),exists=True,negative=False,magnitude='1') for i in range(4)])]
    L=matrix([[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]])
    N=matrix([[3,0,1,0],[0,3,0,1],[1,0,3,0],[0,1,0,3]])
    a=declaration('unfilled-square',7,L,matrix([[1],[1],[0],[0]]),matrix([[1,0,0,0],[0,1,0,0]]),[2],s,1)
    b=declaration('filled-square',8,N,matrix([[0],[0],[-1],[-1]]),matrix([[1,0,1,0],[0,1,0,1]]),[3],ns,1, matrix([[2,'1/3'],['1/3',3]]))
    a['operators'][0].pop('map');a['operators'][0].update(hodge_grade=1,method='Original boundary tower with declared identity metrics')
    b['operators'][0].pop('map');b['operators'][0].update(hodge_grade=1,method='Original boundary tower including the explicitly attached face with identity metrics')
    return a,b,identity(4),identity(1),identity(2)
def connect(binary,root,old_id,new_id,old_field,new_field,J,K,KA,old_carrier,new_carrier,old_source,new_source,old_output,new_output):
    d=json.loads((root/'p1/study.json').read_text());d['id']='response-coordinator';d['native']['path']='p1/'+d['native']['path']
    for document in d['documents']:document['path']='p1/'+document['path']
    evidence=d['operators'][-1]['evidence']
    d['publications']=[]
    for p in ['p0','p1']:
        path=root/p/'study.json';cat=run(binary,'query','--study',path,'--query',save(root/'inspect.json',dict(operation='inspect')))
        d['publications'].append(dict(id=p,study=dict(path=p+'/study.json',sha256=digest(path)),source_digest=cat['source_digest'],meaning='Actual independently admitted publication '+p,evidence=evidence))
    save(root/'study.json',d)
    def comparison(m,a,b,meaning):return dict(old_carrier=str(a),new_carrier=str(b),map=m,meaning=meaning,evidence=evidence)
    a=json.loads((root/'p0/study.json').read_text());b=json.loads((root/'p1/study.json').read_text())
    endpoint=lambda p,x,f:dict(publication=p,resolving='L',realization='I',field=f,source_sequence=x['native']['cutoff'])
    q=dict(operation='declared_response_transition',old=endpoint('p0',a,old_field),new=endpoint('p1',b,new_field),**{'lambda':'1'},comparisons=dict(ambient=comparison(J,old_carrier,new_carrier,'Declared complete map of original ambient coordinates'),source=comparison(K,old_source,new_source,'Declared complete map of original source coordinates')))
    save(root/'response-query.json',q);response=run(binary,'query','--study',root/'study.json','--query',root/'response-query.json');save(root/'response.json',response)
    observed=copy.deepcopy(q);observed['observation']=dict(old=old_id,new=new_id,comparison=comparison(KA,old_output,new_output,'Declared complete map of output coordinates'));observed['pairing']=dict(metric='M')
    save(root/'observed-query.json',observed);save(root/'observed.json',run(binary,'query','--study',root/'study.json','--query',root/'observed-query.json'))
    workflow=dict(schema='gadget.scientific.workflow.v1',id='declared-native-response',source_digest=response['source_digest'],limits=dict(depth=2,nodes=3,model_calls=0,parallel_workers=2),definitions=[],entry=dict(id='joined',type='parallel',workers=2,steps=[dict(id='response',type='query',query=q),dict(id='observed',type='query',query=observed)]))
    save(root/'workflow-request.json',workflow);receipt=run(binary,'workflow','--study',root/'study.json','--request',root/'workflow-request.json');assert receipt['status']=='complete',receipt;save(root/'workflow.json',receipt)
    return q,observed
def authored(binary,producer,root,case):
    root.mkdir(parents=True,exist_ok=False)
    a,b,J,K,KA=case
    for p,spec in [('p0',a),('p1',b)]:
        path=save(root/(p+'-declarations.json'),spec)
        subprocess.run([str(producer),str(path),str(root/p)],check=True,capture_output=True,text=True)
    connect(binary,root,'A','A','u','u',J,K,KA,10,10,20,20,30,30)
def clinical(binary,producer,original,root,offset=0):
    root.mkdir(parents=True,exist_ok=False)
    publications=[]
    for p,t in [('p0',offset),('p1',offset+1)]:
        source_study=original/('p'+str(t))/'study.json';d=json.loads(source_study.read_text())
        catalog=run(binary,'query','--study',source_study,'--query',save(root/'inspect.json',dict(operation='inspect')))['result']
        f=next(x for x in catalog['fields'] if x['id']=='all.time');ambient=f['carrier'];n=f['coordinates'];assert n in (80,160,228)
        spec=dict(id='declared-clinical-'+str(t),namespace=d['native']['source_namespace'],operators=[op('L',910001,'resolving',ambient,ambient,hodge_grade=0,method='Original sex-membership boundary tower with original declared identity metrics'),op('I',910002,'realization',ambient,ambient,identity(n)),op('A',910003,'accession',ambient,990003,boundary_transpose_grade=1,method='Algebraic boundary transpose, original sex-membership relations'),op('M',910004,'metric',990003,990003,identity(2))],carriers=[dict(key='990003',grade=1,size=2)])
        path=save(root/(p+'-declarations.json'),spec);subprocess.run([str(producer),str(path),str(root/p),str(source_study)],check=True,capture_output=True,text=True);publications.append((ambient,n))
    a,n=publications[0];b,m=publications[1];J=dict(rows=m,columns=n,entries=[dict(row=i,column=i,value='1') for i in range(n)])
    connect(binary,root,'A','A','all.time','all.time',J,J,identity(2),a,b,a,b,990003,990003)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--producer',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--clinical',type=Path);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    authored(a.binary.resolve(),a.producer.resolve(),a.out/'square',square())
    rectangular=(declaration('two-records',17,matrix([[2,1],[0,3]]),matrix([[1,0,2],[0,1,-1]]),matrix([[1,2]]),[3,5,7]),declaration('three-records',21,matrix([[1,0,2],[0,2,1],[0,0,3]]),matrix([[1,0],[2,1],[0,1]]),matrix([[1,0,1],[0,1,2]]),[11,13]),matrix([[1,0],[0,1],[1,-1]]),matrix([[1,0,0],[0,1,0]]),matrix([[1],[2]]))
    authored(a.binary.resolve(),a.producer.resolve(),a.out/'rectangular',rectangular)
    if a.clinical:
        clinical(a.binary.resolve(),a.producer.resolve(),a.clinical.resolve(),a.out/'clinical')
        clinical(a.binary.resolve(),a.producer.resolve(),a.clinical.resolve(),a.out/'clinical-late',1)
    from verify import verify
    save(a.out/'verification.json',dict(cases=[verify(x) for x in sorted(a.out.iterdir()) if x.is_dir()]))
