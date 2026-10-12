"""Complete typed and temporal moments on original and rectangular source carriers."""
from pathlib import Path
import json,copy,subprocess,time,sys,shutil
import sympy as S
binary,consumer,fixture,o=(Path(value).resolve() for value in sys.argv[1:]);shutil.rmtree(o,ignore_errors=True);o.mkdir(parents=True);rows=[]
B=[S.Matrix([[-1,0],[1,0],[0,0]]),S.Matrix([[-1,0,0],[1,-1,0],[0,1,0]]),S.Matrix([[-1,0],[1,0],[0,0]])];C=[S.zeros(2,0),S.Matrix([[0],[0],[1]]),S.zeros(2,0)];coordinates=[[0,2],[0,1,2],[0,2]];names=['baseline','followup','final'];xs=[S.Matrix([1,2]),S.Matrix([3,1,6]),S.Matrix([13,-7])];ys=[S.Matrix([1,2]),S.Matrix([2,1,4]),S.Matrix([3,2])];maps=[S.Matrix([[1,0],[0,1],[1,1]]),S.Matrix([[1,0,1],[0,1,-1]])]
for case_name in ['identity','weighted']:
 case=o/case_name;study=case/'study';subprocess.run([str(fixture),str(study)],check=True,capture_output=True,text=True);metadata=json.loads((study/'study.json').read_text());history=json.loads((study/'history-query.json').read_text())
 metric0=S.diag(1,2,4) if case_name=='weighted' else S.eye(3);all_m1=[S.Rational(2),S.Rational(3),S.Rational(5)] if case_name=='weighted' else [S.Rational(1)]*3;m2=S.Rational(7) if case_name=='weighted' else S.Rational(1)
 if case_name=='weighted':
  for context in metadata['contexts']:
   context['metrics']=[dict(kind='diagonal',values=['1','2','4']),dict(kind='diagonal',values=['2','3','5']),dict(kind='diagonal',values=['7'])]
  (study/'study.json').write_text(json.dumps(metadata))
 M=[S.diag(*[all_m1[i] for i in indices]) for indices in coordinates];L=[metric.inv()*boundary.T*metric0*boundary + higher*(S.eye(higher.cols)/m2)*higher.T*metric for metric,boundary,higher in zip(M,B,C)];queries={};independent=0
 def run(name,q):
  path=case/(name+'.json');path.write_text(json.dumps(q));p=subprocess.run([str(binary),'query','--study',str(study/'study.json'),'--query',str(path)],capture_output=True,text=True);assert p.returncode==0,p.stderr;queries[name]=q;return json.loads(p.stdout)['result']
 for endpoint in range(3):
  for order in [0,1,17]:
   q=dict(operation='typed_moment',context=names[endpoint],field=names[endpoint]+'.x',other_field=names[endpoint]+'.weights',order=order);actual=run('typed-'+names[endpoint]+'-'+str(order),q);expected=(xs[endpoint].T*M[endpoint]*(L[endpoint]**order)*ys[endpoint])[0];assert S.Rational(actual['rational'])==expected;independent+=1
 for old in range(2):
  for order in [0,1,17]:
   q=dict(operation='temporal_moment_change',context=names[old],other_context=names[old+1],field=names[old]+'.x',other_field=names[old]+'.weights',new_field=names[old+1]+'.x',new_other_field=names[old+1]+'.weights',order=order,correspondence=history['correspondences'][old][1]);actual=run('change-'+str(old)+'-'+str(order),q);j=maps[old];f=M[old]*(L[old]**order);g=M[old+1]*(L[old+1]**order);a=j*xs[old];b=j*ys[old];ex=xs[old+1]-a;ey=ys[old+1]-b;expected={'pullback_form':(xs[old].T*(j.T*g*j-f)*ys[old])[0],'left_innovation':(ex.T*g*b)[0],'right_innovation':(a.T*g*ey)[0],'innovation_interaction':(ex.T*g*ey)[0]};expected['total']=sum(expected.values());interval=S.Rational(2 if old==0 else 3);expected.update(interval=interval,finite_rate=expected['total']/interval)
   for key,value in expected.items():assert S.Rational(actual[key]['rational'])==value,(case_name,key,value,actual[key]);independent+=1
  zero=copy.deepcopy(q);zero['order']=1;zero['correspondence']['entries']=[];actual=run('zero-comparison-'+str(old),zero);expected=(xs[old+1].T*M[old+1]*L[old+1]*ys[old+1])[0]-(xs[old].T*M[old]*L[old]*ys[old])[0];assert S.Rational(actual['total']['rational'])==expected;independent+=1
 path=case/'queries.json';path.write_text(json.dumps(queries));start=time.monotonic();p=subprocess.run([str(consumer),str(study/'study.json'),str(path),str(case/'native'),'moment'],capture_output=True,text=True,timeout=900);row=dict(case=case_name,query_count=len(queries),independent_checks=independent,native_exit_code=p.returncode,native_stdout=p.stdout,native_stderr=p.stderr,seconds=time.monotonic()-start);rows.append(row);(o/'contracts.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True);assert p.returncode==0,p.stderr

print(json.dumps(dict(status='PASS',cases=rows)))
