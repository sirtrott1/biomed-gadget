"""Independent interval, controlled projection, neighborhood and lineage oracles."""
from fractions import Fraction as Q
import copy
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)
checks = 0
requests = []

def check(value, message):
    global checks
    checks += 1
    assert value, message

def run(*args, fail=False):
    p = subprocess.run([str(binary), *map(str,args)], capture_output=True, text=True)
    if fail:
        check(p.returncode != 0, 'Expected rejection: '+str(args))
        return p.stderr
    if p.returncode: raise RuntimeError(str(args)+': '+p.stderr)
    return json.loads(p.stdout)

def save(name, value):
    p=root/name; p.write_text(json.dumps(value,indent=2)+'\n'); return p

def file(name, content):
    p=root/name; p.write_text(content)
    return dict(path=name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())

frame=dict(assembly='authored',sequence_sha256='sha256:'+'0'*64,coordinate_convention='zero_based_half_open',contigs=[dict(id='chr',length=1000000),dict(id='chrB',length=1000000)])
design=dict(question='Separate declared contrast fields from control structure',population='Authored finite observations',selection='All original records',limitations=['Finite algebraic example'])

def admit(name, resources=(), assays=()):
    spec=dict(schema='gadget.assays.import.v1',id=name,design=design,resources=list(resources),assays=list(assays),links=[])
    run('import-assays','--spec',save(name+'-import.json',spec),'--out',root/name)
    return root/name/'workspace.json'

def bed(name, records):
    body=''.join(f'{c}\t{a}\t{b}\t{id}\t{w}\t{s}\n' for c,a,b,id,w,s in records)
    return dict(id=name,format='bed',source=file(name,body),reference=frame,weight_column=4,method='Authored interval weights')

def query(name, workspace, request, fail=False):
    result=run('assay-query','--workspace',workspace,'--request',save(name+'-request.json',request),fail=fail)
    if not fail:
        requests.append(copy.deepcopy(request));save(name+'-result.json',result);return result['result']
    return result

def sparse(m):
    out=[[Q(0) for _ in range(m['columns'])] for _ in range(m['rows'])]
    for e in m['entries']:out[e['row']][e['column']]=Q(e['value'])
    return out

def rank(rows):
    a=copy.deepcopy(rows)
    if not a:return 0
    n=len(a[0]);r=0
    for c in range(n):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];p=a[r][c];a[r]=[x/p for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                q=a[i][c];a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def solve(a,b):
    a=[list(row)+[value] for row,value in zip(a,b)]
    for c in range(len(a)):
        pivot=next(i for i in range(c,len(a)) if a[i][c])
        a[c],a[pivot]=a[pivot],a[c];p=a[c][c];a[c]=[x/p for x in a[c]]
        for i in range(len(a)):
            if i!=c:
                q=a[i][c];a[i]=[x-q*y for x,y in zip(a[i],a[c])]
    return [row[-1] for row in a]

def weighted(a,b,lengths):return sum((x*y*l for x,y,l in zip(a,b,lengths)),Q(0))

def oracle(result, request, raw):
    records=[(s,*row) for s,rows in enumerate(raw) for row in rows]
    frames=sorted({(r[1],r[-1] if request['strand_policy']=='separate' else 'ignored') for r in records})
    atoms=[]
    for c,s in frames:
        points=sorted({v for r in records if r[1]==c and (request['strand_policy']=='ignore' or r[-1]==s) for v in r[2:4]})
        for a,b in zip(points,points[1:]):
            members=[j for j,r in enumerate(records) if r[1]==c and (request['strand_policy']=='ignore' or r[-1]==s) and r[2]<=a and b<=r[3]]
            if members:atoms.append((c,s,a,b,members))
    check(len(result['segments'])==len(atoms),'Complete atom carrier')
    lengths=[Q(b-a) for c,s,a,b,j in atoms]
    check(list(map(Q,result['metric_diagonal']))==lengths,'Physical length metric')
    bmatrix=[[Q(j in members) for j in range(len(records))] for c,s,a,b,members in atoms]
    fields=[]
    for s,sample in enumerate(request['samples']):
        field=[sum((Q(records[j][5]) for j in members if records[j][0]==s),Q(0))/Q(sample['scale']) for c,st,a,b,members in atoms]
        fields.append(field)
        check(list(map(Q,result['sample_fields'][s]['values']))==field,'Independent coverage action')
        # Independent per-base count integral, including overlapping repeats and zeros.
        integral=Q(0)
        for c,st in frames:
            ends=[r[3] for r in records if r[1]==c and (request['strand_policy']=='ignore' or r[-1]==st)]
            for x in range(max(ends,default=0)):
                integral+=sum((Q(r[5]) for r in records if r[0]==s and r[1]==c and (request['strand_policy']=='ignore' or r[-1]==st) and r[2]<=x<r[3]),Q(0))/Q(sample['scale'])
        check(Q(result['sample_fields'][s]['integral'])==integral,'Independent per-base coverage integral')
    for i,(c,st,a,b,members) in enumerate(atoms):
        out=result['segments'][i]
        check((out['contig'],out['strand'],out['start'],out['end'],out['source_record_indices'])==(c,st,a,b,members),'All physical memberships')
    check(result['native_occurrence_count']==sum(len(a[-1]) for a in atoms),'Every native occurrence')
    check(len(result['source_records'])==len(records),'Original intervals remain distinct')
    check(len(result['relation_declaration']['relations'])==len(records),'Interval-primary source')
    for j,r in enumerate(records):
        relation=result['relation_declaration']['relations'][j]
        check([p['cell'] for p in relation['participants']]==['segment:'+str(i) for i,a in enumerate(atoms) if j in a[-1]],'Complete arity and source record map')
        check(all(p['exists'] and not p['negative'] and p['magnitude']=='1' for p in relation['participants']),'Declared unsigned physical membership')
        source=result['source_records'][j];f=root/source['resource'];span=source['source'][0]
        check(f.read_bytes()[span['start']:span['end']].decode()==span['text'],'Original interval byte span')
        check(source['weight']==str(Q(r[5])),'Original zero/repeated weight')
    ids=[s['id'] for s in request['samples']]
    columns=[fields[ids.index(s)] for s in request['controls']]
    if request['intercept']:columns.append([Q(1)]*len(atoms))
    independent=[]
    for col in columns:
        rows=list(map(list,zip(*(independent+[col]))))
        if rank(rows)>len(independent):independent.append(col)
    gram=[[weighted(a,b,lengths) for b in independent] for a in independent]
    check(result['controls']['rank']==len(independent),'Control span rank')
    check(result['controls']['factorizations']==1,'Shared canonical prepared solve')
    residuals=[]
    for c,contrast in enumerate(request['contrasts']):
        inp=[sum((Q(t['coefficient'])*fields[ids.index(t['sample'])][i] for t in contrast['terms']),Q(0)) for i in range(len(atoms))]
        co=solve(gram,[weighted(col,inp,lengths) for col in independent]) if independent else []
        represented=[sum((a*col[i] for a,col in zip(co,independent)),Q(0)) for i in range(len(atoms))]
        residual=[a-b for a,b in zip(inp,represented)];residuals.append(residual)
        out=result['contrasts'][c]
        for key,expected in [('input',inp),('represented',represented),('residual',residual)]:
            check(list(map(Q,out[key]))==expected,'Independent '+key+' field')
            check(Q(out['quadrances'][key])==weighted(expected,expected,lengths),'Exact '+key+' quadrance')
        co_source=[sum((Q(t['coefficient'])*Q(r[5])/Q(request['samples'][r[0]]['scale']) for t in contrast['terms'] if t['sample']==ids[r[0]]),Q(0)) for r in records]
        check(list(map(Q,out['source_coefficients']))==co_source,'Original interval contrast cochain')
        check(out['source_action_exact'],'Commuting native source action')
        check(weighted(represented,residual,lengths)==0,'Orthogonal decomposition')
        check(all(weighted(col,residual,lengths)==0 for col in columns),'All redundant controls annihilated')
        for cut,component in zip(map(Q,request['thresholds']),out['residual_threshold_components']):
            expected=[]
            for i,(contig,st,a,b,members) in enumerate(atoms):
                if residual[i]<cut:continue
                if expected and expected[-1]['contig']==contig and expected[-1]['strand']==st and expected[-1]['end']==a:expected[-1]['end']=b
                else:expected.append(dict(contig=contig,strand=st,start=a,end=b))
            check(component['regions']==expected,'Exact physical threshold components')
            check(Q(component['retained_length'])==sum(l for l,x in zip(lengths,residual) if x>=cut),'Complete retained length')
    expected=[[weighted(a,b,lengths) for b in residuals] for a in residuals]
    check(sparse(result['residual_gram'])==expected,'Signed residual pairings')
    if request.get('source_kernel'):
        k=sparse(result['source_kernel']['frame'])
        check(len(records)-result['source_kernel']['nullity']==rank(bmatrix),'Independent complete boundary rank')
        check(rank(k)==result['source_kernel']['nullity'],'Independent source kernel frame')
        check(all(sum((a*k[j][c] for j,a in enumerate(row)),Q(0))==0 for row in bmatrix for c in range(result['source_kernel']['nullity'])),'Original-variable boundary kernel')

raw=[ [('chr',0,4,'repeat','2','+'),('chr',0,4,'repeat','1','+'),('chr',2,8,'zero','0','+'),('chrB',1,4,'b','1','-')],
      [('chr',0,8,'c','1','+'),('chrB',0,4,'b','2','-')],
      [('chr',1,6,'t','3/2','+'),('chr',4,8,'minus','1','-')] ]
resources=[bed('s'+str(i)+'.bed',r) for i,r in enumerate(raw)]
workspace=admit('controlled',resources)
samples=[dict(id='s'+str(i),resource=r['id'],scale=str(i+1),scale_provenance='Declared divisor') for i,r in enumerate(resources)]
request=dict(schema='gadget.assays.query.v1',operation='controlled_coverage',samples=samples,contrasts=[dict(id='mix',terms=[dict(sample='s0',coefficient='-1/2'),dict(sample='s2',coefficient='2')],evidence='Declared contrast'),dict(id='control-change',terms=[dict(sample='s1',coefficient='1'),dict(sample='s0',coefficient='-1')],evidence='Declared comparison')],controls=['s1'],intercept=True,metric='reference_length',strand_policy='separate',thresholds=['-1','0','1/3'],absence_policy='measured_zero_within_union',units='declared weight per scale',evidence='Original intervals in one declared reference frame',source_kernel=True)
for name,changes in [('declared',{}),('no-control',dict(controls=[],intercept=False)),('intercept-only',dict(controls=[],intercept=True)),('dependent',dict(controls=['s0','s1','s2'],intercept=True)),('ignore-strand',dict(strand_policy='ignore'))]:
    q={**request,**changes};result=query(name,workspace,q);oracle(result,q,raw)
for key,value in [('controls',['missing']),('controls',['s1','s1']),('metric','identity'),('absence_policy','missing'),('source_kernel','yes')]:
    q={**request,key:value};query('reject-'+str(checks),workspace,q,fail=True)
q=copy.deepcopy(request);q['contrasts'][0]['terms'][0]['sample']='missing';query('unknown-term',workspace,q,fail=True)
q=copy.deepcopy(request);q['samples'][0]['scale']='0';query('zero-scale',workspace,q,fail=True)
different=copy.deepcopy(resources[1]);different['reference']['assembly']='other-frame'
w=admit('different-reference',[resources[0],different,resources[2]])
query('different-reference',w,request,fail=True)
q=copy.deepcopy(request);q['contrasts'][0]['terms'].append(copy.deepcopy(q['contrasts'][0]['terms'][0]));query('duplicate-term',workspace,q,fail=True)

# Random geometries, rational weights, both strand semantics and dependent controls.
rng=random.Random(771302)
for trial in range(8):
    raw2=[]
    for s in range(3):
        rows=[]
        for j in range(3):
            a=rng.randrange(0,8);b=rng.randrange(a+1,12)
            rows.append(('chr',a,b,'r'+str(j),str(Q(rng.randrange(0,6),rng.randrange(1,4))),rng.choice(['+','-'])))
        raw2.append(rows)
    res=[bed(f'random-{trial}-{s}.bed',rows) for s,rows in enumerate(raw2)]
    w=admit('random-'+str(trial),res);q=copy.deepcopy(request)
    for s in range(3):q['samples'][s]['resource']=res[s]['id']
    q['strand_policy']='ignore' if trial%2 else 'separate';q['controls']=['s0','s1'];q['intercept']=bool(trial%3)
    oracle(query('random-'+str(trial),w,q),q,raw2)
# Former sample, contrast and threshold ceilings, with complete outputs.
many=[dict(id='s'+str(i),resource=resources[0]['id'],scale='1',scale_provenance='Identity divisor') for i in range(129)]
legacy=dict(schema='gadget.assays.query.v1',operation='chromatin_contrast',samples=many,contrasts=[dict(id='c'+str(i),treated='s0',control='s'+str(i)) for i in range(129)],strand_policy='separate',thresholds=[str(i) for i in range(129)],absence_policy='measured_zero_within_union')
r=query('many-declarations',workspace,legacy)
check(len(r['sample_fields'])==129 and len(r['contrasts'])==129,'All declared samples and contrasts')
check(all(len(c['threshold_components'])==129 and all(Q(x)==0 for x in c['values']) for c in r['contrasts']),'Complete threshold and zero fields')
# Complete source record admission beyond the former combined interval ceiling.
res=bed('many-records.bed',[('chr',0,1,'repeat','1','+')]*100001)
w=admit('many-records',[res]);q={**legacy,'samples':[dict(id='s',resource=res['id'],scale='1',scale_provenance='Identity divisor')],'contrasts':[dict(id='same',treated='s',control='s')],'thresholds':['0']}
r=query('many-records',w,q)
check(len(r['segments'][0]['memberships'])==100001 and r['sample_fields'][0]['values']==['100001'],'Complete original record carrier')
# Actual legacy interval serialization beyond the former atom ceiling.
res=bed('many-atoms.bed',[('chr',i,i+1,str(i),'1','+') for i in range(20001)])
w=admit('many-atoms',[res]);q=dict(schema='gadget.assays.query.v1',operation='interval_coverage',resource=res['id'],strand_policy='separate',thresholds=[str(i) for i in range(129)])
r=query('many-atoms',w,q)
check(len(r['segments'])==20001 and len(r['threshold_components'])==129 and r['coverage_integral']=='20001','Complete atoms and threshold declarations')

# Spatial independent pairwise oracle, exact axis exclusion and original order.
def spatial(name, positions, diagonal, radius, complete_oracle=True):
    dims=len(diagonal);columns=['x'+str(i) for i in range(dims)]
    axis=file(name+'-axis.tsv','id\t'+'\t'.join(columns)+'\tcoordinate_frame\tcoordinate_units\n'+''.join(str(i)+'\t'+'\t'.join(map(str,p))+'\tsection\tum\n' for i,p in enumerate(positions)))
    features=file(name+'-features.tsv','id\nf\n');matrix=file(name+'-matrix.tsv','id\t'+'\t'.join(map(str,range(len(positions))))+'\nf\t'+'\t'.join(['1']*len(positions))+'\n')
    assay=dict(id='points',kind='spatial_expression',format='dense_tsv',source=matrix,features=features,observations=axis,source_namespace='772001',carrier='observations',quantity='read_count',units='counts',channel='E',method='Declared coordinates',arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA'])
    w=admit(name,assays=[assay]);q=dict(schema='gadget.assays.query.v1',operation='spatial_relations',assay='points',coordinates=columns,frame='section',units='um',metric_diagonal=list(map(str,diagonal)),radius_squared=str(radius),relation_type='physical.neighborhood',evidence='Declared physical coordinate frame',execution_certificate=True)
    r=query(name,w,q);pairs=[]
    if complete_oracle:
        for i,p in enumerate(positions):
            for j in range(i+1,len(positions)):
                d=sum((Q(a)*(Q(x)-Q(y))**2 for a,x,y in zip(diagonal,p,positions[j])),Q(0))
                if d<=Q(radius):pairs.append((i,j,d))
        actual=[(int(rel['participants'][0]['cell']),int(rel['participants'][1]['cell']),Q(dist['squared_distance'])) for rel,dist in zip(r['relations']['relations'],r['distances'])]
        check(actual==pairs,'Complete spatial neighborhoods and native exact distances')
    check(r['execution']['complete'],'Exact index completion')
    return w,q,r
points=[[Q(rng.randrange(-5,6),2) for _ in range(9)] for _ in range(27)]
points[6]=points[3][:]
for radius in ['0','3/2','50']:spatial('spatial-'+radius.replace('/','-'),points,[Q(i+1,3) for i in range(9)],radius)
w,q,r=spatial('large-spatial',[[3*i]+[0]*8 for i in range(4101)],[1]*9,'1',False)
check(not r['distances'] and r['execution']['native_metric_comparisons']==0,'4101 observations and 9 dimensions with exact exclusion')
w,q,r=spatial('dense-spatial',[[0,0]]*145,[1,1],'0')
check(len(r['distances'])==10440,'Complete dense relation output beyond former ceiling')
q['frame']='wrong';query('bad-frame',w,q,fail=True)
q['frame']='section';q['metric_diagonal'][0]='0';query('nonpositive-metric',w,q,fail=True)

# Generic assay metrics are retained on each real storage partition.
axis=file('metric-features.tsv','id\na\nb\nc\nd\n');obs=file('metric-observations.tsv','id\nx\n');data=file('metric-matrix.tsv','id\tx\na\t1\nb\t2\nc\t3\nd\t4\n')
metric_assay=dict(id='metric',kind='declared_expression',format='dense_tsv',source=data,features=axis,observations=obs,source_namespace='772011',carrier='features',quantity='declared_field',units='weights',channel='E',method='Declared metric fixture',arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA'],partition_size=2,metrics=[dict(kind='diagonal',values=['2','3','5','7'])])
w=admit('metric-partitions',assays=[metric_assay])
# Native query execution validates each partition's actual metric dimensions.
query('metric-partitions',w,dict(schema='gadget.assays.query.v1',operation='native_query',assay='metric',query=dict(operation='inspect'),partitions=[0,1]))
for label,values in [('short',['2','3']),('nonpositive',['2','0','5','7'])]:
    bad=copy.deepcopy(metric_assay);bad['metrics'][0]['values']=values
    spec=dict(schema='gadget.assays.import.v1',id='invalid-metric',design=design,resources=[],assays=[bad],links=[])
    run('import-assays','--spec',save('invalid-'+label+'.json',spec),'--out',root/('invalid-'+label),fail=True)

# Materialized contrast fields and interval-primary native carrier are portable.
original_result=query('materialize-oracle',workspace,request)
material=dict(schema='gadget.assays.materialize.v1',id='fields',source_namespace='771992',design=design,recipe=request,channel='E')
run('derive-assays','--workspace',workspace,'--request',save('materialize.json',material),'--out',root/'derived')
w=root/'derived/workspace.rfcxs'
for method in ['inspect','hodge','green_action']:
    inner=dict(operation=method)
    if method!='inspect':inner.update(context='all',field='all.interval-contrast:0' if method=='hodge' else 'all.contrast:0:residual')
    if method=='green_action':inner['lambda']='1'
    native=query('native-'+method,w,dict(schema='gadget.assays.query.v1',operation='native_query',assay='fields',query=inner))
    if method != 'inspect':
        certificate=native['answers'][0]['certificate']
        actual=[Q(v['value']['rational']) for v in certificate['result']]
        seg=original_result['segments'];rels=original_result['relation_declaration']['relations'];index={a['id']:i for i,a in enumerate(seg)}
        boundary=[[Q(0) for _ in rels] for _ in seg]
        for j,rel in enumerate(rels):
            for p in rel['participants']:boundary[index[p['cell']]][j]+=Q(p['magnitude'])
        lengths=list(map(Q,original_result['metric_diagonal']))
        if method=='hodge':
            x=list(map(Q,original_result['contrasts'][0]['source_coefficients']))
            physical=[sum((q*v for q,v in zip(row,x)),Q(0)) for row in boundary]
            expected=[sum((boundary[i][j]*lengths[i]*physical[i] for i in range(len(seg))),Q(0)) for j in range(len(rels))]
            check(actual==expected,'Independent grade-one Hodge with physical grade-zero metric')
        else:
            laplacian=[[sum((boundary[i][k]*boundary[j][k] for k in range(len(rels))),Q(0))*lengths[j] for j in range(len(seg))] for i in range(len(seg))]
            resolvent=[[laplacian[i][j]+Q(i==j) for j in range(len(seg))] for i in range(len(seg))]
            expected=solve(resolvent,list(map(Q,original_result['contrasts'][0]['residual'])))
            check(actual==expected,'Independent grade-zero Green solve with retained physical metric')
            check(certificate['residual_verified'],'Native Green residual')
chain=query('native-chain',w,dict(schema='gadget.assays.query.v1',operation='native_query',assay='fields',query=dict(operation='chain_certificate',context='all',grade=1,coefficient_domain='rational',field='all.interval-contrast:0')))['answers'][0]['certificate']['result']
x=list(map(Q,original_result['contrasts'][0]['source_coefficients']))
columns=[]
for row in boundary:
    if rank(list(map(list,zip(*(columns+[row])))))>len(columns):columns.append(row)
gram=[[sum((a*b for a,b in zip(p,q)),Q(0)) for q in columns] for p in columns]
co=solve(gram,[sum((a*b for a,b in zip(col,x)),Q(0)) for col in columns])
gradient=[sum((a*col[i] for a,col in zip(co,columns)),Q(0)) for i in range(len(x))]
harmonic=[a-b for a,b in zip(x,gradient)]
check(list(map(Q,chain['sectors']['gradient']))==gradient,'Independent native gradient sector')
check(list(map(Q,chain['sectors']['harmonic']))==harmonic,'Independent native harmonic sector')
check(all(Q(v)==0 for v in chain['sectors']['curl']),'No invented grade-two curl')
check(chain['integral_homology']['free_rank']==original_result['source_kernel']['nullity'],'Whole interval-source homology')
check(all(sum((a*b for a,b in zip(row,harmonic)),Q(0))==0 for row in boundary),'Complete harmonic boundary vanishes')
from method_answers import qualify
checks += qualify(binary, root/'typed-answers', workspace,
                  {name: (requests[i], json.loads((root/(name+'-result.json')).read_text()))
                   for i, name in enumerate(['declared','no-control','intercept-only','dependent','ignore-strand'])})
shutil.rmtree(workspace.parent)
for res in resources:(root/res['source']['path']).unlink()
receipt=run('replay-derivation','--workspace',w,'--assay','fields','--out',root/'replay')
check(receipt['status']=='PASS','Source-only complete mathematical replay')
# A public selected ChIP protocol runs the same method, then native field actions.
public=Path(__file__).resolve().parents[1]/'examples/bioinformatics/controlled-coverage'
run('import-assays','--spec',public/'import.json','--out',root/'public')
plan=run('pipeline-prepare','--workspace',root/'public/workspace.json','--request',public/'pipeline.json')
save('public-plan.json',plan)
run('pipeline','--workspace',root/'public/workspace.json','--request',root/'public-plan.json','--out',root/'public-run')
shutil.rmtree(root/'public')
receipt=run('replay-derivation','--workspace',root/'public-run/contrast/workspace/workspace.rfcxs','--assay','coverage-fields','--out',root/'public-replay')
check(receipt['status']=='PASS','Public selected pipeline replay')
save('schema-requests.json',requests)
save('checks.json',dict(status='PASS',independent_checks=checks,actual_query_schemas=len(requests)))
print(f'{checks} independent physical geometry checks PASS; {len(requests)} actual queries')
