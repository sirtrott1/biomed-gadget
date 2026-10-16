"""Independent complete read placement, allocation, coupling and replay oracles."""
from fractions import Fraction as Q
import copy
import hashlib
import itertools
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
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
    if fail:
        check(p.returncode != 0, f'Expected rejection: {args}')
        return p.stderr
    if p.returncode:
        raise RuntimeError(f'{args}: {p.stderr}')
    return json.loads(p.stdout)

def save(name, value):
    p = root / name
    p.write_text(json.dumps(value, indent=2)+'\n')
    return p

def file(name, data):
    p = root/name
    p.write_bytes(data.encode() if isinstance(data, str) else data)
    return dict(path=name, sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())

def resource(name, fmt, data, **kw):
    return dict(id=name, format=fmt, source=file(name,data), method='Authored independent oracle', **kw)

def admit(name, resources):
    spec=dict(schema='gadget.assays.import.v1',id=name,design=dict(question='Complete assignment ambiguity',population='Authored finite experiment',selection='Every supplied source record',limitations=['No population inference']),assays=[],resources=resources,links=[])
    manifest=root/name/'workspace.json'
    run('import-assays','--spec',save(name+'-import.json',spec),'--out',root/name)
    return manifest

def query(name, workspace, request, fail=False):
    requests.append(copy.deepcopy(request))
    result=run('assay-query','--workspace',workspace,'--request',save(name+'-request.json',request),fail=fail)
    if not fail:save(name+'-result.json',result)
    return result['result'] if not fail else result

reference={'A1':'ACCGGTTACC','A2':'ACCAATT','B':'TTTCCGA'}
raw=[['CGG','ACC','ACC','TCC','GGG'],['CGG','CGG','CGG','ACC','GGG']]
fa=''.join('>'+k+'\n'+v+'\n' for k,v in reference.items())
ref=resource('transcripts.fa','fasta',fa)
ref['reference']=dict(assembly='declared-transcripts',sequence_sha256=ref['source']['sha256'],coordinate_convention='zero_based_half_open',contigs=[dict(id=k,length=len(v)) for k,v in reference.items()])
reads=[resource('baseline.fastq','fastq',''.join(f'@b{i}\n{s}\n+\n'+('I'*len(s))+'\n' for i,s in enumerate(raw[0]))),resource('followup.fa','fasta',''.join(f'>f{i}\n{s}\n' for i,s in enumerate(raw[1])))]
workspace=admit('raw', [ref,*reads])
samples=[dict(id=name,resource=r['id'],scale='1',scale_provenance='Identity count scale') for name,r in zip(['baseline','followup'],reads)]
targets=[dict(id=k,candidates=[k]) for k in reference]+[dict(id='geneA',candidates=['A1','A2']),dict(id='all',candidates=[*reference,'__unassigned__']),dict(id='unassigned',candidates=['__unassigned__'])]
request=dict(schema='gadget.assays.query.v1',operation='expression_allocations',samples=samples,evidence_source=dict(mode='transcript_sequences',reference_resource='transcripts.fa',strands=['+'],ambiguity='literal',maximum_mismatches=0,evidence='All supplied transcript sequence windows'),targets=targets,contrasts=[dict(id='change',left='followup',right='baseline')],map_evidence='Declared A1 and A2 isoforms of geneA',kernel_certificate=True)

def compositions(n,k):
    if k == 1:
        yield [n];return
    for i in range(n+1):
        for rest in compositions(n-i,k-1):yield [i,*rest]

def class_assignments(allocation):
    candidates=allocation['candidates']
    outputs=[]
    for sample in allocation['samples']:
        choices=[]
        for c,count in zip(allocation['classes'],sample['counts']):
            tuples=[]
            for x in compositions(int(count),len(c['candidates'])):
                v=[Q(0)]*len(candidates)
                for name,value in zip(c['candidates'],x):v[candidates.index(name)]=value
                tuples.append(v)
            choices.append(tuples)
        totals=[]
        for assignment in itertools.product(*choices):
            totals.append([sum((v[i] for v in assignment),Q(0)) for i in range(len(candidates))])
        outputs.append(totals)
    return list(itertools.product(*outputs))

def oracle(result, declaration):
    allocation=result['allocation_request']; candidates=allocation['candidates']
    sample_ids=[s['id'] for s in declaration['samples']]
    indices=[[candidates.index(c) for c in t['candidates']] for t in declaration['targets']]
    tuples=class_assignments(allocation)
    def count(point,s,t):return sum((point[s][i] for i in indices[t]),Q(0))
    if 'coupling' in declaration:
        retained=[]
        for point in tuples:
            valid=True
            for c in declaration['coupling']['constraints']:
                v=sum((Q(term['coefficient'])*count(point,sample_ids.index(term['sample']),[t['id'] for t in declaration['targets']].index(term['target'])) for term in c['terms']),Q(0));rhs=Q(c['rhs'])
                valid &= v == rhs if c['sense']=='eq' else v <= rhs if c['sense']=='le' else v >= rhs
            if valid:retained.append(point)
        tuples=retained
    check(result['expression']['status'] == ('feasible' if tuples else 'infeasible'),'Feasibility oracle')
    if not tuples:
        check(not result['expression']['samples'],'Infeasible constraints cannot publish expression observations');return
    for s,sample in enumerate(declaration['samples']):
        for t,observed in enumerate(result['expression']['samples'][s]['bounds']):
            all_values=[Q(sample['scale'])*count(point,s,t) for point in tuples]
            lo,hi=min(all_values),max(all_values)
            check((Q(observed['lower']),Q(observed['upper']))==(lo,hi),'Sharp joint sample bound')
            check(observed['identifiable']==(lo==hi),'Identifiability oracle')
    for c,contrast in enumerate(declaration['contrasts']):
        l,r=sample_ids.index(contrast['left']),sample_ids.index(contrast['right'])
        for t,observed in enumerate(result['expression']['contrasts'][c]['bounds']):
            vals=[Q(declaration['samples'][l]['scale'])*count(point,l,t)-Q(declaration['samples'][r]['scale'])*count(point,r,t) for point in tuples]
            lo,hi=min(vals),max(vals)
            check((Q(observed['lower']),Q(observed['upper']))==(lo,hi),'Sharp shared-carrier contrast bound')
            direction='positive_for_every_allocation' if lo>0 else 'negative_for_every_allocation' if hi<0 else 'zero_for_every_allocation' if lo==hi==0 else 'unresolved'
            check(observed['direction']==direction,'Allocation-invariant direction')

result=query('raw',workspace,request);oracle(result,request)
for s,readset in enumerate(raw):
    evidence=result['evidence_records'][s]['evidence'];matches=evidence['search_result']['matches']
    for i,read in enumerate(readset):
        expected=[]
        for name,seq in reference.items():
            for start in range(len(seq)-len(read)+1):
                if seq[start:start+len(read)] == read:expected.append((name,start))
        witnesses=evidence['observations'][i]['match_indices']
        actual=[(matches[j]['reference_id'],matches[j]['blocks'][0][0]) for j in witnesses]
        check(actual==expected,'Every placement and its original read index')
        check(evidence['observations'][i]['source_evidence']['record_index']==i,'Original record identity')
        members=evidence['classes'][evidence['observations'][i]['class_index']]['candidates']
        check(members==[k for k in reference if any(n==k for n,p in expected)] if expected else members==['__unassigned__'],'All transcript candidates or explicit unassigned outcome')
check(result['expression']['samples'][0]['bounds'][3]['identifiable'],'Gene identifiable despite isoform ambiguity')
check(not result['expression']['samples'][0]['bounds'][0]['identifiable'],'Isoform ambiguity remains')
from method_answers import qualify
checks += qualify(binary, root/'typed-answers', workspace, {'read-expression': (request, result)})
check(len(result['evidence_records'][0]['evidence']['observations'][1]['match_indices'])==3,'Repeated transcript windows count once per source read')
scaled=copy.deepcopy(request);scaled['samples'][1]['scale']='1/2';scaled['samples'][1]['scale_provenance']='Authored positive scale; not estimated'
scaled_result=query('scaled',workspace,scaled);oracle(scaled_result,scaled)

coupled=copy.deepcopy(request);coupled['coupling']=dict(witness_representation='factorized',evidence='Authored independent coupling hypotheses',constraints=[dict(id='cross-sample',terms=[dict(sample='followup',target='A1',coefficient='1'),dict(sample='baseline',target='A1',coefficient='-1')],sense='eq',rhs='1',evidence='Declared count difference'),dict(id='baseline-capacity',terms=[dict(sample='baseline',target='A1',coefficient='1')],sense='ge',rhs='2',evidence='Declared source count lower bound'),dict(id='followup-capacity',terms=[dict(sample='followup',target='A1',coefficient='1')],sense='le',rhs='3',evidence='Declared source count upper bound')])
cr=query('coupled',workspace,coupled);oracle(cr,coupled)
check(cr['expression']['samples'][0]['bounds'][0]['lower']=='2' and cr['expression']['samples'][1]['bounds'][0]['lower']=='3','Coupling changes identifiability')
for c in cr['expression']['certificates']:
    check(c['result']['complete'],'Complete existing native solver certificate')
    check(len(c['request']['variables'])==len(cr['expression']['variable_origins']),'All original allocation coordinates')
explicit=copy.deepcopy(coupled);explicit['coupling']['witness_representation']='explicit'
er=query('explicit',workspace,explicit);oracle(er,explicit)
for c in er['expression']['certificates']:
    for w in c['result']['optimizers']:
        v=list(map(Q,w['values']));obj=c['request']['objective']
        check(sum((Q(a)*b for a,b in zip(obj['coefficients'],v)),Q(0))==Q(c['result']['minimum_or_maximum']),'Independent native optimizer objective')
        for row in c['request']['constraints']:
            value=sum((Q(a)*b for a,b in zip(row['coefficients'],v)),Q(0));rhs=Q(row['rhs'])
            check(value==rhs if row['sense']=='eq' else value<=rhs if row['sense']=='le' else value>=rhs,'All original conservation and coupling rows')
infeasible=copy.deepcopy(coupled);infeasible['coupling']['constraints'][0]['rhs']='100'
oracle(query('infeasible',workspace,infeasible),infeasible)
materialize=dict(schema='gadget.assays.materialize.v1',id='infeasible',source_namespace='900822',channel='E',recipe=infeasible,design=dict(question='Check infeasibility publication',population='Authored data',selection='Every supplied record',limitations=[]))
run('derive-assays','--workspace',workspace,'--request',save('infeasible-derive.json',materialize),'--out',root/'infeasible-derived',fail=True)
check(not (root/'infeasible-derived').exists(),'Infeasibility leaves no measured endpoint publication')


# Equivalence class files may use different candidate orders and repeated semantic
# classes: every original row stays in evidence while the declared map sums counts.
def eqtext(ids, classes):return '\n'.join([str(len(ids)),str(len(classes)),*ids,*[str(len(p))+' '+ ' '.join(str(ids.index(i)) for i in p)+' '+str(n) for p,n in classes]])+'\n'
eqs=[resource('one.eq','salmon_eq',eqtext(['A1','A2','B'],[(['A1','A2'],2),(['A1'],1),(['B'],1),(['A1','A2'],1)])),resource('two.eq','salmon_eq',eqtext(['B','A2','A1'],[(['A1'],2),(['A2'],2),(['B'],1)]))]
eq_workspace=admit('classes',eqs)
eq_request=copy.deepcopy(request);eq_request['samples']=[dict(id=n,resource=r['id'],scale=q,scale_provenance='Authored scale') for n,r,q in zip(['baseline','followup'],eqs,['1','2/3'])];eq_request['targets']=targets[:4];eq_request['evidence_source']=dict(mode='equivalence_classes',candidate_search=dict(method='Authored unweighted classes',completeness='Supplied complete declared equivalence classes',reference_sha256=ref['source']['sha256']))
oracle(query('classes',eq_workspace,eq_request),eq_request)
# Exact counts beyond JSON int64 and every long source label remain intact.
long='transcript-'+('x'*5000);big=2**70
big_workspace=admit('large-count',[resource('big.eq','salmon_eq',eqtext([long],[([long],big)]))])
big_request=copy.deepcopy(eq_request);big_request['samples']=[dict(id='sample',resource='big.eq',scale='1',scale_provenance='Identity')];big_request['targets']=[dict(id=long,candidates=[long])];big_request['contrasts']=[]
br=query('large-count',big_workspace,big_request)
check(br['expression']['samples'][0]['bounds'][0]['lower']==str(big),'Exact large count and unrestricted source label')

# Small independent randomized read/candidate and contrast experiments.
rng=random.Random(4815)
for trial in range(8):
    refs={f't{i}':''.join(rng.choice('AC') for _ in range(7)) for i in range(3)}
    data=[[''.join(rng.choice('AC') for _ in range(3)) for i in range(3)] for s in range(2)]
    rr=resource(f'transcripts{trial}.fa','fasta',''.join(f'>{k}\n{v}\n' for k,v in refs.items()))
    rr['reference']=dict(assembly='random-transcripts',sequence_sha256=rr['source']['sha256'],coordinate_convention='zero_based_half_open',contigs=[dict(id=k,length=len(v)) for k,v in refs.items()])
    rs=[resource(f'random{trial}-{s}.fa','fasta',''.join(f'>r{i}\n{v}\n' for i,v in enumerate(data[s]))) for s in range(2)]
    w=admit('random'+str(trial),[rr,*rs]);q=copy.deepcopy(request);q['evidence_source']['reference_resource']=rr['id'];q['samples']=[dict(id=str(s),resource=r['id'],scale='1/2' if s else '2/3',scale_provenance='Authored rational scale') for s,r in enumerate(rs)];q['contrasts']=[dict(id='difference',left='1',right='0')];q['targets']=[dict(id=k,candidates=[k]) for k in refs]+[dict(id='pair',candidates=['t0','t1']),dict(id='all',candidates=[*refs,'__unassigned__'])]
    answer=query('random'+str(trial),w,q);oracle(answer,q)
    for s,reads_ in enumerate(data):
        evidence=answer['evidence_records'][s]['evidence']
        for i,read in enumerate(reads_):
            expected={k for k,v in refs.items() if read in v}
            observed=set(evidence['classes'][evidence['observations'][i]['class_index']]['candidates'])
            check(observed==(expected or {'__unassigned__'}),'Independent random literal candidate universe')

# Junction compatible prealigned evidence through the same expression contract.
if run('bioinformatics')['htslib_enabled']:
    genome=resource('genome.fa','fasta','>chr\n'+'A'*100+'\n')
    frame=dict(assembly='authored-genome',sequence_sha256=genome['source']['sha256'],coordinate_convention='zero_based_half_open',contigs=[dict(id='chr',length=100)])
    gtf=''.join(f'chr\tfixture\texon\t{a+1}\t{b}\t.\t+\t.\ttranscript_id "{t}";\n' for t,exons in [('T1',[(0,4),(8,12)]),('T2',[(0,4),(6,12)]),('T3',[(0,12)])] for a,b in exons)
    ann=resource('exons.gtf','gtf',gtf,reference=frame)
    def row(q,flag,pos,cigar,seq):return f'{q}\t{flag}\tchr\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*len(seq))+'\tRG:Z:lane\n'
    header='@HD\tVN:1.6\n@SQ\tSN:chr\tLN:100\n@RG\tID:lane\tSM:s\n'
    sam=header+row('T1-splice',0,1,'4M4N4M','A'*8)+row('T2-splice',0,1,'4M2N4M','A'*8)+row('inside',0,1,'2M','AA')+row('deletion',0,1,'4M4D4M','A'*8)+row('alternative',0,1,'2M','AA')+row('alternative',256,30,'2M','AA')+row('reverse',16,1,'2M','AA')+row('filtered',1024,1,'2M','AA')
    samres=resource('aligned.sam','sam',sam,reference=frame)
    w=admit('aligned',[genome,ann,samres]);q=copy.deepcopy(request);q['samples']=[dict(id='sample',resource=samres['id'],scale='1',scale_provenance='Identity')];q['targets']=[dict(id=k,candidates=[k]) for k in ['T1','T2','T3','__unassigned__']];q['contrasts']=[];q['evidence_source']=dict(mode='annotation_alignments',annotation=ann['id'],feature_type='exon',group_attribute='transcript_id',overlap='all_aligned_bases',strand_policy='same',filter=dict(exclude_flags=1024,minimum_mapping_quality=0,unknown_mapping_quality='include'),candidate_search=dict(method='Authored retained SAM records',completeness='All supplied source mappings',reference_sha256=frame['sequence_sha256']),junction_policy='declared_exon_chain')
    answer=query('aligned',w,q);oracle(answer,q)
    evidence=answer['evidence_records'][0]['evidence'];allocation=evidence['allocation_request'];by_query={s['query']:allocation['classes'][s['class_index']]['candidates'] for s in evidence['segments']}
    check(by_query['T1-splice']==['T1'] and by_query['T2-splice']==['T2'],'CIGAR junctions require exact declared adjacent exons')
    check(by_query['deletion']==['T3'],'Deletion across an intron is not a splice junction')
    check(by_query['alternative']==['T1','T2','T3','__unassigned__'],'Retained incompatible alternative remains unassigned possibility')
    check(by_query['reverse']==['__unassigned__'] and evidence['excluded_record_indices']==[7],'Strand and source filter evidence')
    ignored=copy.deepcopy(q);ignored['evidence_source']['junction_policy']='ignore';other=query('ignored-junctions',w,ignored)
    check(other['evidence_records'][0]['evidence']['allocation_request']['classes'] != allocation['classes'],'Declared overlap and exact transcript compatibility are different contracts')
    invalid=copy.deepcopy(q);invalid['evidence_source']['overlap']='any_aligned_base';query('invalid-junction-policy',w,invalid,fail=True)
    # Exceed both former 10,000 feature and ten million pair comparison limits.
    count=10003;seq=resource('large-genome.fa','fasta','>chr\n'+'A'*(count*2+1)+'\n');fr=copy.deepcopy(frame);fr['sequence_sha256']=seq['source']['sha256'];fr['contigs'][0]['length']=count*2+1
    large=resource('large-exons.gtf','gtf',''.join(f'chr\tfixture\texon\t{2*i+1}\t{2*i+1}\t.\t+\t.\ttranscript_id "t{i}";\n' for i in range(count)),reference=fr)
    sh='@HD\tVN:1.6\n@SQ\tSN:chr\tLN:'+str(count*2+1)+'\n@RG\tID:lane\tSM:s\n'
    aligned=resource('large-aligned.sam','sam',sh+''.join(row(f'r{i}',0,2*i+1,'1M','A') for i in range(1001)),reference=fr)
    lw=admit('large-annotation',[seq,large,aligned]);fq=dict(schema='gadget.assays.query.v1',operation='feature_counts',resource=aligned['id'],annotation=large['id'],feature_type='exon',group_attribute='transcript_id',unit='sequenced_segment',overlap='all_aligned_bases',strand_policy='same',filter=q['evidence_source']['filter'],candidate_search=dict(method='Complete authored source',completeness='Every supplied mapping',reference_sha256=fr['sequence_sha256']),junction_policy='declared_exon_chain')
    large_result=query('large-annotation',lw,fq)
    check(len(large_result['features'])==count,'All features above former ceiling')
    check(len(large_result['segments'])==1001 and large_result['execution']['candidate_feature_comparisons']==1001,'Complete interval join excludes disjoint geometry without pair ceiling')
    check(len(large_result['allocation_certificate']['result']['samples'][0]['bounds'])==count,'Complete target coordinate output')

# Scientist recipes and fully detached source replay.
settings=dict(schema='gadget.pipeline.default.v1',profile='rna_seq',question='Which gene counts and changes are identified while isoforms remain ambiguous?',goals=[dict(method='expression_allocations',id='expression',parameters={k:v for k,v in request.items() if k not in ['schema','operation']},output=dict(id='expression',source_namespace='900821',channel='E',design=dict(question='Inspect sharp assignment endpoints',population='Authored paired samples',selection='Every supplied read',limitations=['Bounds retain candidate ambiguity'])))])
plan=run('pipeline-prepare','--workspace',workspace,'--request',save('pipeline.json',settings));save('plan.json',plan)
run('pipeline','--workspace',workspace,'--request',root/'plan.json','--out',root/'pipeline-results')
lineage_workspace=root/'pipeline-results'/'expression'/'workspace.rfcxs'
# Read the emitted result index instead of assuming an internal directory layout.
paths=list((root/'pipeline-results').rglob('workspace.rfcxs'))
check(len(paths)==1,'Single materialized expression workspace')
lineage_workspace=paths[0]
run('assay-query','--workspace',lineage_workspace,'--request',save('derived-inspect.json',dict(schema='gadget.assays.query.v1',operation='native_query',assay='expression',query=dict(operation='inspect'))))
shutil.rmtree(workspace.parent)
for r in [ref,*reads]:(root/r['source']['path']).unlink()
replayed=run('replay-derivation','--workspace',lineage_workspace,'--assay','expression','--out',root/'replayed')
check(replayed['status']=='PASS','Replay reconstructs every generating result after original source removal')
# Invalid declared identities and constraints must fail, never guess mappings.
bad=copy.deepcopy(eq_request);bad['targets'][0]['candidates']=['absent'];query('bad-target',eq_workspace,bad,fail=True)
bad=copy.deepcopy(eq_request);bad['samples'][1]['resource']=bad['samples'][0]['resource'];query('duplicate-source',eq_workspace,bad,fail=True)
bad=copy.deepcopy(eq_request);bad['samples'][0]['scale']='0';query('zero-scale',eq_workspace,bad,fail=True)
bad=copy.deepcopy(eq_request);bad['contrasts'][0]['left']='missing';query('bad-contrast',eq_workspace,bad,fail=True)
# Optional Python schema validation is separate from the C++ admission oracle.
try:
    import jsonschema
except ImportError:
    schema_checks=0
else:
    schema_root=Path(__file__).resolve().parents[1]/'schemas'
    schema=json.loads((schema_root/'assays-query.schema.json').read_text())
    resolver=jsonschema.RefResolver(base_uri=schema_root.as_uri()+'/',referrer=schema)
    validator=jsonschema.Draft202012Validator(schema,resolver=resolver)
    for r in requests:validator.validate(r)
    schema_checks=len(requests)
save('acceptance.json',dict(status='PASS',checks=checks,schema_requests=schema_checks,requests=len(requests),complete_assignment_oracles=True,detached_source_replay=True))
print(json.dumps(dict(status='PASS',checks=checks,requests=len(requests),schema_requests=schema_checks)))
