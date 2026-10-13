"""Independent source-incidence, path, rational and integral chain witnesses."""
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from research_cli import columns, dot, mv, projection, rref


def dense(m):
    a=[[Q(0) for _ in range(m['columns'])] for _ in range(m['rows'])]
    for i,j,x in m['entries']:a[i][j]=Q(x)
    return a


def main():
    binary,root=Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve()
    if root.exists():shutil.rmtree(root)
    root.mkdir(parents=True);checks=0;requests={};results={}
    def check(v,msg):
        nonlocal checks
        checks+=1
        if not v:raise AssertionError(msg)
    def save(name,d):
        p=root/name;p.write_text(json.dumps(d,indent=2)+'\n');return p
    def file(name,body):
        p=root/name;p.write_text(body)
        return dict(path=name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())
    def run(*args,fail=False):
        p=subprocess.run([str(binary),*map(str,args)],text=True,capture_output=True)
        if fail:check(p.returncode!=0,'invalid structural declaration rejects');return
        if p.returncode:raise RuntimeError(str(args)+': '+p.stderr)
        return json.loads(p.stdout)
    design=dict(question='Relations, retained incidence cancellation and coefficient-dependent classes',population='Authored finite complexes',selection='Every declared cell and occurrence',limitations=['Source-bound mathematical examples'])
    participant=lambda c,q='1',exists=True:dict(cell=c,exists=exists,negative=Q(q)<0,magnitude=str(abs(Q(q))))
    edges=[dict(id='e'+str(i),meaning='Oriented square edge',participants=[participant('v'+str(i),'-1'),participant('v'+str((i+1)%4))]) for i in range(4)]
    def face(coeff):return dict(id='f',grade=2,meaning='Authored higher relation',participants=[participant('e'+str(i),coeff) for i in range(4)])
    branching=[dict(id='branch',meaning='Primary branching relation',participants=[participant('v0','-1')]+[participant('v'+str(i),'1/3') for i in range(1,4)])]
    branching += [dict(id='pair'+str(i),meaning='Declared comparison relation',participants=[participant('v0','-1'),participant('v'+str(i))]) for i in range(1,4)]
    branching += [dict(id='fill',grade=2,meaning='Authored nested branching relation',participants=[participant('branch','-1')]+[participant('pair'+str(i),'1/3') for i in range(1,4)])]
    inactive=[dict(id='loop',meaning='Retained cancellation with an inactive occurrence',participants=[participant('v0','-1'),participant('v0','1'),participant('v0','7',False)])]
    assays=[]
    for name,rels,vertices in [('square',edges,4),('filled',edges+[face('1')],4),('torsion',edges+[face('2')],4),('branching',branching,4),('inactive',inactive,1)]:
        assays.append(dict(id=name,kind='quantitative',format='dense_tsv',source=file(name+'.tsv','id\tsample\n'+''.join(f'v{i}\t{i+1}\n' for i in range(vertices))),
            features=file(name+'-features.tsv','id\n'+''.join(f'v{i}\n' for i in range(vertices))),observations=file(name+'-observations.tsv','id\nsample\n'),source_namespace=str(920000+len(assays)*1000),carrier='features',quantity='authored_value',units='units',channel='observation',method='Authored contract input',arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA'],relations=file(name+'-relations.json',json.dumps(dict(relations=rels)))))
    resource=dict(id='reads',format='fasta',source=file('reads.fa','>a\nAACG\n>b\nAACT\n>r\nAAAA\n>r2\nAAAA\n'),method='Read-walk test')
    run('import-assays','--spec',save('import.json',dict(schema='gadget.assays.import.v1',id='structure',design=design,assays=assays,resources=[resource],links=[])),'--out',root/'workspace')
    workspace=root/'workspace/workspace.json'
    assembly=dict(schema='gadget.assays.query.v1',operation='sequence_assembly_graph',resource='reads',k=3,ambiguity='literal_tokens',evidence='All read occurrences')
    derive=dict(schema='gadget.assays.materialize.v1',id='assembly',source_namespace='929000',design=design,recipe=assembly,channel='observation')
    run('derive-assays','--workspace',workspace,'--request',save('derive.json',derive),'--out',root/'assembly')
    def execute(name,assay,query):
        ws=root/'assembly/workspace.rfcxs' if assay=='assembly' else workspace
        req=dict(schema='gadget.assays.query.v1',operation='native_query',assay=assay,query=query)
        requests[name]=req;answer=run('assay-query','--workspace',ws,'--request',save(name+'-request.json',req))
        save(name+'-result.json',answer);results[name]=answer
        if len(sys.argv) > 3:
            native = subprocess.run([str(Path(sys.argv[3]).resolve()), str(ws),
                str(save(name+'-native-queries.json', {name: query})), str(root/(name+'-native')), assay],
                text=True, capture_output=True)
            if native.returncode: raise RuntimeError(name+': '+native.stderr)
            check(json.loads(native.stdout)['status']=='PASS', 'actual native structure and portable replay')
        return answer['result']['answers'][0]['certificate']['result']
    def chain(name,assay,values=None,domain='rational',field=None,grade=1):
        req=dict(operation='chain_certificate',context='all',grade=grade,coefficient_domain=domain)
        if field:req['field']=field
        else:req['values']=list(map(str,values))
        out=execute(name,assay,req)
        x=list(map(Q,out['rational_coordinates']));b=dense(out['boundary_matrix']);u=dense(out['upper_boundary_matrix']);f=out['facts']
        check(mv(b,x)==list(map(Q,out['boundary'])),'independent native boundary action')
        check(f['cycle']==not_any(mv(b,x)),'cycle membership from original boundary')
        g=projection(b,x) if b else [Q(0)]*len(x);c=projection(columns(u),x);h=[a-b-d for a,b,d in zip(x,g,c)]
        check(g==list(map(Q,out['sectors']['gradient'])),'independent gradient projection')
        check(c==list(map(Q,out['sectors']['curl'])),'independent curl projection')
        check(h==list(map(Q,out['sectors']['harmonic'])),'independent harmonic projection')
        check(dot(g,c)==dot(g,h)==dot(c,h)==0,'exact sector orthogonality')
        check(Q(out['quadrance']['total'])==dot(x,x)==sum(Q(v) for k,v in out['quadrance'].items() if k!='total'),'sector quadrance identity')
        check(list(map(Q,out['complete_field_coordinates']['reconstructed']))==x,'complete coordinates invert exactly')
        if f['rational_boundary']:check(mv(u,list(map(Q,out['boundary_preimage'])))==x,'full rational higher-chain witness')
        sums=[Q(0)]*len(b)
        for o in out['occurrences']:
            value=Q(o['rational_coefficient'])*x[o['upper']] if o['exists'] else Q(0)
            check(value==Q(o['contribution']),'signed source occurrence contribution')
            sums[o['lower']]+=value
        check(sums==list(map(Q,out['boundary'])),'occurrence gluing keeps every original slot')
        ih=out['integral_homology'];basis=dense(ih['cycle_basis']);boundaries=dense(ih['boundaries_in_cycle_coordinates']);left=dense(ih['smith_left']);right=dense(ih['smith_right']);diagonal=dense(ih['smith_diagonal'])
        def mm(a,b,n):return [[dot(row,col) for col in (columns(b) if b else [[] for _ in range(n)])] for row in a]
        check(mm(mm(left,boundaries,ih['smith_right']['rows']),right,ih['smith_right']['columns'])==diagonal,'independent full Smith factorization')
        return out
    not_any=lambda a:not any(a)
    loop=chain('self-loop','assembly',[0,0,0,1]);check(loop['facts']['cycle'] and not loop['facts']['zero_chain'] and not loop['facts']['rational_boundary'],'a nonzero loop is a cycle, not a point or filled boundary')
    check(loop['facts']['cancelled_lower_coordinates']==1 and loop['facts']['nonzero_occurrence_contributions']==2,'two retained opposite contributions, zero lower chain')
    counts=chain('assembly-count-field','assembly',field='all.kmer_occurrence_count');check(counts['rational_coordinates']==['2','1','1','4'] and counts['sectors']['harmonic']==['0','0','0','4'],'measured self-loop count remains independent of cancellation')
    chain('zero-chain','assembly',[0,0,0,0]);chain('route-chain','assembly',[1,1,0,0]);chain('vertex-chain','square',[1,0,0,0],grade=0)
    square=chain('unfilled-cycle','square',[1]*4);check(square['integral_homology']['free_rank']==1 and square['facts']['nonzero_integral_homology_class'],'unfilled cycle free class')
    filled=chain('filled-cycle','filled',[1]*4);check(filled['sectors']['curl']==['1']*4 and filled['facts']['integral_boundary'],'authored higher relation turns the closed field into curl')
    torsion=chain('torsion-cycle','torsion',[1]*4,domain='retained_lattice');check(torsion['integral_homology']['torsion_invariants']==['2'] and torsion['facts']['rational_boundary'] and not torsion['facts']['integral_boundary'],'Z/2 class erased by rationalization')
    double=chain('twice-torsion-cycle','torsion',[2]*4,domain='retained_lattice');check(double['facts']['integral_boundary'],'twice the torsion generator has an integral filling')
    half=chain('half-cycle','torsion',[Q(1,2)]*4);check(not half['facts']['belongs_to_retained_lattice'] and half['facts']['integral_boundary'] is None,'fractional rational chains do not silently become integral')
    branch=chain('branching-cycle','branching',[-1,Q(1,3),Q(1,3),Q(1,3)]);check(branch['facts']['cycle'] and branch['facts']['rational_boundary'],'multiway and nested branching closure')
    check(dense(branch['boundary_matrix'])==[[-1,-1,-1,-1],[Q(1,3),1,0,0],[Q(1,3),0,1,0],[Q(1,3),0,0,1]],'branching source magnitudes are preserved without pairwise expansion')
    check(dense(branch['upper_boundary_matrix'])==[[-1],[Q(1,3)],[Q(1,3)],[Q(1,3)]],'nested branching source magnitudes are preserved')
    inactive_result=chain('inactive-occurrence','inactive',[1]);check(inactive_result['facts']['retained_occurrences']==3 and inactive_result['facts']['active_occurrences']==2 and inactive_result['occurrences'][2]['contribution']=='0','inactive asserted slots stay present')
    paths=[dict(id='p',steps=[dict(index=i,coefficient='1') for i in [0,1]]),dict(id='q',steps=[dict(index=i,coefficient='-1') for i in [3,2]])]
    for assay in ['square','filled','torsion']:
        p=execute('paths-'+assay,assay,dict(operation='realized_path_comparison',context='all',grade=1,paths=paths))
        check([s['index'] for s in p['paths'][1]['steps']]==[3,2],'signed return path retains its composable traversal order')
        check(p['realization_rank']==2 and p['endpoint_rank']==1 and p['directions_lost_by_endpoint_view']==1,'endpoint view discards a retained independent direction')
        check(dense(p['shared_field_form'])==[[2,0],[0,2]] and dense(p['endpoint_form'])==[[2,2],[2,2]],'complete signed field and endpoint forms differ')
        check(p['pairs'][0]['equal_boundaries'] and not p['pairs'][0]['equal_fields'],'equal endpoints retain distinct paths')
        sector='harmonic' if assay=='square' else 'curl';check(p['pairs'][0]['sectors'][sector]==['1']*4,'closed path difference sector depends on asserted higher structure')
    repeats=[dict(id='three-loops',steps=[dict(index=3,coefficient='1')]*3),dict(id='none',steps=[])]
    p=execute('paths-repeats','assembly',dict(operation='realized_path_comparison',context='all',grade=1,paths=repeats));check(len(p['paths'][0]['steps'])==3 and p['paths'][0]['field']==['0','0','0','3'],'ordered repeats remain separate from their field realization')
    births=[dict(grade=0,index=i,level=0) for i in range(4)]+[dict(grade=1,index=i,level=1) for i in range(4)]+[dict(grade=2,index=0,level=2)]
    req=dict(schema='gadget.assays.query.v1',operation='relational_persistence',assay='torsion',partition=0,levels=['0','1','2'],births=births,filtration_meaning='Vertices, loop, twice-loop filling',evidence='Every authored birth',lattice_snapshots=True)
    requests['lattice-filtration']=req;p=run('assay-query','--workspace',workspace,'--request',save('lattice-filtration-request.json',req));save('lattice-filtration-result.json',p);results['lattice-filtration']=p
    check(p['result']['stages'][2]['grades'][1]['integral_homology']==dict(free_rank=0,torsion_invariants=['2'],retained_cell_scales=['1']*4),'retained torsion survives the rational barcode death')
    for suffix,change in [('fractional-lattice',dict(coefficient_domain='retained_lattice',values=['1/2']*4)),('missing-coordinate',dict(values=['1'])),('extra-member',dict(guess='anything'))]:
        bad=dict(requests['unfilled-cycle']['query'],**change);run('query','--study',root/'workspace/assays/0/partitions/0/study.json','--query',save('invalid-'+suffix+'.json',bad),fail=True)
    invalid=dict(requests['assembly-count-field']['query'],grade=0);run('assay-query','--workspace',root/'assembly/workspace.rfcxs','--request',save('invalid-field-grade.json',dict(schema='gadget.assays.query.v1',operation='native_query',assay='assembly',query=invalid)),fail=True)
    check(all(o.get('evidence') for o in loop['occurrences']),'source occurrence spans available for review')
    def workflow(name,value):
        inspection=run('assay-query','--workspace',root/'assembly/workspace.rfcxs','--request',save(name+'-query.json',dict(schema='gadget.assays.query.v1',operation='native_query',assay='assembly',query=dict(operation='inspect'))))
        binding=inspection['result']['answers'][0]['source_index']['source_digest']
        claim=dict(path=['facts','zero_chain'],value_type='boolean',comparison='equal',value=value)
        entry=dict(id='whole',type='sequence',steps=[dict(id='chain',type='query',query=requests['self-loop']['query']),dict(id='check',type='check_claims',**{'from':'chain'},claims=[claim]),dict(id='decision',type='branch',**{'from':'check'},condition=dict(path=['contradicted'],value_type='count',comparison='equal',value=0),on_true=dict(id='agreed',type='query',query=dict(operation='coverage',field='all.kmer_occurrence_count')),on_false=dict(id='inspect',type='query',query=dict(operation='inspect')))])
        w=dict(schema='gadget.scientific.workflow.v1',id=name,source_digest=binding,limits=dict(depth=8,nodes=8,model_calls=0,parallel_workers=2),definitions=[],entry=entry)
        answer=run('workflow','--workspace',root/'assembly/workspace.rfcxs','--assay','assembly','--request',save(name+'-workflow.json',w));save(name+'-workflow-result.json',answer)
        check(answer['status']=='complete','structural claims drive an owning workflow')
        events=answer['trace'];event=next(e for e in events if e['type']=='branch');check(event['branch_taken']==('on_true' if not value else 'on_false'),'zero boundary does not certify a zero chain')
    workflow('correct-structural-claim',False);workflow('wrong-structural-claim',True)
    from method_answers import qualify
    checks += qualify(binary, root/'typed-answers', workspace,
                      {name: (requests[name], results[name]) for name in ['paths-square','paths-filled','paths-torsion']})
    save('cases.json',dict(requests=requests,checks=checks));save('verification.json',dict(status='PASS',independent_checks=checks,queries=len(requests)))
    print(f'PASS: {checks} independent structural checks across {len(requests)} requests')


if __name__=='__main__':main()
