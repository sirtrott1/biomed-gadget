"""Independent finite search, allocation and filtration checks; retained replay inputs."""
import copy
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys


def main():
    binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    checks = 0

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(name, value):
        p = root / name
        p.write_text(json.dumps(value, indent=2) + '\n')
        return p

    def source(name, value):
        p = root / name
        p.write_text(value)
        return dict(path=name, sha256='sha256:' + hashlib.sha256(p.read_bytes()).hexdigest())

    def run(*args, fail=False):
        p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True)
        if fail:
            check(p.returncode != 0, 'invalid request rejected')
            return
        if p.returncode:
            raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)

    design = dict(question='Finite source and constraint witnesses', population='Authored sequences and cells',
                  selection='All retained inputs', limitations=['Authored verification examples'])
    resources = []
    def resource(name, data, fmt='fasta'):
        resources.append(dict(id=name, format=fmt, source=source(name, data), method='Authored contract input'))
    reference = 'AACCGGTTACGTAACCGGTT'
    queries = dict(exact='ACGT', splice='CCTT', ambiguous='ACGN', reverse='AACC', unmatched='AAAA')
    resource('reference.fa', '>chrA\n'+reference+'\n')
    resource('reads.fa', ''.join('>'+k+'\n'+v+'\n' for k,v in queries.items()))
    assembly_reads = dict(a='AACG', b='AACT', repeat='AAAA', repeat2='AAAA')
    resource('assembly.fa', ''.join('>'+k+'\n'+v+'\n' for k,v in assembly_reads.items()))
    resource('assembly.fq', ''.join('@'+k+'\n'+v+'\n+\n'+'I'*len(v)+'\n' for k,v in assembly_reads.items()), 'fastq')
    participants = lambda a,b: [dict(cell=a,exists=True,negative=True,magnitude='1'),dict(cell=b,exists=True,negative=False,magnitude='1')]
    relations = [dict(id='e'+str(i),meaning='Authored cycle boundary',participants=participants('v'+str(i),'v'+str((i+1)%4))) for i in range(4)]
    for i in range(2):
        relations.append(dict(id='f'+str(i),grade=2,meaning='Independent authored filling',participants=[dict(cell='e'+str(j),exists=True,negative=False,magnitude=str(i+1)) for j in range(4)]))
    topology = source('relations.json', json.dumps(dict(relations=relations)))
    assay = dict(id='cells', kind='quantitative',format='dense_tsv',source=source('matrix.tsv','id\tsample\n'+''.join(f'v{i}\t{i}\n' for i in range(4))),
        features=source('features.tsv','id\n'+''.join(f'v{i}\n' for i in range(4))),observations=source('observations.tsv','id\nsample\n'),
        source_namespace='710000',carrier='features',quantity='authored_value',units='units',channel='observation',method='Authored verification',
        arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA'],relations=topology)
    import_spec=dict(schema='gadget.assays.import.v1',id='finite-services',design=design,assays=[assay],resources=resources,links=[])
    run('import-assays','--spec',save('import.json',import_spec),'--out',root/'workspace')
    workspace = root/'workspace/workspace.json'
    requests, results = {}, {}
    def execute(name, request):
        request = dict(schema='gadget.assays.query.v1', **request)
        requests[name] = request
        result=run('assay-query','--workspace',workspace,'--request',save(name+'-request.json',request))
        save(name+'-result.json',result)
        results[name]=result
        return result['result']

    for ambiguity in ['literal','possible','definite']:
        req=dict(operation='reference_sequence_matches',query_resource='reads.fa',reference_resource='reference.fa',strands=['+','-'],
            ambiguity=ambiguity,maximum_mismatches=1,minimum_overhang=2,
            junctions=[dict(id='j1',contig='chrA',donor=4,acceptor=6,strand='+',evidence='Authored splice interval')],evidence='Complete finite placements')
        result=execute('matches-'+ambiguity,req)
        masks=dict(A={0},C={1},G={2},T={3},N={0,1,2,3})
        def compatible(a,b):
            return a==b if ambiguity=='literal' else (bool(masks[a]&masks[b]) if ambiguity=='possible' else masks[a]<=masks[b])
        expected=[]
        for qid,sequence in queries.items():
            for strand in ['+','-']:
                q=sequence if strand=='+' else sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
                placements=[(list(range(start,start+len(q))),((start,start+len(q)),),None) for start in range(len(reference)-len(q)+1)]
                if strand=='+':
                    for split in range(2,len(q)-1):
                        if 4>=split and 6+len(q)-split<=len(reference):
                            placements.append((list(range(4-split,4))+list(range(6,6+len(q)-split)),((4-split,4),(6,6+len(q)-split)),'j1'))
                for positions,blocks,j in placements:
                    mismatches=sum(not compatible(q[i],reference[p]) for i,p in enumerate(positions))
                    if mismatches<=1:
                        expected.append((qid,strand,blocks,j,str(mismatches)))
        actual=[(m['query_id'],m['strand'],tuple(map(tuple,m['blocks'])),m.get('junction',{}).get('id'),m['mismatches']) for m in result['matches']]
        check(actual==expected, 'every reference/splice placement, strand and mismatch')
        check(result['complete'] and len(result['queries'])==len(queries),'including unsupported query identities')

    assembly_req=dict(operation='sequence_assembly_graph',resource='assembly.fa',k=3,ambiguity='literal_tokens',evidence='Every read k-mer',
                      route=dict(start='AA',end='CG',objective='fewest_steps'))
    assembly=execute('assembly',assembly_req)
    spectrum=Counter(s[i:i+3] for s in assembly_reads.values() for i in range(len(s)-2))
    check({e['sequence']:int(e['count']) for e in assembly['edges']}==spectrum,'all literal k-mer counts')
    check(assembly['route']['sequence']=='AACG','native declared endpoint route')
    check(assembly['cycle_frame_dimension']==1,'self-loop retained as a closed direction')
    check([w['reconstructed_sequence'] for w in assembly['read_walks']]==list(assembly_reads.values()),'every read walk reconstructs')
    check(any(n['sequence']=='AA' and n['branch'] for n in assembly['nodes']),'ambiguity visible at shared prefix')
    fq=execute('assembly-fastq', dict(assembly_req,resource='assembly.fq'))
    check(fq['edges']==assembly['edges'],'FASTA/FASTQ numerical parity')

    base=dict(operation='constrained_allocations',domain='bounded_integer',variables=[dict(id='x',lower='0',upper='5'),dict(id='y',lower='0',upper='5')],
        constraints=[dict(id='total',coefficients=['1','1'],sense='eq',rhs='5'),dict(id='capacity',coefficients=['2','1'],sense='le',rhs='7')],
        objective=dict(sense='minimize',coefficients=['1','0']),units='counts',evidence='Declared coupled assignment constraints')
    result=execute('allocation-integer',base)
    check(result['minimum_or_maximum']=='0' and result['tested_candidates']==6 and result['declared_box_points']==36 and result['feasible_points_or_vertices']==3,'coupled integer feasible set after exact equality elimination')
    check([(b['lower'],b['upper']) for b in result['bounds']]==[('0','2'),('3','5')],'sharp coupled ranges')
    rational=copy.deepcopy(base); rational['domain']='bounded_rational'; rational['constraints'][1]['rhs']='13/2'; rational['objective']['sense']='maximize'
    result=execute('allocation-rational',rational)
    check(result['minimum_or_maximum']=='3/2','exact rational vertex objective')
    check(result['optimizers'][0]['values']==['3/2','7/2'],'exact rational vertex witness')
    infeasible=copy.deepcopy(base); infeasible['constraints'][0]['rhs']='11'
    check(execute('allocation-infeasible',infeasible)['status']=='infeasible','infeasibility is complete enumeration')
    tied=copy.deepcopy(base); tied['objective']['coefficients']=['1','1']
    check(len(execute('allocation-tied',tied)['optimizers'])==3,'all minimizing assignments retained')
    degenerate=copy.deepcopy(rational); degenerate['variables'][0]['upper']='0'
    check(execute('allocation-degenerate',degenerate)['optimizers'][0]['values']==['0','5'],'lower-dimensional rational polytope')
    rng=random.Random(710)
    for trial in range(16):
        req=copy.deepcopy(base)
        req['constraints']=[dict(id='linear',coefficients=[str(rng.randint(-3,3)) for _ in range(2)],sense=rng.choice(['eq','le','ge']),rhs=str(rng.randint(-5,12)))]
        req['objective']=dict(sense=rng.choice(['minimize','maximize']),coefficients=[str(rng.randint(-3,3)) for _ in range(2)])
        req['variables'][0]['lower']='-1';req['variables'][0]['upper']='3'
        result=execute(f'allocation-oracle-{trial}',req)
        c=req['constraints'][0]; sense=c['sense'];rhs=Q(c['rhs'])
        feasible=[]
        for x in itertools.product(range(-1,4),range(6)):
            lhs=sum(Q(a)*v for a,v in zip(c['coefficients'],x))
            if {'eq':lhs==rhs,'le':lhs<=rhs,'ge':lhs>=rhs}[sense]:
                feasible.append((x,sum(Q(a)*v for a,v in zip(req['objective']['coefficients'],x))))
        check(result['feasible_points_or_vertices']==len(feasible),'independent integer feasible count')
        if feasible:
            optimum=(min if req['objective']['sense']=='minimize' else max)(v for _,v in feasible)
            check(Q(result['minimum_or_maximum'])==optimum,'independent integer objective')
            check({tuple(map(Q,w['values'])) for w in result['optimizers']}=={x for x,v in feasible if v==optimum},'all independent optimizers')

    for trial in range(12):
        req=copy.deepcopy(rational)
        req['variables']=[dict(id='x',lower='-1/2',upper='7/2'),dict(id='y',lower='0',upper='5')]
        req['constraints']=[dict(id='c'+str(i),coefficients=[str(rng.randint(-3,3)) for _ in range(2)],
            sense=rng.choice(['eq','le','ge']),rhs=str(Q(rng.randint(-8,12),2))) for i in range(2)]
        req['objective']=dict(sense=rng.choice(['minimize','maximize']),coefficients=[str(rng.randint(-3,3)) for _ in range(2)])
        result=execute(f'rational-oracle-{trial}',req)
        lines=[(tuple(map(Q,c['coefficients'])),Q(c['rhs'])) for c in req['constraints']]+[((Q(1),Q(0)),Q(-1,2)),((Q(1),Q(0)),Q(7,2)),((Q(0),Q(1)),Q(0)),((Q(0),Q(1)),Q(5))]
        vertices=set()
        for (a,b),(c,d) in itertools.combinations(lines,2):
            determinant=a[0]*c[1]-a[1]*c[0]
            if not determinant:
                continue
            x=((b*c[1]-a[1]*d)/determinant,(a[0]*d-b*c[0])/determinant)
            if not (Q(-1,2)<=x[0]<=Q(7,2) and 0<=x[1]<=5):
                continue
            admitted=True
            for constraint in req['constraints']:
                lhs=sum(Q(v)*coordinate for v,coordinate in zip(constraint['coefficients'],x))
                rhs=Q(constraint['rhs'])
                if not {'eq':lhs==rhs,'le':lhs<=rhs,'ge':lhs>=rhs}[constraint['sense']]:
                    admitted=False
            if admitted:
                vertices.add(x)
        check(result['feasible_points_or_vertices']==len(vertices),'independent polygon vertex count')
        if vertices:
            objective=lambda x:sum(Q(c)*v for c,v in zip(req['objective']['coefficients'],x))
            optimum=(min if req['objective']['sense']=='minimize' else max)(map(objective,vertices))
            check(Q(result['minimum_or_maximum'])==optimum,'independent rational objective')
            check({tuple(map(Q,w['values'])) for w in result['optimizers']}=={x for x in vertices if objective(x)==optimum},'complete independently derived optimal face')
            check([(Q(b['lower']),Q(b['upper'])) for b in result['bounds']]==[(min(x[i] for x in vertices),max(x[i] for x in vertices)) for i in range(2)],'sharp rational feasible ranges')
        else:
            check(result['status']=='infeasible','rational infeasibility')

    births=[dict(grade=0,index=i,level=0) for i in range(4)]+[dict(grade=1,index=i,level=1) for i in range(4)]+[dict(grade=2,index=i,level=i+2) for i in range(2)]
    persistence=dict(operation='relational_persistence',assay='cells',partition=0,levels=['0','1','2','3'],births=births,
        filtration_meaning='Vertices, complete cycle, first filling, second independent filling',evidence='Every cell birth declared')
    result=execute('persistence',persistence)
    barcode={(b['grade'],b['birth_index'],b.get('death_index'),b['multiplicity']) for b in result['barcode']}
    check(barcode=={(0,0,1,3),(0,0,None,1),(1,1,2,1),(2,3,None,1)},'cycle death and higher cycle birth barcode')
    check(result['rank_reconstruction_verified'],'every interval reconstructs induced ranks')
    check([[g['betti'] for g in s['grades']] for s in result['stages']]==[[4,0,0],[1,1,0],[1,0,0],[1,0,1]],'all native stage Betti values')
    for inc in result['inclusions']:
        check(inc['chain_map_verified'], 'native inclusion verification')
        rank=sum(b['multiplicity'] for b in result['barcode'] if b['grade']==inc['grade'] and b['birth_index']<=inc['from'] and inc['to']<b.get('death_index',4))
        check(rank==inc['homology_image_rank'],'independent barcode reconstruction')
    empty=copy.deepcopy(persistence);empty['levels']=['-1','0','1','2','3']
    for b in empty['births']:
        b['level']+=1
    empty_result=execute('persistence-empty-first-level',empty)
    check(all(g['betti']==0 for g in empty_result['stages'][0]['grades']),'empty first level preserves all declared grades')
    invalid=copy.deepcopy(persistence);invalid['births'][-1]['level']=0
    run('assay-query','--workspace',workspace,'--request',save('invalid-closure.json',dict(schema='gadget.assays.query.v1',**invalid)),fail=True)
    invalid=copy.deepcopy(persistence);invalid['births'][-1]=invalid['births'][0]
    run('assay-query','--workspace',workspace,'--request',save('invalid-birth.json',dict(schema='gadget.assays.query.v1',**invalid)),fail=True)

    for name in ['assembly','allocation-integer','allocation-rational']:
        req=dict(schema='gadget.assays.materialize.v1',id='derived-'+name,source_namespace=str(810000+len(name)*1000),design=design,recipe=requests[name],channel='observation')
        run('derive-assays','--workspace',workspace,'--request',save(name+'-derive.json',req),'--out',root/('derived-'+name))
        shutil.move(root/'workspace',root/'original-unavailable')
        try:
            receipt=run('replay-derivation','--workspace',root/('derived-'+name)/'workspace.json','--assay','derived-'+name,'--out',root/('replayed-'+name))
            check(receipt['status']=='PASS','generation replay without original source workspace')
            save(name+'-replay.json',receipt)
        finally:
            shutil.move(root/'original-unavailable',root/'workspace')
        if name=='assembly':
            carrier=run('assay-query','--workspace',root/('derived-'+name)/'workspace.json','--request',save('assembly-native-request.json',dict(schema='gadget.assays.query.v1',operation='native_query',assay='derived-assembly',query=dict(operation='hodge',context='all',field='all.kmer_occurrence_count'))))
            save('assembly-native-result.json',carrier)
            check(len(carrier['result']['answers'])==1,'generated assembly accepted for native Hodge fields')
    bad=copy.deepcopy(requests['matches-literal']);bad['junctions'][0]['acceptor']=len(reference)+1
    run('assay-query','--workspace',workspace,'--request',save('invalid-junction.json',bad),fail=True)
    bad=copy.deepcopy(requests['allocation-integer']);bad['variables'][0]['upper']='1/2'
    run('assay-query','--workspace',workspace,'--request',save('invalid-integer.json',bad),fail=True)
    from method_answers import qualify
    checks += qualify(binary, root/'typed-answers', workspace,
                      {name: (request, results[name]) for name, request in requests.items()})
    save('verification.json',dict(schema='gadget.finite.services.contracts.v1',status='PASS',checks=checks,requests=list(requests)))
    print(f'PASS: {checks} independent finite-service checks')


if __name__=='__main__':
    main()
