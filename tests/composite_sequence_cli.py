"""Independent exact cover, exon-chain and observed-carrier qualification."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys


def main():
    binary, root = (Path(x).resolve() for x in sys.argv[1:3])
    shutil.rmtree(root,ignore_errors=True);root.mkdir(parents=True)
    rng=random.Random(41603);checks=0;requests=[]
    def check(value,label):
        nonlocal checks
        checks+=1
        if not value:raise AssertionError(label)
    def save(name,data):
        p=root/name;p.write_text(json.dumps(data,indent=2)+'\n');return p
    def run(*args):
        p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
        if p.returncode:raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)
    def workspace(name,resources):
        specs=[]
        for id,records in resources.items():
            p=root/(name+'-'+id+'.fa');p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in records.items()))
            specs.append(dict(id=id,format='fasta',method='Authored exact source contracts',source=dict(path=p.name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())))
        spec=dict(schema='gadget.assays.import.v1',id=name,design=dict(question='Complete source-sized sequence contracts',population='Authored sequences',selection='All source records',limitations=['Nonclinical verification fixture']),assays=[],resources=specs,links=[])
        run('import-assays','--spec',save(name+'-import.json',spec),'--out',root/name)
        return root/name/'workspace.json'
    def query(name,ws,req):
        req=dict(req,schema='gadget.assays.query.v1');requests.append(req)
        r=run('assay-query','--workspace',ws,'--request',save(name+'-request.json',req));save(name+'-result.json',r);return r['result']
    def cover_request(resource='reads',minimum=2):
        return dict(operation='read_overlap_reconstruction',resource=resource,minimum_overlap=minimum,objective='minimum_suffix_prefix_read_cover',evidence='Every source read and all minimum covers')
    def route(route,carrier,cost):
        cursor=carrier['source_vertex'];counts=Counter(route['ordered_relations']);total=0
        for e in route['ordered_relations']:
            a,b,c=carrier['steps'][e];check(cursor==a,'source ordered route continuity');cursor=b;total+=int(c)
        check(cursor==carrier['sink_vertex'] and total==cost,'exact route endpoint and cost')
        field={v['relation']:int(v['coefficient']) for v in route['realized_field']}
        check(field==counts,'ordered history realizes its complete field')
    def verify_cover(reads,minimum,r,expected=None):
        if expected is None:
            candidates=[]
            for order in itertools.permutations(range(len(reads))):
                options=[[n for n in range(minimum,min(len(reads[a]),len(reads[b]))+1) if reads[a][-n:]==reads[b][:n]] for a,b in zip(order,order[1:])]
                for shared in itertools.product(*options):
                    spelling=reads[order[0]]+''.join(reads[b][n:] for b,n in zip(order[1:],shared))
                    candidates.append((len(spelling),spelling,order,shared))
            best=min((x[0] for x in candidates),default=None)
            expected={(s,order,shared) for length,s,order,shared in candidates if length==best}
        else:best=len(next(iter(expected))[0]) if expected else None
        actual={(l['spelling'],tuple(l['read_order']),tuple(l['overlap_lengths'])) for l in r['layouts']}
        check(actual==expected and len(actual)==len(r['layouts']),'every independently minimized source layout and spelling')
        check((None if r['minimum_length'] is None else int(r['minimum_length']))==best,'independent minimum length')
        check(len(r['source_reads'])==len(reads),'all original distinct source reads')
        check(r['identifiable_sequence']==(len({s for s,_,_ in expected})==1),'identifiability from complete minimum set')
        for l in r['layouts']:
            route(l['route'],r['carrier'],best)
            for p in l['gluing']['placements']:
                check(l['spelling'][p['start']:p['end']]==reads[p['source_read']],'each categorical restriction restores its original read')
            check(l['gluing']['complete_overlap_gluing_verified'] and l['gluing']['every_read_restriction_verified'],'native complete categorical gluing')
        for state in r['carrier']['states']:
            mask=int(state['source_read_mask']);check({i for i in range(len(reads)) if mask>>i&1}==set(state['source_reads']),'exact original subset identities')
            check(state['last_read'] in state['source_reads'],'last read belongs to the original subset')
        representatives=[Counter(r['layouts'][i]['route']['ordered_relations']) for i in r['path_comparison']['representative_layouts']]
        for i,a in enumerate(representatives):
            for j,b in enumerate(representatives):
                check(int(r['path_comparison']['shared_field_gram'][i][j])==sum(v*b.get(e,0) for e,v in a.items()),'independent full shared-field Gram')
        check(r['complete'],'complete declared cover')
    for i in range(10):
        reads=[''.join(rng.choice('ACGT') for _ in range(rng.randint(2,5))) for _ in range(rng.randint(2,5))]
        ws=workspace('random-cover-'+str(i),dict(reads={str(j):v for j,v in enumerate(reads)}))
        req=cover_request(minimum=1);verify_cover(reads,1,query('random-cover-'+str(i),ws,req))
    reads=['AAC','ACG','CGA','GAC','ACT']
    ws=workspace('dead-suffix',dict(reads={str(i):v for i,v in enumerate(reads)}))
    result=query('dead-suffix',ws,cover_request());verify_cover(reads,2,result)
    check(any(s.get('minimum_suffix_length','missing') is None for s in result['carrier']['states']),'unreachable suffix is absent, never a supported zero')
    check(result['prefix_dual_verified'] and result['viable_suffix_dual_verified'],'native complete prefix and viable suffix certificates')
    check(any(not s['eligible_start'] for s in result['start_reachability']),'native exclusion witnesses exercised')
    n=130;genome=''.join(rng.choice('ACGT') for _ in range(n+11));reads=[genome[i:i+12] for i in range(n)]
    check(all([length for length in range(11,13) if reads[a][-length:]==reads[b][:length]]==([11] if b==a+1 else []) for a in range(n) for b in range(n) if a!=b),'independent structured overlap graph is one exact chain')
    ws=workspace('long-chain',dict(reads={str(i):v for i,v in enumerate(reads)}))
    req=cover_request(minimum=11);result=query('long-chain',ws,req)
    verify_cover(reads,11,result,{(genome,tuple(range(n)),(11,)*(n-1))})
    check(len(result['carrier']['states'])==n and result['native_route_calls']==1,'130-source chain retains one complete native route with all source states')
    check(int(result['carrier']['states'][-1]['source_read_mask'])==(1<<n)-1,'source domain exceeds machine-word masks exactly')
    long=['A'*140+'C'*140,'C'*140+'G'*140]
    ws=workspace('long-overlap',dict(reads={str(i):v for i,v in enumerate(long)}))
    result=query('long-overlap',ws,cover_request(minimum=140));verify_cover(long,140,result)
    duplicate=['ACGTAC']*7;ws=workspace('many-layouts',dict(reads={str(i):v for i,v in enumerate(duplicate)}))
    result=query('many-layouts',ws,cover_request())
    expected={(duplicate[0],order,(6,)*6) for order in itertools.permutations(range(7))}
    verify_cover(duplicate,2,result,expected);check(len(result['layouts'])==5040,'all layouts beyond former 4096 policy')
    variable=['A'+c+'A' for c in 'CGTRY'];ws=workspace('many-spellings',dict(reads={str(i):v for i,v in enumerate(variable)}))
    result=query('many-spellings',ws,cover_request(minimum=1));verify_cover(variable,1,result)
    check(len(result['sequences'])==120,'all minimum spellings beyond former 32 policy')

    def splice_request(**kw):
        return dict(dict(operation='sequence_splice_discovery',query_resource='queries',reference_resource='reference',strands=['+'],minimum_exon_bases=2,minimum_intron_bases=4,maximum_intron_bases=6,maximum_junctions=6,motifs=['GT-AG'],evidence='All exact declared exon chains'),**kw)
    exons=[''.join(rng.choice('ACGT') for _ in range(36)) for _ in range(6)]
    q=''.join(exons);ref='GTCCAG'.join(exons)
    ws=workspace('six-exons',dict(queries=dict(q=q),reference=dict(r=ref)))
    req=splice_request(minimum_exon_bases=36,maximum_junctions=5)
    result=query('six-exons',ws,req)
    check(len(result['matches'])==1 and len(result['matches'][0]['junctions'])==5,'complete five-junction source chain beyond former four policy')
    match=result['matches'][0];check(''.join(ref[b['start']:b['end']] for b in match['blocks'])==q,'all original exons reconstruct the query')
    route(match['route'],result['carriers'][match['carrier']],5)
    exons=[''.join(rng.choice('ACGT') for _ in range(129)) for _ in range(2)]
    ws=workspace('large-exons',dict(queries=dict(q=''.join(exons)),reference=dict(r='GTCCAG'.join(exons))))
    req=splice_request(minimum_exon_bases=129,maximum_junctions=1)
    result=query('large-exons',ws,req);check(len(result['matches'])==1 and len(result['matches'][0]['junctions'])==1,'258-base query and 129-base full overhang from intrinsic-sized seeds')
    qs={str(i):'ACG' for i in range(257)};ws=workspace('many-splice-queries',dict(queries=qs,reference=dict(r='CCCCCC')))
    motifs=[x[:2]+'-'+x[2:] for x in map(''.join,itertools.islice(itertools.product('ACGT',repeat=4),17))]
    result=query('many-splice-queries',ws,splice_request(motifs=motifs,maximum_intron_bases=64000001))
    check(len(result['queries'])==257 and not result['matches'],'all queries, 17 declared motifs and source-sized intron bound admitted')
    ws=workspace('many-exons',dict(queries=dict(q='A'*200+'CC'),reference=dict(r='A'*20)))
    result=query('many-exons',ws,splice_request(maximum_intron_bases=64000001))
    check(int(result['work']['exon_candidates'])>20000 and not result['matches'],'complete partial exon carrier beyond former 20000 policy, without a full query match')

    def assembly_request(k):return dict(operation='sequence_assembly_graph',resource='reads',k=k,ambiguity='literal_tokens',evidence='Complete authored occurrence source')
    def verify_assembly(reads,k,r):
        counts=Counter(s[i:i+k] for s in reads for i in range(len(s)-k+1))
        check({e['sequence']:int(e['count']) for e in r['edges']}==counts,'independent complete original k-mer counts')
        check(r['source_occurrences']==sum(counts.values())==r['count_operator_nonzeros'],'native count map retains exactly one original occurrence membership')
        check(len(r['read_walks'])==len(reads),'every original source read walk')
        for i,w in enumerate(r['read_walks']):
            check(w['reconstructed_sequence']==reads[i] and len(w['steps'])==len(reads[i])-k+1,'each source walk reconstructs its full source')
        # Independent component/rank oracle is test-only; production uses native rank.
        parent=list(range(len(r['nodes'])))
        def find(a):
            while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
            return a
        for e in r['edges']:
            a,b=find(e['source']),find(e['target']);parent[a]=b
        rank=len(parent)-len({find(i) for i in range(len(parent))})
        check(r['cycle_frame_dimension']==len(r['edges'])-rank,'independent complete incidence rank and cycle dimension')
        if len(r['nodes'])==1:
            relation=r['relation_declaration']['relations'][0]
            check(len(relation['participants'])==2 and {p['negative'] for p in relation['participants']}=={True,False},'self-loop retains both signed source occurrences')
    reads=['A'*30]*4100;ws=workspace('many-observed-occurrences',dict(reads={str(i):s for i,s in enumerate(reads)}))
    result=query('many-observed-occurrences',ws,assembly_request(2));verify_assembly(reads,2,result)
    check(result['source_occurrences']==118900,'source read and occurrence counts exceed former assembly policies')
    def code(i):
        chars=[]
        for _ in range(7):chars.append('ACGT'[i%4]);i//=4
        return ''.join(chars)
    reads=['AC'+code(i)+'GT' for i in range(10003)]
    ws=workspace('many-kmers',dict(reads={str(i):s for i,s in enumerate(reads)}))
    result=query('many-kmers',ws,assembly_request(11));verify_assembly(reads,11,result)
    check(len(result['edges'])==10003 and len(result['nodes'])==20006,'complete distinct source edges and nodes beyond former policies')
    reads=['A'*129+'C','C'*129+'G'];ws=workspace('wide-k',dict(reads={str(i):s for i,s in enumerate(reads)}))
    result=query('wide-k',ws,assembly_request(129));verify_assembly(reads,129,result)

    reads=['A'*5000+'C'];ws=workspace('long-route-labels',dict(reads=dict(q=reads[0])))
    req=dict(assembly_request(5001),route=dict(start='A'*5000,end='A'*4999+'C',objective='fewest_steps'))
    result=query('long-route-labels',ws,req);verify_assembly(reads,5001,result)
    check(result['route']['found'] and result['route']['sequence']==reads[0],'complete source-sized k-mer endpoint labels beyond former 4096-string policy')

    example=Path(__file__).resolve().parents[1]/'examples/bioinformatics/composite-sequences'
    shutil.copytree(example,root/'portable-source')
    run('import-assays','--spec',root/'portable-source/import.json','--out',root/'portable')
    ws=root/'portable/workspace.json'
    plan=run('pipeline-prepare','--workspace',ws,'--request',example/'pipeline.json')
    result=run('pipeline','--workspace',ws,'--request',save('portable-plan.json',plan),'--out',root/'portable-run')
    check(len(result['results'])==3 and result['dispatcher']=='rfc::runtime::WorkGraph','all three source-selected methods execute through the native work graph')
    recipe=dict(assembly_request(3),schema='gadget.assays.query.v1');requests.append(recipe)
    derive=dict(schema='gadget.assays.materialize.v1',id='observed-sequence-structure',source_namespace='416030',channel='observation',design=dict(question='Observe full literal k-mer structure',population='Authored reads',selection='All source occurrences',limitations=['Method example']),recipe=recipe)
    run('derive-assays','--workspace',ws,'--request',save('derive.json',derive),'--out',root/'derived')
    withheld=root/'source-unavailable';withheld.mkdir()
    for name in ['portable','portable-source']:shutil.move(root/name,withheld/name)
    try:
        r=run('replay-derivation','--workspace',root/'derived/workspace.json','--assay','observed-sequence-structure','--out',root/'replayed')
        check(r['status']=='PASS','complete observed assembly replays with original workspace, inputs and import declaration unavailable')
    finally:
        for name in ['portable','portable-source']:shutil.move(withheld/name,root/name)
        withheld.rmdir()
    caps=run('bioinformatics');check(not caps['limits'],'no advertised sequence policy ceilings remain')
    save('requests.json',requests);save('acceptance.json',dict(status='PASS',checks=checks,requests=len(requests)))
    print(f'{checks} independent composite sequence checks passed')

if __name__=='__main__':main()
