"""Independent full grid dual/history and exhaustive compatibility search contracts."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys


def main():
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    binary, root = map(lambda p: Path(p).resolve(), sys.argv[1:3])
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    checks = 0
    requests = []
    rng = random.Random(41602)

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(name, value):
        p = root/name
        p.write_text(json.dumps(value,indent=2)+'\n')
        return p

    def run(*args):
        p = subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
        if p.returncode:
            raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)

    def source(name, text):
        p=root/name
        p.write_text(text)
        return dict(path=name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())

    design=dict(question='Exact source sequence contracts',population='Authored sequences',selection='Complete declarations',limitations=['Verification fixture'])
    def workspace(name, queries, references):
        resources=[dict(id='queries',format='fasta',source=source(name+'-queries.fa',''.join('>'+k+'\n'+v+'\n' for k,v in queries.items())),method='Authored contract source'),
                   dict(id='reference',format='fasta',source=source(name+'-reference.fa',''.join('>'+k+'\n'+v+'\n' for k,v in references.items())),method='Authored contract source')]
        spec=dict(schema='gadget.assays.import.v1',id=name,design=design,assays=[],resources=resources,links=[])
        run('import-assays','--spec',save(name+'-import.json',spec),'--out',root/name)
        return root/name/'workspace.json'

    def query(name, ws, request):
        request=dict(request,schema='gadget.assays.query.v1')
        requests.append(request)
        result=run('assay-query','--workspace',ws,'--request',save(name+'-request.json',request))
        save(name+'-result.json',result)
        return result['result']

    masks=dict(A=1,C=2,G=4,T=8,U=8,R=5,Y=10,S=6,W=9,K=12,M=3,B=14,D=13,H=11,V=7,N=15)
    complement=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    def compatible(a,b,rule,pair=False):
        if rule=='literal':return a==b
        if rule=='possible':return bool(masks[a]&masks[b])
        return (masks[a]==masks[b] and masks[a] in [1,2,4,8]) if pair else not (masks[a]&~masks[b])

    def alignment_request(q,r,rule='literal',costs=('0','5/6','3/4','2/3')):
        return dict(operation='exact_sequence_alignment',query_resource='queries',query_id=q,reference_resource='reference',reference_id=r,
            ambiguity=rule,mode='global_linear_gap',scoring=dict(zip(['match','mismatch','query_gap','reference_gap'],costs),units='declared_cost',evidence='Authored exact edit costs'))

    def alignment_oracle(a,b,request,result):
        m,n=len(a),len(b)
        c={k:Q(v) for k,v in request['scoring'].items() if k in ['match','mismatch','query_gap','reference_gap']}
        rule=request['ambiguity']
        # Independent Fraction DP is a test oracle; production routing is native.
        prefix={(0,0):Q(0)}
        for i in range(m+1):
            for j in range(n+1):
                if i==j==0:continue
                options=[]
                if i:options.append(prefix[i-1,j]+c['reference_gap'])
                if j:options.append(prefix[i,j-1]+c['query_gap'])
                if i and j:options.append(prefix[i-1,j-1]+c['match' if compatible(a[i-1],b[j-1],rule,True) else 'mismatch'])
                prefix[i,j]=min(options)
        suffix={(m,n):Q(0)}
        for i in reversed(range(m+1)):
            for j in reversed(range(n+1)):
                if (i,j)==(m,n):continue
                options=[]
                if i<m:options.append(suffix[i+1,j]+c['reference_gap'])
                if j<n:options.append(suffix[i,j+1]+c['query_gap'])
                if i<m and j<n:options.append(suffix[i+1,j+1]+c['match' if compatible(a[i],b[j],rule,True) else 'mismatch'])
                suffix[i,j]=min(options)
        steps=[]
        width=n+1
        for i in range(m+1):
            for j in range(n+1):
                if i<m:steps.append(((i,j),(i+1,j),'I',c['reference_gap']))
                if j<n:steps.append(((i,j),(i,j+1),'D',c['query_gap']))
                if i<m and j<n:
                    eq=compatible(a[i],b[j],rule,True)
                    steps.append(((i,j),(i+1,j+1),'=' if eq else 'X',c['match' if eq else 'mismatch']))
        check(Q(result['minimum_cost'])==prefix[m,n],'independent global minimum')
        check(result['vertices']==(m+1)*(n+1) and result['relations']==len(steps),'complete original grid carrier')
        check(result['primal_dual_verified'] and result['suffix_dual_verified'],'both native dual certificates')
        expected={e for e,(u,v,op,cost) in enumerate(steps) if prefix[u]+cost+suffix[v]==prefix[m,n]}
        actual={int(p['relation_identity'])-1000000 for p in result['all_optimal_relations']}
        check(expected==actual,'complete optimal-step union')
        histories=[]
        for history in result['route_witnesses']:
            point=divmod(history['from_vertex'],width)
            points=[point];costs=[Q(0)]
            for e in history['ordered_relations']:
                u,v,op,cost=steps[e]
                check(u==point,'retained native history is continuous')
                point=v;points.append(point);costs.append(costs[-1]+cost)
            check(point==divmod(history['to_vertex'],width) and costs[-1]==Q(history['minimum_cost']),
                  'retained native history reconstructs endpoint and cost')
            histories.append((points,costs))
        check(result['native_route_calls']==len(histories) and len(histories)<=2*(m+1)*(n+1)+1,'every actual native call retained once')
        for witness in result['prefix_suffix_witnesses']:
            point=(witness['query_prefix'],witness['reference_prefix'])
            check(Q(witness['minimum_prefix_cost'])==prefix[point] and Q(witness['minimum_suffix_cost'])==suffix[point],
                  'every independent original-coordinate prefix/suffix cost')
            for key,expected_cost in [('prefix_witness',prefix[point]),('suffix_witness',suffix[point])]:
                origin=witness[key];points,costs=histories[origin['route']];p=origin['position']
                value=costs[p] if key=='prefix_witness' else costs[-1]-costs[p]
                check(points[p]==point and value==expected_cost,'every coordinate reconstructs from its retained native subroute')
        fields=result['relation_declaration']['fields']
        check(all(len(f['values'])==len(steps) for f in fields),'all original relation fields and supported zeros retained')
        check(len(result['relation_declaration']['relations'])==len(steps) and
              all(len(r['participants'])==2 for r in result['relation_declaration']['relations']),
              'every signed pair occurrence retained')
        chosen=[int(p['relation_identity'])-1000000 for p in result['path']]
        check(result['realized_path']['ordered_relations']==chosen,'native ordered realization retains chosen history')
        field={x['relation']:Q(x['coefficient']) for x in result['realized_path']['realized_field']}
        check(field=={e:Q(1) for e in chosen},'native field equals the chosen ordered source path')
        check(result['aligned_query'].replace('-','')==a and result['aligned_reference'].replace('-','')==b,
              'aligned spellings reconstruct original sources')
        return prefix,suffix

    a={f'q{i}':''.join(rng.choice('ACGTRYN') for _ in range(rng.randint(1,7))) for i in range(15)}
    b={f'r{i}':''.join(rng.choice('ACGTRYN') for _ in range(rng.randint(1,7))) for i in range(15)}
    a.update(large=''.join(rng.choice('ACGT') for _ in range(32)),long='A'*257,tied='ACGT',ambiguous='R',wide_coefficient='A')
    b.update(large=''.join(rng.choice('ACGT') for _ in range(37)),long='A',tied='TGCA',ambiguous='R',wide_coefficient='T')
    ws=workspace('alignment-source',a,b)
    for i in range(15):
        for rule in ['literal','possible','definite']:
            req=alignment_request('q'+str(i),'r'+str(i),rule)
            result=query(f'pair-{i}-{rule}',ws,req)
            alignment_oracle(a['q'+str(i)],b['r'+str(i)],req,result)
    for name,costs in [('large',('0','1','1','1')),('long',('0','1','1','1')),('tied',('0','0','0','0')),('ambiguous',('0','1','1','1')),('wide_coefficient',('0','1'+'0'*5000,'1','1'))]:
        req=alignment_request(name,name,'definite' if name=='ambiguous' else 'literal',costs)
        result=query('pair-'+name,ws,req)
        alignment_oracle(a[name],b[name],req,result)
        if name=='large':check(result['vertices']==1254,'full grid exceeds former 256 policy')
        if name=='long':check(result['vertices']==516,'literal source sequence exceeds former 255 policy')
        if name=='wide_coefficient':
            field=next(f for f in result['relation_declaration']['fields'] if f['id']=='step_cost')
            check(result['minimum_cost']=='2' and any(p['value']=='1'+'0'*5000 for p in field['values']),'5001-digit cost remains on its original relation while minimum chooses the complete two-gap route')
        if name=='tied':check(len(result['all_optimal_relations'])==result['relations'],'all zero-cost route alternatives survive')
        if name=='ambiguous':check(result['minimum_cost']=='1' and result['compatibility_rule']=='equal_singleton_IUPAC_sets','established pair definite rule explicitly declared')

    def matching_request(rule='literal',budget=1,junctions=None,overhang=1):
        return dict(operation='reference_sequence_matches',query_resource='queries',reference_resource='reference',strands=['+','-'],
            ambiguity=rule,maximum_mismatches=budget,minimum_overhang=overhang,junctions=junctions or [],evidence='Complete declared compatibility universe')

    def matching_oracle(queries,refs,req,result):
        expected=[];universe=0
        for qid,q in queries.items():
            for strand in req['strands']:
                seq=q if strand=='+' else q.translate(complement)[::-1]
                for rid,r in refs.items():
                    placements=[(list(range(start,start+len(q))),((start,start+len(q)),),None,None) for start in range(max(0,len(r)-len(q)+1))]
                    for j in req['junctions']:
                        if j['contig']!=rid or j['strand']!=strand:continue
                        for split in range(req['minimum_overhang'],len(q)-req['minimum_overhang']+1):
                            if split<=j['donor'] and len(q)-split<=len(r)-j['acceptor']:
                                start=j['donor']-split
                                placements.append((list(range(start,j['donor']))+list(range(j['acceptor'],j['acceptor']+len(q)-split)),
                                    ((start,j['donor']),(j['acceptor'],j['acceptor']+len(q)-split)),j['id'],split))
                    universe+=len(placements)
                    for positions,blocks,junction,split in placements:
                        count=sum(not compatible(seq[i],r[p],req['ambiguity']) for i,p in enumerate(positions))
                        if count<=req['maximum_mismatches']:expected.append((qid,rid,strand,blocks,junction,split,str(count)))
        actual=[(m['query_id'],m['reference_id'],m['strand'],tuple(map(tuple,m['blocks'])),m.get('junction',{}).get('id'),m.get('query_split_in_reference_orientation'),m['mismatches']) for m in result['matches']]
        check(actual==expected,'all independent matching placements in stable source order')
        check(int(result['candidate_universe'])==universe,'complete declared candidate universe')
        check([q['query_id'] for q in result['queries']]==list(queries),'every original query identity including no-match')
        for m in result['matches']:
            check(m['query_evidence']['source_start']<m['query_evidence']['source_end'] and
                  m['reference_evidence']['source_start']<m['reference_evidence']['source_end'],'original source spans retained')
        check(result['complete'],'complete matching result')

    dna={f'dna{i}':''.join(rng.choice('ACGT') for _ in range(rng.randint(1,11))) for i in range(20)}
    reads=dict(dna,ambiguous='RYSN',longer='A'*41,unmapped='GGGGGGGG',one='A')
    refs=dict(chrA=''.join(rng.choice('ACGTN') for _ in range(38)),chrB='ACGTTGCAACGTTGCAACGTTGCAACGTTGCA')
    junctions=[dict(id=f'j{i}',contig=rid,donor=d,acceptor=c,strand=s,evidence='Authored source junction')
               for i,(rid,d,c,s) in enumerate([('chrA',9,13,'+'),('chrA',11,19,'-'),('chrB',8,16,'+'),('chrB',8,16,'+')])]
    mws=workspace('matching-source',reads,refs)
    for rule in ['literal','possible','definite']:
        for budget in [0,1,2,12]:
            req=matching_request(rule,budget,junctions,2)
            matching_oracle(reads,refs,req,query('matching-'+rule+'-'+str(budget),mws,req))
    literal_reads=dict(first='NCGTACGT',second='ACGTNCGT',both='NCGTNCGT',
                       no_seed='NNG',all_unknown='NN',iupac='RY',short='N')
    literal_refs=dict(repeated='ACGTACGTACGT',other='TGCAACGT')
    literal_ws=workspace('literal-query-symbols',literal_reads,literal_refs)
    for budget in [0,1,2,8]:
        req=matching_request('literal',budget)
        matching_oracle(literal_reads,literal_refs,req,
            query('literal-query-symbols-'+str(budget),literal_ws,req))
    mixed_refs=dict(literal='NCGTNCGT',iupac='RYACGT',dna='ACGTACGT')
    mixed_ws=workspace('literal-reference-symbols',literal_reads,mixed_refs)
    for budget in [0,1,2]:
        req=matching_request('literal',budget)
        matching_oracle(literal_reads,mixed_refs,req,
            query('literal-reference-symbols-'+str(budget),mixed_ws,req))
    boundary_reads=dict(joined='AATT',contiguous='CCTT',short='AA',longer='AACCTTA')
    boundary_refs=dict(joined='AACCTT',contiguous='AATT')
    boundary_junctions=[dict(id='boundary-'+strand,contig='joined',donor=2,acceptor=4,
        strand=strand,evidence='Both exon ends meet the source boundary') for strand in ['+','-']]
    bws=workspace('matching-exon-boundaries',boundary_reads,boundary_refs)
    for overhang in [1,2,3]:
        req=matching_request('literal',0,boundary_junctions,overhang)
        matching_oracle(boundary_reads,boundary_refs,req,
            query('matching-exon-boundaries-'+str(overhang),bws,req))
    req=matching_request('definite',0);req['strands']=['+']
    aws=workspace('contained-mask',dict(q='R'),dict(r='R'))
    result=query('contained-mask',aws,req)
    check(len(result['matches'])==1 and result['compatibility_rule']=='query_IUPAC_set_contained_in_reference',
          'established directed subset rule remains separate from pair singleton equality')

    # Former source/policy bounds, with independently determined complete outputs.
    q=''.join(rng.choice('ACGT') for _ in range(1050));r='T'*25000+q+'T'*35000
    lws=workspace('long-read',dict(q=q),dict(r=r))
    req=matching_request('literal',0);req['strands']=['+']
    result=query('long-read',lws,req)
    check(len(result['matches'])==1 and result['matches'][0]['blocks']==[[25000,26050]],'long complete source read exact placement')
    check(int(result['candidate_universe'])*len(q)>50000000 and int(result['base_comparisons'])<50000000,
          'complete index proof covers universe beyond former comparison limit')
    check(result['index']['reference_base_visits']==len(r),'entire retained reference indexed without a sampled position set')
    for offset in [5,len(q)//2+5]:
        unknown=q[:offset]+'N'+q[offset+1:]
        nws=workspace('long-read-symbol-'+str(offset),dict(q=unknown),dict(r=r))
        req=matching_request('literal',1);req['strands']=['+']
        result=query('long-read-symbol-'+str(offset),nws,req)
        expected=[m.start() for m in re.finditer('(?='+re.escape(unknown).replace('N','[ACGT]')+')',r)]
        check([m['blocks'][0][0] for m in result['matches']]==expected and expected==[25000],
              'complete independent literal-symbol placement family')
        check(all(m['mismatches']=='1' for m in result['matches']),
              'literal unknown retains its native mismatch count')
        check(int(result['base_comparisons'])<int(result['candidate_universe'])*len(q),
              'literal unknown retains indexed candidate filtering')

    many={f'q{i}':'ACGTCAGT' for i in range(1025)}
    ws_many=workspace('many-queries',many,dict(r='TTTTTTTTTT'))
    req=matching_request('literal',0);req['strands']=['+']
    result=query('many-queries',ws_many,req)
    check(len(result['queries'])==1025 and not result['matches'],'query count beyond former 1024 policy with no-match identities')
    repetitive=workspace('many-matches',dict(q='A'),dict(r='A'*10007))
    result=query('many-matches',repetitive,req)
    check(len(result['matches'])==10007 and result['matches'][-1]['blocks']==[[10006,10007]],'all accepted matches beyond former 10000 policy')

    junc=[dict(id=f'j{i}',contig='r',donor=2,acceptor=4,strand='+',evidence='Distinct declared source assertion') for i in range(10001)]
    jws=workspace('many-junctions',dict(q='AATT'),dict(r='AACCTT'))
    req=matching_request('literal',0,junc,2);req['strands']=['+']
    result=query('many-junctions',jws,req)
    check(len(result['matches'])==10001 and [m['junction']['id'] for m in result['matches']]==[j['id'] for j in junc],
          'every repeated-location junction assertion remains distinct')

    ows=workspace('wide-overhang',dict(q='A'*1025+'C'*1025),dict(r='A'*1025+'T'*10+'C'*1025))
    req=matching_request('literal',0,[dict(id='junction',contig='r',donor=1025,acceptor=1035,strand='+',evidence='Authored exon join')],1025);req['strands']=['+']
    result=query('wide-overhang',ows,req)
    check(len(result['matches'])==1 and result['matches'][0]['blocks']==[[0,1025],[1035,2060]],'source-sized overhang beyond former 1024 policy')
    bws=workspace('wide-budget',dict(q='A'*1050),dict(r='T'*1060))
    req=matching_request('literal',1050);req['strands']=['+']
    result=query('wide-budget',bws,req)
    check(len(result['matches'])==11 and all(m['mismatches']=='1050' for m in result['matches']),
          'budget at least query length uses the complete fallback without clipping')

    big='T'*10000023
    huge_ws=workspace('large-reference',dict(q='ACGTCAGTACGTCAGTACGTCAGTACGTCAGTACGT'),dict(r=big))
    req=matching_request('literal',0);req['strands']=['+']
    result=query('large-reference',huge_ws,req)
    check(not result['matches'] and result['index']['reference_base_visits']==len(big),
          'complete reference beyond former ten-million-base policy')

    # Derived alignment fields replay without original source files/workspace.
    recipe=alignment_request('large','large',costs=('0','1','1','1'))
    derive=dict(schema='gadget.assays.materialize.v1',id='alignment-fields',source_namespace='416020',design=design,recipe=dict(recipe,schema='gadget.assays.query.v1'),channel='observation')
    run('derive-assays','--workspace',ws,'--request',save('derive.json',derive),'--out',root/'derived')
    withheld=root/'original-unavailable'
    withheld.mkdir()
    originals=['alignment-source','alignment-source-queries.fa','alignment-source-reference.fa','alignment-source-import.json']
    for name in originals:shutil.move(root/name,withheld/name)
    try:
        replay=run('replay-derivation','--workspace',root/'derived/workspace.rfcxs','--assay','alignment-fields','--out',root/'replayed')
        check(replay['status']=='PASS','complete larger alignment replays with original workspace, source files and import declaration unavailable')
        save('replay.json',replay)
    finally:
        for name in originals:shutil.move(withheld/name,root/name)
        withheld.rmdir()

    example=Path(__file__).resolve().parents[1]/'examples/bioinformatics/sequence-witnesses'
    run('import-assays','--spec',example/'import.json','--out',root/'portable')
    pws=root/'portable/workspace.json'
    defaults=json.loads((example/'pipeline.json').read_text())
    plan=run('pipeline-prepare','--workspace',pws,'--request',save('portable-plan-request.json',defaults))
    index=run('pipeline','--workspace',pws,'--request',save('portable-plan.json',plan),'--out',root/'portable-run')
    check(len(index['results'])==2 and index['dispatcher']=='rfc::runtime::WorkGraph','both sequence profiles execute through native default pipeline')
    save('requests.json',requests)
    save('acceptance.json',dict(status='PASS',checks=checks,requests=len(requests)))
    print(f'{checks} independent sequence witness checks passed')


if __name__=='__main__':main()
