"""Independent allocation factor and complete indexed sequence contracts."""
import argparse
import hashlib
import itertools
import json
import os
from fractions import Fraction as Q
from pathlib import Path
import random
import shutil
import subprocess
import time


def rank(rows):
    rows = [[Q(x) for x in row] for row in rows]
    if not rows:
        return 0
    k = 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        divisor = rows[k][j]
        rows[k] = [x / divisor for x in rows[k]]
        for i in range(len(rows)):
            if i != k and rows[i][j]:
                factor = rows[i][j]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[k])]
        k += 1
        if k == len(rows):
            break
    return k


def compositions(n, k):
    if k == 1:
        return [(n,)]
    return [(i, *rest) for i in range(n + 1) for rest in compositions(n - i, k - 1)]


def edit(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        row = [i]
        for j, y in enumerate(b, 1):
            row.append(min(previous[j] + 1, row[-1] + 1, previous[j - 1] + (x != y)))
        previous = row
    return previous[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=Path)
    parser.add_argument('root', type=Path)
    parser.add_argument('--baseline', type=Path, default=os.environ.get('GADGET_TEST_BASELINE_BINARY'))
    args = parser.parse_args()
    binary, root = args.binary.resolve(), args.root.resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    checks = 0
    timings = {}

    def check(v, message):
        nonlocal checks
        checks += 1
        if not v:
            raise AssertionError(message)

    def save(name, value):
        p = root / name
        p.write_text(json.dumps(value, indent=2) + '\n')
        return p

    def run(*argv, fail=False, executable=binary):
        start = time.perf_counter()
        p = subprocess.run([str(executable), *map(str, argv)], capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        if fail:
            check(p.returncode != 0, 'invalid source must reject')
            return None, elapsed
        if p.returncode:
            raise RuntimeError(f'{argv}: {p.stderr}')
        return json.loads(p.stdout), elapsed

    def allocation(name, candidates, classes, counts, targets, factors=True, **extra):
        req = dict(schema='gadget.assay.allocations.v1', operation='allocation_bounds', candidates=candidates,
                   classes=classes, samples=[dict(id='sample', counts=list(map(str, counts)), scale='1',
                   scale_provenance='Declared raw counts')], targets=targets, contrasts=[],
                   class_partition='disjoint_observations', candidate_search=dict(method='Authored candidate set',
                   completeness='All candidates in this fixture', reference_sha256='sha256:'+'0'*64),
                   evidence=['Authored qualification source'], units='assigned_fragment_count')
        req.update(extra)
        if factors:
            req['kernel_certificate'] = True
        path = save(name+'-request.json', req)
        answer, elapsed = run('allocation-bounds', '--request', path)
        save(name+'-result.json', answer)
        timings[name] = elapsed
        return req, answer['result']

    def verify_factors(req, result, full_rank=True, allocation_occurrence_count=None):
        cert = result['kernel_certificate']
        check(cert['target_coordinate_order'] == req['targets'], 'common target coordinate order')
        offset = column = 0
        images = []
        for c, f in zip(req['classes'], cert['factors']):
            check(f['class_id'] == c['id'] and f['candidates'] == c['candidates'], 'retained class/candidate order')
            check(f['native_relation_identity'] == 1000000 + len(images), 'source native relation identity')
            check(f['occurrence_offset'] == offset and f['kernel_offset'] == column, 'independent occurrence and kernel offsets')
            k = len(c['candidates']) - 1
            check(f['kernel_dimension'] == k, 'independent class rank-nullity')
            columns = [dict() for _ in range(k)]
            for i, j, q in f['kernel_frame']:
                check(i < len(c['candidates']) and j < k and Q(q) != 0, 'sparse kernel coordinate')
                check(i not in columns[j], 'unique kernel coordinate')
                columns[j][i] = Q(q)
            for j, col in enumerate(columns):
                check(sum(col.values()) == 0, 'class sum of kernel column is exactly zero')
                check({i:q for i,q in col.items() if i > 0} == {j+1:Q(1)}, 'free-coordinate identity proves complete independent kernel')
            observed = {(t,j):Q(q) for t,j,q in f['target_image']}
            expected = {}
            for t, target in enumerate(req['targets']):
                selected = set(target['candidates'])
                for j, col in enumerate(columns):
                    value = sum(q for i,q in col.items() if c['candidates'][i] in selected)
                    if value:
                        expected[t,j] = value
            check(observed == expected, 'native target image equals independently aggregated kernel')
            images.append((k, observed))
            offset += len(c['candidates']); column += k
        expected_occurrences=result['native_occurrence_count'] if allocation_occurrence_count is None else allocation_occurrence_count
        check(offset == expected_occurrences and column == result['allocation_nullity'], 'complete source factors')
        check(result['constraint_rank'] == len(req['classes']) and result['kernel_residual_zero'], 'full retained constraint rank/residual')
        if full_rank:
            rows = [[] for _ in req['targets']]
            for k, image in images:
                for t, row in enumerate(rows):
                    row.extend(image.get((t,j), Q(0)) for j in range(k))
            check(result['target_ambiguity_rank'] == rank(rows), 'common-target image rank, not sum of local ranks')

    example_path = Path(__file__).resolve().parents[1] / 'examples/bioinformatics/allocation-factors.json'
    example = json.loads(example_path.read_text())
    answer,_ = run('allocation-bounds','--request',example_path)
    save('example-result.json',answer)
    verify_factors(example,answer['result'])
    check([(Q(b['lower']),Q(b['upper'])) for b in answer['result']['contrasts'][0]['bounds']]==[(-12,-9),(-3,0),(-10,-2)],
          'public illustration exact assignment-robust contrasts')

    rng = random.Random(20261004)
    for case in range(24):
        names = ['a', 'b', 'c', 'd']
        classes = [dict(id=str(i), candidates=rng.sample(names, rng.randint(1,3))) for i in range(3)]
        counts = [rng.randint(0,2) for _ in classes]
        targets = [dict(id=s, candidates=[s]) for s in names]+[dict(id='ab', candidates=['a','b']),dict(id='all', candidates=names)]
        req, result = allocation('finite-'+str(case), names, classes, counts, targets)
        verify_factors(req, result)
        possibilities = []
        for joint in itertools.product(*(compositions(n,len(c['candidates'])) for n,c in zip(counts,classes))):
            values = {s:0 for s in names}
            for c, xs in zip(classes,joint):
                for s,x in zip(c['candidates'],xs):
                    values[s] += x
            possibilities.append([sum(values[s] for s in t['candidates']) for t in targets])
        for t,b in enumerate(result['samples'][0]['bounds']):
            lo,hi = min(x[t] for x in possibilities), max(x[t] for x in possibilities)
            check((Q(b['lower']),Q(b['upper']),b['identifiable']) == (lo,hi,lo==hi), 'independent exhaustive marginal bounds')
        if args.baseline and case < 3:
            plain = {k:v for k,v in req.items() if k != 'kernel_certificate'}
            path = save('baseline-request.json',plain)
            old,_ = run('allocation-bounds','--request',path,executable=args.baseline.resolve())
            fresh,_ = run('allocation-bounds','--request',path)
            check(old['result'] == fresh['result'], 'small allocation mathematical payload preserved')

    # Both classes vary along the same target direction: aggregate rank is one.
    req,result=allocation('shared-target', ['a','b'], [dict(id=str(i),candidates=['a','b']) for i in range(3)], [1,2,0],
                          [dict(id='a',candidates=['a']),dict(id='b',candidates=['b'])])
    verify_factors(req,result)
    check(result['allocation_nullity']==3 and result['target_ambiguity_rank']==1,'shared target dependence retained')

    req,result=allocation('zero-count', ['a','b'], [dict(id='retained',candidates=['a','b'])], [0],
                          [dict(id='a',candidates=['a']),dict(id='b',candidates=['b'])])
    verify_factors(req,result)
    check(result['native_occurrence_count']==2 and result['allocation_nullity']==1 and result['target_ambiguity_rank']==1,
          'zero measurement preserves the homogeneous source kernel')
    check(all(Q(b['lower'])==Q(b['upper'])==0 and b['identifiable'] for b in result['samples'][0]['bounds']),
          'nonnegative zero-count feasible bounds remain distinct from structural ambiguity')

    n=4001
    names=['c'+str(i) for i in range(n)]
    req,result=allocation('sparse-targets',names,[dict(id=str(i),candidates=[s]) for i,s in enumerate(names)],list(range(n)),
                          [dict(id=s,candidates=[s]) for s in names])
    check(n*n>10000000 and result['kernel_certificate']['bound_map_nonzeros']==2*n,'sparse maps above old logical-entry ceiling')
    check(result['constraint_rank']==n and result['allocation_nullity']==result['target_ambiguity_rank']==0,'singleton factors, all coordinates retained')
    for i,b in enumerate(result['samples'][0]['bounds']):
        check((b['id'],Q(b['lower']),Q(b['upper']))==(names[i],i,i),'all sparse target bounds retained')

    names=['t'+str(i) for i in range(5001)]
    targets=[dict(id='all',candidates=names),dict(id='even',candidates=names[::2])]
    wide_request,wide_result=allocation('wide-class',names,[dict(id='wide',candidates=names)],[17],targets)
    req,result=wide_request,wide_result
    verify_factors(req,result,full_rank=False)
    check(result['native_occurrence_count']==sum(len(c['candidates']) for c in wide_request['classes']),
          'standalone class source retains every original allocation occurrence')
    check(result['allocation_nullity']==5000 and result['target_ambiguity_rank']==1,'large class native kernel and target rank')
    check([(Q(b['lower']),Q(b['upper'])) for b in result['samples'][0]['bounds']]==[(17,17),(0,17)],'wide-class sharp bounds')
    samples=[dict(id='s'+str(i),counts=[str(i)],scale='1/2',scale_provenance='Declared half-unit scale') for i in range(131)]
    contrasts=[dict(id='contrast'+str(i),left='s'+str(i+1),right='s'+str(i)) for i in range(130)]
    req,result=allocation('samples-contrasts',['a'],[dict(id='a',candidates=['a'])],[0],[dict(id='a',candidates=['a'])],samples=samples,contrasts=contrasts)
    check(len(result['samples'])==131 and len(result['contrasts'])==130,'all samples and contrasts retained')
    for c in result['contrasts']:
        check(Q(c['bounds'][0]['lower'])==Q(c['bounds'][0]['upper'])==Q(1,2),'exact scaled contrast')

    resources=[]
    def resource(name,records):
        p=root/name;p.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in records.items()))
        resources.append(dict(id=name,format='fasta',method='Authored exact source',source=dict(path=name,sha256='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())))
    reference={'r':'AAAAAAAAAA','varied':'ACGTACATGTNNACGT'}
    queries={'q':'AAAAAAAAAA'}
    resource('reference.fa',reference);resource('queries.fa',queries)
    tiny_reads={f'q{i}':s for i,s in enumerate(map(''.join,itertools.product('AC',repeat=2)))}
    tiny_reference={f'r{i}':s for i,s in enumerate(map(''.join,itertools.product('AC',repeat=3)))}
    tiny_reference.update(ambiguity='ANCA',left='A',right='C')
    resource('partition-queries.fa',tiny_reads);resource('partition-reference.fa',tiny_reference)
    repeated_reads=['AC','AC','CA']
    repeated=root/'repeated.fastq'
    repeated.write_bytes(b'@same first lane\r\nac\r\n+same first lane\r\n!!\r\n'
                         b'@same second lane\nAC\n+\n55\n@same\nCA\n+\n~~\n')
    resources.append(dict(id=repeated.name,format='fastq',method='Distinct original read occurrences with repeated labels',
                          source=dict(path=repeated.name,sha256='sha256:'+hashlib.sha256(repeated.read_bytes()).hexdigest())))
    long_size=20065
    resource('long-reference.fa',{'long':'A'*long_size});resource('long-query.fa',{'long':'A'*long_size})
    # Many source reads, and enough placements to cross the old 10000 result cap.
    resource('many-query.fa',{'read'+str(i):'A' for i in range(1025)})
    salmon=root/'classes.eq'; salmon.write_text('5001\n1\n'+'\n'.join(names)+'\n5001 '+ ' '.join(map(str,range(5001)))+' 17\n')
    resources.append(dict(id='classes',format='salmon_eq',method='Unweighted declared class source',source=dict(path=salmon.name,sha256='sha256:'+hashlib.sha256(salmon.read_bytes()).hexdigest())))
    spec=dict(schema='gadget.assays.import.v1',id='sequence-assignment',assays=[],resources=resources,links=[],design=dict(question='Complete exact allocation factors and sequence placements',population='Authored sequences/classes',selection='All source records',limitations=['Qualification fixtures']))
    run('import-assays','--spec',save('import.json',spec),'--out',root/'workspace')
    workspace=root/'workspace/workspace.json'
    def query(name,**kw):
        req=dict(schema='gadget.assays.query.v1',**kw); answer,elapsed=run('assay-query','--workspace',workspace,'--request',save(name+'-request.json',req));save(name+'-result.json',answer);timings[name]=elapsed;return answer['result']
    base=dict(operation='indexed_sequence_alignment',query_resource='queries.fa',reference_resource='reference.fa',strands=['+','-'],scoring='unit_edit_cost',maximum_edits=9,evidence='Every admitted placement')
    result=query('nine-edits',**base)
    reverse=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    expected=set()
    for strand in base['strands']:
        for rid,r in reference.items():
            q=queries['q'] if strand=='+' else reverse(queries['q'])
            for start in range(len(r)):
                for end in range(start+1,min(len(r),start+len(q)+9)+1):
                    cost=edit(q,r[start:end])
                    if cost<=9:expected.add((rid,strand,start,end,cost))
    check({(m['reference_id'],m['strand'],m['start'],m['end'],int(m['minimum_edits'])) for m in result['matches']}==expected,'complete placements beyond former edit budget')
    for m in result['matches']:
        c=result['carriers'][m['carrier']]; current=c['prefix_coordinates'].index([0,0]);cost=0;field={}
        for e in m['route']['ordered_relations']:
            a,b,w=c['steps'][e];check(a==current,'contiguous native ordered route');current=b;cost+=Q(w);field[e]=field.get(e,0)+1
        check(c['prefix_coordinates'][current]==[10,m['end']-m['start']] and cost==Q(m['minimum_edits']),'independent endpoint and cost')
        check(field=={v['relation']:int(v['coefficient']) for v in m['route']['realized_field']},'complete realization field')
    result=query('long-exact',**dict(base,query_resource='long-query.fa',reference_resource='long-reference.fa',strands=['+'],maximum_edits=0))
    check(result['work']['candidate_starts']>20000 and result['work']['route_vertices']==long_size+1,'source-derived candidate and native carrier sizes')
    check(len(result['matches'])==1 and result['matches'][0]['aligned_query']=='A'*long_size,'complete long-read result')
    c=result['carriers'][0];m=result['matches'][0]
    check(c['vertices']==long_size+1 and len(c['steps'])==long_size and len(m['route']['ordered_relations'])==long_size,'complete long source carrier/history')
    for i,(coord,step) in enumerate(zip(c['prefix_coordinates'],c['steps'])):
        check(coord==[i,i] and step==[i,i+1,'0'],'all long-read source coordinates/steps retained')
    result=query('many-placements',**dict(base,query_resource='many-query.fa',reference_resource='reference.fa',strands=['+'],maximum_edits=0))
    check(len(result['queries'])==1025 and len(result['matches'])>10000,'all source reads and placements beyond former output cap')
    expected_count=sum(v.count('A') for v in reference.values())
    check(all(q['matches']==expected_count for q in result['queries']),'every source read has every exact placement')
    check(result['work']['read_partitions']==5 and result['work']['peak_partition_reads']==256,
          'all 1025 original reads processed through complete bounded batches')
    for budget in [0,1]:
        oracle=set()
        for qid,q in tiny_reads.items():
            for strand in ['+','-']:
                oriented=q if strand=='+' else reverse(q)
                for rid,r in tiny_reference.items():
                    for start in range(len(r)):
                        for end in range(start+1,len(r)+1):
                            cost=edit(oriented,r[start:end])
                            if cost<=budget:oracle.add((qid,rid,strand,start,end,cost))
        common=None
        for batch,width in [(1,1),(3,2),(20,5)]:
            control=dict(read_partition_records=batch,reference_block_bases=width)
            result=query(f'partition-{budget}-{batch}-{width}',**dict(base,
                query_resource='partition-queries.fa',reference_resource='partition-reference.fa',
                maximum_edits=budget,execution=control))
            actual={(m['query_id'],m['reference_id'],m['strand'],m['start'],m['end'],
                     int(m['minimum_edits'])) for m in result['matches']}
            check(actual==oracle,'exhaustive placements across read batches, seams and repeats')
            payload={k:v for k,v in result.items() if k!='work'}
            if common is None:common=payload
            else:check(payload==common,'complete routes, fields, carriers and source spans invariant under partition widths')
            check(result['work']['peak_partition_reads']<=batch and
                  result['work']['peak_normalized_reference_block_bases']<=width,
                  'execution widths bound normalized batches and blocks, without truncation')
    for i,execution in enumerate([
        dict(read_partition_records=0,reference_block_bases=2),
        dict(read_partition_records=1,reference_block_bases=0),
        dict(read_partition_records=-1,reference_block_bases=2),
        dict(read_partition_records=1.5,reference_block_bases=2),
        dict(read_partition_records=1),
        dict(read_partition_records=1,reference_block_bases=2,ignored=True)]):
        run('assay-query','--workspace',workspace,'--request',save(f'invalid-partition-{i}.json',
            dict(schema='gadget.assays.query.v1',**dict(base,execution=execution))),fail=True)
    for budget in [0,1]:
        expected=set()
        for occurrence,read in enumerate(repeated_reads):
            for strand in ['+','-']:
                oriented=read if strand=='+' else reverse(read)
                for rid,reference_bases in tiny_reference.items():
                    for start in range(len(reference_bases)):
                        for end in range(start+1,len(reference_bases)+1):
                            cost=edit(oriented,reference_bases[start:end])
                            if cost<=budget:expected.add((occurrence,rid,strand,start,end,cost))
        for batch,width in [(1,1),(3,2)]:
            result=query(f'repeated-edits-{budget}-{batch}',**dict(base,
                query_resource=repeated.name,reference_resource='partition-reference.fa',maximum_edits=budget,
                execution=dict(read_partition_records=batch,reference_block_bases=width)))
            actual={(m['query_record_index'],m['reference_id'],m['strand'],m['start'],m['end'],int(m['minimum_edits']))
                    for m in result['matches']}
            check(actual==expected and len(result['matches'])==len(expected),
                  'every repeated original read has every independently enumerated placement')
            check([q['query_record_index'] for q in result['queries']]==[0,1,2] and
                  all(q['query_id']=='same' for q in result['queries']),
                  'repeated read labels retain distinct ordered summary occurrences')
            for m in result['matches']:
                check(m['query_id']=='same' and m['query_record_index']==m['query_evidence']['record_index'] and
                      m['query_evidence']['source_sha256']=='sha256:'+hashlib.sha256(repeated.read_bytes()).hexdigest(),
                      'indexed placement source occurrence differs from retained original bytes')
            for q in result['queries']:
                check(q['matches']==sum(row[0]==q['query_record_index'] for row in expected),
                      'repeated read summary merged occurrence counts')
        mismatch=query(f'repeated-mismatches-{budget}',operation='reference_sequence_matches',
            query_resource=repeated.name,reference_resource='partition-reference.fa',strands=['+','-'],
            ambiguity='literal',maximum_mismatches=budget,junctions=[],minimum_overhang=1,
            evidence='Every fixed length original read occurrence')
        mismatch_expected=set()
        for occurrence,read in enumerate(repeated_reads):
            for strand in ['+','-']:
                oriented=read if strand=='+' else reverse(read)
                for rid,reference_bases in tiny_reference.items():
                    for start in range(len(reference_bases)-len(read)+1):
                        cost=sum(a!=b for a,b in zip(oriented,reference_bases[start:start+len(read)]))
                        if cost<=budget:mismatch_expected.add((occurrence,rid,strand,start,cost))
        actual={(m['query_record_index'],m['reference_id'],m['strand'],m['blocks'][0][0],int(m['mismatches']))
                for m in mismatch['matches']}
        check(actual==mismatch_expected and len(mismatch['matches'])==len(mismatch_expected),
              'complete mismatch placements retain repeated original occurrences')
        check([q['query_record_index'] for q in mismatch['queries']]==[0,1,2],
              'mismatch summaries retain original occurrence order')
    spliced=query('repeated-splice',operation='sequence_splice_discovery',query_resource=repeated.name,
        reference_resource='partition-reference.fa',strands=['+','-'],minimum_exon_bases=2,
        minimum_intron_bases=4,maximum_intron_bases=4,maximum_junctions=1,motifs=['GT-AG'],
        evidence='Every original short read; exact contiguous placements in this source')
    check([q['query_record_index'] for q in spliced['queries']]==[0,1,2] and
          all(q['query_id']=='same' for q in spliced['queries']),
          'splice summaries retain repeated original read occurrences')
    for m in spliced['matches']:
        check(m['query_id']=='same' and m['query_record_index']==m['query_evidence']['record_index'],
              'splice placement preserves original read occurrence')
    covered=query('repeated-cover',operation='read_overlap_reconstruction',resource=repeated.name,
        minimum_overlap=1,objective='minimum_suffix_prefix_read_cover',evidence='Every original repeated record')
    expected_layouts=set()
    best=None
    for order in itertools.permutations(range(3)):
        spelling=repeated_reads[order[0]]
        for a,b in zip(order,order[1:]):
            overlap=max(n for n in range(1,3) if repeated_reads[a][-n:]==repeated_reads[b][:n])
            spelling+=repeated_reads[b][overlap:]
        if best is None or len(spelling)<best:best=len(spelling);expected_layouts=set()
        if len(spelling)==best:expected_layouts.add((tuple(order),spelling))
    actual_layouts={(tuple(p['source_read'] for p in layout['gluing']['placements']),layout['spelling'])
                    for layout in covered['layouts']}
    check(actual_layouts==expected_layouts and len(covered['layouts'])==len(expected_layouts),
          'every minimum literal read cover retains all repeated original occurrences')
    query('inspect-no-execution',operation='inspect')
    run('assay-query','--workspace',workspace,'--request',save('foreign-execution.json',
        dict(schema='gadget.assays.query.v1',operation='inspect',
             execution=dict(read_partition_records=1,reference_block_bases=1))),fail=True)
    result=query('wide-salmon',operation='allocation_bounds',resource='classes',targets=targets,candidate_search=req['candidate_search'],scale='1',scale_provenance='raw',kernel_certificate=True)
    class_occurrences=sum(len(c['candidates']) for c in wide_request['classes'])
    aggregate_targets=[t for t in wide_request['targets'] if len(t['candidates'])!=1]
    target_occurrences=sum(1+len(t['candidates']) for t in aggregate_targets)
    check(result['native_occurrence_count']==class_occurrences+target_occurrences and
          result['native_relation_count']==len(wide_request['classes'])+len(aggregate_targets),
          'native Salmon diagram retains class and complete declared aggregate target occurrences')
    check(result['allocation_nullity']==5000 and result['target_ambiguity_rank']==1,
          'workspace Salmon adapter preserves the original exact allocation family')
    salmon_request={**wide_request,'classes':[dict(id='class:0',candidates=wide_request['classes'][0]['candidates'])]}
    verify_factors(salmon_request,result,full_rank=False,allocation_occurrence_count=class_occurrences)
    check(result['samples'][0]['bounds']==wide_result['samples'][0]['bounds'],
          'native Salmon bounds agree completely with the original unweighted class source')
    # The source hypothesis e<m remains necessary; arbitrary budget ceilings do not.
    run('assay-query','--workspace',workspace,'--request',save('invalid-edit.json',dict(schema='gadget.assays.query.v1',**dict(base,maximum_edits=10))),fail=True)
    invalid=dict(req);invalid['kernel_certificate']='true'
    run('allocation-bounds','--request',save('invalid-factor.json',invalid),fail=True)
    caps,_=run('bioinformatics')
    check(not any(k in caps['limits'] for k in ['sequence_resource_bytes','indexed_edit_budget','indexed_seed_hits','indexed_candidate_starts','indexed_total_route_vertices']),'capabilities remove replaced/stale limits')
    check(caps['indexed_sequence_execution']['routes'].endswith('dijkstra_exact endpoints'),'canonical ownership advertised')
    save('acceptance.json',dict(status='PASS',independent_checks=checks,runtime_seconds=timings,qualification='CPU exact fixtures; runtime excludes build'))
    print(json.dumps(dict(status='PASS',independent_checks=checks,runtime_seconds=timings)))


if __name__=='__main__':
    main()
