"""Independent exhaustive edit, exon-chain and read-cover expectations."""
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
    requests, results = {}, {}

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(name, value):
        p = root / name
        p.write_text(json.dumps(value, indent=2) + '\n')
        return p

    def run(*args, fail=False):
        p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
        if fail:
            check(p.returncode != 0, 'invalid declaration must reject')
            if isinstance(fail, str):
                check(fail in p.stderr, 'unconverted publisher refuses before native publication')
            return
        if p.returncode:
            raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)

    resources = []
    def resource(name, records):
        p = root / name
        p.write_text(''.join('>' + k + '\n' + v + '\n' for k, v in records.items()))
        resources.append(dict(id=name, format='fasta', method='Authored exhaustive qualification fixture',
                              source=dict(path=name, sha256='sha256:' + hashlib.sha256(p.read_bytes()).hexdigest())))

    reverse = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    rng = random.Random(20261003)
    reference = dict(first='GACGTACCTTAGCGATACGGTANACGT', second=''.join(rng.choice('ACGT') for _ in range(30)))
    queries = dict(exact='ACGTA', inserted='ACGGTA', deleted='ACTA', reverse='TAGCG', no_hit='TTTTTTTTT')
    queries.update({'random' + str(i): ''.join(rng.choice('ACGT') for _ in range(3 + i)) for i in range(4)})
    resource('reference.fa', reference); resource('queries.fa', queries)
    splice_ref = dict(plus='TTAAAGTGGAGCCCGTTCAGTTTAA')
    splice_ref['minus'] = reverse(splice_ref['plus'])
    splice_ref['unresolved-boundaries'] = 'AAAUTAGCCCGTANTTT'
    splice_queries = dict(two_junctions='AAACCCTTT', single='AAACCC', contiguous='CCCGTT', absent='ACACAC')
    resource('splice-reference.fa', splice_ref); resource('splice-queries.fa', splice_queries)
    edited_refs = dict(plus='AAGTAGCC', minus=reverse('AAGTAGCC'), repeated='AAAGTAGCCC')
    edited_queries = dict(exact='AACC', inserted='AATCC')
    resource('edited-reference.fa', edited_refs); resource('edited-queries.fa', edited_queries)
    resource('edited-multi-reference.fa', dict(two='AAGTAGCCGTAGTT'))
    resource('edited-multi-query.fa', dict(two='AACATT'))
    resource('edited-wide-reference.fa', dict(tie='CA'))
    resource('edited-wide-query.fa', dict(tie='AC'))
    covers = {
        'unique': dict(a='AACG', b='CGTT', c='TTGC'),
        'ambiguous': dict(a='ATG', b='TGA', c='GAT'),
        'duplicates': dict(a='AACG', a2='AACG', b='CGTT'),
        'disconnected': dict(a='AAAA', b='CCCC'),
        'contained': dict(a='AACGTA', b='CGT'),
        'singleton': dict(a='ACGT'),
    }
    for i in range(5):
        s = ''.join(rng.choice('ACGT') for _ in range(14))
        covers['random' + str(i)] = {str(j): s[j * 2:j * 2 + 6] for j in range(4)}
    for name, reads in covers.items():
        resource('cover-' + name + '.fa', reads)
    resource('invalid-query.fa', dict(ambiguous='ACNT'))
    resource('too-many.fa', {str(i): 'ACGT' for i in range(9)})
    spec = dict(schema='gadget.assays.import.v1', id='sequence-reconstruction', assays=[], resources=resources, links=[],
                design=dict(question='Complete finite edit placements, exon chains and read covers',
                            population='Authored sequences', selection='All source records',
                            limitations=['Deterministic source and exhaustive mathematical qualification examples']))
    run('import-assays', '--spec', save('import.json', spec), '--out', root / 'workspace')
    workspace = root / 'workspace/workspace.json'

    def execute(name, req):
        req = dict(schema='gadget.assays.query.v1', **req)
        requests[name] = req
        answer = run('assay-query', '--workspace', workspace, '--request', save(name + '-request.json', req))
        save(name + '-result.json', answer); results[name] = answer
        return answer['result']

    def verify_route(route, carrier, expected_cost, source, target):
        field = dict((e['relation'], int(e['coefficient'])) for e in route['realized_field'])
        check(field == dict((e, route['ordered_relations'].count(e)) for e in set(route['ordered_relations'])), 'ordered path realizes complete field')
        cost = 0; boundary = [0] * carrier['vertices']; current = source
        for e in route['ordered_relations']:
            a, b, c = carrier['steps'][e]
            check(a == current, 'ordered route is contiguous')
            current = b; cost += int(c)
        for e, coefficient in field.items():
            a, b, _ = carrier['steps'][e]
            boundary[a] -= coefficient; boundary[b] += coefficient
        demand = [0] * carrier['vertices']; demand[source] -= 1; demand[target] += 1
        check(current == target and boundary == demand and cost == expected_cost, 'independent boundary and exact route cost')
        check(route['boundary_verified'] and route['cost_verified'], 'native route checks present')

    def edit(a, b):
        previous = list(range(len(b) + 1))
        for i, x in enumerate(a, 1):
            row = [i]
            for j, y in enumerate(b, 1):
                row.append(min(previous[j] + 1, row[-1] + 1, previous[j - 1] + (x != y)))
            previous = row
        return previous[-1]

    indexed = dict(operation='indexed_sequence_alignment', query_resource='queries.fa', reference_resource='reference.fa',
                   strands=['+', '-'], scoring='unit_edit_cost', maximum_edits=0, evidence='All unit-edit placements')
    for errors in range(3):
        req = dict(indexed, maximum_edits=errors)
        answer = execute('indexed-' + str(errors), req)
        expected = set()
        for qid, original in queries.items():
            for strand in req['strands']:
                q = original if strand == '+' else reverse(original)
                for rid, r in reference.items():
                    for start in range(len(r)):
                        for length in range(len(q) - errors, min(len(q) + errors, len(r) - start) + 1):
                            cost = edit(q, r[start:start + length])
                            if cost <= errors:
                                expected.add((qid, rid, strand, start, start + length, cost))
        observed = set()
        for m in answer['matches']:
            key = (m['query_id'], m['reference_id'], m['strand'], m['start'], m['end'], int(m['minimum_edits']))
            observed.add(key)
            q = queries[m['query_id']] if m['strand'] == '+' else reverse(queries[m['query_id']])
            check(m['aligned_query'].replace('-', '') == q and m['aligned_reference'].replace('-', '') == reference[m['reference_id']][m['start']:m['end']], 'alignment preserves both complete source sequences')
            check(sum(x != y for x, y in zip(m['aligned_query'], m['aligned_reference'])) == key[-1], 'reported edit history cost')
            c = answer['carriers'][m['carrier']]
            positions = c['prefix_coordinates']; target = positions.index([len(q), m['end'] - m['start']])
            verify_route(m['route'], c, key[-1], positions.index([0, 0]), target)
            check(all(abs(i - j) <= errors for i, j in positions), 'only provably admissible band vertices')
        check(observed == expected and len(observed) == len(answer['matches']), 'every exhaustive gapped placement retained exactly once')
        check(answer['complete'], 'indexed completeness')
        for summary in answer['queries']:
            check(summary['matches'] == sum(m[0] == summary['query_id'] for m in expected), 'unmatched and matched source identities')

    def minimum_histories(a, b):
        distance = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
        for i in range(len(a) + 1):
            distance[i][0] = i
        for j in range(len(b) + 1):
            distance[0][j] = j
        for i, x in enumerate(a, 1):
            for j, y in enumerate(b, 1):
                distance[i][j] = min(distance[i-1][j] + 1, distance[i][j-1] + 1,
                                     distance[i-1][j-1] + (x != y))
        histories = []
        pending = [(len(a), len(b), ())]
        while pending:
            i, j, reversed_steps = pending.pop()
            if i == j == 0:
                histories.append(tuple(reversed(reversed_steps)))
            if i and distance[i-1][j] + 1 == distance[i][j]:
                pending.append((i-1, j, reversed_steps + ((i-1, j, i, j),)))
            if j and distance[i][j-1] + 1 == distance[i][j]:
                pending.append((i, j-1, reversed_steps + ((i, j-1, i, j),)))
            if i and j and distance[i-1][j-1] + (a[i-1] != b[j-1]) == distance[i][j]:
                pending.append((i-1, j-1, reversed_steps + ((i-1, j-1, i, j),)))
        return histories

    family_cases = {}
    for errors in range(3):
        req = dict(schema='gadget.assays.query.v1', **dict(indexed, maximum_edits=errors,
                                                       edit_histories='all_minimum'))
        answer = run('assay-query', '--workspace', workspace, '--request',
                     save('minimum-families-' + str(errors) + '-request.json', req))
        save('minimum-families-' + str(errors) + '-result.json', answer)
        family_cases[str(errors)] = answer
        total = 0
        for m in answer['result']['matches']:
            q = queries[m['query_id']] if m['strand'] == '+' else reverse(queries[m['query_id']])
            r = reference[m['reference_id']][m['start']:m['end']]
            expected = minimum_histories(q, r)
            carrier = answer['result']['carriers'][m['carrier']]
            positions = carrier['prefix_coordinates']
            observed = []
            for history in m['minimum_edit_histories']:
                route = history['route']
                observed.append(tuple(tuple(positions[carrier['steps'][e][0]] +
                                            positions[carrier['steps'][e][1]])
                                      for e in route['ordered_relations']))
                verify_route(route, carrier, int(m['minimum_edits']), positions.index([0, 0]),
                             positions.index([len(q), len(r)]))
                check(history['aligned_query'].replace('-', '') == q and
                      history['aligned_reference'].replace('-', '') == r,
                      'minimum history preserves complete oriented source words')
            check(sorted(observed) == sorted(expected) and len(observed) == int(m['minimum_edit_history_count']),
                  'independent complete ordered minimum history family')
            total += len(observed)
        check(total == int(answer['result']['work']['minimum_edit_histories']), 'exact total minimum family count')
    save('minimum-history-families.json', family_cases)

    splice = dict(operation='sequence_splice_discovery', query_resource='splice-queries.fa', reference_resource='splice-reference.fa',
                  strands=['+', '-'], minimum_exon_bases=3, minimum_intron_bases=4, maximum_intron_bases=10,
                  maximum_junctions=2, motifs=['GT-AG'], evidence='Discover all exact exon chains from reference sequence')
    for maximum in range(3):
        req = dict(splice, maximum_junctions=maximum)
        answer = execute('splice-' + str(maximum), req)
        expected = set()
        for qid, original in splice_queries.items():
            for strand in req['strands']:
                q = original if strand == '+' else reverse(original)
                for rid, r in splice_ref.items():
                    def extend(i, previous, blocks):
                        if i == len(q):
                            expected.add((qid, rid, strand, tuple(blocks)))
                            return
                        if len(blocks) > maximum:
                            return
                        for end in range(i + 3, len(q) + 1):
                            for start in range(len(r) - (end - i) + 1):
                                if r[start:start + end - i] != q[i:end]: continue
                                if previous is not None:
                                    gap = start - previous
                                    if not 4 <= gap <= 10: continue
                                    motif = r[previous:previous + 2] + '-' + r[start - 2:start] if strand == '+' else reverse(r[start - 2:start]) + '-' + reverse(r[previous:previous + 2])
                                    if motif not in req['motifs']: continue
                                extend(end, start + end - i, blocks + [(i, end, start, start + end - i)])
                    extend(0, None, [])
        observed = set()
        for m in answer['matches']:
            blocks = tuple((b['query_begin'], b['query_end'], b['start'], b['end']) for b in m['blocks'])
            observed.add((m['query_id'], m['reference_id'], m['strand'], blocks))
            c = answer['carriers'][m['carrier']]
            verify_route(m['route'], c, len(m['junctions']), c['source_vertex'], c['sink_vertex'])
            check(m['query_reconstruction_verified'], 'full exon reconstruction')
            check(all(j['origin'] == 'reference_sequence_and_exact_exon_match' for j in m['junctions']), 'no annotation supplied for discovered junctions')
        check(expected == observed and len(observed) == len(answer['matches']), 'every independently enumerated exon chain')
        check(answer['complete'], 'splice completeness')
    check(any(len(m['junctions']) == 2 for m in results['splice-2']['result']['matches']), 'two-junction example exercised')
    check(any(m['strand'] == '-' and m['junctions'] for m in results['splice-2']['result']['matches']), 'negative transcript motif example exercised')
    check(all(m['reference_id'] != 'unresolved-boundaries' for m in results['splice-2']['result']['matches']), 'unresolved motif candidates preserve valid chains elsewhere')

    def edited_splice_expectations(q, r, strand, req):
        # Enumerate independent exon partitions and complete substring edit
        # tables. This does not construct the native biological prefix graph.
        candidates = {}
        h, budget = req['minimum_exon_bases'], req['maximum_edits']
        def extend(qb, previous, blocks, route, cost):
            if qb == len(q):
                key = (blocks[0][2], blocks[-1][3], len(blocks) - 1)
                if key not in candidates or cost < candidates[key][0]:
                    candidates[key] = (cost, [])
                if cost == candidates[key][0]:
                    candidates[key][1].append((tuple(blocks), route))
                return
            if len(blocks) > req['maximum_junctions']:
                return
            for qe in range(qb + h, len(q) + 1):
                for rb in range(len(r)):
                    jump = ()
                    if previous is not None:
                        gap = rb - previous
                        if not req['minimum_intron_bases'] <= gap <= req['maximum_intron_bases']:
                            continue
                        motif = r[previous:previous+2] + '-' + r[rb-2:rb] if strand == '+' else reverse(r[rb-2:rb]) + '-' + reverse(r[previous:previous+2])
                        if motif not in req['motifs']:
                            continue
                        jump = ((qb, previous, qb, rb, 'N'),)
                    for re in range(rb + h, min(len(r), rb + qe - qb + budget - cost) + 1):
                        local_cost = edit(q[qb:qe], r[rb:re])
                        if local_cost + cost > budget:
                            continue
                        for local in minimum_histories(q[qb:qe], r[rb:re]):
                            lifted = []
                            for a, b, c, d in local:
                                op = 'D' if a == c else 'I' if b == d else '=' if q[qb+a] == r[rb+b] else 'X'
                                lifted.append((qb+a, rb+b, qb+c, rb+d, op))
                            extend(qe, re, blocks + [(qb, qe, rb, re)], route + jump + tuple(lifted), cost + local_cost)
        extend(0, None, [], (), 0)
        return candidates

    edited_cases = {}
    cases = [(0, 1, 'edited-queries.fa', 'edited-reference.fa', edited_queries, edited_refs),
             (1, 1, 'edited-queries.fa', 'edited-reference.fa', edited_queries, edited_refs),
             (1, 2, 'edited-multi-query.fa', 'edited-multi-reference.fa', dict(two='AACATT'), dict(two='AAGTAGCCGTAGTT')),
             (4, 65536, 'edited-wide-query.fa', 'edited-wide-reference.fa', dict(tie='AC'), dict(tie='CA'))]
    exercised = set()
    for case, (budget, maximum, query_resource, reference_resource, reads, references) in enumerate(cases):
        req = dict(schema='gadget.assays.query.v1', **dict(splice, query_resource=query_resource,
            reference_resource=reference_resource, minimum_exon_bases=2, maximum_junctions=maximum,
            maximum_edits=budget, scoring='unit_edit_cost'))
        answer = run('assay-query', '--workspace', workspace, '--request', save('edited-splice-' + str(case) + '-request.json', req))
        save('edited-splice-' + str(case) + '-result.json', answer); result = answer['result']
        edited_cases[str(case)] = answer
        expected, observed = {}, {}
        for qid, original in reads.items():
            for strand in req['strands']:
                q = original if strand == '+' else reverse(original)
                for rid, r in references.items():
                    for key, (cost, family) in edited_splice_expectations(q, r, strand, req).items():
                        expected[(qid, rid, strand, *key)] = (cost, sorted(family))
        for m in result['matches']:
            q = reads[m['query_id']] if m['strand'] == '+' else reverse(reads[m['query_id']])
            r = references[m['reference_id']]
            blocks = tuple((b['query_begin'], b['query_end'], b['start'], b['end']) for b in m['blocks'])
            c = result['carriers'][m['carrier']]
            route, qp, rp, transitions = m['route'], 0, m['start'], []
            for e in route['ordered_relations']:
                op = c['operations'][e]
                if op == 'T': continue
                if op == 'N':
                    junction = m['junctions'][sum(t[-1] == 'N' for t in transitions)]
                    check(junction['donor'] == rp, 'intron begins at actual route prefix')
                    nq, nr = qp, junction['acceptor']
                else:
                    nq, nr = qp + (op != 'D'), rp + (op != 'I')
                transitions.append((qp, rp, nq, nr, op)); qp, rp = nq, nr
                exercised.add(op)
            check(qp == len(q) and rp == m['end'], 'edited splice preserves complete endpoint prefixes')
            cost = int(m['minimum_edits'])
            verify_route(route, c, cost, c['source_vertex'], result['minimum_history_groups'][m['minimum_history_group']]['sink_vertex'])
            check(sum(op in 'XID' for *_, op in transitions) == cost, 'complete unit edit cost excludes introns')
            for b in m['blocks']:
                check(b['aligned_query'].replace('-', '') == q[b['query_begin']:b['query_end']] and
                      b['aligned_reference'].replace('-', '') == r[b['start']:b['end']], 'each edited exon preserves both complete source words')
            key = (m['query_id'], m['reference_id'], m['strand'], m['start'], m['end'], len(m['junctions']))
            if key not in observed: observed[key] = (cost, [])
            check(observed[key][0] == cost, 'one native minimum per declared endpoint and junction count')
            observed[key][1].append((blocks, tuple(transitions)))
        observed = {key: (cost, sorted(family)) for key, (cost, family) in observed.items()}
        check(observed == expected, 'all independent edited exon partitions and ordered minimizing histories')
        check(len(result['matches']) == int(result['work']['minimum_edit_histories']), 'exact edited splice total history count')
        for i, group in enumerate(result['minimum_history_groups']):
            members = [m for m in result['matches'] if m['minimum_history_group'] == i]
            check(len(members) == int(group['histories']) and
                  [int(m['history_index']) for m in members] == list(range(len(members))), 'complete group size and retained occurrence order')
    check(set('=XIDN') <= exercised, 'match substitution insertion deletion and intron witnesses exercised')
    check(any(m['strand'] == '-' and m['junctions'] for m in edited_cases['1']['result']['matches']), 'edited negative-strand transcript motifs')
    check(any(len(m['junctions']) == 2 for m in edited_cases['2']['result']['matches']), 'edited two-junction witness')
    save('edited-splice-families.json', edited_cases)

    for name, reads in covers.items():
        req = dict(operation='read_overlap_reconstruction', resource='cover-' + name + '.fa', minimum_overlap=2,
                   objective='minimum_suffix_prefix_read_cover', evidence='Complete source read covers')
        answer = execute('cover-' + name, req)
        values = list(reads.values()); candidates = []
        for order in itertools.permutations(range(len(values))):
            options = []
            for a, b in zip(order, order[1:]):
                options.append([n for n in range(2, min(len(values[a]), len(values[b])) + 1) if values[a][-n:] == values[b][:n]])
            for overlaps in itertools.product(*options):
                s = values[order[0]] + ''.join(values[b][n:] for b, n in zip(order[1:], overlaps))
                candidates.append((len(s), s, order, overlaps))
        minimum = min((c[0] for c in candidates), default=None)
        expected = {(s, order, overlaps) for length, s, order, overlaps in candidates if length == minimum}
        observed = {(l['spelling'], tuple(l['read_order']), tuple(l['overlap_lengths'])) for l in answer['layouts']}
        check(expected == observed and len(observed) == len(answer['layouts']), 'all overlap lengths and permutations independently exhaustively minimized')
        check((None if answer['minimum_length'] is None else int(answer['minimum_length'])) == minimum, 'native global cover minimum')
        check(answer['identifiable_sequence'] == (len({s for s, _, _ in expected}) == 1), 'identifiability derives from every minimum sequence')
        for layout in answer['layouts']:
            c = answer['carrier']; verify_route(layout['route'], c, minimum, c['source_vertex'], c['sink_vertex'])
            g = layout['gluing']; check(g['complete_overlap_gluing_verified'] and g['every_read_restriction_verified'], 'native categorical gluing')
            for p in g['placements']:
                check(layout['spelling'][p['start']:p['end']] == values[p['source_read']], 'each original read reconstructs at its declared coordinates')
        representatives = [answer['layouts'][i] for i in answer['path_comparison']['representative_layouts']]
        for i, a in enumerate(representatives):
            af = {e['relation']: int(e['coefficient']) for e in a['route']['realized_field']}
            for j, b in enumerate(representatives):
                bf = {e['relation']: int(e['coefficient']) for e in b['route']['realized_field']}
                check(int(answer['path_comparison']['shared_field_gram'][i][j]) == sum(v * bf.get(e, 0) for e, v in af.items()), 'complete transition-carrier Gram retains shared route fields')
        check(answer['complete'], 'cover completeness')
    check(len(results['cover-ambiguous']['result']['sequences']) == 3, 'distinct cycle-derived minimum reconstructions retained')
    check(len(results['cover-duplicates']['result']['layouts']) == 2, 'duplicate source read identities survive')

    batch_request = dict(operation='query_batch', concurrent_jobs=2,
                         queries=[dict(id=name, query=requests[name])
                                  for name in ['indexed-1', 'splice-2', 'cover-ambiguous']])
    caps = run('bioinformatics')
    if all(job['query']['operation'] in caps['native_publication_methods'] for job in batch_request['queries']):
        batch = execute('parallel-source-queries', batch_request)
        for item in batch['results']:
            check(item['answer']['result'] == results[item['id']]['result'], 'joined native parallel source queries equal serial complete results')
            check(item['admitted_workers'] == 1, 'native assigned worker admission retained')
        check(batch['dispatcher'] == 'rfc::runtime::WorkGraph', 'canonical dispatcher owns and joins all jobs')
    else:
        run('assay-query', '--workspace', workspace, '--request',
            save('parallel-source-queries-request.json', dict(schema='gadget.assays.query.v1', **batch_request)),
            fail='source query child has no native publisher')

    for i, req in enumerate([
        dict(indexed, scoring='approximate'), dict(indexed, maximum_edits=3),
        dict(indexed, query_resource='invalid-query.fa'), dict(indexed, extra=True),
        dict(splice, minimum_exon_bases=1), dict(splice, maximum_junctions=-1),
        dict(splice, motifs=['G?-AG']), dict(splice, minimum_intron_bases=11),
        dict(splice, maximum_edits=1), dict(splice, scoring='unit_edit_cost'),
        dict(splice, maximum_edits=-1, scoring='unit_edit_cost'),
        dict(splice, maximum_edits=1, scoring='affine'),
        dict(operation='query_batch', concurrent_jobs=2, queries=[dict(id='nested', query=dict(schema='gadget.assays.query.v1', operation='query_batch'))]),
        dict(operation='query_batch', concurrent_jobs=0, queries=[dict(id='x', query=dict(schema='gadget.assays.query.v1', operation='inspect'))]),
        dict(operation='read_overlap_reconstruction', resource='too-many.fa', minimum_overlap=0,
             objective='minimum_suffix_prefix_read_cover', evidence='Reject invalid overlap declaration'),
    ]):
        run('assay-query', '--workspace', workspace, '--request', save('reject-' + str(i) + '.json', dict(schema='gadget.assays.query.v1', **req)), fail=True)
    caps = run('bioinformatics')
    check(all(op in caps['workspace_methods'] for op in ['indexed_sequence_alignment', 'sequence_splice_discovery', 'read_overlap_reconstruction']), 'all three generic services discoverable')
    from method_answers import qualify
    checks += qualify(binary, root/'typed-answers', workspace,
                      {name: (request, results[name]) for name, request in requests.items()})
    save('cases.json', dict(requests=requests, results=results))
    save('verification.json', dict(status='pass', checks=checks, queries=len(requests),
         independent_expectations=['Exhaustive unit edit distances at every substring placement',
                                   'Exhaustive exon partitions and reference placements with transcript motifs',
                                   'All source-read permutations and all compatible overlap lengths',
                                   'Independent signed boundary, cost, gluing and shared-field checks']))
    print(json.dumps(dict(status='pass', checks=checks, queries=len(requests))))


if __name__ == '__main__':
    main()
