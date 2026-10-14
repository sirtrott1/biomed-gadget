"""Exhaustive joint substring alignment expectations for short original reads."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess


def qualify(binary, root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    checks = 0
    reverse = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(path, value):
        path.write_text(json.dumps(value, indent=2) + '\n')
        return path

    def run(*args):
        p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True)
        if p.returncode:
            raise RuntimeError(p.stderr)
        return json.loads(p.stdout)

    def alignments(q, s, budget, paired):
        results = []
        for begin in range(len(s)):
            for end in range(begin + 1, len(s) + 1):
                pending = [(0, begin, 0, 0, (), frozenset())]
                while pending:
                    i, j, cost, pairs, steps, support = pending.pop()
                    if i == len(q) and j == end:
                        if pairs >= paired:
                            results.append((cost, begin, end, steps, support))
                        continue
                    if cost < budget:
                        if i < len(q):
                            pending.append((i+1, j, cost+1, pairs, steps+((i, j, i+1, j, 'I'),), support))
                        if j < end:
                            pending.append((i, j+1, cost+1, pairs, steps+((i, j, i, j+1, 'D'),), support))
                    if i < len(q) and j < end:
                        equal = q[i] == s[j]
                        if equal or cost < budget:
                            pending.append((i+1, j+1, cost+(not equal), pairs+1,
                                            steps+((i, j, i+1, j+1, '=' if equal else 'X'),),
                                            support | {j} if equal else support))
        return results

    def expectations(reads, req):
        budget, paired = req['maximum_edits'], req['minimum_paired_bases']
        if any(len(s) < paired for s in reads):
            return set()
        for length in range(1, sum(map(len, reads)) + 1):
            candidates = []
            for bases in itertools.product('ACGT', repeat=length):
                spelling = ''.join(bases)
                for strands in itertools.product(req['strands'], repeat=len(reads)):
                    options = [alignments(q if strand == '+' else reverse(q), spelling, budget, paired)
                               for q, strand in zip(reads, strands)]
                    for family in itertools.product(*options):
                        cost = sum(a[0] for a in family)
                        if cost > budget or set().union(*(a[4] for a in family)) != set(range(length)):
                            continue
                        signature = tuple((strand, a[1], a[2], a[3]) for strand, a in zip(strands, family))
                        candidates.append((cost, spelling, signature))
            if candidates:
                cost = min(a[0] for a in candidates)
                return {(spelling, signature) for e, spelling, signature in candidates if e == cost}
        return set()

    cases = [(['AC', 'AG'], 1, 2, ['+']), (['AA', 'AA'], 2, 1, ['+']),
             (['AA', 'A'], 0, 1, ['+']), (['A', 'C'], 0, 1, ['+']),
             (['AT', 'AT'], 0, 2, ['+', '-']), (['A', 'A'], 0, 2, ['+']),
             (['AT', 'ACT', 'ACT'], 1, 2, ['+'])]
    reports = []
    exercised = set()
    for index, (reads, budget, paired, strands) in enumerate(cases):
        case = root / str(index); case.mkdir(exist_ok=True)
        # Repeated labels, mixed line framing and quality are original evidence.
        raw = ''.join('@same' + ('\r\n' if r == 0 else '\n') + q + '\n+\n' + 'I'*len(q) + '\n'
                      for r, q in enumerate(reads)).encode()
        (case/'reads.fastq').write_bytes(raw)
        spec = dict(schema='gadget.assays.import.v1', id='consensus-' + str(index), assays=[], links=[],
                    resources=[dict(id='reads', format='fastq', method='Complete original consensus fixture',
                                    source=dict(path='reads.fastq', sha256='sha256:' + hashlib.sha256(raw).hexdigest()))],
                    design=dict(question='Every joint consensus under declared unit edits', population='Original short reads',
                                selection='All supplied occurrences', limitations=['Exhaustive finite verification fixture']))
        run('import-assays', '--spec', save(case/'import.json', spec), '--out', case/'workspace')
        request = dict(schema='gadget.assays.query.v1', operation='read_edit_consensus', resource='reads',
                       strands=strands, maximum_edits=budget, minimum_paired_bases=paired,
                       objective='minimum_length_then_total_unit_edits', ambiguity='literal_ACGT',
                       evidence='Complete joint substring alignment family')
        answer = run('assay-query', '--workspace', case/'workspace/workspace.rfcxs',
                     '--request', save(case/'request.json', request))['result']
        save(case/'result.json', answer)
        defaults = dict(schema='gadget.pipeline.default.v1', profile='sequence',
                        question=spec['design']['question'],
                        goals=[dict(id='consensus', method='read_edit_consensus',
                                    parameters={k: v for k, v in request.items()
                                                if k not in ('schema', 'operation')})])
        plan = run('pipeline-prepare', '--workspace', case/'workspace/workspace.rfcxs',
                   '--request', save(case/'default.json', defaults))
        check(len(plan['steps']) == 1 and plan['steps'][0]['request'] == request,
              'default consensus preserves every declared scientific parameter')
        index_result = run('pipeline', '--workspace', case/'workspace/workspace.rfcxs',
                           '--request', save(case/'plan.json', plan), '--out', case/'pipeline')
        check(index_result['status'] == 'complete' and len(index_result['results']) == 1,
              'default consensus executes through the existing pipeline dispatcher')
        receipt = index_result['results'][0]
        for member in ('request', 'result'):
            archive = case/'pipeline'/receipt[member]
            check(receipt[member + '_sha256'] == 'sha256:' + hashlib.sha256(archive.read_bytes()).hexdigest(),
                  'pipeline receipt identifies the complete native archive')
        retained = run('record-view', '--input', case/'pipeline'/receipt['result'],
                       '--member', 'result')
        check(retained['result'] == answer,
              'default pipeline retains every independently verified consensus family')
        expected = expectations(reads, request)
        observed = []
        for family in answer['families']:
            s = family['spelling']; signature = []; support = set(); cost = 0
            check(len(family['read_alignments']) == len(reads), 'original repeated source occurrences remain distinct')
            for r, alignment in enumerate(family['read_alignments']):
                q = reads[r] if alignment['strand'] == '+' else reverse(reads[r])
                i, j = 0, alignment['start']; steps = []; paired_count = 0
                for step in alignment['steps']:
                    op = step['operation']; exercised.add(op)
                    a, b, c, d = step['query_begin'], step['output_begin'], step['query_end'], step['output_end']
                    check((a, b) == (i, j) and (c-a, d-b) == ((1, 0) if op == 'I' else (0, 1) if op == 'D' else (1, 1)), 'complete monotone original-read alignment')
                    if op in '=X':
                        check((q[a] == s[b]) == (op == '='), 'declared match and substitution use actual source tokens')
                        paired_count += 1
                    if op == '=': support.add(b)
                    if op != 'D':
                        original = a if alignment['strand'] == '+' else len(q)-1-a
                        check(step['original_offset'] == original, 'signed strand map addresses original source token')
                    cost += op != '='; steps.append((a, b, c, d, op)); i, j = c, d
                check(i == len(q) and j == alignment['end'] and paired_count == alignment['paired_bases'] and paired_count >= paired, 'no query clipping and every original token consumed')
                signature.append((alignment['strand'], alignment['start'], alignment['end'], tuple(steps)))
            check(support == set(range(len(s))) and cost == family['edits'] <= budget, 'every emitted token has exact original support and explicit edit cost')
            route, carrier = family['route'], answer['carrier']; current = carrier['source_vertex']; total = 0
            for e in route['ordered_relations']:
                a, b, c = carrier['steps'][e]
                check(a == current, 'native joint optimization route continuity'); current = b; total += int(c)
            check(current == carrier['sink_vertex'] and total == (budget+1)*len(s)+cost == int(answer['minimum_cost']), 'native exact lexicographic objective witness')
            observed.append((s, tuple(signature)))
        check(set(observed) == expected and len(observed) == len(expected) == int(answer['family_count']), 'independent complete joint consensus family without administrative duplicates')
        check(answer['complete'] and answer['identifiable_sequence'] == (len({s for s, _ in expected}) == 1), 'declared ambiguity derives from every retained minimizing spelling')
        check(answer['status'] == 'no_admissible_consensus' if not expected else answer['minimum_length'] == len(next(iter(expected))[0]), 'complete no-family or length result')
        reports.append(dict(index=index, families=len(observed), request=request))
        if index == 0:
            for changes in [dict(maximum_edits=-1), dict(minimum_paired_bases=0),
                            dict(objective='majority'), dict(ambiguity='unknown'),
                            dict(strands=['+', '+']), dict(extra=True)]:
                bad = dict(request, **changes)
                p = subprocess.run([str(binary), 'assay-query', '--workspace', str(case/'workspace/workspace.rfcxs'),
                                    '--request', str(save(case/('reject-' + str(checks) + '.json'), bad))],
                                   text=True, capture_output=True)
                check(p.returncode != 0, 'invalid declared consensus profile refuses before publication')
    check({'=', 'X', 'I', 'D'} <= exercised, 'exact substitution insertion and deletion witnesses exercised')
    save(root/'verification.json', dict(status='PASS', checks=checks, cases=reports,
                                      oracle='Every short output word, orientation, substring and ordered unit-edit alignment; exact source support and declared joint minimum'))
    return checks


if __name__ == '__main__':
    import sys
    print(qualify(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()))
