"""Independent bounded-domain oracles and affine-factor reconstruction contracts."""
import copy
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
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    binary, root = map(lambda p: Path(p).resolve(), sys.argv[1:3])
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    checks = 0
    requests = []

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(name, data):
        path = root / name
        path.write_text(json.dumps(data, indent=2) + '\n')
        return path

    def source(name, data):
        path = root / name
        path.write_text(data)
        return dict(path=name, sha256='sha256:' + hashlib.sha256(path.read_bytes()).hexdigest())

    def run(*args, fail=False):
        p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
        if fail:
            check(p.returncode != 0, 'invalid declaration rejected')
            return
        if p.returncode:
            raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)

    assay = dict(id='fixture', kind='quantitative', format='dense_tsv',
        source=source('matrix.tsv', 'id\tsample\ngene\t1\n'),
        features=source('features.tsv', 'id\ngene\n'),
        observations=source('observations.tsv', 'id\nsample\n'),
        source_namespace='845100', carrier='features', quantity='authored_value',
        units='units', channel='observation', method='Authored verification',
        arithmetic='exact_rational', implicit_entries='unobserved', missing_tokens=['NA'])
    spec = dict(schema='gadget.assays.import.v1', id='coupled-allocation-contracts',
        design=dict(question='Exact coupled allocation contracts', population='Authored domains',
                    selection='Complete declared domain', limitations=['Verification fixture']),
        assays=[assay], resources=[], links=[])
    run('import-assays', '--spec', save('import.json', spec), '--out', root/'workspace')
    workspace = root/'workspace/workspace.json'

    def execute(name, request):
        request = dict(request, schema='gadget.assays.query.v1')
        requests.append(request)
        result = run('assay-query', '--workspace', workspace,
                     '--request', save(name+'-request.json', request))
        save(name+'-result.json', result)
        return result['result']

    def declaration(n, domain='bounded_integer', upper='2'):
        return dict(operation='constrained_allocations', domain=domain,
            variables=[dict(id=f'x{i}', lower='0', upper=upper) for i in range(n)],
            constraints=[], objective=dict(sense='minimize', coefficients=['0']*n),
            units='declared_units', evidence='Authored exact contract domain', factor_certificate=True)

    def linear(name, coeff, sense, rhs):
        return dict(id=name, coefficients=list(map(str, coeff)), sense=sense, rhs=str(rhs))

    def action(matrix, vector):
        result = [Q(0)]*matrix['rows']
        check(matrix['columns'] == len(vector), 'certificate matrix action dimensions')
        for i, j, q in matrix['entries']:
            result[i] += Q(q)*vector[j]
        return result

    def columns(matrix):
        return [[next((Q(q) for i,j,q in matrix['entries'] if i == row and j == col), Q(0))
                 for row in range(matrix['rows'])] for col in range(matrix['columns'])]

    def reconstruct(result):
        cert = result['factor_certificate']
        if result['status'] == 'infeasible':
            return set()
        tuples = [b['optimal_parameter_tuples'] for b in result['block_results']]
        witnesses = set()
        for product in itertools.product(*tuples):
            z = [Q(0)]*cert['kernel']['columns']
            for block, values in zip(cert['blocks'], product):
                for j, value in zip(block['parameter_indices'], values):
                    z[j] = Q(value)
            delta = action(cert['kernel'], z)
            witnesses.add(tuple(Q(x)+v for x,v in zip(cert['offset'], delta)))
        return witnesses

    def certificate(result, request):
        cert = result['factor_certificate']
        carrier = result['constraint_carrier']
        n, m = len(request['variables']), len(request['constraints'])
        check(carrier['variable_count'] == n and carrier['relation_count'] == m and
              carrier['occurrence_count'] == n*m, 'every declared term including zero remains an occurrence')
        check(cert['declaration']['variables'] == request['variables'] and
              cert['declaration']['constraints'] == request['constraints'] and
              cert['declaration']['objective'] == request['objective'], 'original declaration and ordering retained')
        if not cert['consistent']:
            if result['contradiction']['kind'] == 'inconsistent_equalities':
                lam = list(map(Q, result['contradiction']['left_null_vector']))
                e = cert['equality_operator']
                left = [Q(0)]*e['columns']
                for i,j,q in e['entries']:
                    left[j] += lam[i]*Q(q)
                value = sum(a*Q(b) for a,b in zip(lam, cert['equality_rhs']))
                check(not any(left) and value != 0 and value == Q(result['contradiction']['rhs_action']),
                      'independent left-null contradiction')
            return
        k, e = cert['kernel'], cert['equality_operator']
        x0 = list(map(Q, cert['offset']))
        check(action(e, x0) == list(map(Q, cert['equality_rhs'])), 'exact particular equality solution')
        for col in columns(k):
            check(not any(action(e, col)), 'independent equality-kernel residual')
        free = cert['free_variable_indices']
        check(len(free) == k['columns'] and n-k['columns'] == result['equality_rank'], 'full affine dimension')
        km = {(i,j):Q(q) for i,j,q in k['entries']}
        for a, i in enumerate(free):
            check(x0[i] == 0 and all(km.get((i,j),Q(0)) == int(a==j) for j in range(len(free))),
                  'parameters retain original source coordinates')
        groups = cert['group_indicators']
        grouped = {i:j for i,j,q in groups['entries'] if Q(q) == 1}
        check(len(grouped) == k['columns'] and len(groups['entries']) == len(grouped),
              'complete disjoint component indicators')
        for col in columns(groups):
            check(not any(action(cert['grouping_contrasts'], col)), 'independent grouping-kernel residual')
        for matrix in [k, cert['reduced_inequalities']]:
            supports = {}
            for i,j,q in matrix['entries']:
                supports.setdefault(i,set()).add(grouped[j])
            check(all(len(s) == 1 for s in supports.values()), 'each source variable and inequality belongs to one factor')

    def feasible(request, x):
        if any(not Q(v['lower']) <= a <= Q(v['upper']) for v,a in zip(request['variables'],x)):
            return False
        for row in request['constraints']:
            value = sum(Q(c)*a for c,a in zip(row['coefficients'],x))
            rhs = Q(row['rhs'])
            if not {'eq':value == rhs, 'le':value <= rhs, 'ge':value >= rhs}[row['sense']]:
                return False
        return True

    def verify(request, result, points):
        check(int(result['feasible_points_or_vertices']) == len(points), 'complete independent feasible count')
        check(result['complete'] and (result['status'] == 'feasible') == bool(points), 'complete feasibility status')
        certificate(result, request)
        if not points:
            check(not result['bounds'] and not result['optimizers'] and
                  int(result['optimal_points_or_vertices']) == 0, 'empty feasible set has no optimum')
            return
        values = {x:sum(Q(c)*a for c,a in zip(request['objective']['coefficients'],x)) for x in points}
        opt = (min if request['objective']['sense'] == 'minimize' else max)(values.values())
        optimal = {x for x,v in values.items() if v == opt}
        check(Q(result['minimum_or_maximum']) == opt, 'independent exact objective')
        check(int(result['optimal_points_or_vertices']) == len(optimal), 'complete optimal count')
        actual = {tuple(map(Q,w['values'])) for w in result['optimizers']}
        if result['witness_representation'] == 'explicit':
            check(actual == optimal and len(actual) == len(result['optimizers']), 'every explicit optimizer exactly once')
            for w in result['optimizers']:
                check(list(map(Q,w['constraint_values'])) == [sum(Q(c)*a for c,a in zip(row['coefficients'],map(Q,w['values'])))
                      for row in request['constraints']], 'original row action witness')
        else:
            check(not actual, 'factor output contract remains explicitly labeled')
        check(reconstruct(result) == optimal, 'complete independent Cartesian reconstruction')
        for i,b in enumerate(result['bounds']):
            check(tuple(map(Q,[b['lower'],b['upper'],b['optimal_lower'],b['optimal_upper']])) ==
                  (min(x[i] for x in points),max(x[i] for x in points),
                   min(x[i] for x in optimal),max(x[i] for x in optimal)), 'sharp original-coordinate feasible/optimal bounds')

    rng = random.Random(8451)
    for trial in range(55):
        n = rng.randint(1,4)
        req = declaration(n)
        req['variables'] = [dict(id=f'x{i}',lower=str(rng.randint(-2,0)),upper=str(rng.randint(0,3))) for i in range(n)]
        req['constraints'] = [linear(f'c{i}',[Q(rng.randint(-3,3),rng.randint(1,3)) for _ in range(n)],
                                     rng.choice(['eq','le','ge']),Q(rng.randint(-4,6),rng.randint(1,3)))
                              for i in range(rng.randint(0,4))]
        req['objective'] = dict(sense=rng.choice(['minimize','maximize']),
                               coefficients=[str(Q(rng.randint(-2,2),rng.randint(1,3))) for _ in range(n)])
        if trial % 2:
            req['witness_representation'] = 'factorized'
        points = {tuple(map(Q,x)) for x in itertools.product(*[range(int(v['lower']),int(v['upper'])+1) for v in req['variables']]) if feasible(req,x)}
        verify(req, execute(f'integer-oracle-{trial}',req), points)

    def solve(rows, rhs):
        # Independent Fraction elimination is a test oracle, not production code.
        n = len(rhs)
        a = [list(row)+[value] for row,value in zip(rows,rhs)]
        for j in range(n):
            pivot = next((i for i in range(j,n) if a[i][j]),None)
            if pivot is None:
                return None
            a[j],a[pivot] = a[pivot],a[j]
            v = a[j][j]
            a[j] = [q/v for q in a[j]]
            for i in range(n):
                if i != j:
                    v = a[i][j]
                    a[i] = [q-v*r for q,r in zip(a[i],a[j])]
        return tuple(row[-1] for row in a)

    for trial in range(35):
        n = rng.randint(1,3)
        req = declaration(n,'bounded_rational')
        req['variables'] = [dict(id=f'x{i}',lower=str(Q(rng.randint(-3,0),2)),upper=str(Q(rng.randint(1,6),2))) for i in range(n)]
        req['constraints'] = [linear(f'c{i}',[Q(rng.randint(-2,2),rng.randint(1,2)) for _ in range(n)],
                                     rng.choice(['eq','le','ge']),Q(rng.randint(-3,5),2))
                              for i in range(rng.randint(0,3))]
        req['objective'] = dict(sense=rng.choice(['minimize','maximize']), coefficients=[str(rng.randint(-2,2)) for _ in range(n)])
        if trial % 2:
            req['witness_representation'] = 'factorized'
        planes = [(list(map(Q,c['coefficients'])),Q(c['rhs'])) for c in req['constraints']]
        for i,v in enumerate(req['variables']):
            row = [Q(int(i==j)) for j in range(n)]
            planes.extend([(row,Q(v['lower'])),(row,Q(v['upper']))])
        vertices = set()
        for active in itertools.combinations(planes,n):
            x = solve([a for a,b in active],[b for a,b in active])
            if x is not None and feasible(req,x):
                vertices.add(x)
        verify(req,execute(f'rational-oracle-{trial}',req),vertices)

    # Equalities eliminate an otherwise enormous original box and preserve rows.
    req = declaration(128, upper='1000000000000')
    req['constraints'] = [linear(f'fixed-{i}',[int(i==j) for j in range(128)],'eq',i) for i in range(128)]
    result = execute('full-rank-large-box',req)
    verify(req,result,{tuple(map(Q,range(128)))})
    check(result['tested_candidates'] == 0 and result['affine_dimension'] == 0 and
          int(result['declared_box_points']) == 1000000000001**128, 'fixed equality coordinates require no box search')

    req = declaration(40, upper='1000000000000')
    req['variables'][-1]['upper'] = '7'
    req['constraints'] = [linear(f'linked-{i}',[int(j==i)-int(j==39) for j in range(40)],'eq',0) for i in range(39)]
    result = execute('one-free-coordinate',req)
    verify(req,result,{tuple([Q(i)]*40) for i in range(8)})
    check(result['tested_candidates'] == 8 and result['affine_dimension'] == 1,
          'source box replaced by exact one-dimensional affine fiber')

    req = declaration(130, upper='1')
    req['witness_representation'] = 'factorized'
    result = execute('binary-product-130',req)
    certificate(result,req)
    check(result['independent_blocks'] == 130 and result['tested_candidates'] == 260 and
          int(result['feasible_points_or_vertices']) == 2**130 and
          int(result['optimal_points_or_vertices']) == 2**130, 'complete 130-way product with exact large counts')
    check(all(b['optimal_parameter_tuples'] == [['0'],['1']] for b in result['block_results']),
          'every independent binary ambiguity retained')

    req = declaration(2000,upper='1')
    req['witness_representation'] = 'factorized'
    result = execute('sparse-2000-factor-product',req)
    check(result['independent_blocks'] == 2000 and result['tested_candidates'] == 4000 and
          int(result['feasible_points_or_vertices']) == 2**2000 and
          int(result['optimal_points_or_vertices']) == 2**2000,
          'complete 2000-factor domain prepared from supported entries')
    cert = result['factor_certificate']
    check(cert['kernel']['entries'] == [[i,i,'1'] for i in range(2000)] and
          cert['group_indicators']['entries'] == [[i,i,'1'] for i in range(2000)] and
          all(b['optimal_parameter_tuples'] == [['0'],['1']] for b in result['block_results']),
          'source identity frames and every local tuple survive large factor preparation')

    req = declaration(13,upper='1')
    result = execute('complete-8192-witnesses',req)
    check(len(result['optimizers']) == 8192 and int(result['optimal_points_or_vertices']) == 8192,
          'explicit output exceeds former 4096 witness policy without truncation')
    check({tuple(map(int,w['values'])) for w in result['optimizers']} == set(itertools.product(range(2),repeat=13)),
          'all 8192 witnesses in original coordinate order')

    req = declaration(130,upper='1')
    req['objective']['coefficients'] = ['1']*130
    result = execute('unique-optimum-large-product',req)
    check(result['optimizers'][0]['values'] == ['0']*130 and result['tested_candidates'] == 260 and
          int(result['feasible_points_or_vertices']) == 2**130, 'explicit unique optimum in a huge independent domain')

    req = declaration(14,'bounded_rational',upper='10')
    req['constraints'] = [linear(f'eq-{i}',[int(j==i)-int(j==13) for j in range(14)],'eq',0) for i in range(13)]
    result = execute('rational-reduced-dimension',req)
    verify(req,result,{tuple([Q(i)]*14) for i in [0,10]})
    check(result['tested_candidates'] == 2 and result['affine_dimension'] == 1,
          'rational bases use affine dimension and deduplicated bounds')

    req = declaration(3,'bounded_rational',upper='10')
    req['constraints'] = [linear(f'scaled-{i}',[i+1,i+1,i+1],'le',5*(i+1)) for i in range(150)]
    result = execute('scaled-active-hyperplanes',req)
    verify(req,result,{(Q(0),Q(0),Q(0)),(Q(5),Q(0),Q(0)),(Q(0),Q(5),Q(0)),(Q(0),Q(0),Q(5))})
    check(result['tested_candidates'] == 35 and result['constraint_carrier']['occurrence_count'] == 450,
          'source rows preserved while scaled active hyperplanes share search work')

    req = declaration(2,upper='5')
    req['constraints'] = [linear('congruence',[2,3],'eq',9)]
    result = execute('integer-congruence',req)
    verify(req,result,{(Q(0),Q(3)),(Q(3),Q(1))})
    check(result['tested_candidates'] == 6, 'fractional affine reconstruction never rounded')

    req = declaration(2)
    req['constraints'] = [linear('first',[1,1],'eq',1),linear('second',[2,2],'eq',3)]
    result = execute('equality-contradiction',req)
    verify(req,result,set())
    check(result['contradiction']['kind'] == 'inconsistent_equalities', 'native inconsistent equality witness')
    req['constraints'] = [linear('fractional',[2,0],'eq',1)]
    result = execute('fixed-integrality-contradiction',req)
    verify(req,result,set())
    check(result['contradiction']['kind'] == 'nonintegral_fixed_coordinate', 'constant rational coordinate excludes integer domain')

    # A constraint canceled by equality substitution must not couple factors.
    req = declaration(4)
    req['constraints'] = [linear('a',[1,-1,0,0],'eq',0),linear('b',[0,0,1,-1],'eq',0),
                          linear('canceled',[1,-1,1,-1],'le',0)]
    result = execute('canceled-coupling',req)
    verify(req,result,{(Q(a),Q(a),Q(b),Q(b)) for a,b in itertools.product(range(3),repeat=2)})
    check(result['independent_blocks'] == 2 and result['tested_candidates'] == 6,
          'actual reduced support establishes independence')
    req['constraints'].append(linear('actual-coupling',[1,0,1,0],'le',2))
    result = execute('actual-coupling',req)
    verify(req,result,{(Q(a),Q(a),Q(b),Q(b)) for a,b in itertools.product(range(3),repeat=2) if a+b <= 2})
    check(result['independent_blocks'] == 1 and result['tested_candidates'] == 9,
          'shared reduced normal keeps coupled parameters together')

    req = declaration(1)
    huge = '1'+'0'*5000
    req['constraints'] = [linear('large-exact-coefficient',[huge],'eq',huge)]
    result = execute('large-exact-coefficient',req)
    verify(req,result,{(Q(1),)})
    check(result['optimizers'][0]['values'] == ['1'], 'rational input follows source representation beyond former digit ceiling')

    example = json.loads((Path(__file__).resolve().parents[1]/'examples/bioinformatics/coupled-allocation-factors.json').read_text())
    result = execute('portable-example',example)
    points = {tuple(map(Q,[a,12-a,c,9-c,e,5-e])) for a,c,e in itertools.product(range(13),range(10),range(6)) if a+c<=13}
    verify(example,result,points)
    check(result['independent_blocks'] == 2 and result['tested_candidates'] == 136,
          'portable coupled capacity and independent ambiguity fiber')

    parameters = {k:v for k,v in example.items() if k not in ['schema','operation']}
    defaults = dict(schema='gadget.pipeline.default.v1',profile='custom',
                    question='Resolve declared coupled allocation ambiguity',
                    goals=[dict(method='constrained_allocations',parameters=parameters)])
    plan = run('pipeline-prepare','--workspace',workspace,'--request',save('pipeline-default.json',defaults))
    index = run('pipeline','--workspace',workspace,'--request',save('pipeline-plan.json',plan),
                '--out',root/'pipeline-run')
    check(len(index['results']) == 1 and index['dispatcher'] == 'rfc::runtime::WorkGraph',
          'factorized query executes through native pipeline orchestration')
    design = spec['design']
    derived = dict(schema='gadget.assays.materialize.v1',id='coupled-bounds',
                   source_namespace='845200',design=design,recipe=example,channel='observation')
    run('derive-assays','--workspace',workspace,'--request',save('derive.json',derived),
        '--out',root/'derived')
    shutil.move(root/'workspace',root/'original-unavailable')
    try:
        replay = run('replay-derivation','--workspace',root/'derived/workspace.rfcxs',
                     '--assay','coupled-bounds','--out',root/'replayed')
        check(replay['status'] == 'PASS', 'factorized producer replays from retained lineage without original workspace')
        save('replay.json',replay)
    finally:
        shutil.move(root/'original-unavailable',root/'workspace')

    for field,value in [('witness_representation','truncated'),('factor_certificate','yes')]:
        bad = dict(example,**{field:value})
        run('assay-query','--workspace',workspace,'--request',save('bad-'+field+'.json',bad),fail=True)
    save('requests.json',requests)
    save('acceptance.json',dict(checks=checks,requests=len(requests),status='PASS'))
    print(f'{checks} independent coupled allocation checks passed')


if __name__ == '__main__':
    main()
