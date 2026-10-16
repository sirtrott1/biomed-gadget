"""Independent rational and source support checks for declared multivariable comparisons."""
import copy
import csv
import importlib.util
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

from research_cli import rref, solve


def main(binary, root):
    shutil.rmtree(root, ignore_errors=True)
    module_path = Path(__file__).resolve().parents[1]/'examples/bioinformatics/linear-models/replay.py'
    spec = importlib.util.spec_from_file_location('linear_protocol', module_path)
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    queries = demo.stage(root/'inputs')
    checks, rejected = 0, 0

    def check(value):
        nonlocal checks
        assert value
        checks += 1

    def command(*args, fail=False):
        nonlocal rejected
        p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True, timeout=180)
        check(bool(p.returncode) == fail)
        if fail:
            rejected += 1
            return None
        return json.loads(p.stdout)

    command('import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    workspace = root/'workspace/workspace.json'
    def execute(query, fail=False):
        request = dict(schema='gadget.assays.query.v1', operation='native_query', assay='measurements', query=query)
        demo.save(root/'request.json', request)
        answer = command('assay-query', '--workspace', workspace, '--request', root/'request.json', fail=fail)
        return answer['result']['answers'][0]['certificate'] if answer else None

    source = list(csv.DictReader((root/'inputs/observations.tsv').open(), delimiter='\t'))
    response = next(csv.reader((root/'inputs/response.tsv').read_text().splitlines()[1:], delimiter='\t'))[1:]
    declaration = json.loads((root/'inputs/import.json').read_text())['assays'][0]
    fields = {'all.response': list(map(Q, response))}
    for covariate in declaration['observation_fields']:
        values = []
        for record in source:
            token = record[covariate['column']]
            values.append(None if token in covariate['missing_tokens'] else
                          Q(token == covariate['level']) if covariate['kind'] == 'level_indicator' else Q(token))
        fields['all.' + covariate['id']] = values

    def dot(a, b, weights):
        return sum((x*y*w for x,y,w in zip(a,b,weights)), Q(0))

    def represented(design, y, weights):
        _, pivots = rref(design)
        basis = [[row[j] for row in design] for j in pivots]
        if not basis:
            return [Q(0)] * len(y)
        coefficients = solve([[dot(a,b,weights) for b in basis] for a in basis],
                             [dot(a,y,weights) for a in basis])
        return [sum((basis[j][i]*coefficients[j] for j in range(len(basis))), Q(0)) for i in range(len(y))]

    def product(values):
        out = Q(1)
        for value in values:
            out *= value
        return out

    def verify(receipt, query):
        check(receipt['query'] == query)
        result = receipt['result']
        complete = set(query['responses'] + query['controls'])
        for term in query['terms']:
            complete.update(term['factors'])
        if isinstance(query['metric'], dict):
            complete.add(query['metric']['field'])
        selections = {'all': query['selection']}
        for s in query.get('strata', []):
            selections[s['id']] = list(dict.fromkeys(query['selection'] + s['selection']))
        for stratum in result['strata']:
            expected = [i for i in range(len(source)) if all(fields[f][i] is not None for f in complete)
                        and all(fields[f][i] == 1 for f in selections[stratum['id']])]
            check(stratum['coordinates'] == expected)
            check(stratum['entities'] == [source[i]['id'] for i in expected])
            check(stratum['n'] == len(expected))
            for excluded in stratum['excluded']:
                i = excluded['coordinate']
                check(i not in expected and excluded['entity'] == source[i]['id'])
                check(set(excluded['unobserved_fields']) == {f for f in complete if fields[f][i] is None})
                check(set(excluded['unselected_indicators']) == {f for f in selections[stratum['id']] if fields[f][i] != 1})
            if not expected:
                check(stratum['status'] == 'no_joint_observations' and 'models' not in stratum)
                continue
            weights = [fields[query['metric']['field']][i] if isinstance(query['metric'], dict) else Q(1) for i in expected]
            check(list(map(Q, stratum['metric_diagonal'])) == weights)
            fits = {}
            for model in stratum['models']:
                design = demo.dense(model['design'])
                meanings = {c['id']: c for c in result['columns']}
                truth = [[product(fields[f][i] for f in meanings[c]['factors']) for c in model['columns']] for i in expected]
                check(design == truth)
                _, pivots = rref(design)
                check(model['rank'] == len(pivots))
                check(model['coefficient_nullity'] == model['design']['columns'] - len(pivots))
                kernel = demo.dense(model['coefficient_kernel'])
                check(model['coefficient_kernel']['rows'] == model['design']['columns'])
                for k in range(model['coefficient_nullity']):
                    check(all(sum((r[j]*kernel[j][k] for j in range(len(r))), Q(0)) == 0 for r in design))
                check(model['coefficient_identified'] == [all(v == 0 for v in row) for row in kernel])
                fits[model['id']] = {}
                for reading in model['responses']:
                    y = [fields[reading['field']][i] for i in expected]
                    check(list(map(Q, reading['input'])) == y)
                    fitted = represented(design, y, weights)
                    residual = [a-b for a,b in zip(y,fitted)]
                    check(list(map(Q, reading['represented'])) == fitted)
                    check(list(map(Q, reading['residual'])) == residual)
                    check(Q(reading['residual_quadrance']) == dot(residual,residual,weights))
                    coefficient = list(map(Q, reading['particular_coefficients']))
                    check([sum((a*b for a,b in zip(r,coefficient)), Q(0)) for r in design] == fitted)
                    check(all(Q(v) == 0 for v in reading['normal_equations']))
                    for contrast in reading.get('coefficient_contrasts', []):
                        w = list(map(Q, contrast['weights']))
                        ambiguity = [sum((w[j]*kernel[j][k] for j in range(len(w))), Q(0)) for k in range(model['coefficient_nullity'])]
                        check(list(map(Q, contrast['kernel_pairings'])) == ambiguity)
                        check(contrast['identified'] == all(v == 0 for v in ambiguity))
                        check(contrast['value'] is None if not contrast['identified'] else Q(contrast['value']) == sum((a*b for a,b in zip(w,coefficient)),Q(0)))
                    fits[model['id']][reading['field']] = (fitted,residual)
            for comparison in stratum['comparisons']:
                for reading in comparison['responses']:
                    full = fits['full'][reading['field']]
                    reduced = fits[comparison['reduced']][reading['field']]
                    contribution = [a-b for a,b in zip(full[0],reduced[0])]
                    gain = dot(reduced[1],reduced[1],weights) - dot(full[1],full[1],weights)
                    check(list(map(Q,reading['added_representation'])) == contribution)
                    check(Q(reading['added_quadrance']) == gain == dot(contribution,contribution,weights))
                    check(gain >= 0 and dot(contribution,full[1],weights) == 0)
        return result

    receipts = {}
    for name, query in queries.items():
        receipts[name] = execute(query)
        if name != 'relation_fields':
            verify(receipts[name],query)
        else:
            result = receipts[name]['result']
            model = result['strata'][0]['models'][1]
            check(result['grade'] == 1 and result['strata'][0]['n'] == 3)
            check(list(map(Q,model['responses'][0]['particular_coefficients'])) == [Q(3),Q(2)])
            check(model['responses'][0]['represented'] == ['3','5','7'])
            check(model['responses'][0]['residual_quadrance'] == '0')
        demo.save(root/f'{name}.certificate.json', receipts[name])
    continuous = receipts['continuous']['result']['strata'][0]
    full = continuous['models'][1]
    check(list(map(Q,full['responses'][0]['particular_coefficients'])) == list(map(Q,[5,2,2,-3,4,5])))
    check(Q(full['responses'][0]['residual_quadrance']) == Q(16,3))
    gains = {r['reduced']: Q(r['responses'][0]['added_quadrance']) for r in continuous['comparisons']}
    check(gains['without_interaction'] == 100 and gains['without_arm'] == 196 and gains['without_z'] == 216)
    redundant = receipts['redundant']['result']['strata'][0]['models'][1]
    check(redundant['coefficient_nullity'] == 2)
    contrasts = redundant['responses'][0]['coefficient_contrasts']
    check(contrasts[0]['value'] == '2' and contrasts[1]['value'] is None and contrasts[2]['value'] == '4')
    check(receipts['near_collinear']['result']['strata'][0]['models'][1]['rank'] == 3)
    check(receipts['near_collinear']['result']['strata'][0]['models'][1]['responses'][0]['residual_quadrance'] == '0')
    check(all(r['residual_quadrance'] == '0' for r in receipts['category_quadratic']['result']['strata'][0]['models'][1]['responses']))
    def without_execution(value):
        if isinstance(value,dict): return {k:without_execution(v) for k,v in value.items() if k!='execution'}
        if isinstance(value,list): return [without_execution(v) for v in value]
        return value
    serial = {**queries['categories'],'workers':1}
    check(without_execution(execute(serial)['result']) == without_execution(receipts['categories']['result']))
    multi = receipts['categories']['result']['strata'][0]['models'][1]['execution']
    check(multi['gram_preparations'] == 1 and multi['applications'] == 3 and multi['worker_slots'] == 2)
    partial = copy.deepcopy(queries['continuous'])
    partial['terms'][0]['factors'] = ['all.partial_x']
    check(verify(execute(partial),partial)['strata'][0]['n'] == 23)
    zero = {**queries['categories'], 'terms':[], 'intercept':False, 'ablations':[]}
    check(verify(execute(zero),zero)['strata'][0]['models'][1]['rank'] == 0)
    constant = {**queries['categories'], 'selection':['all.zero']}
    constant_result = verify(execute(constant),constant)
    check(constant_result['strata'][0]['models'][1]['coefficient_nullity'] == 1)
    for changes in [dict(responses=[]), dict(responses=['all.x']), dict(selection=['all.x']),
        dict(metric='invented'), dict(metric=dict(kind='diagonal_field',field='all.x',meaning='Not positive')),
        dict(terms=[dict(id='intercept',factors=['all.x'])]), dict(terms=[dict(id='x',factors=[])]),
        dict(terms=[dict(id='x',factors=['all.x']),dict(id='x',factors=['all.z'])]),
        dict(controls=['all.x','all.x']), dict(support='impute'), dict(workers=0),
        dict(ablations=[dict(id='full',remove=['x'])]), dict(ablations=[dict(id='wrong',remove=['unknown'])]),
        dict(coefficient_contrasts=[dict(id='wrong',weights=[dict(column='unknown',value='1')])]),
        dict(coefficient_contrasts=[dict(id='wrong',weights=[dict(column='x',value='1/0')])]),
        dict(strata=[dict(id='all',selection=[])]), dict(extra=True)]:
        execute({**queries['categories'],**changes},fail=True)
    execute({**queries['relation_fields'],'responses':['all.response']},fail=True)
    check(execute(queries['categories']) == receipts['categories'])
    print(json.dumps(dict(status='PASS', checks=checks, rejections=rejected, authored_cases=len(queries))))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
