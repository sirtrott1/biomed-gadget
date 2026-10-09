#!/usr/bin/env python3
"""Independently check physical fields and global optimum certificates with rational arithmetic."""
import argparse
import csv
import json
from fractions import Fraction as Q
from pathlib import Path


def matrix(record):
    out = [[Q(0)] * record['columns'] for _ in range(record['rows'])]
    for entry in record['entries']:
        out[entry['row']][entry['column']] = Q(entry['value'])
    return out


def apply(a, x):
    return [sum((v*w for v, w in zip(row, x)), Q(0)) for row in a]


def verify(root):
    certificate = json.loads((root/'certificate.json').read_text())
    result, query = certificate['result'], certificate['query']
    physical = result['physical_shape']
    study = json.loads(next((root/'workspace').rglob('study.json')).read_text())
    metadata = list(csv.DictReader((root/'inputs/observations.tsv').open(), delimiter='\t'))
    n = len(metadata)
    metric = study['contexts'][0]['metrics'][0]
    weights = [Q(1)]*n if metric['kind'] == 'identity' else list(map(Q, metric['values']))
    fields = {}
    with (root/'inputs/measurements.tsv').open() as handle:
        rows = csv.reader(handle, delimiter='\t')
        assert next(rows)[1:] == [r['id'] for r in metadata]
        for row in rows:
            fields['all.'+row[0]] = list(map(Q, row[1:]))
    spec = json.loads((root/'inputs/import.json').read_text())['assays'][0]
    for field in spec['observation_fields']:
        fields['all.'+field['id']] = [Q(int(r[field['column']] == field['level']))
            if field['kind'] == 'level_indicator' else Q(r[field['column']]) for r in metadata]
    responses = {r['id']: r for r in query['responses']}
    states = {s['id']: s for s in result['states']}
    checks = 0

    def check(value):
        nonlocal checks
        if not value:
            raise AssertionError('Independent physical response check '+str(checks))
        checks += 1

    for fit in physical['results']:
        p = physical['preparations'][fit['preparation']]
        p = dict(physical['families'][p['family']], **p)
        check((fit['response'], fit['model']) == (p['response'], p['model']))
        response = responses[fit['response']]
        original = fields[states[fit['state']]['response_fields'][fit['response']]]
        if response['kind'] == 'ordinal':
            codes = [Q(v['code']) for v in response['levels']]
            inputs = [[Q(int(codes.index(v) >= cut)) for v in original] for cut in range(1, len(codes))]
        else:
            inputs = [original]
        f, c, a, gram = map(matrix, (p['frame'], p['source_normals'], p['normals'], p['gram']))
        expected_rows, expected_rhs = [], []
        def row(entries, bound):
            values = [Q(0)]*(n*len(inputs))
            for i,value in entries.items():
                values[i] = Q(value)
            expected_rows.append(values);expected_rhs.append(Q(bound))
        for i in range(n):
            if response['kind'] in ('count','binomial_count'):
                row({i:-1},0)
                if response['kind'] == 'binomial_count':
                    row({i:1},fields[response['denominator']][i])
            elif response['kind'] == 'ordinal':
                row({i:1},1)
                for k in range(1,len(inputs)):
                    row({(k-1)*n+i:-1,k*n+i:1},0)
                row({(len(inputs)-1)*n+i:-1},0)
        for declaration in query['physical_shape']['dose_constraints']:
            if declaration['response'] != response['id']:
                continue
            dose = next(d for d in query['dose_comparisons'] if d['id'] == declaration['comparison'])
            groups = [[i for i in range(n) if fields[g['indicator']][i] == 1] for g in dose['groups']]
            levels = [fields[dose['field']][indices[0]] for indices in groups]
            denominator = fields[response['denominator']] if 'denominator' in response else [Q(1)]*n
            sign = -1 if dose['direction'] == 'nondecreasing' else 1
            for axis in range(len(inputs)):
                for k in range(len(groups)-1):
                    entries = {}
                    spacing = levels[k+1]-levels[k]
                    for direction,indices in ((-1,groups[k]),(1,groups[k+1])):
                        for i in indices:
                            entries[axis*n+i] = sign*direction/(spacing*len(indices)*denominator[i])
                    row(entries,0)
        check(c == expected_rows and list(map(Q,p['rhs'])) == expected_rhs)
        v, v0, rhs = map(lambda x: list(map(Q, x)), (fit['coordinates'], fit['unconstrained_coordinates'], p['rhs']))
        predicted = apply(f, v)
        for i, row in enumerate(a):
            for j in range(len(v)):
                check(row[j] == sum((c[i][k]*f[k][j] for k in range(len(f))), Q(0)))
        for i in range(len(v)):
            for j in range(len(v)):
                check(gram[i][j] == sum((weights[k % n]*f[k][i]*f[k][j] for k in range(len(f))), Q(0)))
        slack = [b-value for b, value in zip(rhs, apply(c, predicted))]
        check(slack == list(map(Q, fit['slack'])) and all(s >= 0 for s in slack))
        duals, active = list(map(Q, fit['multipliers'])), fit['active_rows']
        check(len(duals) == len(active) and len(active) == len(set(active)))
        check(all(v >= 0 for v in duals) and all(slack[i] == 0 for i in active))
        gradient = apply(gram, [x-y for x,y in zip(v,v0)])
        stationarity = [gradient[j]+sum((duals[k]*a[i][j] for k,i in enumerate(active)), Q(0)) for j in range(len(v))]
        check(stationarity == list(map(Q, fit['half_gradient_stationarity'])) and all(x == 0 for x in stationarity))
        objective = Q(0)
        for k, part in enumerate(fit['fits']):
            represented = predicted[k*n:(k+1)*n]
            residual = [x-y for x,y in zip(inputs[k], represented)]
            check(represented == list(map(Q, part['objective_represented'])))
            check(residual == list(map(Q, part['objective_residual'])))
            q = sum((w*x*x for w,x in zip(weights,residual)), Q(0))
            check(q == Q(part['residual_quadrance']))
            objective += q
            denominator = fields[response['denominator']] if 'denominator' in response else [Q(1)]*n
            readout = [x/e for x,e in zip(represented,denominator)]
            check(readout == list(map(Q, part['represented']['sectors']['full'])))
            if response['kind'] in ('count','binomial_count','ordinal'):
                check(all(x >= 0 for x in readout))
                if response['kind'] != 'count':
                    check(all(x <= 1 for x in readout))
            for field in ('represented','residual'):
                sectors = {name:list(map(Q, values)) for name,values in part[field]['sectors'].items()}
                check(all(sectors['full'][i] == sum(sectors[name][i] for name in ('gradient','curl','harmonic')) for i in range(n)))
            for d in query['physical_shape']['dose_constraints']:
                if d['response'] != response['id']:
                    continue
                dose = next(x for x in query['dose_comparisons'] if x['id'] == d['comparison'])
                means = []
                for group in dose['groups']:
                    values = [readout[i] for i in range(n) if fields[group['indicator']][i] == 1]
                    means.append(sum(values,Q(0))/len(values))
                sign = 1 if dose['direction'] == 'nondecreasing' else -1
                check(all(sign*(b-a) >= 0 for a,b in zip(means,means[1:])))
        check(objective == Q(fit['joint_residual_quadrance']))
        correction = [x-y for x,y in zip(v,v0)]
        cost = sum((x*y for x,y in zip(correction,apply(gram,correction))),Q(0))
        check(objective == Q(fit['unconstrained_joint_quadrance'])+cost)
        check(cost == Q(fit['constraint_cost']))
        if response['kind'] == 'ordinal':
            categories = [list(map(Q,x)) for x in fit['ordinal_level_fields']]
            check(all(0 <= value <= 1 for row in categories for value in row))
            check(all(sum(row[i] for row in categories) == 1 for i in range(n)))
            for k, category in enumerate(categories):
                check(category == [(predicted[(k-1)*n+i] if k else Q(1))-
                                    (predicted[k*n+i] if k < len(inputs) else Q(0)) for i in range(n)])
    for c in physical['nested_comparisons']:
        check(Q(c['joint_gain']) >= 0)
        check(Q(c['joint_gain']) == Q(c['displacement_quadrance'])+Q(c['residual_cross_term']))
    return dict(status='PASS',checks=checks,reference='Independent Fraction arithmetic, original measurements and primal/dual global optimum proof')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.run.resolve())))
