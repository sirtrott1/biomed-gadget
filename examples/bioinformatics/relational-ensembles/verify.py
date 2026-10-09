#!/usr/bin/env python3
"""Independent exact reference for retained source fields and ensemble certificates."""
import argparse
import csv
import json
from pathlib import Path

import sympy as s


def sparse(record):
    out=s.zeros(record['rows'],record['columns'])
    for entry in record['entries']:out[entry['row'],entry['column']]=s.Rational(entry['value'])
    return out


def vector(values):return s.Matrix([s.Rational(v) for v in values]) if values else s.zeros(0,1)


def frame(columns, n):return s.Matrix.hstack(*columns) if columns else s.zeros(n,0)


def projector(f, metric):return f*(f.T*metric*f).inv()*f.T*metric if f.cols else s.zeros(f.rows)


def verify(root, certificate=None):
    root=Path(root)
    cert=certificate or json.loads((root/'certificate.json').read_text())
    q=cert['query'];out=cert['result'];grade=out['grade'];coordinates=out['coordinate_carriers']
    relations=json.loads((root/'inputs/relations.json').read_text())
    n=[len(c) for c in coordinates]
    positions=[{c['entity']:i for i,c in enumerate(cells)} for cells in coordinates]
    boundaries=[s.zeros(0,n[0])]+[s.zeros(n[k-1],n[k]) for k in range(1,len(n))]
    for relation in relations['relations']:
        k=relation['grade'];j=positions[k]['relation:'+relation['id']]
        for occurrence in relation['participants']:
            if occurrence['exists']:
                entity=occurrence['cell'] if k==1 else 'relation:'+occurrence['cell']
                boundaries[k][positions[k-1][entity],j]+=s.Rational(occurrence['magnitude'])*(-1 if occurrence['negative'] else 1)
    metrics=[]
    for k,m in enumerate(cert['context']['metrics']):
        metrics.append(s.eye(n[k]) if m['kind']=='identity' else s.diag(*[s.Rational(m['values'][cell['ambient_coordinate']]) for cell in coordinates[k]]))
    M=metrics[grade];checks=0
    def check(condition):
        nonlocal checks
        assert condition
        checks+=1
    check(sparse(out['metric'])==M)
    for k in range(1,len(boundaries)-1):check(boundaries[k]*boundaries[k+1]==s.zeros(n[k-1],n[k+1]))
    lower=boundaries[grade]
    upper=boundaries[grade+1] if grade+1<len(n) else s.zeros(n[grade],0)
    lower_adj=M.inv()*lower.T*metrics[grade-1] if grade else s.zeros(n[grade],0)
    upper_adj=metrics[grade+1].inv()*upper.T*M if grade+1<len(n) else s.zeros(0,n[grade])
    g=frame(lower_adj.columnspace(),n[grade]);c=frame(upper.columnspace(),n[grade])
    h=frame(lower.col_join(upper_adj).nullspace(),n[grade])
    P=[s.eye(n[grade]),projector(g,M),projector(c,M),projector(h,M)]
    Ldown=lower_adj*lower;Lup=upper*upper_adj;L=Ldown+Lup
    check(sum(P[1:],s.zeros(n[grade]))==s.eye(n[grade]))
    for name,F in zip(('gradient','curl'),(g,c)):
        actual=sparse(out['frames'][name]['frame'])
        check(actual.cols==F.cols and projector(actual,M)==projector(F,M))
        check(sparse(out['frames'][name]['gram'])==actual.T*M*actual)
    check(out['frames']['harmonic']['rank']==h.cols and not out['frames']['harmonic']['basis_materialized'])
    original={}
    if grade:
        original={f['id']:s.Matrix([s.Rational(next(v['value'] for v in f['values'] if 'relation:'+v['relation']==cell['entity'])) for cell in coordinates[grade]]) for f in relations['fields'] if f['grade']==grade and all(any('relation:'+v['relation']==cell['entity'] for v in f['values']) for cell in coordinates[grade])}
    else:
        with (root/'inputs/measurements.tsv').open() as handle:
            rows=list(csv.reader(handle,delimiter='\t'));header=rows[0];table={r[0]:r[1:] for r in rows[1:]}
        original={name:s.Matrix([s.Rational(table[cell['entity']][j]) for cell in coordinates[0]]) for j,name in enumerate(header[1:])}
    weights=vector(out['weights']);sample_fields={}
    for item in out['samples']:
        x=original[item['field'].removeprefix('all.')]
        sample_fields[item['field']]=[p*x for p in P]
        for name,v in zip(('full','gradient','curl','harmonic'),sample_fields[item['field']]):check(vector(item[name])==v)
        check(vector(item['hodge_down'])==Ldown*x and vector(item['hodge_up'])==Lup*x and vector(item['hodge'])==L*x)
        expected=[s.zeros(size,1) for size in n]
        if grade:expected[grade-1]=lower*x
        if grade+1<len(n):expected[grade+1]=upper_adj*x
        check([vector(v) for v in item['dirac']]==expected)
        check(vector(item['boundary'])==lower*x)
        check(vector(item['upper_adjoint'])==upper_adj*x)
        for green in item['green']:
            check(vector(green['field'])==(s.eye(n[grade])+s.Rational(green['lambda'])*L).inv()*x)
            check(green['residual_exact'])
        check(item['reconstruction_exact'] and item['orthogonality_exact'] and item['dirac_square_exact'])
    means={}
    for group in out['groups']:
        ids=group['fields'];w=vector(group['conditional_weights']);D=s.diag(*w)
        expected_mass=sum(weights[q['fields'].index(name)] for name in ids)
        check(s.Rational(group['original_mass'])==expected_mass)
        check(w==s.Matrix([weights[q['fields'].index(name)]/expected_mass for name in ids]))
        means[group['id']]=[];centered=[];variances=[]
        for channel,name in enumerate(('full','gradient','curl','harmonic')):
            values=[sample_fields[id][channel] for id in ids]
            mean=sum((weight*v for weight,v in zip(w,values)),s.zeros(n[grade],1))
            C=s.Matrix.hstack(*[v-mean for v in values]).T
            Sigma=C.T*D*C;K=C*M*C.T
            V=s.trace(M*Sigma);S2=s.trace(M*Sigma*M*Sigma)
            reading=group['ensembles'][name]
            check(vector(reading['mean'])==mean and sparse(reading['centered_fields'])==C)
            check(reading['rank']==C.rank())
            check(s.Rational(reading['variance'])==V and s.Rational(reading['second_contraction'])==S2)
            check(V==s.trace(D*K) and S2==s.trace(D*K*D*K))
            check(reading['contraction_space']==('samples' if len(ids)<n[grade] else 'source_coordinates'))
            check(sparse(reading['second_moment_or_gram'])==(K if len(ids)<n[grade] else Sigma))
            if V:
                check(s.Rational(reading['effective_dimension_2'])==V**2/S2)
                check(s.Rational(reading['collision_probability'])==S2/V**2)
            else:check(reading['effective_dimension_2'] is None and reading['collision_probability'] is None)
            centered.append(C);variances.append(V);means[group['id']].append(mean)
        check(variances[0]==sum(variances[1:]))
        for coupling in group['sector_coupling']:
            a=('full','gradient','curl','harmonic').index(coupling['left']);b=('full','gradient','curl','harmonic').index(coupling['right'])
            if variances[a] and variances[b]:
                cross=centered[a].T*D*centered[b]
                gamma=s.trace(M*cross*M*cross.T)/(variances[a]*variances[b])
                check(s.Rational(coupling['value'])==gamma and coupling['status']=='completed')
            else:check(coupling['value'] is None and coupling['status']=='zero_sector_variance')
    for comparison in out['comparisons']:
        d=comparison['declaration'];left=means[d['left']];right=means[d['right']]
        for channel,name in enumerate(('full','gradient','curl','harmonic')):
            delta=right[channel]-left[channel]
            check(vector(comparison['right_minus_left'][name]['field'])==delta)
            check(s.Rational(comparison['right_minus_left'][name]['quadrance'])==(delta.T*M*delta)[0])
    check(out['execution']['projector_preparations']==2)
    check(out['execution']['applications_per_projector']==len(q['fields']))
    check(out['execution']['completed_tiles']==len(q['fields']))
    return dict(status='PASS',checks=checks,fields=len(q['fields']),source_dimensions=n)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    a=p.parse_args();print(json.dumps(verify(a.run)))
