#!/usr/bin/env python3
"""Independent rational reference for complete retained exchange experiments."""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import sympy as s


def verify(root,certificate="certificate.json",output="independent-verification.json"):
    root=root.resolve();cert=json.loads((root/certificate).read_text());query=cert['query'];out=cert['result']
    spec=json.loads((root/'inputs/import.json').read_text())['assays'][0]
    relations=json.loads((root/'inputs/relations.json').read_text())
    with (root/'inputs'/spec['source']['path']).open() as handle:
        rows=list(csv.reader(handle,delimiter='\t'))
    coordinates=[r[0] for r in rows[1:]] if spec['carrier']=='features' else rows[0][1:]
    grades=[coordinates]
    for k in range(1,max(r.get('grade',1) for r in relations['relations'])+1):
        grades.append(['relation:'+r['id'] for r in relations['relations'] if r.get('grade',1)==k])
    B=[None]*len(grades)
    for k in range(1,len(grades)):
        B[k]=s.zeros(len(grades[k-1]),len(grades[k]));positions={x:i for i,x in enumerate(grades[k-1])}
        for j,r in enumerate(r for r in relations['relations'] if r.get('grade',1)==k):
            for p in r['participants']:
                name=p['cell'] if k==1 else 'relation:'+p['cell']
                if p['exists']:B[k][positions[name],j]+=(-1 if p['negative'] else 1)*s.Rational(p['magnitude'])
    metrics=[]
    for k,ids in enumerate(grades):
        m=spec.get('metrics',[dict(kind='identity')]*len(grades))[k]
        metrics.append(s.eye(len(ids)) if m['kind']=='identity' else s.diag(*map(s.Rational,m['values'])))
    fields={}
    if spec['carrier']=='features':
        for j,name in enumerate(rows[0][1:]):fields['all.'+name]=s.Matrix([s.Rational(r[j+1]) for r in rows[1:]])
    else:
        for r in rows[1:]:fields['all.'+r[0]]=s.Matrix(list(map(s.Rational,r[1:])))
        metadata={r['id']:r for r in csv.DictReader((root/'inputs/observations.tsv').open(),delimiter='\t')}
        for f in spec.get('observation_fields',[]):
            fields['all.'+f['id']]=s.Matrix([int(metadata[name][f['column']]==f['level']) if f['kind']=='level_indicator' else s.Rational(metadata[name][f['column']]) for name in coordinates])
    for f in relations.get('fields',[]):
        positions={name:i for i,name in enumerate(grades[f['grade']])};x=s.zeros(len(positions),1)
        for v in f['values']:x[positions['relation:'+v['relation']]]=s.Rational(v['value'])
        fields['all.'+f['id']]=x
    k=out['grade'];M=metrics[k];n=M.rows
    adj=M.inv()*B[k].T*metrics[k-1] if k else s.zeros(n,0)
    upper=B[k+1] if k+1<len(grades) else s.zeros(n,0)
    def projector(matrix):
        columns=matrix.columnspace()
        F=s.Matrix.hstack(*columns) if columns else s.zeros(n,0)
        return F*(F.T*M*F).inv()*F.T*M if columns else s.zeros(n)
    Pg,Pc=projector(adj),projector(upper);Ps=dict(full=s.eye(n),gradient=Pg,curl=Pc,harmonic=s.eye(n)-Pg-Pc)
    L=adj*B[k] if k else s.zeros(n)
    if k+1<len(grades):L+=upper*metrics[k+1].inv()*upper.T*M
    Q=lambda x:(x.T*M*x)[0]
    checks=0
    def check(c):
        nonlocal checks
        if not c:raise AssertionError('Independent reference mismatch')
        checks+=1
    for j in range(2,len(grades)):check(B[j-1]*B[j]==s.zeros(len(grades[j-2]),len(grades[j])))
    blocks=[b['coordinates'] for b in out['exchange_blocks']]
    # Match the native mixed-radix order: the first block changes fastest.
    assignments=[list(reversed(v)) for v in itertools.product(*(list(itertools.permutations(b)) for b in reversed(blocks)))]
    check(len(assignments)==math.prod(math.factorial(len(b)) for b in blocks))
    scenarios=query.get('uncertainty',{}).get('scenarios',[dict(id='observed',weight='1',fields=query['fields'])])
    independent=[]
    for decl,scenario in zip(scenarios,out['scenarios']):
        scores=[]
        for state,assigned in enumerate(assignments):
            check(out['states'][state]['block_source_coordinates']==[list(b) for b in assigned])
            P=s.eye(n)
            for target,origin in zip(blocks,assigned):
                for i,j in zip(target,origin):P[i,i]=0;P[i,j]=1
            readings=[P*fields[f] for f in decl['fields']];row=[]
            for h in query['hypotheses']:
                x=readings[query['fields'].index(h['field'])];den=Q(x)
                if h['statistic']=='sector_fraction':num=Q(Ps[h['sector']]*x)
                elif h['statistic']=='probe_fraction':
                    y=Ps[h['sector']]*x;den=Q(y)
                    F=s.Matrix.hstack(*(Ps[h['sector']]*fields[f] for f in h['probes']))
                    captured=projector(F)*y;num=Q(captured)
                    if 'readings' in scenario['statistics'][state]:
                        detail=scenario['statistics'][state]['readings']['hypotheses'][len(row)]
                        check(s.Matrix(list(map(s.Rational,detail['captured'])))==captured)
                        check(Q(captured)+Q(y-captured)==den)
                else:
                    green=(s.eye(n)+s.Rational(h['lambda'])*L).inv()*x;num=Q(x-green)
                    check(green+s.Rational(h['lambda'])*L*green==x)
                row.append(num/den if den else None)
            check([str(v) if v is not None else None for v in row]==scenario['statistics'][state]['statistics']);scores.append(row)
            if 'readings' in scenario['statistics'][state]:
                for a,reading in enumerate(scenario['statistics'][state]['readings']['fields']):
                    for key in Ps:check(s.Matrix(list(map(s.Rational,reading[key])))==Ps[key]*readings[a])
        complete=all(v is not None for row in scores for v in row)
        maxima=[max(row) for row in scores] if complete else None
        weights=[s.Rational(state['mass']) for state in out['states']]
        for h,hyp in enumerate(scenario['hypotheses']):
            values=[row[h] for row in scores]
            if any(v is None for v in values):check(hyp['reference_tail_mass'] is None);continue
            obs=values[0];tail=sum(w for w,v in zip(weights,values) if v>=obs)
            check(s.Rational(hyp['reference_tail_mass'])==tail)
            check(s.Rational(hyp['tie_mass'])==sum(w for w,v in zip(weights,values) if v==obs))
            if maxima is not None:check(s.Rational(hyp['max_statistic_tail_mass'])==sum(w for w,v in zip(weights,maxima) if v>=obs))
            for threshold in set(values):
                level=sum(w for w,v in zip(weights,values) if v>=threshold)
                check(sum(w for w,v in zip(weights,values) if sum(q for q,x in zip(weights,values) if x>=v)<=level)<=level)
        if maxima is not None:check(s.Rational(scenario['joint_family']['global_tail_mass'])==sum(w for w,v in zip(weights,maxima) if v>=maxima[0]))
        independent.append(scores)
    if 'scenario_sensitivity' in out:
        for h,summary in enumerate(out['scenario_sensitivity']):
            values=[scores[0][h] for scores in independent];weights=[s.Rational(v['weight']) for v in scenarios]
            mean=sum(w*v for w,v in zip(weights,values))
            check(s.Rational(summary['weighted_mean'])==mean)
            check(s.Rational(summary['weighted_variance'])==sum(w*(v-mean)**2 for w,v in zip(weights,values)))
    result=dict(status='PASS',checks=checks,reference='Independent exact SymPy source matrices, metric frames, Green solves and complete weighted enumeration')
    (root/output).write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    result=verify(a.run)
    if (a.run/'alignment-certificate.json').exists():
        alignment=verify(a.run,'alignment-certificate.json','alignment-independent-verification.json')
        result=dict(status='PASS',checks=result['checks']+alignment['checks'],families=[result,alignment])
    print(json.dumps(result))
