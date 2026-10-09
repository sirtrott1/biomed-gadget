#!/usr/bin/env python3
"""Independent rational checks from original tables and relation declarations."""
import argparse,csv,json
from pathlib import Path
import sympy as S


def matrix(d):
    m=S.zeros(d['rows'],d['columns'])
    for e in d['entries']:m[e['row'],e['column']]=S.Rational(e['value'])
    return m

def vector(a):return S.Matrix([S.Rational(v) for v in a])

def projector(x,m):
    columns=x.columnspace();f=S.Matrix.hstack(*columns) if columns else S.zeros(x.rows,0)
    return f*(f.T*m*f).inv()*f.T*m if f.cols else S.zeros(x.rows)


def verify(root):
    cert=json.loads((root/'certificate.json').read_text());r=cert['result'];q=cert['query']
    study=json.loads(next((root/'workspace').rglob('study.json')).read_text())
    spec=json.loads((root/'inputs/import.json').read_text())['assays'][0]
    metadata=list(csv.DictReader((root/'inputs/observations.tsv').open(),delimiter='\t'))
    ids=[p['id'] for p in metadata];n=len(ids);fields={}
    with (root/'inputs/measurements.tsv').open() as f:
        rows=csv.reader(f,delimiter='\t');assert next(rows)[1:]==ids
        for row in rows:fields['all.'+row[0]]=vector(row[1:])
    for c in spec['observation_fields']:
        fields['all.'+c['id']]=vector([str(int(row[c['column']]==c['level'])) if c['kind']=='level_indicator' else row[c['column']] for row in metadata])
    relations=json.loads((root/'inputs/relations.json').read_text())['relations']
    cell_ids=[[c['id'].removeprefix('relation:') for c in study['cells'] if c['grade']==k] for k in range(len(study['contexts'][0]['metrics']))]
    boundaries=[S.zeros(0,n)]
    for k in range(1,len(cell_ids)):
        b=S.zeros(len(cell_ids[k-1]),len(cell_ids[k]))
        for rel in relations:
            if rel['grade']!=k:continue
            j=cell_ids[k].index(rel['id'])
            for p in rel['participants']:
                if p['exists']:b[cell_ids[k-1].index(p['cell']),j]+=(-1 if p['negative'] else 1)*S.Rational(p['magnitude'])
        boundaries.append(b)
    boundaries.append(S.zeros(len(cell_ids[-1]),0))
    metrics=[S.eye(len(cell_ids[k])) if m['kind']=='identity' else S.diag(*map(S.Rational,m['values'])) for k,m in enumerate(study['contexts'][0]['metrics'])]
    m=metrics[0];curl=projector(boundaries[1],m);harmonic=S.eye(n)-curl
    l=boundaries[1]*metrics[1].inv()*boundaries[1].T*m;parts=dict(full=S.eye(n),gradient=S.zeros(n),curl=curl,harmonic=harmonic)
    g=(S.eye(n)+S.Rational(q['green_parameter'])*l).inv() if 'green_parameter' in q else None
    checks=0
    def check(value):
        nonlocal checks
        assert value
        checks+=1
    def equal(a,b):check(a==b)
    if (root/'source/source.json').exists():
        receipt=json.loads((root/'source/source.json').read_text())
        with (root/'source/counts.tsv').open() as handle:
            rows=csv.reader(handle,delimiter='\t');samples=next(rows)[1:];totals=[0]*len(samples)
            for row in rows:totals=[a+int(b) for a,b in zip(totals,row[1:])]
        equal(totals,[receipt['library_totals'][sample] for sample in samples])
        equal(vector(totals),fields['all.retained_total'])
    equal(r['hodge_dimensions'],dict(gradient=0,curl=boundaries[1].rank(),harmonic=n-boundaries[1].rank()))
    for k,b in enumerate(r['source_boundaries']):equal(matrix(b),boundaries[k+1])
    responses={d['id']:d for d in q['responses']};models={};axis_values={};scores=[]
    terms=q['model']['terms'];columns=[('intercept',[])] if q['model']['intercept'] else []
    columns += [(d['id'],d['factors']) for d in terms]
    removal=dict(baseline={d['id'] for d in terms},full=set(),**{a['id']:set(a['remove']) for a in q['model'].get('ablations',[])})
    for p in r['preparations']:
        e=S.diag(*fields[p['denominator_field']]) if p['denominator_field'] else S.eye(n)
        equal(matrix(p['readout']),e.inv());equal(matrix(p['readout_metric_pullback_map']),e)
        for d in p['models']:
            xs=[]
            for name,factors in columns:
                if name in removal[d['id']]:continue
                xs.append(S.Matrix([S.prod(fields[f][i] for f in factors) for i in range(n)]))
            a=S.Matrix.hstack(*xs) if xs else S.zeros(n,0);x=e*a
            equal(matrix(d['design']),x);equal(d['rank'],x.rank())
            kernel=matrix(d['coefficient_kernel']);equal(x*kernel,S.zeros(n,kernel.cols));equal(kernel.cols,x.cols-x.rank());equal(kernel.rank(),kernel.cols)
            frame=matrix(d['independent_frame']);equal(frame.rank(),x.rank());equal(S.Matrix.hstack(x,frame).rank(),x.rank())
            models[(p['denominator_field'],d['id'])]=(x,projector(x,m),e.inv())
    def reading(data,x):
        for name,op in parts.items():equal(vector(data['sectors'][name]),op*x);equal(S.Rational(data['quadrances'][name]),(x.T*op.T*m*op*x)[0])
        equal(vector(data['hodge_action']),l*x);equal(vector(data['upper_adjoint']),metrics[1].inv()*boundaries[1].T*m*x)
        if g is not None:equal(vector(data['green']),g*x)
    for state in r['states']:
        objectives={};state_scores={}
        for f in state['fits']:
            response=responses[f['response']];x0=fields[state['response_fields'][response['id']]]
            if response['kind']=='ordinal':
                cut=int(f['axis'].split(':')[-1]);codes=[S.Rational(k['code']) for k in response['levels']]
                x0=S.Matrix([int(codes.index(v)>=cut) for v in x0])
            x,p,readout=models[(response.get('denominator'),f['model'])];represented=p*x0;residual=x0-represented;cost=(residual.T*m*residual)[0]
            equal(vector(f['objective_input']),x0);equal(vector(f['objective_represented']),represented);equal(vector(f['objective_residual']),residual)
            equal(x*vector(f['particular_coefficients']),represented);equal(S.Rational(f['residual_quadrance']),cost);equal(vector(f['normal_equations']),S.zeros(x.cols,1))
            reading(f['observed'],readout*x0);reading(f['represented'],readout*represented);reading(f['residual'],readout*residual)
            objectives[(f['axis'],f['model'])]=(represented,cost)
            if f['model']=='full':axis_values[(state['id'],f['axis'])]=readout*x0
        for c in state['nested_comparisons']:
            full,cf=objectives[(c['axis'],'full')];reduced,cr=objectives[(c['axis'],c['reduced'])]
            equal(vector(c['represented_change']),full-reduced);equal(S.Rational(c['gain']),cr-cf)
            equal(None if c['fraction'] is None else S.Rational(c['fraction']),None if cr==0 else (cr-cf)/cr)
            if c['reduced']=='baseline':state_scores[c['axis']]=None if cr==0 else (cr-cf)/cr
        scores.append(state_scores)
        for c in state['ordinal_level_readings']:
            fs=[f for f in state['fits'] if f['response']==c['response'] and f['model']==c['model']]
            cuts=[vector(f['represented']['sectors']['full']) for f in fs];levels=[S.ones(n,1)]+cuts+[S.zeros(n,1)]
            for k,actual in enumerate(c['level_fields']):equal(vector(actual),levels[k]-levels[k+1])
            equal(c['admissible'],all(v>=0 and v<=1 for a in c['level_fields'] for v in vector(a)))
        for dose in state['dose_comparisons']:
            d=next(d for d in q['dose_comparisons'] if d['id']==dose['id']);means=S.zeros(len(d['groups']),n);lift=S.zeros(n,len(d['groups']));levels=[]
            for k,group in enumerate(d['groups']):
                rows=[i for i,v in enumerate(fields[group['indicator']]) if v==1];levels.append(fields[d['field']][rows[0]])
                for i in rows:means[k,i]=S.Rational(1,len(rows));lift[i,k]=1
            equal(matrix(dose['mean_map']),means);equal(matrix(dose['lift_map']),lift);equal(vector(dose['dose_levels']),S.Matrix(levels))
            differences=S.zeros(len(levels)-1,len(levels))
            for k in range(len(levels)-1):differences[k,k]=-1/(levels[k+1]-levels[k]);differences[k,k+1]=-differences[k,k]
            equal(matrix(dose['divided_difference_map']),differences)
            for f in dose['responses']:
                x=axis_values[(state['id'],f['axis'])];equal(vector(f['means']),means*x);equal(vector(f['divided_differences']),differences*means*x);equal(vector(f['within_group_field']),x-lift*means*x)
    if 'reference_law' in q:
        weights=[S.Rational(d['mass']) for d in q['reference_law']['states']];observed=r['observed_state_index']
        calibration=r['reference_calibration'];complete=all(all(v is not None for v in score.values()) for score in scores)
        for h in calibration['hypotheses']:
            values=[score[h['id']] for score in scores]
            equal(h['status'],'completed' if all(v is not None for v in values) else 'undefined_statistic')
            if all(v is not None for v in values):
                equal(S.Rational(h['reference_tail_mass']),sum(w for w,v in zip(weights,values) if v>=values[observed]))
                if complete:equal(S.Rational(h['max_statistic_tail_mass']),sum(w for w,score in zip(weights,scores) if max(score.values())>=values[observed]))
    return dict(status='PASS',checks=checks,reference='SymPy exact rational projection; original tables and source declarations')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    result=verify(a.run.resolve());(a.run/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
