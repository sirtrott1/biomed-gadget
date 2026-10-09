#!/usr/bin/env python3
"""Independent exact arithmetic for temporal and censored certificates."""
import csv,itertools,json
from fractions import Fraction as Q
from pathlib import Path

def vec(x):return list(map(Q,x))
def matrix(x):
    out=[[Q(0)]*x['columns'] for _ in range(x['rows'])]
    if 'values' in x:return [vec(row) for row in x['values']]
    for e in x['entries']:out[e['row']][e['column']]=Q(e['value'])
    return out

def add(a,b):return [x+y for x,y in zip(a,b)]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def apply(a,x):return [dot(row,x) for row in a]
def mul(a,b):return [[dot(row,col) for col in zip(*b)] for row in a]
def transpose(a):return list(map(list,zip(*a)))
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def solve(a,b):
    a=[list(row)+[value] for row,value in zip(a,b)];n=len(a)
    for k in range(n):
        pivot=next(i for i in range(k,n) if a[i][k]);a[k],a[pivot]=a[pivot],a[k];d=a[k][k];a[k]=[v/d for v in a[k]]
        for i in range(n):
            if i!=k:
                d=a[i][k];a[i]=[v-d*w for v,w in zip(a[i],a[k])]
    return [row[-1] for row in a]

def verify_history(root):
    certificate=json.loads((root/'certificate.json').read_text());r=certificate['result'];q=certificate['query'];checks=0
    def check(value):
        nonlocal checks
        assert value;checks+=1
    B=[[[Q(-1),Q(0)],[Q(1),Q(0)],[Q(0),Q(0)]],[[Q(-1),Q(0),Q(0)],[Q(1),Q(-1),Q(0)],[Q(0),Q(1),Q(0)]],[[Q(-1),Q(0)],[Q(1),Q(0)],[Q(0),Q(0)]]]
    L=[[[Q(2),Q(0)],[Q(0),Q(0)]],[[Q(2),Q(-1),Q(0)],[Q(-1),Q(2),Q(0)],[Q(0),Q(0),Q(1)]],[[Q(2),Q(0)],[Q(0),Q(0)]]]
    D=[[[Q(1),Q(0)],[Q(0),Q(2)]],[[Q(2),Q(0),Q(0)],[Q(0),Q(1),Q(0)],[Q(0),Q(0),Q(4)]],[[Q(3),Q(0)],[Q(0),Q(2)]]]
    xs=[vec([1,2]),vec([3,1,6]),vec([13,-7])];js=[matrix(step['source_maps'][1]) for step in r['steps']]
    if q['word']==[dict(kind='boundary')]:words=[[b] for b in B]
    else:words=[[d,l,b] for d,l,b in zip(D,L,B)]
    composed=[]
    for t,(e,x,factors) in enumerate(zip(r['endpoints'],xs,words)):
        check(vec(e['field'])==x);p=eye(len(x))
        for factor in factors:p=mul(factor,p)
        composed.append(p);check(vec(e['word_field'])==apply(p,x))
    channels=[apply(mul(js[1],js[0]),xs[0])]
    for t,step in enumerate(r['steps']):
        j=js[t];innovation=sub(xs[t+1],apply(j,xs[t]));check(vec(step['structure']['field_innovation'])==innovation)
        sector_fields=step['structure']['sector_transport'];check([sum(Q(sector_fields[a][b][i]) for a in range(3) for b in range(3)) for i in range(len(xs[t+1]))]==apply(j,xs[t]))
        sector_maps=[matrix(m) for row in step['structure']['sector_maps'] for m in row]
        check([[sum(m[i][k] for m in sector_maps) for k in range(len(xs[t]))] for i in range(len(xs[t+1]))]==j)
        oldg=solve([[v+int(i==k) for k,v in enumerate(row)] for i,row in enumerate(L[t])],xs[t]);newg=solve([[v+int(i==k) for k,v in enumerate(row)] for i,row in enumerate(L[t+1])],xs[t+1]);check(vec(step['structure']['old']['green'])==oldg);check(vec(step['structure']['new']['green'])==newg);check(vec(step['structure']['total'])==sub(newg,apply(j,oldg)))
        delta=sub(apply(composed[t+1],xs[t+1]),apply(composed[t],xs[t]));check(vec(step['field_delta'])==delta)
        check(Q(step['field_delta_quadrance'])==dot(delta,delta));gram=matrix(step['field_delta_channel_gram']);check(sum(map(sum,gram))==dot(delta,delta))
        fs=[vec(f['field']) for f in step['word_factors']]+[apply(composed[t+1],innovation)]
        check(gram==[[dot(a,b) for b in fs] for a in fs])
        oldform=mul(transpose(composed[t]),composed[t]);newform=mul(transpose(composed[t+1]),composed[t+1]);jx=apply(j,xs[t]);terms=step['endpoint_quadrance_change'];expected=[dot(jx,apply(newform,jx))-dot(xs[t],apply(oldform,xs[t])),dot(innovation,apply(newform,jx)),dot(jx,apply(newform,innovation)),dot(innovation,apply(newform,innovation))]
        for key,value in zip(('pullback_form','left_innovation','right_innovation','innovation_interaction'),expected):check(Q(terms[key])==value)
        check(Q(terms['total'])==sum(expected));check(Q(terms['total'])==dot(apply(composed[t+1],xs[t+1]),apply(composed[t+1],xs[t+1]))-dot(apply(composed[t],xs[t]),apply(composed[t],xs[t])))
        check(Q(terms['finite_rate'])==Q(terms['total'])/Q(terms['interval']))
        # Native factor attribution is compared to independent matrix products,
        # retaining order. Reordering D and L changes the middle endpoint.
        if len(words[t])>1 and t==1:check(mul(D[t],L[t])!=mul(L[t],D[t]))
        channels.append(apply(js[1],innovation) if t==0 else innovation)
    check([vec(c['field']) for c in r['terminal_channels']]==channels);check([sum(c[i] for c in channels) for i in range(2)]==xs[-1]);check(matrix(r['terminal_channel_gram'])==[[dot(a,b) for b in channels] for a in channels]);check(Q(r['terminal_quadrance'])==Q(218))
    return dict(status='PASS',checks=checks)

def verify_endpoints(root):
    certificate=json.loads((root/'certificate.json').read_text());r=certificate['result'];q=certificate['query'];rows=list(csv.DictReader((root/'inputs/observations.csv').open()));checks=0
    def check(value):
        nonlocal checks
        assert value;checks+=1
    check(len(rows)==len(r['observations']));groups=[list(range(len(rows)))]+[[i for i,row in enumerate(rows) if row[id.removeprefix('all.')]=='1'] for id in q['strata']];readout=matrix(r['readout']);check(readout==[[Q(int(i in group),len(group)) if group else Q(0) for i in range(len(rows))] for group in groups])
    for reading in r['readings']:
        s=Q(reading['time']);lower=[];upper=[];rmin=[];rmax=[]
        for row in rows:
            l=Q(row['lower']);u=None if row['upper_finite']=='0' else Q(row['upper']);lc=row['lower_closed']=='1';uc=row['upper_closed']=='1'
            lower.append(Q(u is not None and u<=s));upper.append(Q(l<s or (l==s and lc)));rmin.append(min(l,s));rmax.append(min(u,s) if u is not None else s)
        for name,expected in [('event_field_lower',lower),('event_field_upper',upper),('restricted_time_lower_field',rmin),('restricted_time_upper_field',rmax)]:check(vec(reading[name])==expected)
        for name,x in [('cumulative_event_fraction_lower',lower),('cumulative_event_fraction_upper',upper),('restricted_event_free_time_infimum',rmin),('restricted_event_free_time_supremum',rmax)]:check(vec(reading[name])==apply(readout,x))
        contrast=matrix(r['contrast_readout']);check(contrast==[[a-b for a,b in zip(readout[1],readout[2])]])
        for name,lo,hi,minimize in [('event_fraction_contrast_lower',lower,upper,True),('event_fraction_contrast_upper',lower,upper,False),('restricted_time_contrast_infimum',rmin,rmax,True),('restricted_time_contrast_supremum',rmin,rmax,False)]:
            expected=[sum((coefficient*(lo[i] if (coefficient>0)==minimize else hi[i]) for i,coefficient in enumerate(row)),Q(0)) for row in contrast];check(vec(reading[name])==expected)
        family=reading['admissible_event_family'];frame=matrix(family['frame']);ambiguous=[i for i in range(len(rows)) if lower[i]!=upper[i]];check(family['frame']['columns']==len(ambiguous));check(vec(family['offset'])==lower)
        for j,i in enumerate(ambiguous):check([frame[k][j] for k in range(len(rows))]==[Q(k==i) for k in range(len(rows))])
        # The importer retains category membership as positive source columns.
        category='sex' if 'sex' in rows[0] else 'group'
        labels=[row[category] for row in rows];sizes={label:labels.count(label) for label in labels}
        gradient=[[Q(int(a==b),sizes[a]) for b in labels] for a in labels]
        expected_gradient=apply(gradient,lower)
        check(vec(family['sectors']['curl']['offset'])==expected_gradient)
        check(vec(family['sectors']['gradient']['offset'])==[Q(0)]*len(rows))
        check(vec(family['sectors']['harmonic']['offset'])==sub(lower,expected_gradient))
        check(matrix(family['sectors']['curl']['frame'])==mul(gradient,frame) if ambiguous else family['sectors']['curl']['frame']['columns']==0)
        sectors=family['sectors'];check([sum(Q(sectors[k]['offset'][i]) for k in ('gradient','curl','harmonic')) for i in range(len(rows))]==lower)
        sector_frames={k:matrix(v['frame']) for k,v in sectors.items()};check([[sum(sector_frames[k][i][j] for k in ('gradient','curl','harmonic')) for j in range(len(ambiguous))] for i in range(len(rows))]==frame)
        # Independently enumerate every binary completion in the small authored
        # example. Larger clinical domains use the exact separable signed sum.
        for k,v in sectors.items():
            offset=vec(v['offset']);f=sector_frames[k]
            lo=[a+sum(min(x,Q(0)) for x in line) for a,line in zip(offset,f)];hi=[a+sum(max(x,Q(0)) for x in line) for a,line in zip(offset,f)];check(vec(v['coordinate_lower'])==lo);check(vec(v['coordinate_upper'])==hi)
            if len(rows)==8:
                outputs=[add(offset,apply(f,list(map(Q,bits)))) for bits in itertools.product((0,1),repeat=len(ambiguous))]
                check([min(v[i] for v in outputs) for i in range(8)]==lo);check([max(v[i] for v in outputs) for i in range(8)]==hi)
    return dict(status='PASS',checks=checks)

if __name__=='__main__':
    import sys
    root=Path(sys.argv[1]);c=json.loads((root/'certificate.json').read_text());out=verify_history(root) if c['query']['operation']=='temporal_field_history' else verify_endpoints(root);print(json.dumps(out))
