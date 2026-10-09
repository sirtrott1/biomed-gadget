#!/usr/bin/env python3
"""Independent rational checks of publication-indexed certificates."""
import argparse,csv,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
checks=0

def check(value):
    global checks
    assert value;checks+=1


def sparse(m):return {(e['row'],e['column']):Q(e['value']) for e in m['entries']}
def apply(m,x):
    y=[Q(0)]*m['rows']
    for (i,j),v in sparse(m).items():y[i]+=v*x[j]
    return y

def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def vector(x):return list(map(Q,x))
def sha(path):return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()


def common(root,name='certificate'):
    cert=json.loads((root/(name+'.json')).read_text());result=cert['result'];query=cert['query'];check(result['identity_verified'])
    for e in result['endpoints']:
        binding=e['declaration']['publication_binding'];path=root/e['declaration']['publication']/'study.json';metadata=json.loads(path.read_text())
        check(binding['declaration_sha256']==sha(path));check(binding['native_sha256']==sha(path.parent/metadata['native']['path']))
        check(binding['source_namespace']==metadata['native']['source_namespace'])
        check([Q(v['execution_value']['rational']) for v in e['declaration']['observations']]==vector(e['field']))
    for t,s in enumerate(result['steps']):
        a,b=result['endpoints'][t:t+2];j=s['source_maps'][a['grade']];x,y=vector(a['field']),vector(b['field']);jx=apply(j,x);eta=[v-u for u,v in zip(jx,y)]
        check(eta==vector(s['structure']['field_innovation']))
        total=[Q(0)]*len(s['field_delta'])
        for c in s['field_delta_channels']:total=[u+v for u,v in zip(total,vector(c['field']))]
        check(total==vector(s['field_delta']));check(dot(total,total)==Q(s['field_delta_quadrance']))
        gram=sparse(s['field_delta_channel_gram']);cs=[vector(c['field']) for c in s['field_delta_channels']]
        for i,u in enumerate(cs):
            for k,v in enumerate(cs):check(gram.get((i,k),Q(0))==dot(u,v))
        terms=s['endpoint_quadrance_change'];check(sum(Q(terms[k]) for k in ('pullback_form','left_innovation','right_innovation','innovation_interaction'))==Q(terms['total']))
        check(Q(terms['total'])==dot(vector(b['word_field']),vector(b['word_field']))-dot(vector(a['word_field']),vector(a['word_field'])))
        check(Q(terms['finite_rate'])*Q(terms['interval'])==Q(terms['total']))
        maps=s['structure']['sector_maps'];summed={}
        for row in maps:
            for m in row:
                for ij,v in sparse(m).items():summed[ij]=summed.get(ij,Q(0))+v
        check({ij:v for ij,v in summed.items() if v}==sparse(j))
        green=s['structure'];check(vector(green['total'])==[sum(vs,Q(0)) for vs in zip(*[vector(green[k]) for k in ('realization_term','innovation_term','operator_term')])])
        for g in s['source_changes']:
            for c in g['occurrence_changes']:
                check(Q(c['new_coefficient'])-Q(c['old_coefficient'])==Q(c['coefficient_change']))
                check(sum(Q(c[k]) for k in ('existence_term','orientation_term','magnitude_term'))==Q(c['coefficient_change']))
    channels=[vector(c['field']) for c in result['terminal_channels']];terminal=vector(result['endpoints'][-1]['field']);check([sum(v,Q(0)) for v in zip(*channels)]==terminal)
    g=sparse(result['terminal_channel_gram'])
    for i,u in enumerate(channels):
        for j,v in enumerate(channels):check(g.get((i,j),Q(0))==dot(u,v))
    check(sum(g.values(),Q(0))==dot(terminal,terminal)==Q(result['terminal_quadrance']))
    return cert


def authored(root):
    import sympy as S
    for name in ('certificate','relations','completion'):
        cert=common(root,name);result=cert['result'];ends=result['endpoints'];query=cert['query'];operators=[]
        for e in ends:
            d=e['declaration'];grades=d['source_grades'];bs=[S.zeros(0,len(grades[0]))]
            for g in range(1,len(grades)):
                lo={c['key']:i for i,c in enumerate(grades[g-1])};up={c['key']:i for i,c in enumerate(grades[g])};b=S.zeros(len(lo),len(up))
                for o in d['source_occurrences']:
                    if o['grade']==g:b[lo[o['lower_key']],up[o['parent_key']]]+=S.Rational(o['coefficient'])
                bs.append(b)
            check(bs[1]*bs[2]==S.zeros(bs[1].rows,bs[2].cols))
            k=e['grade'];op=S.eye(len(grades[k]));h=bs[k].T*bs[k]+(bs[k+1]*bs[k+1].T if k+1<len(bs) else S.zeros(len(grades[k])))
            for f in query['word']:
                if f['kind']=='hodge':op=h*op
                elif f['kind']=='green':op=(S.eye(h.rows)+S.Rational(f['lambda'])*h).inv()*op
                elif f['kind']=='boundary':op=bs[k]*op;k-=1
            x=S.Matrix([S.Rational(v) for v in e['field']]);check(list(op*x)==[S.Rational(v) for v in e['word_field']]);operators.append(op)
            # Full independent projection for each actually admitted source.
            k=e['grade'];lower=bs[k].T;upper=bs[k+1] if k+1<len(bs) else S.zeros(len(grades[k]),0)
            def projector(m):
                basis=m.columnspace()
                if not basis:return S.zeros(m.rows)
                f=S.Matrix.hstack(*basis);return f*(f.T*f).inv()*f.T
            pg,pc=projector(lower),projector(upper);ph=S.eye(pg.rows)-pg-pc
            if e is not ends[-1]:
                old=result['steps'][ends.index(e)]['structure']['old']
            else:old=result['steps'][-1]['structure']['new']
            for label,p in [('gradient',pg),('curl',pc),('harmonic',ph)]:check(list(p*x)==[S.Rational(v) for v in old[label]])
            check(list((S.eye(h.rows)+h).inv()*x)==[S.Rational(v) for v in old['green']])
            # Other independent checks below use full operators even where
            # source sectors have zero output for the current particular x.
            check(pg*pc==S.zeros(pg.rows));check(ph*ph==ph)
        for t,s in enumerate(result['steps']):
            def matrix(m):
                x=S.zeros(m['rows'],m['columns'])
                for ij,v in sparse(m).items():x[ij[0],ij[1]]=S.Rational(v.numerator,v.denominator)
                return x
            j=matrix(s['source_maps'][ends[t]['grade']]);ko=operators[t].T*operators[t];kn=operators[t+1].T*operators[t+1];x=S.Matrix(ends[t]['field']);y=S.Matrix(ends[t+1]['field']);eta=y-j*x
            expected=[(x.T*(j.T*kn*j-ko)*x)[0],(eta.T*kn*j*x)[0],(x.T*j.T*kn*eta)[0],(eta.T*kn*eta)[0]]
            for k,v in zip(('pullback_form','left_innovation','right_innovation','innovation_interaction'),expected):check(S.Rational(s['endpoint_quadrance_change'][k])==v)
    check([len(e['field']) for e in result['endpoints']]==[3,3,2])
    complete=json.loads((root/'certificate.json').read_text())['result'];check(complete['steps'][0]['source_changes'][0]['born_cell_keys']==['40']);check(complete['steps'][0]['source_changes'][0]['removed_cell_keys']==['10'])
    conditional=json.loads((root/'completion.json').read_text())['result'];o=conditional['endpoints'][1]['declaration']['observations'][1];check(not o['reported'] and o['observed_value'] is None and Q(o['execution_value']['rational'])==4)
    check(not conditional['steps'][0]['retained_observation_alignment'][1]['new']['reported'])
    return dict(status='PASS',checks=checks)


def clinical(root):
    cert=common(root);r=cert['result'];rows=list(csv.DictReader((root/'original/lung.csv').open()));check(len(rows)==228);sizes=[80,160,228]
    word_fields=[]
    for t,e in enumerate(r['endpoints']):
        original=rows[:sizes[t]];x=[Q(row['time']) for row in original];check(vector(e['field'])==x)
        staged=list(csv.DictReader((root/('input'+str(t))/'observations.csv').open()));check(staged==original)
        by={sex:[i for i,row in enumerate(original) if Q(row['sex'])==sex] for sex in (Q(1),Q(2))};sums={sex:sum((x[i] for i in ids),Q(0)) for sex,ids in by.items()};y=[sums[Q(row['sex'])]/(len(by[Q(row['sex'])])+1) for row in original];word_fields.append(y);check(vector(e['word_field'])==y)
        # For the authored original sex-membership relation, independent
        # projection gives its mean on each group and the within-group field.
        for g in e['declaration']['source_grades'][0]:check(g['label']==original[int(g['key'])-1]['source_row'])
        structure=r['steps'][t]['structure']['old'] if t<2 else r['steps'][-1]['structure']['new']
        coexact=[sums[Q(row['sex'])]/len(by[Q(row['sex'])]) for row in original]
        check(vector(structure['curl'])==coexact)
        check(vector(structure['harmonic'])==[u-v for u,v in zip(x,coexact)])
        check(vector(structure['gradient'])==[Q(0)]*len(x))
        check(vector(structure['green'])==[u-v for u,v in zip(x,y)])
        # The complete original variables remain admitted, including missingness.
        metadata=json.loads((root/('p'+str(t))/'study.json').read_text());fieldids={f['id'] for f in metadata['fields']};check(all('all.'+c in fieldids for c in rows[0] if c!='source_row'))
    for t,step in enumerate(r['steps']):
        nx=sizes[t+1]-sizes[t];check(step['source_maps'][0]['rows']==sizes[t+1] and step['source_maps'][0]['columns']==sizes[t]);check(len(step['source_changes'][0]['born_cell_keys'])==nx);check(not step['source_changes'][0]['removed_cell_keys'])
        innovation=vector(step['structure']['field_innovation']);check(innovation[:sizes[t]]==[Q(0)]*sizes[t]);check(innovation[sizes[t]:]==[Q(row['time']) for row in rows[sizes[t]:sizes[t+1]]])
        check(Q(step['endpoint_quadrance_change']['total'])==dot(word_fields[t+1],word_fields[t+1])-dot(word_fields[t],word_fields[t]))
    check(vector(r['terminal_channels'][0]['field'])==[Q(row['time']) for row in rows[:80]]+[Q(0)]*148)
    return dict(status='PASS',checks=checks)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args();receipt=clinical(a.root) if (a.root/'original').exists() else authored(a.root);save=a.root/'independent-checks.json';save.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
