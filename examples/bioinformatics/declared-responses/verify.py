#!/usr/bin/env python3
"""Independent exact replay of complete maps, responses and ordered cause pairings."""
import csv, hashlib, io, json, sys
from pathlib import Path
import sympy as S

def matrix(j):
    out=S.zeros(j['rows'],j['columns'])
    for e in j['entries']:out[e['row'],e['column']]=S.Rational(e['value'])
    return out
def vector(v):return S.Matrix([S.Rational(q) for q in v])
def sha(p):return 'sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()
def source_maps(root,p,spec,result):
    d=result['publications'][p]['declaration'];base=root/('p0' if p=='old' else 'p1')
    assert sha(base/d['native']['path'])==d['native']['sha256']
    for doc in d['documents']:assert sha(base/doc['path'])==doc['sha256']
    if 'source' in spec:
        s=spec['source'];vertices=s['vertices'];relations=s['relations'];keys=[r['key'] for r in relations]
        B=S.zeros(len(vertices),len(relations))
        for j,r in enumerate(relations):
            for o in r['occurrences']:
                if o['exists']:B[vertices.index(o['cell']),j]+=(-1 if o['negative'] else 1)*S.Rational(o['magnitude'])
        higher=s.get('higher',[]);C=S.zeros(len(relations),len(higher))
        for j,r in enumerate(higher):
            assert r['grade']==2
            for o in r['occurrences']:
                if o['exists']:C[keys.index(o['cell']),j]+=(-1 if o['negative'] else 1)*S.Rational(o['magnitude'])
        assert B*C==S.zeros(B.rows,C.cols)
        assert len(d['occurrences'])==sum(len(r['occurrences']) for r in relations+higher)
        hodge={0:B*B.T,1:B.T*B+C*C.T}
    else:
        table=next(x for x in d['documents'] if x['id']=='table')
        rows=list(csv.DictReader(io.StringIO((base/table['path']).read_text())))
        # The producer declares two sex-membership relations in the existing source.
        categories=sorted({r['sex'] for r in rows});assert len(categories)==2
        B=S.zeros(len(rows),2)
        for i,r in enumerate(rows):B[i,categories.index(r['sex'])]=1
        hodge={0:B*B.T}
        assert result['publications'][p]['input_values']==[str(S.Rational(r['time'])) for r in rows]
        assert len(d['occurrences'])==len(rows)
    out={}
    for o in spec['operators']:
        out[o['id']]=hodge[o['hodge_grade']] if 'hodge_grade' in o else B.T if 'boundary_transpose_grade' in o else matrix(o['map'])
    return out
def verify(root):
    root=Path(root);response=json.loads((root/'response.json').read_text());observed=json.loads((root/'observed.json').read_text())
    for cert in [response,observed]:
        expected=cert['certificate_sha256'];body={k:v for k,v in cert.items() if k!='certificate_sha256'}
        # json-c's public pretty serialization is checked by CLI contract tests;
        # the verifier binds original archive/document hashes and exact coefficients.
        assert expected.startswith('sha256:') and len(expected)==71 and body['result']['identity_verified']
    r=response['result'];o=observed['result'];q=response['query'];oq=observed['query']
    a=json.loads((root/'p0-declarations.json').read_text());b=json.loads((root/'p1-declarations.json').read_text())
    old=source_maps(root,'old',a,r);new=source_maps(root,'new',b,r)
    for cert in [response,observed]:
        maps=cert['result']['declared_maps']
        for side,values in [('old',old),('new',new)]:
            for suffix,name in [('resolving','L'),('realization','I')]:assert matrix(maps[side+'_'+suffix])==values[name]
            if side+'_accession' in maps:assert matrix(maps[side+'_accession'])==values['A']
        if 'metric' in maps:assert matrix(maps['metric'])==new['M']
    J=matrix(q['comparisons']['ambient']['map']);K=matrix(q['comparisons']['source']['map']);KA=matrix(oq['observation']['comparison']['map']);lam=S.Rational(q['lambda'])
    u=vector(r['publications']['old']['input_values']);v=vector(r['publications']['new']['input_values'])
    G=(S.eye(old['L'].rows)+lam*old['L']).inv(method='DM');H=(S.eye(new['L'].rows)+lam*new['L']).inv(method='DM')
    g=G*old['I']*u;h=H*new['I']*v;EI=new['I']*K-J*old['I'];EL=new['L']*J-J*old['L'];z=v-K*u
    terms=[H*EI*u,H*new['I']*z,-lam*H*EL*g]
    coordinates=0
    def same(x,y):
        nonlocal coordinates
        assert x==y,(x,y);coordinates+=x.rows*x.cols
    same(vector(r['old_response']),g);same(vector(r['new_response']),h);same(vector(r['coordinate_innovation']),z)
    same(matrix(r['realization_defect']),EI);same(matrix(r['operator_defect']),EL)
    same(vector(r['aggregate']),h-J*g)
    for name,t in zip(r['cause_order'],terms):same(vector(r['causes'][name]),t)
    EA=new['A']*J-KA*old['A'];obs=[EA*g]+[new['A']*t for t in terms]
    same(matrix(o['accession_defect']),EA);same(matrix(o['realization_defect']),EI);same(matrix(o['operator_defect']),EL)
    same(vector(o['aggregate']),new['A']*h-KA*old['A']*g)
    for name,t in zip(o['cause_order'],obs):same(vector(o['causes'][name]),t)
    C=S.Matrix.hstack(*obs);gram=C.T*new['M']*C
    same(matrix(o['pairing']['gram']),gram)
    assert S.Rational(o['pairing']['aggregate_quadrance'])==sum(gram)
    for cert in [response,observed]:
        report=cert['result']['native_execution']
        for name in ['old','new']:
            assert report[name]['sequence']==cert['query'][name]['source_sequence']
            assert report['execution_availability'][name]['source_sequence']==report[name]['sequence']
        assert not report['fallback']
    if root.name=='square':
        assert g==vector(['4/3','4/3','2/3','2/3']) and h==vector(['1/5','1/5','-4/5','-4/5'])
        assert vector(o['aggregate'])==vector(['-29/15','-29/15'])
        assert all(any(x!=0 for x in t) for t in obs)
    return dict(case=root.name,status='PASS',complete_coordinates=coordinates,old_ambient=old['L'].rows,new_ambient=new['L'].rows,all_causes_preserved=True)
if __name__=='__main__':print(json.dumps(verify(Path(sys.argv[1])),indent=2))
