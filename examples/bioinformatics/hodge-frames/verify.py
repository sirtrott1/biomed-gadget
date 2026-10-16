#!/usr/bin/env python3
"""Independent source, frame and full field checks using rational matrices."""
import csv, ctypes, ctypes.util, hashlib, json, sys
from pathlib import Path
import sympy as S

def matrix(d):
    result=S.zeros(int(d['rows']),int(d['columns']))
    for e in d['entries']: result[int(e['row']),int(e['column'])]=S.Rational(e['value'])
    return result
def vector(values):return S.Matrix(len(values),1,[S.Rational(x) for x in values])
def digest(p):return 'sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()
def certificate_digest(document):
    # Certificates use json-c's existing pretty byte contract, including empty
    # arrays and objects. Reuse that serializer solely for byte/hash checking.
    codec=ctypes.CDLL(ctypes.util.find_library('json-c'))
    codec.json_tokener_parse.argtypes=[ctypes.c_char_p];codec.json_tokener_parse.restype=ctypes.c_void_p
    codec.json_object_to_json_string_ext.argtypes=[ctypes.c_void_p,ctypes.c_int]
    codec.json_object_to_json_string_ext.restype=ctypes.c_char_p
    codec.json_object_put.argtypes=[ctypes.c_void_p]
    owner=codec.json_tokener_parse(json.dumps(document,ensure_ascii=False).encode())
    assert owner
    try:return 'sha256:'+hashlib.sha256(codec.json_object_to_json_string_ext(owner,2|16)).hexdigest()
    finally:codec.json_object_put(owner)
def verify(root):
    c=json.loads((root/'result.json').read_text());r=c['result'];d=r['publication']['declaration']
    assert digest(root/'study'/d['native']['path'])==d['native']['sha256']
    for doc in d['documents']:assert digest(root/'study'/doc['path'])==doc['sha256']
    # The certificate seal binds all fields, frames, provenance and declarations.
    seal=c.pop('certificate_sha256');assert certificate_digest(c)==seal
    sizes=[len(x) for x in r['publication']['grade_coordinates']];boundaries=[S.zeros(0,sizes[0])]
    spec=json.loads((root/'spec.json').read_text()) if (root/'spec.json').exists() else None
    if spec:
        source=spec['source'];keys=[source['vertices'],[x['key'] for x in source['relations']]]
        for k in range(2,len(sizes)):keys.append([x['key'] for x in source.get('higher',[]) if x['grade']==k])
        for k in range(1,len(sizes)):
            cells=source['relations'] if k==1 else [x for x in source['higher'] if x['grade']==k]
            B=S.zeros(sizes[k-1],sizes[k])
            for j,cell in enumerate(cells):
                for o in cell['occurrences']:
                    if o['exists']:B[keys[k-1].index(o['cell']),j]+=(-1 if o['negative'] else 1)*S.Rational(o['magnitude'])
            boundaries.append(B)
        assert len(d['occurrences'])==sum(len(x['occurrences']) for x in source['relations']+source.get('higher',[]))
        metrics=[S.eye(n) for n in sizes] if r['metrics']['kind']=='unit' else [matrix(o['map']) for o in spec['operators']]
        source_fields={x['id']:x['values'] for x in spec['fields']}
    else:
        doc=next(x for x in d['documents'] if x['id']=='table')
        rows=list(csv.DictReader((root/'study'/doc['path']).open()))
        categories=sorted({row['sex'] for row in rows});B=S.zeros(len(rows),len(categories))
        for i,row in enumerate(rows):B[i,categories.index(row['sex'])]=1
        boundaries.append(B);metrics=[S.eye(n) for n in sizes]
        source_fields={'all.time':[row['time'] for row in rows],'all.age':[row['age'] for row in rows]}
    for k in range(2,len(sizes)):assert boundaries[k-1]*boundaries[k]==S.zeros(sizes[k-2],sizes[k])
    k=r['grade'];n=sizes[k];B=boundaries[k];C=boundaries[k+1] if k+1<len(sizes) else S.zeros(n,0);M=metrics[k]
    originals={'gradient':B.T,'curl':C,'harmonic':B.col_join(C.T*M)}
    frames={name:matrix(value) for name,value in r['frames'].items()};coordinates=0
    def same(a,b):
        nonlocal coordinates
        assert a.shape==b.shape and a==b,(root.name,a.shape,b.shape)
        coordinates+=a.rows*a.cols
    for name,original in originals.items():
        reduction=r['reductions'][name];A=matrix(reduction['original']);F=matrix(reduction['image']);K=matrix(reduction['kernel']);R=matrix(reduction['coefficients'])
        same(A,original);same(F*R,A);same(A*K,S.zeros(A.rows,K.cols))
        pivots=list(A.rref()[1]);assert reduction['pivot_columns']==pivots
        assert reduction['free_columns']==[j for j in range(A.cols) if j not in pivots]
        same(F,A[:,pivots]);assert K.cols==A.cols-len(pivots)
        assert K.rank()==K.cols and F.rank()==F.cols
        if name=='gradient':same(M*frames[name],F)
        elif name=='curl':same(frames[name],F)
        else:same(frames[name],K)
    assert sum(F.cols for F in frames.values())==n
    for name,F in frames.items():
        for other,G in frames.items():
            if other!=name:same(F.T*M*G,S.zeros(F.cols,G.cols))
    projectors={name:F*(F.T*M*F).inv()*F.T*M if F.cols else S.zeros(n,n) for name,F in frames.items()}
    for f in r['fields']:
        x=vector(f['input']);same(x,vector(source_fields[f['id']]))
        parts=[]
        for name,F in frames.items():
            value=vector(f['sectors'][name]);a=vector(f['frame_coordinates'][name])
            same(value,projectors[name]*x);same(F*a,value);parts.append(value)
        same(sum(parts,S.zeros(n,1)),x);same(vector(f['reconstructed']),x)
        report=f['native_execution'];assert report['source_sequence']==c['query']['source_sequence']
        assert report['carrier']==c['query']['grade_carriers'][k] and not report['fallback']
    if root.name=='unfilled-square':assert frames['harmonic'].cols==1 and frames['curl'].cols==0
    if root.name=='filled-square':assert frames['harmonic'].cols==0 and frames['curl'].cols==1
    return dict(case=root.name,status='PASS',coordinates=coordinates,grade=k,source_size=n,
                dimensions={name:F.cols for name,F in frames.items()})
if __name__=='__main__':print(json.dumps(verify(Path(sys.argv[1])),indent=2))
