#!/usr/bin/env python3
"""Replay complete metric relative changes and check them with independent Fractions."""
import argparse
import copy
import importlib.util
import json
import subprocess
from fractions import Fraction as F
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / 'declared-responses/replay.py'
_spec = importlib.util.spec_from_file_location('novelty_sources', _path)
sources = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sources)
save, run, identity, matrix = sources.save, sources.run, sources.identity, sources.matrix

_span_path = Path(__file__).resolve().parents[1] / 'source-admission/replay.py'
_span_spec = importlib.util.spec_from_file_location('novelty_evidence', _span_path)
shared_evidence = importlib.util.module_from_spec(_span_spec)
_span_spec.loader.exec_module(shared_evidence)


def shared_spans(path):
    declaration=json.loads(path.read_text())
    documents={x['id']:(path.parent/x['path']).read_bytes() for x in declaration['documents'] if x['format']=='utf8'}
    for span in shared_evidence.evidence_records(declaration):
        original=documents[span['document']][span['start']:span['end']]
        if 'text' in span:
            assert span['text'].encode()==original
            span.pop('text')
    declaration['schema']='gadget.scientific.study.v2'
    save(path,declaration)



def declaration(name, sequence, values, metric, channels=(), source=None, grade=0):
    fields = [sources.field(name, 100, 10, values), sources.field('zero', 101, 10, [0]*len(values)),
              sources.field('partial', 102, 10, [None]*len(values))]
    for i, (label, value) in enumerate(channels):
        component=sources.field(label, 200+i, 10, value)
        component['quantity']='authored_input_change'
        fields.append(component)
    category = sources.field('category', 999, 10, [1]*len(values))
    category.update(value_type='categorical', codebook={'1': 'declared category'})
    fields.append(category)
    return dict(id='novelty-'+name, namespace='9401', sequence=str(sequence),
                source=source or sources.source(len(values)), carriers=[dict(key='10',grade=grade,size=len(values))],
                operators=[sources.op('M',1,'metric',10,10,metric)], fields=fields)


def connect(binary, root, comparison, old_field='old', new_field='new', causes=None,
            old_carrier='10', new_carrier='10', metric='M'):
    d = json.loads((root/'p1/study.json').read_text())
    d['id'] = 'novelty-coordinator'; d['native']['path'] = 'p1/'+d['native']['path']
    for document in d['documents']:
        document['path'] = 'p1/'+document['path']
    evidence = copy.deepcopy(d['operators'][-1]['evidence'])
    span=evidence[0]
    document=next(x for x in d['documents'] if x['id']==span['document'])
    raw=(root/document['path']).read_bytes()
    if b'"comparison"' in raw:
        span['start']=raw.index(b'"comparison"')
        span['end']=len(raw.rstrip())-1
    if d['schema']=='gadget.scientific.study.v2':
        span.pop('text',None)
    else:
        span['text']=raw[span['start']:span['end']].decode()
    d['publications'] = []
    for name in ['p0','p1']:
        path = root/name/'study.json'
        digest = run(binary,'query','--study',path,'--query',save(root/'inspect.json',dict(operation='inspect')))['source_digest']
        d['publications'].append(dict(id=name,study=dict(path=name+'/study.json',sha256=sources.digest(path)),
            source_digest=digest,meaning='Actual admitted source publication '+name,evidence=evidence))
    save(root/'study.json', d)
    old = json.loads((root/'p0/study.json').read_text()); new = json.loads((root/'p1/study.json').read_text())
    query = dict(operation='metric_relative_novelty',old=dict(publication='p0',source_sequence=old['native']['cutoff'],field=old_field),
        new=dict(publication='p1',source_sequence=new['native']['cutoff'],field=new_field,metric=metric),
        comparison=dict(old_carrier=old_carrier,new_carrier=new_carrier,map=comparison,
            meaning='Complete explicitly declared correspondence between original source coordinates',evidence=evidence))
    if causes is not None:
        query['causes'] = causes
    save(root/'query.json', query)
    receipt = run(binary,'query','--study',root/'study.json','--query',root/'query.json')
    save(root/'result.json', receipt)
    return receipt


def execute(binary, producer, root, a, b, comparison, causes=None):
    root.mkdir(parents=True,exist_ok=False)
    # The complete correspondence is an authored declaration, not fitted from data.
    b = copy.deepcopy(b); b['comparison'] = comparison
    for name, spec in [('p0',a),('p1',b)]:
        path = save(root/(name+'-spec.json'),spec)
        subprocess.run([str(producer),str(path),str(root/name)],check=True,capture_output=True,text=True)
        shared_spans(root/name/'study.json')
    receipt = connect(binary,root,comparison,causes=causes)
    return verify(root,receipt)


def dense(m):
    out = [[F(0) for _ in range(m['columns'])] for _ in range(m['rows'])]
    for e in m['entries']:
        out[e['row']][e['column']] = F(e['value'])
    return out


def mv(a,x): return [sum((v*w for v,w in zip(row,x) if v and w),F(0)) for row in a]
def sub(a,b): return [x-y for x,y in zip(a,b)]
def dot(a,b): return sum((x*y for x,y in zip(a,b) if x and y),F(0))


def solve(a,b):
    """Separate exact Gauss Jordan elimination; verifier only, never production."""
    rows = [list(row)+[value] for row,value in zip(a,b)]
    for k in range(len(b)):
        pivot = next(i for i in range(k,len(b)) if rows[i][k])
        rows[k],rows[pivot] = rows[pivot],rows[k]
        value = rows[k][k]; rows[k] = [x/value for x in rows[k]]
        for i in range(len(b)):
            if i != k:
                value = rows[i][k]
                if value:
                    rows[i] = [x-value*y for x,y in zip(rows[i],rows[k])]
    return [row[-1] for row in rows]


def values(binary, study, field, root):
    cat = run(binary,'query','--study',study,'--query',save(root/'source-inspect.json',dict(operation='inspect')))
    digest,sequence = cat['source_digest'], cat['source_sequence']
    source = run(binary,'source-read','--study',study,'--request',save(root/'source-read.json',
        dict(schema='gadget.scientific.source.read.v1',source_digest=digest,publication='self',source_sequence=sequence,
             stream=dict(kind='field',id=field),offset=0,limit=max(1,next(x for x in cat['result']['fields'] if x['id']==field)['coordinates']))))
    return [F(x['value']) for x in source['records']]


def verify(root, receipt):
    q = receipt['query']; r = receipt['result']; j = dense(q['comparison']['map']); n = q['comparison']['map']['columns']; m = len(j)
    if (root/'p0-spec.json').exists():
        a = json.loads((root/'p0-spec.json').read_text()); b = json.loads((root/'p1-spec.json').read_text())
        fields = {x['id']:[F(v) for v in x['values']] for x in b['fields'] if all(v is not None for v in x['values'])}
        x = [F(v) for v in next(v for v in a['fields'] if v['id']==q['old']['field'])['values']]
        xp = fields[q['new']['field']]; metric = dense(b['operators'][0]['map'])
    else:
        data = json.loads((root/'original-fields.json').read_text()); x = list(map(F,data['old'])); xp = list(map(F,data['new'])); metric = dense(data['metric']); fields = {}
    dual = [[sum((j[t][i]*metric[t][k] for t in range(m) if j[t][i] and metric[t][k]),F(0)) for k in range(m)] for i in range(n)]
    gram = [[sum((dual[i][t]*j[t][k] for t in range(m) if dual[i][t] and j[t][k]),F(0)) for k in range(n)] for i in range(n)]
    lifted = solve(gram,mv(dual,xp)); transported = mv(j,x); innovation = sub(xp,transported)
    coordinates = sub(lifted,x); persistent = mv(j,coordinates); complement = sub(xp,mv(j,lifted))
    expected = dict(transported_old=transported, innovation=innovation,persistent_coordinates=coordinates,
        persistent=persistent,complement=complement,reconstructed=[a+b for a,b in zip(persistent,complement)])
    checks = 0
    for name, vector in expected.items():
        assert list(map(F,r['split'][name])) == vector,(name,r['split'][name],vector); checks += len(vector)
    assert dot(persistent,mv(metric,complement)) == 0
    split_gram = [[dot(a,mv(metric,b)) for b in [persistent,complement]] for a in [persistent,complement]]
    assert dense(r['split']['pairings']) == split_gram
    quadrance = dot(innovation,mv(metric,innovation)); assert F(r['split']['quadrance']) == quadrance; checks += 6
    assert len(r['publications']['old']['source_coordinates']) == n and len(r['publications']['new']['source_coordinates']) == m
    if 'causes' in q:
        c = r['causes']; channels = [fields[name] for name in q['causes']]
        restored = [list(map(F,c['fields'][slot['dictionary_index']]['values'])) for slot in c['word']]
        assert restored == channels and [v['field'] for v in c['word']] == q['causes']
        aggregate = [sum((v[i] for v in channels),F(0)) for i in range(m)]
        measured = [mv(metric,b) for b in channels]
        cause_gram = [[dot(a,b) for b in measured] for a in channels]
        assert dense(c['pairings']) == cause_gram
        assert aggregate == innovation == list(map(F,c['aggregate']))
        assert sum((sum(row,F(0)) for row in cause_gram),F(0)) == quadrance == F(c['quadrance'])
        assert len(c['fields']) == len({next(f['key'] for f in b['fields'] if f['id']==name) for name in q['causes']})
        checks += len(channels)**2+len(channels)*m+2
    else:
        assert 'causes' not in r
    out = dict(status='PASS',old_coordinates=n,new_coordinates=m,grade=r['native_execution']['new']['grade'],
               independent_fraction_values=checks,quadrance=str(quadrance),cause_slots=len(q.get('causes',[])))
    save(root/'verification.json',out);return out


def examples(binary, producer, root):
    root.mkdir(parents=True,exist_ok=False); records={}
    j=matrix([[1,0],[0,1],[1,1]]); m=matrix([[2,0,0],[0,3,0],[0,0,4]])
    old=declaration('old',7,[2,3],identity(2)); new=declaration('new',8,[3,2,8],m,
        [('part',[1,0,0]),('negative',[-2,0,0]),('whole',[1,-1,3])])
    records['algebra']=execute(binary,producer,root/'algebra',old,new,j,['part','part','negative','zero','whole'])
    for name,J,M,x,xp in [
        ('coupled',matrix([[1,2],[0,1],[2,-1]]),matrix([[3,1,0],[1,4,1],[0,1,2]]),[F(1,2),-2],[3,0,-1]),
        ('empty-old',dict(rows=3,columns=0,entries=[]),m,[],[3,2,8]),
        ('empty-both',dict(rows=0,columns=0,entries=[]),identity(0),[],[]),
        ('zero-innovation',identity(2),matrix([[2,1],[1,3]]),[2,3],[2,3])]:
        records[name]=execute(binary,producer,root/name,declaration('old',7,x,identity(len(x))),
            declaration('new',8,xp,M),J,[] if name in ['empty-both','zero-innovation'] else None)
    a=declaration('old',7,[1,0],identity(2));b=declaration('new',8,[0,1],identity(2),[('delta',[-1,1])])
    a['fields'][0]['value_type']=b['fields'][0]['value_type']='indicator'
    records['indicator-change']=execute(binary,producer,root/'indicator-change',a,b,identity(2),['delta','zero'])
    a,b,*_=sources.square()
    old=declaration('old',7,[2,2,0,0],identity(4),source=a['source'],grade=1)
    new=declaration('new',8,[0,0,-3,-3],identity(4),[('representation',[-2,-2,-2,-2]),('amplitude',[0,0,-1,-1])],source=b['source'],grade=1)
    records['square']=execute(binary,producer,root/'square',old,new,identity(4),['representation','amplitude'])
    records['grade-two']=execute(binary,producer,root/'grade-two',declaration('old',7,[2],identity(1),source=b['source'],grade=2),
        declaration('new',8,[3],identity(1),source=b['source'],grade=2),identity(1))
    # A sparse source larger than earlier dimensional ceilings remains complete.
    dimension=257
    records['sparse-257']=execute(binary,producer,root/'sparse-257',declaration('old',7,[1]+[0]*(dimension-1),identity(dimension)),
        declaration('new',8,[2]+[0]*(dimension-1),identity(dimension),[('whole',[1]+[0]*(dimension-1))]),identity(dimension),['zero']*127+['whole'])
    save(root/'verification.json',dict(status='PASS',cases=records));return records


def clinical(binary,producer,original,root):
    root.mkdir(parents=True,exist_ok=False); sizes=[]; carriers=[]; fields=[]
    for p in ['p0','p1']:
        study=original/p/'study.json';d=json.loads(study.read_text())
        cat=run(binary,'query','--study',study,'--query',save(root/'inspect.json',dict(operation='inspect')))['result']
        field=next(x for x in cat['fields'] if x['id']=='all.time');n=field['coordinates'];carrier=field['carrier']
        spec=dict(id='original-clinical-'+p,namespace=d['native']['source_namespace'],operators=[sources.op('M',940001,'metric',carrier,carrier,identity(n),method='Declared unit metric on original follow-up-day coordinates')],fields=[])
        if p=='p1':
            spec['comparison']=dict(rows=n,columns=sizes[0],entries=[dict(row=i,column=i,value='1') for i in range(sizes[0])])
        subprocess.run([str(producer),str(save(root/(p+'-extension.json'),spec)),str(root/p),str(study)],check=True,capture_output=True,text=True)
        shared_spans(root/p/'study.json')
        sizes.append(n);carriers.append(carrier);fields.append(values(binary,root/p/'study.json','all.time',root))
    n,m=sizes;comparison=dict(rows=m,columns=n,entries=[dict(row=i,column=i,value='1') for i in range(n)])
    save(root/'original-fields.json',dict(old=list(map(str,fields[0])),new=list(map(str,fields[1])),metric=identity(m)))
    receipt=connect(binary,root,comparison,old_field='all.time',new_field='all.time',old_carrier=carriers[0],new_carrier=carriers[1])
    return verify(root,receipt)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--binary',required=True,type=Path);p.add_argument('--producer',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--clinical',type=Path);args=p.parse_args()
    result=clinical(args.binary.resolve(),args.producer.resolve(),args.clinical.resolve(),args.out) if args.clinical else examples(args.binary.resolve(),args.producer.resolve(),args.out)
    print(json.dumps(result,indent=2))
