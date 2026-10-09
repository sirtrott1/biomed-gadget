#!/usr/bin/env python3
"""Publish declared sources and execute complete native Hodge decompositions."""
import argparse, copy, importlib.util, json, shutil, subprocess
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / 'declared-responses/replay.py'
_spec = importlib.util.spec_from_file_location('response_sources', _path)
sources = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sources)
save, run, matrix, identity = sources.save, sources.run, sources.matrix, sources.identity

def specification(name, source, grade, positive=False):
    sizes = [len(source['vertices'])]
    if source['relations']:
        sizes.append(len(source['relations']))
    for r in source.get('higher', []):
        while len(sizes) <= r['grade']:
            sizes.append(0)
        sizes[r['grade']] += 1
    carriers = [dict(key=str(10+k), grade=k, size=n) for k, n in enumerate(sizes)]
    operators = []
    for k, n in enumerate(sizes):
        metric = identity(n)
        if positive:
            metric = matrix([[str(2+i) if i == j else ('1/5' if abs(i-j) == 1 else '0')
                              for j in range(n)] for i in range(n)]) if n else identity(0)
        operators.append(sources.op('M'+str(k), 200+k, 'metric', 10+k, 10+k, metric))
    n = sizes[grade]
    fields = [sources.field('x', 100, 10+grade, [str(i+1) for i in range(n)]),
              sources.field('x2', 101, 10+grade, [str((-1)**i*(i+2)) for i in range(n)]),
              sources.field('zero', 102, 10+grade, ['0']*n),
              sources.field('partial', 103, 10+grade, [None]*n)]
    category = sources.field('category', 104, 10+grade, ['1']*n)
    category['value_type'] = 'categorical'; category['codebook'] = {'1': 'declared_category'}
    fields.append(category)
    return dict(id=name, namespace='9201', sequence='7', source=source,
                carriers=carriers, operators=operators, fields=fields)

def authored(binary, producer, root, spec, grade, positive=False):
    root.mkdir(parents=True, exist_ok=False)
    save(root/'spec.json', spec)
    subprocess.run([str(producer), str(root/'spec.json'), str(root/'study')], check=True,
                   capture_output=True, text=True)
    query = dict(operation='hodge_frame_decomposition', publication='self', source_sequence='7',
                 grade=grade, grade_carriers=[c['key'] for c in spec['carriers']],
                 fields=['x', 'x2', 'zero'], metrics=dict(kind='positive', tensors=[o['id'] for o in spec['operators']]) if positive else dict(kind='unit'))
    save(root/'query.json', query)
    certificate = run(binary, 'query', '--study', root/'study/study.json', '--query', root/'query.json')
    save(root/'result.json', certificate)
    workflow = dict(schema='gadget.scientific.workflow.v1', id='native-Hodge-fields',
        source_digest=certificate['source_digest'], limits=dict(depth=2, nodes=3, model_calls=0, parallel_workers=2),
        definitions=[], entry=dict(id='joined', type='parallel', workers=2,
            steps=[dict(id='fields', type='query', query=query),
                   dict(id='repeat', type='query', query=query)]))
    save(root/'workflow-request.json', workflow)
    receipt = run(binary, 'workflow', '--study', root/'study/study.json', '--request', root/'workflow-request.json')
    assert receipt['status'] == 'complete', receipt
    save(root/'workflow.json', receipt)
    return query

def examples(binary, producer, out, clinical=None):
    old, new, *_ = sources.square()
    authored(binary, producer, out/'unfilled-square', specification('unfilled-square', old['source'], 1), 1)
    authored(binary, producer, out/'filled-square', specification('filled-square', new['source'], 1), 1)
    authored(binary, producer, out/'positive-metric', specification('positive-metric', new['source'], 1, True), 1, True)
    branch = dict(vertices=['1','2','3','4'], relations=[dict(key=str(101+i), occurrences=[
        dict(key=str(1001+4*i+j), cell=str(j+1), exists=True, negative=j==0,
             magnitude='1' if j==0 else '1/3') for j in range(4)]) for i in range(4)])
    authored(binary, producer, out/'branching', specification('branching', branch, 1, True), 1, True)
    higher = copy.deepcopy(new['source'])
    higher['higher'] += [dict(grade=2, key='202', occurrences=[dict(key=str(2101+i), cell=str(101+i), exists=True, negative=False, magnitude='1') for i in range(4)]),
        dict(grade=3, key='301', occurrences=[dict(key='3001', cell='201', exists=True, negative=True, magnitude='1'),dict(key='3002', cell='202', exists=True, negative=False, magnitude='1')])]
    authored(binary, producer, out/'grade-three', specification('grade-three', higher, 2, True), 2, True)
    authored(binary, producer, out/'beyond-native-default', specification('65-original-coordinates', sources.source(65), 0), 0)
    authored(binary, producer, out/'empty', specification('empty', sources.source(0), 0), 0)
    if clinical:
        base = out/'recorded-followup'; base.mkdir(parents=True, exist_ok=False)
        original=clinical/'p0/study.json'
        d = json.loads(original.read_text())
        inspection=run(binary,'query','--study',original,
                       '--query',save(base/'inspect-query.json',dict(operation='inspect')))
        admission=json.loads((original.parent/'admission.json').read_text())
        native_carriers=sorted((c for c in admission.get('carriers',[])
                                if c['kind']=='source_grade'),key=lambda c:c['grade'])
        if [c['grade'] for c in native_carriers]==[0,1]:
            # Current table imports already own the complete declared family.
            shutil.copytree(original.parent,base/'study')
            carrier_keys=[c['key'] for c in native_carriers]
        else:
            # A retained earlier publication can lack its relation carrier.
            # Declare it explicitly in a new example owner; preserve the archive.
            native_spec=dict(id='recorded-followup-complete-carriers',namespace=d['native']['source_namespace'],
                sequence=inspection['source_sequence'],carriers=[dict(key='12',grade=1,
                size=sum(c['grade']==1 for c in d['cells']))],operators=[],fields=[])
            save(base/'carrier-declaration.json',native_spec)
            subprocess.run([str(producer),str(base/'carrier-declaration.json'),str(base/'study'),str(original)],
                           capture_output=True,text=True,check=True)
            carrier_keys=['1','12']
        # The source publication retains the original durations, censoring and
        # category-membership relations. These are recorded fields, not event estimates.
        query = dict(operation='hodge_frame_decomposition', publication='self',
            source_sequence=inspection['source_sequence'], grade=0, grade_carriers=carrier_keys,
            fields=['all.time','all.age'], metrics=dict(kind='unit'))
        save(base/'query.json', query)
        save(base/'result.json', run(binary,'query','--study',base/'study/study.json','--query',base/'query.json'))
    return out

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--producer',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--clinical',type=Path); a=p.parse_args()
    examples(a.binary.resolve(), a.producer.resolve(), a.out.resolve(), a.clinical)
    from verify import verify
    records=[verify(x) for x in sorted(a.out.iterdir()) if (x/'result.json').is_file()]
    save(a.out/'verification.json', dict(status='PASS', cases=records,
                                      complete_coordinates=sum(r['coordinates'] for r in records)))
    print(json.dumps(dict(status='PASS',cases=len(records))))
