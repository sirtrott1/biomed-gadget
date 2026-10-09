"""Independent exact large-carrier and complete-certificate regression."""
import hashlib
import json
import subprocess
import time
from fractions import Fraction as Q
from pathlib import Path


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def load(path):
    return json.loads(path.read_text())


def invoke(binary, args, output=None):
    start = time.perf_counter()
    p = subprocess.run([str(binary), *map(str, args)], text=True,
                       capture_output=True, timeout=420)
    assert p.returncode == 0, (args, p.stderr[-4000:])
    if output:
        output.write_text(p.stdout)
        assert output.stat().st_size > 16*1024*1024, 'Fixture did not exercise the former output cap'
    return json.loads(p.stdout), time.perf_counter() - start


def digest(path):
    return 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def control_case(binary, root, workers):
    """Independent 2x2 rational oracle, large carrier and redundant controls."""
    root.mkdir(parents=True)
    n, tested, controls, selections = 4209, 129, 65, 33
    header = ['patient', 'group', 'endpoint'] + [f'x{k}' for k in range(tested)]
    header += [f'control{k}' for k in range(controls)]
    header += [f'selection{k}' for k in range(selections)]
    def measurement(i, k):
        return 17 + 3 * (i % 23) + (k + 1) * (i % 2) + (i % 5 - 2) * (k % 3)
    with (root/'table.tsv').open('w') as stream:
        stream.write('\t'.join(header) + '\n')
        for i in range(n):
            values = [f'person-{i}', 'declared', str(i % 2)]
            values += [str(measurement(i, k)) for k in range(tested)]
            values += [str(i % 23)] * controls
            values += [str(int(i < n-4))] + ['1'] * (selections-1)
            stream.write('\t'.join(values) + '\n')
    dictionary = 'Authored nonclinical counts, indicators and explicit identities.\n'
    (root/'dictionary.txt').write_text(dictionary)
    spec = dict(schema='gadget.table.import.v1', id='large-control-family',
        source_namespace='903100', source=dict(path='table.tsv', sha256=digest(root/'table.tsv'),
        delimiter='\t', header=True), dictionary=dict(path='dictionary.txt',
        sha256=digest(root/'dictionary.txt')), columns=[dict(id=name, channel='observation',
        quantity=name, units='dimensionless', value_type='categorical' if name in ('patient','group')
        else 'indicator' if name == 'endpoint' or name.startswith('selection') else 'quantitative',
        dictionary_start=0, dictionary_end=len(dictionary.encode())) for name in header],
        identity_column='patient', missing_tokens=[], relation_columns=['group'],
        question='Full controlled fields on an explicitly measured large carrier',
        population=f'{n} authored records', selection='declared indicator intersection',
        limitations=['synthetic fixture'], source_url='synthetic://large-control-family',
        attribution='GADGET qualification', license='synthetic')
    save(root/'import.json', spec)
    _, imported = invoke(binary, ['import-table','--spec',root/'import.json','--out',root/'study'])
    query = dict(operation='controlled_pairings', fields=[f'all.x{k}' for k in range(tested)],
        endpoint='all.endpoint', controls=[f'all.control{k}' for k in range(controls)],
        selection=[f'all.selection{k}' for k in range(selections)], support='joint_reported',
        metric='identity', intercept=True, workers=workers, gram_fields=['all.x2','all.x0'])
    save(root/'query.json', query)
    receipt, seconds = invoke(binary, ['query','--study',root/'study/study.json',
        '--query',root/'query.json'], root/'certificate.json')
    result = receipt['result']; count = n-4
    assert result['n'] == count and result['excluded_rows'] == 4
    assert result['control_columns'] == controls+1 and result['control_rank'] == 2
    assert result['entities'] == [f'person-{i}' for i in range(count)]
    assert result['coordinates'] == list(range(count))
    t = [Q(i % 23) for i in range(count)]; y = [Q(i % 2) for i in range(count)]
    s, ss, det = sum(t), dot(t,t), count * dot(t,t) - sum(t)**2
    def represented(x):
        total, tx = sum(x), dot(t,x)
        a, b = (ss * total - s * tx) / det, (count * tx - s * total) / det
        return [a+b*v for v in t]
    py = represented(y); ry = [a-b for a,b in zip(y,py)]
    assert list(map(Q,result['endpoint_input'])) == y
    assert list(map(Q,result['endpoint_represented'])) == py
    assert list(map(Q,result['endpoint_residual'])) == ry
    assert Q(result['endpoint_residual_quadrance']) == dot(ry,ry)
    comparisons = 3*count + 1; selected_residuals = {}
    for k, field in enumerate(result['fields']):
        x = [Q(measurement(i,k)) for i in range(count)]
        p = represented(x); r = [a-b for a,b in zip(x,p)]
        assert field['field'] == f'all.x{k}'
        assert list(map(Q,field['input'])) == x
        assert list(map(Q,field['represented'])) == p
        assert list(map(Q,field['residual'])) == r
        assert Q(field['pairing']) == dot(r,ry)
        assert Q(field['residual_quadrance']) == dot(r,r)
        assert Q(field['endpoint_coefficient']) == dot(r,ry) / dot(ry,ry)
        assert sum(r) == dot(t,r) == dot(p,r) == 0
        comparisons += 3*count+6
        if k in (2,0):
            selected_residuals[k] = r
    expected = [[dot(selected_residuals[a],selected_residuals[b]) for b in (2,0)] for a in (2,0)]
    assert result['gram_fields'] == ['all.x2','all.x0']
    assert [list(map(Q,row)) for row in result['residual_gram']['values']] == expected
    comparisons += 4
    execution = result['execution']
    assert execution['owner'] == 'rfc::numeric::PreparedExactSolve'
    assert execution['gram_preparations'] == 1 and execution['applications'] == tested+1
    assert execution['completed_tiles'] == tested and execution['worker_slots'] == workers
    return dict(status='PASS', source_records=n, selected_records=count, tested_fields=tested,
        declared_controls=controls, selection_indicators=selections, union_fields=tested+controls+selections+1,
        control_rank=2, comparisons=comparisons, input_sha256=digest(root/'table.tsv'),
        source_digest=receipt['source_digest'], certificate_file_sha256=digest(root/'certificate.json'),
        certificate_bytes=(root/'certificate.json').stat().st_size,
        import_seconds=imported, query_seconds=seconds, execution=execution)


def run(binary, root):
    workers=min(4,invoke(binary,["resources"])[0]["worker_capacity"])
    report=control_case(binary,root,workers)
    save(root/"verification.json",report)
    return report


if __name__=="__main__":
    import sys
    print(json.dumps(run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())))
