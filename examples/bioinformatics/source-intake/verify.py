"""Check all source identities, supports, native sectors and evidence selections."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path


def load(path):
    return json.loads(path.read_text())


def rational(value):
    return Fraction(value['rational'] if isinstance(value, dict) else value)


def table(root):
    spec = load(root / 'import.json'); study = load(root / 'study/study.json')
    catalog = load(root / 'study/admission.json'); receipt = load(root / 'study/import-receipt.json')
    for original, retained in [('table.tsv', 'table.txt'), ('dictionary.txt', 'dictionary.txt'),
                               ('import.json', 'import-spec.json')]:
        assert (root / original).read_bytes() == (root / 'study' / retained).read_bytes()
    assert receipt['key_policy'] == 'positive_source_order_per_namespace'
    assert receipt['producer_version'] in (2, 3)
    if receipt['producer_version'] == 3:
        intake = receipt['intake']
        assert receipt['evidence_profile'] == 'quoted_v1'
        assert intake['source_bytes'] == (root / 'table.tsv').stat().st_size
        assert intake['largest_record_bytes'] == max(map(len, (root / 'table.tsv').read_bytes().splitlines()))
        assert intake['token_offsets'] == receipt['rows'] * receipt['columns']
        assert intake['byte_backend'] == 'rexgraph_read_only_mapped_file'
    assert receipt['native_sha256'] == study['native']['sha256']
    assert receipt['native_sha256'] == 'sha256:' + hashlib.sha256((root / 'study' / study['native']['path']).read_bytes()).hexdigest()
    rows = (root / 'table.tsv').read_text().splitlines()[1:]
    cells = study['cells']; grouped = {k: [c for c in cells if c['grade'] == k] for k in (0, 1)}
    assert len(grouped[0]) == len(rows) == receipt['rows']
    assert receipt['columns'] == len(spec['columns'])
    assert receipt['categories'] == len(grouped[1])
    for grade, group in grouped.items():
        assert [int(c['key']) for c in group] == list(range(1, len(group) + 1))
    values = [row.split('\t') for row in rows]
    indices = {column['id']: i for i, column in enumerate(spec['columns'])}
    identity = spec.get('identity_column')
    for i, cell in enumerate(grouped[0]):
        expected = values[i][indices[identity]] if identity else 'record:' + spec['id'] + ':' + str(i + 1)
        assert cell['id'] == cell['label'] == expected
        evidence = cell['evidence'][0]
        assert evidence['document'] == 'table' and evidence['text'] == rows[i]
    groups = []
    for column in spec['relation_columns']:
        membership = {}
        for i, row in enumerate(values):
            value = row[indices[column]]
            if value not in spec['missing_tokens']:
                membership.setdefault(value, []).append(i)
        groups.extend((column, value, membership[value]) for value in sorted(membership))
    assert [cell['id'] for cell in grouped[1]] == [
        'category:' + column + ':' + value for column, value, _ in groups]
    expected_rows = [i for _, _, members in groups for i in members]
    assert len(study['occurrences']) == len(expected_rows)
    source_bytes = (root / 'table.tsv').read_bytes()
    for key, (occurrence, i) in enumerate(zip(study['occurrences'], expected_rows), 1):
        assert occurrence['grade'] == 1 and int(occurrence['key']) == key
        evidence = occurrence['evidence'][0]
        assert evidence['document'] == 'table' and evidence['text'] == rows[i]
        assert source_bytes[evidence['start']:evidence['end']].decode() == rows[i]
    assert catalog['identity_catalog_complete'] and len(catalog['identities']) == len(cells)
    assert {(c['grade'], c['key'], c['id']) for c in catalog['identities']} == {
        (c['grade'], c['key'], c['id']) for c in cells}
    carriers = catalog['carriers']; keys = [c['key'] for c in carriers]
    assert len(keys) == len(set(keys)) == receipt['carriers']
    native = sorted((c for c in carriers if c['kind'] == 'source_grade'), key=lambda c: c['grade'])
    assert len(native) == 2
    for carrier in native:
        assert carrier['extents'] == [len(grouped[carrier['grade']])]
    stalks = [c for c in carriers if c['kind'] == 'cell_stalk']
    assert len(stalks) == len(cells)
    assert {(c['cell']['grade'], c['cell']['key']) for c in stalks} == {(c['grade'], c['key']) for c in cells}
    for carrier in stalks:
        assert carrier['extents'] == [len(spec['relation_columns'])]
    tensors = {}
    for field in study['fields']:
        if field['key'] in tensors:
            assert tensors[field['key']] == field['id'].split('.', 1)[-1]
        tensors[field['key']] = field['id'].split('.', 1)[-1]
    assert sorted(map(int, tensors)) == list(range(1, receipt['tensors'] + 1))
    return dict(rows=len(rows), columns=len(spec['columns']), categories=len(grouped[1]),
                carriers=len(carriers), tensors=len(tensors), complete_identities=len(cells))


def verify(root):
    records = {name: table(root / name) for name in ['records', 'columns', 'hodge']}
    n = records['records']['rows']; summary = load(root / 'records/result.json')['result']
    assert summary['reported'] == n and summary['missing'] == 0
    assert summary['coordinates'] == list(range(n))
    assert rational(summary['mean']) == Fraction(n - 1, 2)
    assert rational(summary['sum_of_squares']) == Fraction((n - 1) * n * (2 * n - 1), 6)
    membership = load(root / 'records/membership-result.json')
    assert len(membership['result']) == n
    identity_index = [c['id'] for c in load(root / 'records/import.json')['columns']].index('identity')
    source_identities = [row.split('\t')[identity_index] for row in
                         (root / 'records/table.tsv').read_text().splitlines()[1:]]
    assert [row['entity'] for row in membership['result']] == source_identities
    assert all([rational(v) for v in row['values']] == [Fraction(1)] for row in membership['result'])
    assert membership['execution']['selected_workers'] <= membership['execution']['cpu_resources']['worker_capacity']
    coverage = load(root / 'columns/result.json')['result']
    assert coverage['reported'] == 2 and coverage['missing'] == 1
    assert [rational(x['value']) if x['value'] else None for x in coverage['observations']] == [Fraction(1, 6), Fraction(0), None]
    hodge = load(root / 'hodge/result.json')['result']
    field = hodge['fields'][0]
    expected = [Fraction(5 if i % 2 == 0 else 6) for i in range(12)]
    assert [rational(v) for v in field['sectors']['gradient']] == [Fraction(0)] * 12
    assert [rational(v) for v in field['sectors']['curl']] == expected
    assert [rational(v) for v in field['sectors']['harmonic']] == [Fraction(i) - expected[i] for i in range(12)]
    assert [rational(v) for v in field['reconstructed']] == list(map(Fraction, range(12)))
    manifest = load(root / 'evidence/snapshot/snapshot.json')
    source_rows = []
    for receipt in manifest['records']:
        body = root / 'evidence/snapshot' / receipt['path']
        assert 'sha256:' + hashlib.sha256(body.read_bytes()).hexdigest() == receipt['sha256']
        source_rows.extend(load(body)['getLindsByUniprotsResponse']['affinities'])
    cov = load(root / 'evidence/coverage-result.json')['result']
    assert cov['binding_records'] == cov['reported'] == len(source_rows) and cov['missing'] == 0
    eq = sum(row['affinity'] == '0' for row in source_rows)
    assert cov['equality_relation_records'] == eq
    selected = load(root / 'evidence/select-result.json')['result']
    assert selected['selected'] == len(selected['records']) == eq
    assert len({r['native_assertion_key'] for r in selected['records']}) == eq
    assert all(r['value_token'] == '0' and r['relation'] == '=' for r in selected['records'])
    affinity = load(root / 'evidence/affinity-result.json')['result']
    assert affinity['selected'] == eq and affinity['contradicted'] == 0
    assert affinity['indeterminate'] == len(source_rows) - eq
    assert len(affinity['records']) == len(source_rows)
    request = load(root / 'database-request.json')
    assert 'offset=1000001' in request['requests'][0]['url']
    assert 'cutoff=1000001' in request['requests'][1]['url']
    assert request['transport'] == load(root / 'database-fetch.json')['transport']
    return dict(status='PASS', tables=records, evidence_responses=len(manifest['records']),
                evidence_observations=len(source_rows), exact_selected_observations=eq,
                native_hodge_coordinates=12)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('directory', type=Path)
    a = p.parse_args(); print(json.dumps(verify(a.directory), indent=2))
