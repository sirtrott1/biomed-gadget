"""Independent Fraction solve and original coordinate observation checks."""
import json
import csv
import hashlib
from fractions import Fraction as F


def dense(m):
    out = [[F(0) for _ in range(m['columns'])] for _ in range(m['rows'])]
    for e in m['entries']:
        out[e['row']][e['column']] = F(e['value'])
    return out


def solve(a, b):
    # This independent verifier is deliberately separate from production kernels.
    n = len(a)
    rows = [{j: q for j, q in enumerate(row) if q} for row in a]
    rhs = list(map(F, b))
    for k in range(n):
        pivot = next(i for i in range(k, n) if rows[i].get(k, 0))
        rows[k], rows[pivot] = rows[pivot], rows[k]
        rhs[k], rhs[pivot] = rhs[pivot], rhs[k]
        for i in range(k+1, n):
            if k not in rows[i]:
                continue
            factor = rows[i][k] / rows[k][k]
            for j, q in rows[k].items():
                value = rows[i].get(j, 0) - factor*q
                if value:
                    rows[i][j] = value
                else:
                    rows[i].pop(j, None)
            rhs[i] -= factor*rhs[k]
    x = [F(0)] * n
    for i in reversed(range(n)):
        x[i] = (rhs[i] - sum(q*x[j] for j, q in rows[i].items() if j > i)) / rows[i][i]
    assert all(sum(q*x[j] for j, q in enumerate(row)) == F(b[i]) for i, row in enumerate(a))
    return x


def verify(root):
    certificate = json.loads((root/'result.json').read_text())
    result = certificate['result']
    a, s = dense(result['operator_map']), dense(result['accession_map'])
    publication = result['publication']
    if publication.get('evidence_profile') == 'document-bound-byte-spans-v1':
        digest = lambda b: 'sha256:' + hashlib.sha256(b).hexdigest()
        assert publication['declaration']['sha256'] == digest((root/'study/study.json').read_bytes())
        documents = {d['id']: d for d in publication['declaration']['documents']}
        def references(value):
            if isinstance(value, dict):
                if 'document_sha256' in value:
                    document = documents[value['document']]
                    raw = (root/'study'/document['path']).read_bytes()
                    assert digest(raw) == document['sha256'] == value['document_sha256']
                    assert digest(raw[value['start']:value['end']]) == value['text_sha256']
                for child in value.values():
                    references(child)
            elif isinstance(value, list):
                for child in value:
                    references(child)
        references(publication)
    spec = json.loads((root/'spec.json').read_text())
    if 'source' in spec:
        declaration = spec['operators'][0]
        if 'map' in declaration:
            assert a == dense(declaration['map'])
        else:
            grade = declaration['hodge_grade']
            source = spec['source']
            ids = [source['vertices'], [r['key'] for r in source['relations']]]
            for relation in source.get('higher', []):
                while len(ids) <= relation['grade']:
                    ids.append([])
                ids[relation['grade']].append(relation['key'])
            def boundary(k):
                relations = source['relations'] if k == 1 else [r for r in source.get('higher', []) if r['grade'] == k]
                b = [[F(0) for _ in ids[k]] for _ in ids[k-1]]
                for j, relation in enumerate(relations):
                    for occurrence in relation['occurrences']:
                        if occurrence['exists']:
                            i = ids[k-1].index(occurrence['cell'])
                            b[i][j] += (-1 if occurrence['negative'] else 1) * F(occurrence['magnitude'])
                return b
            n = len(a)
            expected = [[F(i == j) for j in range(n)] for i in range(n)]
            if grade:
                b = boundary(grade)
                for i in range(n):
                    for j in range(n):
                        expected[i][j] += sum(row[i]*row[j] for row in b)
            if grade+1 < len(ids) and ids[grade+1]:
                c = boundary(grade+1)
                for i in range(n):
                    for j in range(n):
                        expected[i][j] += sum(x*y for x, y in zip(c[i], c[j]))
            assert a == expected, 'operator differs from retained boundary tower'
        assert s == dense(spec['operators'][1]['map'])
        field_values = {f['id']: f['values'] for f in spec['fields']}
        assert all(field['input'] == field_values[field['id']] for field in result['fields'])
    else:
        table = root/'study/original/document-0'
        # Resolve the actual original table by its retained document declaration.
        retained = json.loads((root/'study/original/study.json').read_text())
        table_id = next(d['id'] for d in retained['documents'] if d['path'] == 'table.txt')
        table_record = next(d for d in result['publication']['declaration']['documents'] if d['id'] == table_id)
        with (root/'study'/table_record['path']).open(newline='') as stream:
            observations = list(csv.DictReader(stream))
        assert len(observations) == len(a)
        assert a == [[F(i == j) + F(observations[i]['sex'] == observations[j]['sex'])
                      for j in range(len(a))] for i in range(len(a))]
        for field in result['fields']:
            column = field['id'].split('.', 1)[1]
            assert list(map(F, field['input'])) == [F(row[column]) for row in observations]
    rows = result['fields']
    for field in rows:
        x = solve(a, field['input'])
        observed = [sum(q*x[j] for j, q in enumerate(row)) for row in s]
        assert observed == list(map(F, field['values'])), root
    if result['observation']['kind'] == 'coordinates':
        ordered = result['observation']['coordinates']
        actual = result['publication']['source_coordinates']
        assert result['selected_source_cells'] == [dict(grade=actual[i]['grade'], key=actual[i]['key']) for i in ordered]
    else:
        exact = list(dict.fromkeys(tuple(row) for row in s))
        native = result['native_preparation']
        assert native['distinct_coefficient_rows'] == len(exact)
        assert native['row_reuse_word'] == [exact.index(tuple(row)) for row in s]
        nonzero = sum(any(row) for row in exact)
        assert native['uses_full_solve'] == (nonzero >= len(a))
        assert native['certified_nonzero_rows'] == (0 if native['uses_full_solve'] else nonzero)
    if (root/'declared-result.json').exists():
        declared = json.loads((root/'declared-result.json').read_text())['result']
        assert [f['values'] for f in declared['fields']] == [f['values'] for f in rows]
        assert declared['ordering']['rows'] != declared['ordering']['columns']
    return dict(case=root.name, source_coordinates=len(a), output_rows=len(s), fields=len(rows),
                output_coordinates=len(s)*len(rows), uses_full_solve=result['native_preparation']['uses_full_solve'])
