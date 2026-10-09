"""Compare quoted and shared-span admission of complete declared input tables."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def digest(data):
    return 'sha256:' + hashlib.sha256(data).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    return path


def run(binary, *args):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def evidence_records(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == 'evidence':
                yield from item
            else:
                yield from evidence_records(item)
    elif isinstance(value, list):
        for item in value:
            yield from evidence_records(item)


def semantic_metadata(value):
    if isinstance(value, dict):
        return {key: semantic_metadata(item) for key, item in value.items() if key != 'evidence'}
    if isinstance(value, list):
        return [semantic_metadata(item) for item in value]
    return value


def compare(binary, specification, out):
    """Both imports retain the same original table/dictionary and native archive."""
    spec_path = Path(specification)
    original = json.loads(spec_path.read_text())
    out.mkdir(parents=True, exist_ok=False)
    outputs = {}
    for name, profile in [('quoted', 'quoted_v1'), ('spans', 'document_spans_v2')]:
        spec = {**original, 'evidence_profile': profile}
        # Intake sources are staged beside their specification without changing
        # any original bytes or asserting a correspondence between publications.
        root = out / name
        root.mkdir()
        for kind in ['source', 'dictionary']:
            relative = Path(original[kind]['path'])
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('input artifacts must stay relative to their specification')
            source = spec_path.parent / relative
            target = root / original[kind]['path']
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
        source_spec = save(root / 'spec.json', spec)
        study = root / 'study'
        run(binary, 'import-table', '--spec', source_spec, '--out', study)
        declaration = json.loads((study / 'study.json').read_text())
        receipts = json.loads((study / 'import-receipt.json').read_text())
        records = list(evidence_records(declaration))
        assert all(('text' in record) == (name == 'quoted') for record in records)
        documents = {doc['id']: (study / doc['path']).read_bytes() for doc in declaration['documents']}
        for record in records:
            byte_slice = documents[record['document']][record['start']:record['end']]
            assert byte_slice and 0 <= record['start'] < record['end'] <= len(documents[record['document']])
            if 'text' in record:
                assert record['text'].encode() == byte_slice
        mathematical = {}
        for field in declaration['fields']:
            if field['context'] == 'all' and field['value_type'] != 'categorical' and field['quantity'] != 'declared_category_membership_basis':
                query = save(root / 'query.json', dict(operation='field_summary', field=field['id'], support='reported_only'))
                mathematical[field['id']] = run(binary, 'query', '--study', study / 'study.json', '--query', query)['result']
        query = save(root / 'boundary-query.json', dict(operation='native_boundary', field_set='category.memberships', backend='cpu_scheduled', tile_rows=1))
        mathematical['native_category_boundary'] = run(binary, 'query', '--study', study / 'study.json', '--query', query)['result']
        index = run(binary, 'source-index', '--study', study / 'study.json')
        req = dict(schema='gadget.scientific.source.read.v1',source_digest=index['source_digest'],publication='self',source_sequence=index['source_sequence'],stream=dict(kind='catalog',section='cells'),offset=0,limit=len(declaration['cells']))
        page = run(binary, 'source-read', '--study', study / 'study.json', '--request', save(root / 'read.json', req))
        assert len(page['records']) == len(declaration['cells']) and page['complete']
        for item in page['records']:
            for record in item['evidence']:
                assert record['document_sha256'] == digest(documents[record['document']])
                assert record['text_sha256'] == digest(documents[record['document']][record['start']:record['end']])
        save(root / 'source-index.json', index)
        save(root / 'source-page.json', page)
        save(root / 'mathematical-results.json', mathematical)
        outputs[name] = dict(declaration=declaration, receipts=receipts, mathematical=mathematical, native=digest((study/'source.rfcxs').read_bytes()), bytes=(study/'study.json').stat().st_size, evidence=len(records))
    q, s = outputs['quoted'], outputs['spans']
    assert q['native'] == s['native']
    assert q['mathematical'] == s['mathematical']
    for section in ['cells','occurrences','fields','contexts','field_sets','design']:
        assert semantic_metadata(q['declaration'][section]) == semantic_metadata(s['declaration'][section])
    result = dict(status='PASS', rows=s['receipts']['rows'], columns=s['receipts']['columns'], evidence_spans=s['evidence'], native_archive_identical=True, complete_mathematical_results_identical=True, semantic_identities_identical=True, quoted_declaration_bytes=q['bytes'], span_declaration_bytes=s['bytes'], byte_ratio=s['bytes']/q['bytes'], intake=s['receipts']['intake'], fields=len(q['mathematical']))
    save(out / 'verification.json', result)
    return result


def authored(root, rows):
    root.mkdir(parents=True, exist_ok=False)
    raw = 'identity\tgroup\tx\ty\tz\n' + ''.join(f'sample:{i}\t{"α" if i%2 else "β"}\t{i}/3\t{"?" if i%3==0 else str(i-2)}\t0\n' for i in range(rows))
    dictionary = 'Original authored identities, lexical groups, exact rational x, reported y with ? absence, and reported zero z.\n'
    (root/'table.tsv').write_text(raw)
    (root/'dictionary.txt').write_text(dictionary)
    spec = dict(schema='gadget.table.import.v1',id='source-admission-authored',source_namespace='940701',source=dict(path='table.tsv',sha256=digest(raw.encode()),delimiter='\t',header=True),dictionary=dict(path='dictionary.txt',sha256=digest(dictionary.encode())),columns=[dict(id=n,channel='covariate' if n in ['identity','group'] else 'E',quantity='authored_'+n,units='unit',value_type='categorical' if n in ['identity','group'] else 'quantitative',dictionary_start=0,dictionary_end=len(dictionary.encode())) for n in ['identity','group','x','y','z']],missing_tokens=['?'],relation_columns=['group'],identity_column='identity',arithmetic='exact_rational',contrasts=[dict(id='x_minus_y',left='x',right='y')],question='Do shared source spans preserve the complete exact table and complex?',population='All authored records',selection='Every original input row',limitations=['Authored implementation qualification; no biological interpretation.'],source_url='authored://source-admission',attribution='GADGET contract example',license='Authored fixture')
    return save(root/'spec.json',spec)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--rows',type=int,default=24)
    p.add_argument('--spec',type=Path,help='An existing strict table specification, beside its original files')
    a=p.parse_args()
    if a.rows <= 0:
        p.error('--rows must be positive')
    a.out.mkdir(parents=True,exist_ok=False)
    spec=a.spec or authored(a.out/'original',a.rows)
    result=compare(a.binary.resolve(),spec.resolve(),a.out/'comparison')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
