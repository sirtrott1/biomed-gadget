"""Reproduce complete table admission and offline biomedical evidence queries."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')
    return path


def digest(path):
    return 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary, *args, success=True):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True,
                            text=True, timeout=600)
    if (result.returncode == 0) != success:
        raise AssertionError((args, result.stderr))
    return json.loads(result.stdout) if success else result.stderr


def table(binary, root, names, rows, categorical, relations, *, arithmetic='exact_decimal',
          contrasts=None, identity=None, dictionary_bytes=0):
    root.mkdir(parents=True)
    raw = '\t'.join(names) + '\n' + ''.join('\t'.join(row) + '\n' for row in rows)
    (root / 'table.tsv').write_text(raw)
    definitions = [name + ': authored source variable.\n' for name in names]
    dictionary = ''.join(definitions)
    if dictionary_bytes:
        # Exercise full dictionary retention separately from source dimensions.
        dictionary += 'Authored codec byte contract.\n' * (
            dictionary_bytes // len('Authored codec byte contract.\n') + 1)
    (root / 'dictionary.txt').write_text(dictionary)
    columns = []
    offset = 0
    for name, definition in zip(names, definitions):
        end = offset + len(definition.encode())
        columns.append(dict(
            id=name, channel='observation' if name == identity else 'covariate' if name in categorical else 'E',
            quantity='authored:' + name, units='category' if name in categorical else 'counts',
            value_type='categorical' if name in categorical else 'quantitative',
            dictionary_start=offset, dictionary_end=end))
        offset = end
    spec = dict(schema='gadget.table.import.v1', id=root.name, source_namespace='921701',
                source=dict(path='table.tsv', sha256=digest(root / 'table.tsv'), delimiter='\t', header=True),
                dictionary=dict(path='dictionary.txt', sha256=digest(root / 'dictionary.txt')),
                columns=columns, missing_tokens=['?' + str(i) for i in range(33)],
                relation_columns=relations, arithmetic=arithmetic,
                question='Retain the complete declared source and its exact fields.',
                population='Authored admission experiment', selection='Every supplied record',
                limitations=['Authored software experiment without clinical observations'],
                source_url='authored://source-intake', attribution='GADGET source admission example',
                license='Authored example')
    if contrasts is not None:
        spec['contrasts'] = contrasts
    if identity:
        spec['identity_column'] = identity
    run(binary, 'import-table', '--spec', save(root / 'import.json', spec), '--out', root / 'study')
    return root / 'study'


def query(binary, root, value):
    result = run(binary, 'query', '--study', root / 'study/study.rfcxs',
                 '--query', save(root / 'query.json', value))
    save(root / 'result.json', result)
    return result


def experiments(binary, root, *, records=100001, columns=902, categories=1025,
                responses=33, evidence_rows=10001):
    root.mkdir(parents=True, exist_ok=False)
    names = ['identity', 'category', 'value']
    rows = [['record-' + str(i) + '-' + 'x' * 64, 'category-' + str(i % categories), str(i)]
            for i in range(records)]
    # A lexical identity can exceed the former 4096 byte text profile.
    rows[-1][0] += 'x' * 4097
    table(binary, root / 'records', names, rows, ['identity', 'category'], ['category'],
          identity='identity', dictionary_bytes=(1 << 20) + 1)
    query(binary, root / 'records', dict(operation='field_summary', field='all.value',
                                        support='reported_only'))
    membership = run(binary, 'query', '--study', root / 'records/study/study.rfcxs',
                     '--query', save(root / 'records/membership-query.json', dict(
                         operation='native_boundary', field_set='category.memberships',
                         backend='cpu_scheduled', tile_rows=257)))
    save(root / 'records/membership-result.json', membership)
    relation_names = ['category' + str(i) for i in range(9)]
    quantitative = ['value' + str(i) for i in range(columns - len(relation_names))]
    wide_names = relation_names + quantitative
    wide_rows = [['a'] * len(relation_names) + ['1/3', '1/6'] + ['2'] * (len(quantitative) - 2),
                 ['b'] * len(relation_names) + ['0'] * len(quantitative),
                 ['a'] * len(relation_names) + ['?32', '1/6'] + ['3'] * (len(quantitative) - 2)]
    contrasts = [dict(id='difference' + str(i), left=quantitative[0], right=quantitative[1])
                 for i in range(129)]
    table(binary, root / 'columns', wide_names, wide_rows, relation_names, relation_names,
          arithmetic='exact_rational', contrasts=contrasts)
    query(binary, root / 'columns', dict(operation='coverage', field='all.difference128'))
    # Native Hodge runs directly on the imported owner. No carrier augmentation.
    table(binary, root / 'hodge', ['group', 'value'],
          [['a' if i % 2 == 0 else 'b', str(i)] for i in range(12)], ['group'], ['group'])
    catalog = json.loads((root / 'hodge/study/admission.json').read_text())
    carriers = sorted((c for c in catalog['carriers'] if c['kind'] == 'source_grade'),
                      key=lambda c: c['grade'])
    query(binary, root / 'hodge', dict(operation='hodge_frame_decomposition',
          publication='self', source_sequence=catalog['source_sequence'], grade=0,
          grade_carriers=[c['key'] for c in carriers], fields=['all.value'],
          metrics=dict(kind='unit'), workers=2))
    snapshot = root / 'evidence/snapshot'; snapshot.mkdir(parents=True)
    receipts = []
    for i in range(responses):
        request = dict(id='binding' + str(i), provider='bindingdb', resource='uniprot',
                       identifier='AUTHORED_PROTEIN', cutoff_nm=1000001)
        n = evidence_rows if i == 0 else 1
        data = dict(getLindsByUniprotsResponse=dict(affinities=[
            dict(monomerid=str(r + 1), affinity='0' if r % 2 == 0 else '<1/2',
                 affinity_type='Ki', doi='10.authored/source',
                 source_note='Authored source-byte preservation. ' * (20 if i == 0 else 1))
            for r in range(n)]))
        # BindingDB numerical tokens are declared decimals, not rational source strings.
        for row in data['getLindsByUniprotsResponse']['affinities']:
            if row['affinity'] == '<1/2':
                row['affinity'] = '<0.5'
        body = save(snapshot / (request['id'] + '.json'), data)
        headers = snapshot / (request['id'] + '.headers.txt')
        headers.write_text('Offline authored transport receipt.\n' * (2049 if i == 0 else 1))
        receipts.append(dict(request=request, url='offline:authored', path=body.name,
                             sha256=digest(body), headers_path=headers.name,
                             headers_sha256=digest(headers)))
    manifest = dict(schema='gadget.database.snapshot.v1', records=receipts)
    # Match the application's canonical json-c seal, retaining the exact bytes.
    seal_json(manifest, 'snapshot_sha256')
    save(snapshot / 'snapshot.json', manifest)
    run(binary, 'import-evidence', '--snapshot', snapshot, '--out', root / 'evidence/store')
    for name, q in [('coverage', dict(operation='coverage')),
                    ('select', dict(operation='select', field='relation', comparison='equal', value='=')),
                    ('affinity', dict(operation='affinity_evidence', endpoint='Ki', units='nM', maximum='0.25'))]:
        result = run(binary, 'evidence-query', '--evidence', root / 'evidence/store',
                     '--query', save(root / 'evidence' / (name + '-query.json'), q))
        save(root / 'evidence' / (name + '-result.json'), result)
    fetch = dict(schema='gadget.database.fetch.v1', id='explicit-public-requests', requests=[
        dict(id='chembl', provider='chembl', resource='activity', identifier='CHEMBL941',
             limit=1000, offset=1000001), receipts[0]['request']],
        transport=dict(max_body_bytes=5 << 20, max_header_bytes=128 << 10,
                       max_total_body_bytes=80 << 20, connect_timeout_seconds=11, request_timeout_seconds=31))
    preview = run(binary, 'database-request', '--request', save(root / 'database-fetch.json', fetch))
    save(root / 'database-request.json', preview)


def seal_json(value, member):
    # json-c is the transport codec, not an independent mathematical verifier.
    import ctypes
    import ctypes.util
    lib = ctypes.CDLL(ctypes.util.find_library('json-c'))
    lib.json_tokener_parse.argtypes = [ctypes.c_char_p]; lib.json_tokener_parse.restype = ctypes.c_void_p
    lib.json_object_to_json_string_ext.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.json_object_to_json_string_ext.restype = ctypes.c_char_p
    lib.json_object_put.argtypes = [ctypes.c_void_p]
    obj = lib.json_tokener_parse(json.dumps(value).encode())
    try:
        value[member] = 'sha256:' + hashlib.sha256(lib.json_object_to_json_string_ext(obj, 2 | 64)).hexdigest()
    finally:
        lib.json_object_put(obj)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--records', type=int, default=100001); p.add_argument('--columns', type=int, default=902)
    p.add_argument('--categories', type=int, default=1025); p.add_argument('--responses', type=int, default=33)
    p.add_argument('--evidence-rows', type=int, default=10001)
    a = p.parse_args()
    experiments(a.binary.resolve(), a.out.resolve(), records=a.records, columns=a.columns,
                categories=a.categories, responses=a.responses, evidence_rows=a.evidence_rows)
    from verify import verify
    report = verify(a.out.resolve())
    save(a.out / 'verification.json', report)
    print(json.dumps(report))
