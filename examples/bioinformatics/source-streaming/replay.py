#!/usr/bin/env python3
"""Reassemble original source pages and independently verify every returned byte."""
import argparse
import ctypes
import ctypes.util
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path

previous = Path(__file__).resolve().parents[1] / 'linear-observations/replay.py'
spec = importlib.util.spec_from_file_location('observation_example', previous)
observations = importlib.util.module_from_spec(spec)
spec.loader.exec_module(observations)
run = observations.run


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')
    return path


# The source producer consumes the actual UTF-8 declaration bytes.
observations.save = save
sha = lambda b: 'sha256:' + hashlib.sha256(b).hexdigest()
codec = ctypes.CDLL(ctypes.util.find_library('json-c'))
codec.json_tokener_parse.argtypes = [ctypes.c_char_p]
codec.json_tokener_parse.restype = ctypes.c_void_p
codec.json_object_to_json_string_ext.argtypes = [ctypes.c_void_p, ctypes.c_int]
codec.json_object_to_json_string_ext.restype = ctypes.c_char_p
codec.json_object_put.argtypes = [ctypes.c_void_p]


def encoded(value):
    pointer = codec.json_tokener_parse(json.dumps(value, ensure_ascii=False).encode())
    assert pointer
    try:
        return bytes(codec.json_object_to_json_string_ext(pointer, 2 | 16))
    finally:
        codec.json_object_put(pointer)


def seal(value, member):
    assert value[member] == sha(encoded({k: v for k, v in value.items() if k != member}))


def capture(binary, study, out, field_ids=None, page_records=3, page_bytes=512):
    out.mkdir(parents=True, exist_ok=False)
    index = run(binary, 'source-index', '--study', study)
    seal(index, 'index_sha256')
    save(out/'index.json', index)
    declaration = json.loads(study.read_text())
    documents = {d['id']: (study.parent/d['path']).read_bytes() for d in declaration['documents']}
    for d in declaration['documents']:
        assert sha(documents[d['id']]) == d['sha256']
    all_pages, assembled, totals = [], {}, {}

    def read(stream, limit):
        request = dict(schema='gadget.scientific.source.read.v1', source_digest=index['source_digest'],
                       publication='self', source_sequence=index['source_sequence'], stream=stream,
                       offset=0, limit=limit)
        pages, records, raw, end = [], [], b'', 0
        stream_digest = None
        while request is not None:
            page = run(binary, 'source-read', '--study', study,
                       '--request', save(out/'current-request.json', request))
            seal(page, 'page_sha256')
            assert page['offset'] == end and page['end'] >= end
            assert sha(encoded(page['binding'])) == page['stream_sha256']
            stream_digest = stream_digest or page['stream_sha256']
            assert page['stream_sha256'] == stream_digest
            if 'text' in page:
                content = page['text'].encode()
                assert len(content) == page['end']-end
                assert sha(content) == page['text_sha256']
                raw += content
            else:
                assert len(page['records']) == page['end']-end
                records.extend(page['records'])
            end = page['end']
            assert page['complete'] == (end == page['total'])
            request = page['next_request']
            assert (request is None) == page['complete']
            if request:
                assert request['offset'] == end
            pages.append(page)
        assert end == page['total']
        name = stream.get('section', stream.get('id'))
        name = stream['kind']+'/'+name
        assembled[name] = raw.decode() if stream['kind'] == 'document' else records
        totals[name] = dict(pages=len(pages), items_or_bytes=end)
        all_pages.extend(pages)
        return records, raw

    def evidence(value):
        if isinstance(value, dict):
            if 'text_sha256' in value and 'document' in value:
                assert 'text' not in value
                raw = documents[value['document']]
                assert sha(raw) == value['document_sha256']
                assert sha(raw[value['start']:value['end']]) == value['text_sha256']
            for child in value.values(): evidence(child)
        elif isinstance(value, list):
            for child in value: evidence(child)

    for section, count in index['catalog_sections'].items():
        records, _ = read(dict(kind='catalog', section=section), page_records)
        assert len(records) == count
        evidence(records)
        if section in declaration and section not in ['documents', 'publications']:
            key = (lambda x: (x['grade'], x['key'])) if section in ['cells', 'occurrences'] else (lambda x: (x.get('id', ''), x.get('key', '')))
            assert sorted(map(key, records)) == sorted(map(key, declaration[section]))
    for doc in declaration['documents']:
        if doc['format'] == 'utf8':
            _, raw = read(dict(kind='document', id=doc['id']), page_bytes)
            assert raw == documents[doc['id']]
    names = field_ids if field_ids is not None else [f['id'] for f in declaration['fields']]
    for name in names:
        rows, _ = read(dict(kind='field', id=name), page_records)
        assert [r['source_coordinate'] for r in rows] == list(range(len(rows)))
        assert all((r['value'] is not None) == r['reported'] for r in rows)
    for operation in declaration.get('operators', []):
        rows, _ = read(dict(kind='operator', id=operation['id']), page_records)
        assert [r['source_coordinate'] for r in rows] == list(range(len(rows)))
    save(out/'pages.json', all_pages)
    save(out/'assembled.json', assembled)
    save(out/'verification.json', dict(status='PASS', source_digest=index['source_digest'],
         pages=len(all_pages), streams=totals, all_original_UTF8_bytes_equal=True,
         all_catalogue_identities_equal=True, evidence_spans_verified=True))
    (out/'current-request.json').unlink()
    return assembled


def examples(binary, producer, out, clinical=None):
    out.mkdir(parents=True, exist_ok=False)
    old = observations.sources.square()[0]
    case = out/'authored'
    specification = observations.specification('source-pages', old['source'], 1, 4,
                                               observations.matrix([[1, 0, 0, 0], [0, 1, 0, 0]]))
    specification['research_note'] = ('Original β-catenin and HER2 declarations.\n' * 75)
    observations.execute(binary, producer, case, specification)
    assembled = capture(binary, case/'study/study.json', out/'authored-pages')
    for field in specification['fields']:
        rows = assembled['field/'+field['id']]
        assert [r['value'] for r in rows] == field['values']
    result = json.loads((case/'result.json').read_text())['result']
    for name, matrix in [('R', result['operator_map']), ('S', result['accession_map'])]:
        expected = {(e['row'], e['column']): Fraction(e['value']) for e in matrix['entries']}
        for row in assembled['operator/'+name]:
            assert Fraction(row['value']) == expected.get(tuple(row['coordinates']), Fraction(0))
    original = {o['key']: o for r in old['source']['relations'] for o in r['occurrences']}
    for retained in assembled['catalog/occurrences']:
        assertion = original[retained['key']]
        expected = Fraction(assertion['magnitude']) * (-1 if assertion['negative'] else 1) if assertion['exists'] else 0
        assert Fraction(retained['coefficient']) == expected
    # Retrieval never changes native observation preparations or results.
    repeated = run(binary, 'query', '--study', case/'study/study.json', '--query', case/'query.json')
    assert repeated == json.loads((case/'result.json').read_text())
    if clinical:
        capture(binary, clinical, out/'clinical-pages', ['all.age', 'all.time'], 31, 4096)
    return dict(status='PASS', native_response_unchanged=True, clinical_included=bool(clinical))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['binary', 'producer', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--clinical', type=Path)
    args = parser.parse_args()
    print(json.dumps(examples(args.binary.resolve(), args.producer.resolve(), args.out.resolve(), args.clinical)))
