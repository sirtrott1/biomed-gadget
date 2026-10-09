"""Independent Fraction verification of the complete constructed gluing example."""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path


def verify_sources(result, workspace):
    """Resolve every evidence address against the complete retained publications."""
    manifest_path = Path(workspace)
    manifest = json.loads(manifest_path.read_text())
    assays = {r['declaration']['id']: r for r in manifest['assays']}
    spans = 0
    for binding in result['source_bindings'].values():
        partition = assays[binding['assay']]['partitions'][binding['partition']]
        study_path = manifest_path.parent / partition['path']
        study = json.loads(study_path.read_text())
        assert binding['source_namespace'] == partition['source_namespace']
        assert binding['source_digest'] == partition['source_digest']
        if 'source_pages' not in binding:
            # Compatibility rendering binds original byte artifacts, while the
            # independent structural/field oracle below checks native results.
            assert binding['evidence_profile'] == 'original-source-bytes-v1'
            assert binding['declaration_sha256'] == 'sha256:' + hashlib.sha256(study_path.read_bytes()).hexdigest()
            assert binding['native'] == study['native']
            assert binding['native']['sha256'] == 'sha256:' + hashlib.sha256(
                (study_path.parent / binding['native']['path']).read_bytes()).hexdigest()
            documents = {}
            for document in study['documents']:
                data = (study_path.parent / document['path']).read_bytes()
                assert document['sha256'] == 'sha256:' + hashlib.sha256(data).hexdigest()
                documents[document['id']] = data
            for section in ['cells', 'occurrences', 'fields']:
                for record in study[section]:
                    for span in record.get('evidence', []):
                        data = documents[span['document']]
                        assert 0 <= span['start'] <= span['end'] <= len(data)
                        if 'text' in span:
                            assert data[span['start']:span['end']].decode('utf8') == span['text']
                        spans += 1
            fields = {f['id']: f for f in study['fields']}
            for section in result['sections']:
                for role in ['left', 'right']:
                    if result['source_bindings'][role] is binding:
                        expected = fields[section['declaration'][role + '_field']]
                        assert section[role + '_meaning'] == {k: expected[k] for k in
                            ['id', 'key', 'value_type', 'quantity', 'units', 'channel']}
            continue
        index, pages = binding['source_index'], binding['source_pages']
        assert index['source_digest'] == binding['source_digest']
        assert index['source_sequence'] == binding['source_sequence']
        assert index['evidence_profile'] == binding['evidence_profile'] == 'document-bound-byte-spans-v1'
        documents = {r['id']: (study_path.parent / r['path']).read_bytes() for r in study['documents']}
        digests = {id: 'sha256:' + hashlib.sha256(value).hexdigest() for id, value in documents.items()}
        for section, page in pages.items():
            assert page['complete'] and page['next_request'] is None
            assert page['offset'] == 0 and page['end'] == page['total'] == len(study[section])
            assert page['total'] == index['catalog_sections'][section]
            assert page['binding']['source_digest'] == page['binding']['publication_digest'] == binding['source_digest']
            assert page['binding']['publication'] == 'self'
            assert page['binding']['source_sequence'] == binding['source_sequence']
            identities = lambda r: (r['grade'], r['key']) if section in ['cells', 'occurrences'] else r['id']
            original = {identities(r): r for r in study[section]}
            observed = {identities(r): r for r in page['records']}
            assert len(observed) == len(page['records']) and set(observed) == set(original)
            for identity, value in observed.items():
                source = original[identity]
                if section == 'documents':
                    assert value['sha256'] == digests[value['id']]
                    assert value['bytes'] == len(documents[value['id']])
                target = value['definition'] if section == 'fields' else value
                for key, entry in source.items():
                    if key != 'evidence': assert target[key] == entry
                evidence = target.get('evidence', [])
                assert len(evidence) == len(source.get('evidence', []))
                for ref, span in zip(evidence, source.get('evidence', [])):
                    assert 'text' not in ref
                    assert all(ref[k] == span[k] for k in ['document', 'start', 'end'])
                    data = documents[ref['document']][ref['start']:ref['end']]
                    if 'text' in span: assert data.decode('utf8') == span['text']
                    assert ref['document_sha256'] == digests[ref['document']]
                    assert ref['text_sha256'] == 'sha256:' + hashlib.sha256(data).hexdigest()
                    spans += 1
        meanings = {r['id']: r['definition'] for r in pages['fields']['records']}
        for section in result['sections']:
            for role in ['left', 'right']:
                if result['source_bindings'][role] is binding:
                    assert section[role + '_meaning'] == meanings[section['declaration'][role + '_field']]
    return spans


def verify(result, reference, workspace=None):
    if 'result' in result: result = result['result']
    roles = ['common', 'left', 'right']
    cell_maps = {(int(r['source']), r['from']['grade'], r['from']['key']):
                 (r['to']['grade'], r['to']['key']) for r in result['cell_maps']}
    occ_maps = {(int(r['source']), r['from']['grade'], r['from']['key']):
                (r['to']['grade'], r['to']['key']) for r in result['occurrence_maps']}
    assert len(cell_maps) == len(result['cell_maps'])
    assert len(occ_maps) == len(result['occurrence_maps'])
    original_cells, expected_occurrences, original_occurrences = {}, {}, set()
    for sid, role in enumerate(roles, 1):
        for i, id in enumerate(reference['orders'][role]): original_cells[sid, 0, id] = str(i + 1)
        records = reference['relations'][role]
        for i, r in enumerate(records): original_cells[sid, r['grade'], r['id']] = str(1000000 + i)
        slot = 0
        for r in records:
            g = r['grade']; parent = cell_maps[sid, g, original_cells[sid, g, r['id']]]
            for position, leg in enumerate(r['participants']):
                slot += 1; original_occurrences.add((sid, g, str(slot)))
                destination = occ_maps[sid, g, str(slot)]
                value = dict(parent=parent,
                             participant=cell_maps[sid, g - 1, original_cells[sid, g - 1, leg['cell']]],
                             position=position, exists=leg['exists'], orientation=leg['negative'],
                             coefficient=((-1 if leg['negative'] else 1) * Q(leg['magnitude'])
                                          if leg['exists'] else Q(0)))
                if destination in expected_occurrences: assert expected_occurrences[destination] == value
                expected_occurrences[destination] = value
    assert set(cell_maps) == {(sid, g, key) for (sid, g, _), key in original_cells.items()}
    assert set(occ_maps) == original_occurrences
    assert [r['cells'] for r in result['source_inventory']] == reference['expected_cells']
    assert [r['occurrences'] for r in result['source_inventory']] == reference['expected_occurrences']
    assert len({(r['grade'], r['key']) for r in result['retained_cells']}) == sum(reference['expected_cells'])
    scales = {(r['grade'], r['key']): Q(r['scale']) for r in result['retained_cells']}
    assert len(expected_occurrences) == len(result['retained_occurrences'])
    for r in result['retained_occurrences']:
        expected = expected_occurrences[r['grade'], r['key']]
        assert (r['parent']['grade'], r['parent']['key']) == expected['parent']
        assert (r['participant']['grade'], r['participant']['key']) == expected['participant']
        assert r['position'] == expected['position']
        assert r['aspects']['exists'] == expected['exists']
        assert r['aspects']['orientation'] == expected['orientation']
        assert Q(r['coefficient']) == expected['coefficient']
        # Original coordinated integer coefficients and cell scales realize
        # the rational coefficient; the lattice integer is not that rational.
        assert Q(r['lattice_coefficient']) * scales[expected['participant']] / scales[expected['parent']] == expected['coefficient']
    verified = 0
    for section in result['sections']:
        g = section['lineage'][0]['grade']; expected = {}
        for sid, role in [(2, 'left'), (3, 'right')]:
            if 'source_fields' in reference:
                field = reference['source_fields'][role][section['declaration'][role + '_field']]
                assert field['grade'] == g
                values = {id: Q(v) for id, v in field['values'].items()}
            elif g == 0:
                values = {id: Q(i * 10) for i, id in enumerate(f'{role}-shared-{i}' for i in range(reference['vertices']))}
                values[role + '-unique'] = Q(30 if role == 'left' else 40)
            elif g == 1:
                values = {r['id']: Q(v) for r, v in zip(reference['relations'][role][:4], [2, 3, 5, 7])}
            else: values = {role + '-face': Q(7 if role == 'left' else 11)}
            for id, value in values.items():
                target = cell_maps[sid, g, original_cells[sid, g, id]]
                if target in expected: assert expected[target] == value
                expected[target] = value
        assert len(section['lineage']) == (g + 1 if result['boundary_policy'] == 'to_grade_zero' else 1)
        for level in section['lineage']:
            actual = {(r['grade'], r['key']): Q(r['value']) for r in level['values']}
            assert actual == expected, (actual, expected)
            verified += len(actual)
            boundary = {}
            for occurrence in expected_occurrences.values():
                parent = occurrence['parent']
                if parent not in expected: continue
                # The native boundary retains visited inactive/zero coordinates.
                lower = occurrence['participant']; boundary.setdefault(lower, Q(0))
                if occurrence['exists']: boundary[lower] += occurrence['coefficient'] * expected[parent]
            expected = boundary
    assert result['execution']['participant_publications_reused']
    assert result['execution']['completed_tiles'] == len(result['sections'])
    assert result['constructed_source_status'] == 'prospective'
    spans = verify_sources(result, workspace) if workspace is not None else None
    return dict(status='PASS', verified_coordinates=verified, verified_evidence_spans=spans,
                cells=sum(reference['expected_cells']), occurrences=len(expected_occurrences),
                sections=len(result['sections']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--workspace', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(json.loads(args.result.read_text()), json.loads(args.reference.read_text()), args.workspace)))
