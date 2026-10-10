"""Source annotations linked to a reproducible RNA answer protocol."""
import argparse
import copy
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

QUESTION = (
    'Which gene-level expression changes are fixed across all retained read allocations, '
    'what baseline quality scores are recorded, and which RNA measurement coordinates '
    'are reported or missing? Relate these requests to the declared transcript, read, '
    'field, context and incidence records.'
)


def digest(data):
    return 'sha256:' + hashlib.sha256(data).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')
    return path


def run(binary, *args):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True,
                            text=True, timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr)
    return json.loads(result.stdout)


def prepare(binary, root, fixtures):
    """Stage exactly used inputs; the generated capsule can move independently."""
    if root.exists():
        raise ValueError('The output directory must be new')
    root.mkdir(parents=True)
    inputs = root / 'inputs'
    inputs.mkdir()
    rna = fixtures / 'rna-evidence'
    matrix = fixtures / 'staged-assays'
    spec = json.loads((rna / 'import.json').read_text())
    assay = json.loads((matrix / 'import.json').read_text())['assays'][0]
    spec['id'] = 'scientific-source-context'
    spec['design']['question'] = QUESTION
    spec['assays'] = [assay]
    for item, source in [(r['source'], rna) for r in spec['resources']] + [
            (assay[k], matrix) for k in ['source', 'features', 'observations']]:
        original = source / item['path']
        content = original.read_bytes()
        if digest(content) != item['sha256']:
            raise ValueError('Fixture source hash mismatch: ' + str(original))
        destination = inputs / original.name
        destination.write_bytes(content)
        item['path'] = original.name
    relation = Path(__file__).with_name('relations.json').read_bytes()
    (inputs / 'relations.json').write_bytes(relation)
    assay['relations'] = {'path': 'relations.json', 'sha256': digest(relation)}
    import_path = save(inputs / 'import.json', spec)
    catalog = run(binary, 'import-assays', '--spec', import_path, '--out', root / 'workspace')
    manifest = root / 'workspace/workspace.json'
    owner = catalog['assays'][0]['partitions'][0]
    study_root = root / 'workspace' / Path(owner['path']).parent
    study = json.loads((root / 'workspace' / owner['path']).read_text())
    entries = []

    def native(id, kind, selector, evidence):
        target = dict(kind=kind, assay='rna', partition=0,
                      source_digest=owner['source_digest'], **selector)
        entries.append(dict(id=id, target=target, evidence=evidence))

    def excerpt(document, quote):
        source = next(d for d in study['documents'] if d['id'] == document)
        data = (study_root / source['path']).read_bytes()
        word = quote.encode()
        start = data.index(word)
        return [dict(document=document, start=start, end=start + len(word), text=quote)]

    cell = next(c for c in study['cells'] if c['id'] == 'f1')
    native('feature-A', 'cell', {'id': 'f1'}, cell['evidence'])
    native('reported-RNA', 'field', {'id': 'all.rna-before-1'},
           excerpt('assay_declaration', 'RNA molecules'))
    native('missing-RNA', 'field', {'id': 'all.rna-unobserved-1'},
           excerpt('assay_declaration', 'unobserved'))
    native('RNA-context', 'context', {'id': 'all'},
           excerpt('assay_declaration', 'Authored assay-section example'))
    for index, occurrence in enumerate(study['occurrences']):
        native('incidence-' + str(index), 'occurrence',
               {'grade': occurrence['grade'], 'key': occurrence['key']}, occurrence['evidence'])
    inspections = {}
    for resource_id, record_index, id in [('transcripts.fa', 0, 'transcript-A1'),
                                          ('transcripts.fa', 1, 'transcript-A2'),
                                          ('baseline.fastq', 0, 'baseline-read')]:
        resource = next(r for r in catalog['resources'] if r['declaration']['id'] == resource_id)
        if resource_id not in inspections:
            declaration = resource['declaration']
            specification = dict(schema='gadget.bioinformatics.file.v1',
                                 format=declaration['format'],
                                 source=dict(path=str(Path('workspace') / resource['original']['path']),
                                             sha256=resource['original']['sha256']))
            if 'reference' in declaration:
                specification['reference'] = declaration['reference']
            inspections[resource_id] = run(binary, 'inspect-file', '--spec',
                                           save(root / 'inspection.json', specification))
        record = inspections[resource_id]['records'][record_index]
        data = (root / 'workspace' / resource['original']['path']).read_bytes()
        start, end = record['source_start'], record['source_end']
        entries.append(dict(id=id,
                            target=dict(kind='resource', resource=resource_id, artifact='source',
                                        sha256=resource['original']['sha256'], record=record_index),
                            evidence=[dict(document='source', start=start, end=end,
                                           text=data[start:end].decode())]))
    context = dict(schema='gadget.assays.source.context.v1', source_digest=catalog['manifest_sha256'],
                   question=QUESTION, targets=entries)
    context_path = save(root / 'context.json', context)
    view = run(binary, 'source-context', '--workspace', manifest, '--request', context_path)
    save(root / 'context-catalog.json', view)
    meanings = {
        'feature-A': 'Retained RNA feature labelled gene-A in the feature dictionary.',
        'reported-RNA': 'RNA molecule counts in a reported baseline observation.',
        'missing-RNA': 'Unobserved RNA coordinates remain outside the reported support.',
        'RNA-context': 'The declared whole-assay context retains its own metric.',
        'incidence-0': 'The negative occurrence of feature f1 in the authored relation.',
        'incidence-1': 'The positive occurrence of the same feature f1; its identity is distinct.',
        'incidence-2': 'The inactive f2 occurrence is retained although its realized coefficient is zero.',
        'transcript-A1': 'Transcript A1 is a supplied assignment candidate within declared geneA.',
        'transcript-A2': 'Transcript A2 is another assignment candidate; equality of gene labels supplies no native identity map.',
        'baseline-read': 'The first baseline FASTQ record retains sequence and measured quality symbols.'
    }
    proposal = dict(schema='gadget.assays.annotation.proposal.v1', action='annotate',
                    source_digest=view['source_digest'], context_sha256=view['context_sha256'],
                    question_sha256=view['question_sha256'],
                    mentions=[dict(id=e['id'], target=e['id'], polarity='asserted',
                                   interpretation=meanings[e['id']]) for e in entries])
    annotation_path = save(root / 'annotations.json', proposal)
    response = run(binary, 'annotation-proposal', '--workspace', manifest,
                   '--context', context_path, '--request', annotation_path)
    save(root / 'annotation-receipt.json', response)
    intent = json.loads((rna / 'intent.json').read_text())
    intent['answer_contract']['question'] = QUESTION
    intent['recipes'].append(dict(id='measurement-support', purpose='Reported and missing RNA support',
                                  goal=dict(method='coverage', assay='rna', parameters=dict(
                                      fields=['all.rna-before-1', 'all.rna-unobserved-1']))))
    intent['answer_contract']['requirements'].append(dict(id='measurement-support',
                                                         observable='field_support',
                                                         bindings=['measurement-support']))
    intent['source_context'] = dict(schema='gadget.pipeline.source.context.v1', context=context,
                                    annotations=proposal, requirements=[
                                        dict(id='expression', mentions=['transcript-A1', 'transcript-A2']),
                                        dict(id='quality', mentions=['baseline-read']),
                                        dict(id='measurement-support', mentions=[e['id'] for e in entries
                                             if e['target']['kind'] != 'resource'])])
    intent_path = save(root / 'intent.json', intent)
    recipes = run(binary, 'recipe-catalog', '--workspace', manifest, '--intent', intent_path)
    save(root / 'recipe-catalog.json', recipes)
    selection = dict(schema='gadget.pipeline.recipe.proposal.v1', action='recipes',
                     source_digest=recipes['source_digest'], intent_sha256=recipes['intent_sha256'],
                     goals=[r['id'] for r in intent['recipes']])
    save(root / 'selection.json', selection)
    return manifest, intent_path, root / 'selection.json'


def execute(binary, root, manifest, intent_path, selection_path):
    intent = json.loads(intent_path.read_text())
    planning = run(binary, 'recipe-proposal', '--workspace', manifest, '--intent', intent_path,
                   '--request', selection_path, '--prompt', QUESTION, '--out', root / 'pipeline')
    save(root / 'execution.json', planning)
    answers = run(binary,'record-view','--input',root/'pipeline/answers.rfcxs','--member','result')
    if answers['source_context'] != planning['source_context']:
        raise ValueError('Completed answers do not retain the admitted source context')
    for answer in answers['answers']:
        mapping = next(r for r in intent['source_context']['requirements'] if r['id'] == answer['id'])
        if answer['source_mentions'] != mapping['mentions']:
            raise ValueError('Answer mention mapping mismatch')
    # Relocation removes all runtime access to the initially admitted paths.
    relocated = root / 'relocated-workspace'
    (root / 'workspace').rename(relocated)
    before = json.loads((root / 'context-catalog.json').read_text())
    after = run(binary, 'source-context', '--workspace', relocated / 'workspace.json',
                '--request', root / 'context.json')
    if before != after:
        raise ValueError('Source context changed after relocation')
    (root / 'pipeline').rename(root / 'relocated-pipeline')
    replay = run(binary, 'replay-derivation', '--workspace',
                 root / 'relocated-pipeline/expression/workspace/workspace.json',
                 '--assay', 'expression', '--out', root / 'replayed-expression')
    if replay['status'] != 'PASS':
        raise ValueError('Detached expression replay failed')
    save(root / 'replay-expression.json', replay)
    report = dict(status='PASS', answer_requirements=len(answers['answers']),
                  selected_bindings=len(intent['recipes']),
                  source_targets=len(before['targets']),
                  evidence_identity_verified=True, interpretation_verified=False,
                  relocated_source_context_equal=True, detached_expression_replay='PASS')
    save(root / 'verification.json', report)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    binary, root = args.gadget.resolve(), args.out.resolve()
    manifest, intent, selection = prepare(binary, root, args.fixtures.resolve())
    print(json.dumps(execute(binary, root, manifest, intent, selection)))


if __name__ == '__main__':
    main()
