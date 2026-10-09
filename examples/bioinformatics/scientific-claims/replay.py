"""Check selected source facts and completed analytical results, then relocate and replay."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path


def load(path):
    spec = importlib.util.spec_from_file_location('source_context_protocol', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(binary, root, fixtures):
    demo = load(fixtures / 'source-context/replay.py')
    manifest, intent, selection = demo.prepare(binary, root, fixtures)
    demo.execute(binary, root, manifest, intent, selection)
    manifest = root / 'relocated-workspace/workspace.json'
    context = json.loads((root / 'context.json').read_text())
    answer_path = root / 'relocated-pipeline/answers.rfcxs'
    answers = demo.run(binary,'record-view','--input',answer_path,'--member','result')
    facts = []

    def source(id, target, path, type, meaning):
        facts.append(dict(id=id, meaning=meaning, **{'from': dict(kind='source', target=target)},
                          path=path, value_type=type))

    def result(id, answer, step, path, type, meaning):
        a = next(a for a in answers['answers'] if a['id'] == answer)
        artifact = next(a for a in a['artifacts'] if a['step'] == step)
        facts.append(dict(id=id, meaning=meaning,
                          **{'from': dict(kind='answer', answer=answer, step=step,
                                          binding=artifact['binding'], result=0)},
                          path=path, value_type=type))

    source('reported-coordinates', 'reported-RNA', ['reported'], 'count',
           'Number of source-supported RNA coordinates in the selected field')
    source('unobserved-coordinates', 'missing-RNA', ['reported'], 'count',
           'Number of source-supported coordinates in the unobserved field')
    source('RNA-units', 'reported-RNA', ['units'], 'text', 'Declared units of the RNA field')
    for index in range(3):
        source('incidence-' + str(index) + '-coefficient', 'incidence-' + str(index),
               ['coefficient'], 'rational', 'Realized structural coefficient of this individual incidence')
    source('inactive-exists', 'incidence-2', ['aspects', 'exists'], 'boolean',
           'Retained existence aspect of the inactive occurrence')
    source('transcript-A1-length', 'transcript-A1', ['record', 'length'], 'count',
           'Decoded source length of the supplied A1 sequence')
    result('read-b0-quality-sum', 'quality', 'quality', ['reads', 0, 'quality_score_sum'], 'rational',
           'Sum of stored Phred+33 scores for baseline read b0')
    result('geneA-change-lower', 'expression', 'expression',
           ['expression', 'contrasts', 0, 'bounds', 3, 'lower'], 'rational',
           'Sharp lower bound of followup-minus-baseline assigned read counts for declared geneA')
    result('geneA-change-upper', 'expression', 'expression',
           ['expression', 'contrasts', 0, 'bounds', 3, 'upper'], 'rational',
           'Sharp upper bound of followup-minus-baseline assigned read counts for declared geneA')
    result('A1-change-direction', 'expression', 'expression',
           ['expression', 'contrasts', 0, 'bounds', 0, 'direction'], 'text',
           'Direction classification for the individual A1 allocation contrast')
    result('reported-f1-support', 'measurement-support', 'measurement-support-0',
           ['observations', 0, 'reported'], 'boolean', 'Reported support for f1 in the measured field')
    result('reported-f1-value', 'measurement-support', 'measurement-support-0',
           ['observations', 0, 'value'], 'rational', 'Measured RNA coordinate for f1')
    result('missing-f1-support', 'measurement-support', 'measurement-support-1',
           ['observations', 0, 'reported'], 'boolean', 'Reported support for f1 in the unobserved field')
    request = dict(schema='gadget.scientific.facts.v1', source_digest=context['source_digest'],
                   question=context['question'], context=context,
                   answers=dict(path='relocated-pipeline/answers.rfcxs', sha256=demo.digest(answer_path.read_bytes())),
                   facts=facts)
    request_path = demo.save(root / 'facts.json', request)
    catalog = demo.run(binary, 'fact-catalog', '--workspace', manifest,
                       '--request', request_path, '--artifacts', root)
    demo.save(root / 'facts-catalog.json', catalog)
    view = demo.run(binary, 'fact-catalog', '--workspace', manifest,
                    '--request', request_path, '--artifacts', root, '--view', 'model')
    if [f['id'] for f in view['facts']] != [f['id'] for f in catalog['facts']]:
        raise ValueError('Reading view omitted selected facts')
    demo.save(root / 'facts-model-view.json', view)
    proposal = dict(schema='gadget.scientific.claims.proposal.v1', action='claims',
                    source_digest=catalog['source_digest'], facts_sha256=catalog['facts_sha256'],
                    question_sha256=catalog['question_sha256'],
                    claims=[dict(id=f['id'], fact=f['id'], comparison='equal', value=f['value'])
                            for f in catalog['facts']])
    proposal_path = demo.save(root / 'claims.json', proposal)
    response = demo.run(binary, 'claim-proposal', '--workspace', manifest, '--facts', request_path,
                        '--request', proposal_path, '--artifacts', root)
    if response['status'] != 'matched' or response['certificate']['missing']:
        raise ValueError('Authored literal claims failed')
    demo.save(root / 'claims-receipt.json', response)
    false = copy.deepcopy(proposal)
    substitutions = {'read-b0-quality-sum': '0', 'A1-change-direction': 'positive_for_every_allocation',
                     'reported-coordinates': 0, 'inactive-exists': True}
    for claim in false['claims']:
        if claim['fact'] in substitutions:
            claim['value'] = substitutions[claim['fact']]
    false_path = demo.save(root / 'contradicted-claims.json', false)
    contradicted = demo.run(binary, 'claim-proposal', '--workspace', manifest, '--facts', request_path,
                            '--request', false_path, '--artifacts', root)
    if contradicted['status'] != 'contradicted' or contradicted['certificate']['contradicted'] != 4:
        raise ValueError('False literal claims were not contradicted')
    demo.save(root / 'contradicted-receipt.json', contradicted)
    return manifest, request_path, proposal_path, catalog


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    binary, root, fixtures = args.gadget.resolve(), args.out.resolve(), args.fixtures.resolve()
    _, _, _, before = prepare(binary, root, fixtures)
    moved = root.with_name(root.name + '-relocated')
    root.rename(moved)
    demo = load(fixtures / 'source-context/replay.py')
    after = demo.run(binary, 'fact-catalog', '--workspace', moved / 'relocated-workspace/workspace.json',
                     '--request', moved / 'facts.json', '--artifacts', moved)
    if before != after:
        raise ValueError('Captured facts changed after moving sources and analytical outputs')
    view_after = demo.run(binary, 'fact-catalog', '--workspace', moved / 'relocated-workspace/workspace.json',
                          '--request', moved / 'facts.json', '--artifacts', moved, '--view', 'model')
    if view_after != json.loads((moved / 'facts-model-view.json').read_text()):
        raise ValueError('Reading view changed after relocation')
    response = demo.run(binary, 'claim-proposal', '--workspace', moved / 'relocated-workspace/workspace.json',
                        '--facts', moved / 'facts.json', '--request', moved / 'claims.json', '--artifacts', moved)
    result = dict(status='PASS', selected_facts=len(after['facts']), matched=response['certificate']['matched'],
                  contradicted_examples=4, relocated_catalog_equal=True, interpretation_verified=False)
    demo.save(moved / 'facts-verification.json', result)
    print(json.dumps(dict(result, output=str(moved))))


if __name__ == '__main__':
    main()
