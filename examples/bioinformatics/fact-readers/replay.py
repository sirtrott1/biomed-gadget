"""Reproduce full and focused reading selections over RNA and chromatin results."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def load(path):
    spec = importlib.util.spec_from_file_location('fact_reader_protocol', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')
    return path


def run(binary, *args):
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr)
    return json.loads(p.stdout)


def capture(binary, manifest, root, path):
    catalog = run(binary, 'fact-catalog', '--workspace', manifest, '--request', path,
                  '--artifacts', root)
    view = run(binary, 'fact-catalog', '--workspace', manifest, '--request', path,
               '--artifacts', root, '--view', 'model')
    if [f['id'] for f in view['facts']] != [f['id'] for f in catalog['facts']]:
        raise ValueError('Reading selection changed')
    prefix = path.stem
    save(root / (prefix + '-catalog.json'), catalog)
    save(root / (prefix + '-view.json'), view)
    proposal = dict(schema='gadget.scientific.claims.proposal.v1', action='claims',
                    **{k: catalog[k] for k in ['source_digest', 'facts_sha256', 'question_sha256']},
                    claims=[dict(id=f['id'], fact=f['id'], comparison='equal', value=f['value'])
                            for f in catalog['facts']])
    proposal_path = save(root / (prefix + '-claims.json'), proposal)
    response = run(binary, 'claim-proposal', '--workspace', manifest, '--facts', path,
                   '--request', proposal_path, '--artifacts', root)
    if response['status'] != 'matched' or response['certificate']['missing']:
        raise ValueError('Selected literal values failed admission')
    save(root / (prefix + '-receipt.json'), response)
    return catalog, view


def prepare(binary, root, fixtures):
    root.mkdir(parents=True, exist_ok=False)
    rna = load(fixtures / 'scientific-claims/replay.py')
    manifest, facts_path, _, _ = rna.prepare(binary, root / 'rna', fixtures)
    rna_request = json.loads(facts_path.read_text())
    focused = copy.deepcopy(rna_request)
    selected = {'geneA-change-lower', 'geneA-change-upper', 'A1-change-direction'}
    focused['facts'] = [f for f in focused['facts'] if f['id'] in selected]
    focused_path = save(root / 'rna/expression-facts.json', focused)
    capture(binary, manifest, root / 'rna', focused_path)

    chromatin = root / 'chromatin'
    p = subprocess.run([sys.executable, str(fixtures / 'replay-answers.py'),
                        '--gadget', str(binary), '--protocol', str(fixtures / 'controlled-coverage'),
                        '--out', str(chromatin)], capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr)
    answers_path = chromatin / 'relocated-pipeline/answers.json'
    answers = json.loads(answers_path.read_text())
    intent = json.loads((chromatin / 'inputs/intent.json').read_text())
    selections = []

    def select(id, answer_id, path, value_type, meaning):
        answer = next(a for a in answers['answers'] if a['id'] == answer_id)
        artifact = next(a for a in answer['artifacts'] if a['step'] == answer_id)
        selections.append(dict(id=id, path=path, value_type=value_type, meaning=meaning,
                               **{'from': dict(kind='answer', answer=answer_id,
                                               binding=artifact['binding'], step=artifact['step'], result=0)}))

    select('treated-residual-first', 'contrast', ['contrasts', 0, 'residual', 0], 'rational',
           'Treated-minus-baseline residual on physical segment [0,4), outside the declared control span')
    select('treated-residual-second', 'contrast', ['contrasts', 0, 'residual', 1], 'rational',
           'Treated-minus-baseline residual on physical segment [4,8), outside the declared control span')
    select('treated-residual-energy', 'contrast', ['contrasts', 0, 'quadrances', 'residual'], 'rational',
           'Residual quadrance under the declared physical-length metric')
    select('coverage-kernel-nullity', 'contrast', ['source_kernel', 'nullity'], 'count',
           'Dimension of the exact zero-coverage sector on original interval variables')
    select('control-orthogonality', 'contrast', ['contrasts', 0, 'control_orthogonal'], 'boolean',
           'Checked metric orthogonality of the residual to the declared controls')
    select('field-reconstruction', 'contrast', ['contrasts', 0, 'reconstruction_exact'], 'boolean',
           'Checked exact reconstruction from represented and residual physical fields')
    for coordinate, position in enumerate(['first', 'second', 'third']):
        select('green-' + position, 'green', [coordinate, 'value'], 'rational',
               'Native Green action on physical segment ' + position + ' under the declared operator and metric')
    for sector in ['gradient', 'curl', 'harmonic']:
        select(sector + '-quadrance', 'structure', ['quadrance', sector], 'rational',
               'Native ' + sector + ' sector quadrance of the declared original-interval field')
    request = dict(schema='gadget.scientific.facts.v1', source_digest=answers['source_digest'],
                   question=intent['answer_contract']['question'],
                   answers=dict(path='relocated-pipeline/answers.json',
                                sha256='sha256:' + hashlib.sha256(answers_path.read_bytes()).hexdigest()),
                   facts=selections)
    chrom_manifest = chromatin / 'relocated-source/workspace.json'
    chrom_path = save(chromatin / 'structure-facts.json', request)
    catalog, _ = capture(binary, chrom_manifest, chromatin, chrom_path)
    expected = {'treated-residual-first': '3/2', 'treated-residual-second': '-3/2',
                'treated-residual-energy': '18', 'coverage-kernel-nullity': 8,
                'control-orthogonality': True, 'field-reconstruction': True,
                'green-first': '3/26', 'green-second': '-3/26', 'green-third': '0',
                'gradient-quadrance': '8/5', 'curl-quadrance': '0', 'harmonic-quadrance': '207/5'}
    if {f['id']: f['value'] for f in catalog['facts']} != expected:
        raise ValueError('Chromatin protocol differs from its declared finite result')
    readers = [dict(id='source-and-answers', workspace=str(manifest.relative_to(root)),
                    facts=str(facts_path.relative_to(root)), artifacts='rna', count=15),
               dict(id='expression-bounds', workspace=str(manifest.relative_to(root)),
                    facts=str(focused_path.relative_to(root)), artifacts='rna', count=3),
               dict(id='chromatin-structure', workspace=str(chrom_manifest.relative_to(root)),
                    facts=str(chrom_path.relative_to(root)), artifacts='chromatin', count=12)]
    save(root / 'readers.json', readers)
    return readers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    binary, root, fixtures = args.gadget.resolve(), args.out.resolve(), args.fixtures.resolve()
    readers = prepare(binary, root, fixtures)
    before = {}
    for reader in readers:
        before[reader['id']] = run(binary, 'fact-catalog', '--workspace', root / reader['workspace'],
                                  '--request', root / reader['facts'], '--artifacts', root / reader['artifacts'],
                                  '--view', 'model')
    moved = root.with_name(root.name + '-relocated')
    if moved.exists():
        raise ValueError('The relocated output directory must also be new')
    root.rename(moved)
    for reader in readers:
        after = run(binary, 'fact-catalog', '--workspace', moved / reader['workspace'],
                    '--request', moved / reader['facts'], '--artifacts', moved / reader['artifacts'],
                    '--view', 'model')
        if before[reader['id']] != after:
            raise ValueError('Reading view changed after moving the experiment')
    result = dict(status='PASS', selections={r['id']: r['count'] for r in readers},
                  relocated_views_equal=True, output=str(moved))
    save(moved / 'reader-verification.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
