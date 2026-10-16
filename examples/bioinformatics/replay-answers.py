"""Replay a portable, researcher authored import and answer contract protocol."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    binary, protocol, root = args.gadget.resolve(), args.protocol.resolve(), args.out.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / 'inputs'
    inputs.mkdir()

    def save(name, value):
        path = root / name
        path.write_text(json.dumps(value, indent=2) + '\n')
        return path

    def run(*params):
        p = subprocess.run([str(binary), *map(str, params)], text=True, capture_output=True)
        if p.returncode:
            raise RuntimeError(p.stderr)
        return json.loads(p.stdout)

    for name in ['import.json', 'intent.json']:
        shutil.copyfile(protocol / name, inputs / name)
    declaration = json.loads((inputs / 'import.json').read_text())
    sources = {}

    def stage(value):
        if isinstance(value, list):
            for v in value:
                stage(v)
        elif isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                path = Path(value['path'])
                if path.is_absolute() or '..' in path.parts:
                    raise ValueError('Portable protocol source must be within its input directory')
                data = (protocol / path).read_bytes()
                if value['sha256'] != 'sha256:' + hashlib.sha256(data).hexdigest():
                    raise ValueError('Protocol input hash mismatch: ' + str(path))
                destination = inputs / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(data)
                sources[str(path)] = value['sha256']
            else:
                for v in value.values():
                    stage(v)

    stage(declaration)
    run('import-assays', '--spec', inputs / 'import.json', '--out', root / 'source')
    manifest = root / 'source/workspace.rfcxs'
    intent = json.loads((inputs / 'intent.json').read_text())
    question = intent['answer_contract']['question']
    catalog = run('recipe-catalog', '--workspace', manifest, '--intent', inputs / 'intent.json')
    if catalog['scientific_admission']['status'] != 'satisfied':
        raise ValueError('The authored protocol has an unavailable answer requirement')
    save('catalog.json', catalog)
    proposal = dict(schema='gadget.pipeline.recipe.proposal.v1', source_digest=catalog['source_digest'],
                    intent_sha256=catalog['intent_sha256'], action='recipes',
                    goals=[r['id'] for r in intent['recipes']])
    admission = run('recipe-proposal', '--workspace', manifest, '--intent', inputs / 'intent.json',
                    '--request', save('proposal.json', proposal), '--prompt', question,
                    '--out', root / 'pipeline')
    if admission['status'] != 'admitted':
        raise ValueError('The authored selection did not satisfy its answer contract')
    save('admission.json', admission)
    shutil.move(root / 'source', root / 'relocated-source')
    shutil.move(root / 'pipeline', root / 'relocated-pipeline')
    directory = root / 'relocated-pipeline'
    index = admission['execution']
    answer_bytes = (directory / index['answers']).read_bytes()
    if index['answers_sha256'] != 'sha256:' + hashlib.sha256(answer_bytes).hexdigest():
        raise ValueError('Answer receipt hash mismatch')
    answers = run('record-view','--input',directory/index['answers'],'--member','result')
    certificates = 0
    for answer in answers['answers']:
        for artifact in answer['artifacts']:
            for ref in artifact['results']:
                path = Path(ref['path'])
                if path.is_absolute() or '..' in path.parts:
                    raise ValueError('Answer artifact escapes the pipeline')
                data = (directory / path).read_bytes()
                selection = ref.get('native_selection')
                selectors = ['--assay', selection['assay']] if selection else []
                certificate = run('record-view','--input',directory/path, *selectors)
                if certificate.get('schema')=='gadget.operation.native.v1':certificate=certificate['result']
                if (ref['sha256'] != 'sha256:' + hashlib.sha256(data).hexdigest()
                    or ref['source_digest'] != certificate['source_digest']
                    or ref['query'] != certificate['query']
                    or ref['certificate_sha256'] != certificate['certificate_sha256']
                    or ref['pointer'] != '/result'):
                    raise ValueError('Complete answer certificate mismatch')
                certificates += 1
    replays = []
    for binding in intent['recipes']:
        goal = binding['goal']
        if 'output' not in goal:
            continue
        name, assay = binding['id'], goal['output']['id']
        replay = run('replay-derivation', '--workspace', directory / name / 'workspace/workspace.rfcxs',
                     '--assay', assay, '--out', root / ('replayed-' + name))
        if replay['status'] != 'PASS':
            raise ValueError('Derivation replay failed: ' + name)
        save('replay-' + name + '.json', replay)
        replays.append(name)
    verification = dict(status='PASS', source_digest=catalog['source_digest'],
                        intent_sha256=catalog['intent_sha256'], input_files=sources,
                        answers=len(answers['answers']), complete_certificates=certificates,
                        regenerated_derivations=replays)
    save('verification.json', verification)
    print(json.dumps(verification))


if __name__ == '__main__':
    main()
