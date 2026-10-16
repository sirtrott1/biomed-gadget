"""Actual publications, complete occurrence maps and native glued boundary lineage."""
import copy
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True); root.mkdir(parents=True)
example = Path(__file__).resolve().parents[1] / 'examples/bioinformatics/relational-gluing'
def module(name):
    spec = importlib.util.spec_from_file_location(name, example / (name + '.py'))
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result); return result
make, verify = module('make_example').make_example, module('verify').verify
checks = 0
def run(*args, fail=False):
    global checks
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    assert bool(result.returncode) == fail, (args, result.stdout, result.stderr)
    checks += 1
    return None if fail else json.loads(result.stdout)
def save(name, data):
    path = root / name; path.write_text(json.dumps(data, indent=2) + '\n'); return path
make(root / 'inputs')
run('import-assays', '--spec', root / 'inputs/import.json', '--out', root / 'workspace')
manifest = root / 'workspace/workspace.json'; query = json.loads((root / 'inputs/query.json').read_text())
before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (root / 'workspace').rglob('*') if p.is_file()}
def execute(q, name='query.json', fail=False):
    return run('assay-query', '--workspace', manifest, '--request', save(name, q), fail=fail)
certificate = execute(query); result = certificate['result']
reference = json.loads((root / 'inputs/reference.json').read_text())
verified = verify(result, reference, manifest); assert verified['verified_evidence_spans'] > 0
if len(sys.argv) > 3:
    subprocess.run([sys.argv[3], str(manifest), str(root / 'inputs/query.json')], check=True)
# Complete source catalogues remain available independently; executing a native
# gluing action does not rebuild, seal or embed those presentation pages.
assert all('source_pages' not in b and 'source_index' not in b
           for b in result['source_bindings'].values())
for corrupt in ['source_digest', 'native_sha256', 'declaration_sha256']:
    forged = copy.deepcopy(result)
    binding = forged['source_bindings']['left']
    if corrupt == 'native_sha256': binding['native']['sha256'] = 'sha256:' + '0'*64
    else: binding[corrupt] = 'sha256:' + '0'*64
    try: verify(forged, reference, manifest)
    except AssertionError: pass
    else: raise AssertionError('Forged gluing evidence accepted: ' + corrupt)
serial = copy.deepcopy(query); serial['workers'] = 1
serial_result = execute(serial)['result']; serial_result.pop('execution'); parallel = copy.deepcopy(result); parallel.pop('execution')
assert serial_result == parallel
empty = copy.deepcopy(query); empty['sections'] = []
assert execute(empty)['result']['execution']['completed_tiles'] == 0
none = copy.deepcopy(query); none['boundary_policy'] = 'none'
verify(execute(none), json.loads((root / 'inputs/reference.json').read_text()))
for name, change in [
    ('missing-cell', lambda q: q['embeddings']['common_left']['cells'].pop()),
    ('missing-inactive', lambda q: q['embeddings']['common_left']['occurrences'].pop(4)),
    ('missing-zero', lambda q: q['embeddings']['common_right']['occurrences'].pop(5)),
    ('wrong-position', lambda q: q['embeddings']['common_right']['occurrences'][0].update(to='2')),
    ('duplicate-source', lambda q: q['sources']['left'].update(key='1')),
    ('wrong-key', lambda q: q['sources']['common'].update(key='01')),
    ('wrong-sequence', lambda q: q['sources']['right'].update(available_sequence='00')),
    ('missing-field', lambda q: q['sections'][0].update(right_field='all.absent')),
    ('wrong-units', lambda q: q['sections'][0].update(units='other')),
    ('wrong-grade', lambda q: q['sections'][0].update(right_field='all.edges-0')),
    ('duplicate-section', lambda q: q['sections'].append(q['sections'][0])),
    ('bad-policy', lambda q: q.update(boundary_policy='average')),
    ('zero-workers', lambda q: q.update(workers=0)),
]:
    q = copy.deepcopy(query); change(q); execute(q, name + '.json', fail=True)
# Public default recipe, custom plan and model proposed structured query routes.
prepared = run('pipeline-prepare', '--workspace', manifest, '--request', root / 'inputs/pipeline.json')
plan = save('prepared.json', prepared)
index = run('pipeline', '--workspace', manifest, '--request', plan, '--out', root / 'pipeline')
assert len(index['results']) == 1
verify(json.loads((root / 'pipeline/glue/result.json').read_text()), json.loads((root / 'inputs/reference.json').read_text()))
proposal = dict(schema='gadget.pipeline.proposal.v1', action='pipeline', source_digest=prepared['source_digest'],
                profile='custom', question='Execute the declared gluing diagram',
                steps=[dict(id='glued', type='query', input='source', depends_on=[], request=query)])
admitted = run('pipeline-proposal', '--workspace', manifest, '--request', save('proposal.json', proposal), '--out', root / 'proposed')
assert admitted['status'] == 'admitted'
# Typed scientific answer admission must link to this actual producing method.
intent = dict(schema='gadget.pipeline.intent.v1', profile='multi_omics',
              recipes=[dict(id='glue', purpose='Combine the authored source diagram',
                            goal=dict(method='relational_gluing', parameters={k: v for k, v in query.items() if k not in ['schema', 'operation']}))],
              answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',
                                   question='Which declared sections glue?',
                                   requirements=[dict(id='structure', observable='relational_gluing', bindings=['glue'])]))
intent_path = save('intent.json', intent)
catalog = run('recipe-catalog', '--workspace', manifest, '--intent', intent_path)
answer = dict(schema='gadget.pipeline.recipe.proposal.v1', action='recipes',
              source_digest=catalog['source_digest'], intent_sha256=catalog['intent_sha256'],
              goals=['glue'])
typed = run('recipe-proposal', '--workspace', manifest, '--intent', intent_path,
            '--request', save('typed-proposal.json', answer), '--prompt', intent['answer_contract']['question'],
            '--out', root / 'typed')
assert typed['status'] == 'admitted'
assert {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (root / 'workspace').rglob('*') if p.is_file()} == before
# A whole source carrier larger than both historical 64 and 1024 defaults.
make(root / 'large-inputs', vertices=2049, repetitions=2)
run('import-assays', '--spec', root / 'large-inputs/import.json', '--out', root / 'large-workspace')
large = run('assay-query', '--workspace', root / 'large-workspace/workspace.json', '--request', root / 'large-inputs/query.json')
verify(large, json.loads((root / 'large-inputs/reference.json').read_text()))
assert large['result']['execution']['dimension_limit'] == 4100
# Scale the original incidence fibres, higher fillings and arity three cells,
# rather than only adding isolated vertices to a fixed tiny topology.
make(root / 'incidence-inputs', vertices=34, repetitions=2, modules=17)
run('import-assays', '--spec', root / 'incidence-inputs/import.json', '--out', root / 'incidence-workspace')
incidence_manifest = root / 'incidence-workspace/workspace.json'
incidence = run('assay-query', '--workspace', incidence_manifest, '--request', root / 'incidence-inputs/query.json')
verified_incidence = verify(incidence, json.loads((root / 'incidence-inputs/reference.json').read_text()), incidence_manifest)
assert verified_incidence['cells'] == 221 and verified_incidence['occurrences'] == 442
assert incidence['result']['execution']['dimension_limit'] == 136
# Actual measurement disagreement and incomplete support refuse native gluing.
for name, word in [('disagreement', '999'), ('unsupported', 'NA')]:
    inputs = root / (name + '-inputs'); shutil.copytree(root / 'inputs', inputs)
    path = inputs / 'right-signal.tsv'; rows = path.read_text().splitlines()
    for i, row in enumerate(rows):
        if row.startswith('right-shared-0\t'): rows[i] = 'right-shared-0\t' + word
    path.write_text('\n'.join(rows) + '\n')
    spec = json.loads((inputs / 'import.json').read_text())
    spec['assays'][2]['source']['sha256'] = 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()
    (inputs / 'import.json').write_text(json.dumps(spec))
    run('import-assays', '--spec', inputs / 'import.json', '--out', root / (name + '-workspace'))
    run('assay-query', '--workspace', root / (name + '-workspace/workspace.json'), '--request', inputs / 'query.json', fail=True)
print(json.dumps(dict(status='PASS', process_checks=checks, large_carrier=2051,
                      incidence_modules=17, verified_evidence_spans=verified['verified_evidence_spans'])))
