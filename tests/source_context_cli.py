"""Source owner, annotation, typed mapping and detached execution contracts."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

binary, root, helper = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve()
repo = Path(__file__).resolve().parent.parent
shutil.rmtree(root, ignore_errors=True)
module = importlib.util.spec_from_file_location('source_context_example', repo / 'examples/bioinformatics/source-context/replay.py')
demo = importlib.util.module_from_spec(module)
module.loader.exec_module(demo)
manifest, intent_path, selection_path = demo.prepare(binary, root, repo / 'examples/bioinformatics')
checks = 0
rejections = 0


def check(condition):
    global checks
    assert condition
    checks += 1


def save(name, value):
    return demo.save(root / name, value)


def command(*args, fail=False):
    global rejections
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    check(bool(result.returncode) == fail)
    if fail:
        rejections += 1
        return result.stderr
    return json.loads(result.stdout)


def context_admit(request, fail=False, owner=manifest):
    return command('source-context', '--workspace', owner, '--request', save('active-context.json', request), fail=fail)


def annotate(proposal, fail=False):
    return command('annotation-proposal', '--workspace', manifest, '--context', root / 'context.json',
                   '--request', save('active-annotations.json', proposal), fail=fail)


context = json.loads((root / 'context.json').read_text())
catalog = json.loads((root / 'context-catalog.json').read_text())
proposal = json.loads((root / 'annotations.json').read_text())
intent = json.loads(intent_path.read_text())
check(context_admit(context) == catalog)
targets = {t['id']: t for t in catalog['targets']}
for target in targets.values():
    check(target['evidence_verified'])
    for evidence in target['resolved']['evidence']:
        check(evidence['document_sha256'].startswith('sha256:'))
check(targets['reported-RNA']['resolved']['descriptor']['reported'] == 2)
check(targets['missing-RNA']['resolved']['descriptor']['reported'] == 0)
check(targets['reported-RNA']['resolved']['descriptor']['channel'] == 'E')
check(targets['reported-RNA']['resolved']['descriptor']['units'] == 'RNA molecules')
check(targets['incidence-0']['resolved']['descriptor']['coefficient'] == '-1')
check(targets['incidence-1']['resolved']['descriptor']['coefficient'] == '1')
check(targets['incidence-2']['resolved']['descriptor']['coefficient'] == '0')
check(not targets['incidence-2']['resolved']['descriptor']['aspects']['exists'])
check(targets['incidence-0']['resolved']['descriptor']['lower'] == targets['incidence-1']['resolved']['descriptor']['lower'])
check(targets['incidence-0']['target']['key'] != targets['incidence-1']['target']['key'])

# Native API ownership, mutable-view separation, every retained occurrence,
# UTF-8 split rejection and the removed legacy mention ceiling.
result = subprocess.run([str(helper), 'native', str(root / 'native-fixture')], capture_output=True, text=True, timeout=180)
check(result.returncode == 0)
checks += json.loads(result.stdout)['checks']

for member in ['source_digest', 'schema']:
    value = copy.deepcopy(context)
    value[member] = 'invalid'
    context_admit(value, True)
for field, value in [('targets', []), ('targets', {}), ('question', ''), ('question', 'bad\0question')]:
    invalid = copy.deepcopy(context)
    invalid[field] = value
    context_admit(invalid, True)
invalid = copy.deepcopy(context)
invalid['targets'].append(invalid['targets'][0])
context_admit(invalid, True)
for name in ['feature-A', 'reported-RNA', 'RNA-context']:
    invalid = copy.deepcopy(context)
    next(t for t in invalid['targets'] if t['id'] == name)['target']['id'] = 'invented'
    context_admit(invalid, True)
for member, value in [('assay', 'missing'), ('partition', 1), ('source_digest', 'sha256:' + '0' * 64),
                      ('partition', -1), ('partition', 0.5), ('map', {'invented': True})]:
    invalid = copy.deepcopy(context)
    invalid['targets'][0]['target'][member] = value
    context_admit(invalid, True)
for member, value in [('key', '999'), ('grade', 0), ('grade', 65536)]:
    invalid = copy.deepcopy(context)
    next(t for t in invalid['targets'] if t['id'] == 'incidence-0')['target'][member] = value
    context_admit(invalid, True)
for member, value in [('document', 'missing'), ('text', 'invented'), ('start', -1), ('end', 0), ('end', 100000), ('start', 0.5)]:
    invalid = copy.deepcopy(context)
    invalid['targets'][0]['evidence'][0][member] = value
    context_admit(invalid, True)
invalid = copy.deepcopy(context)
invalid['targets'][0]['evidence'] = []
context_admit(invalid, True)
for member, value in [('sha256', 'sha256:' + '0' * 64), ('record', 1), ('record', 900), ('artifact', 'reference')]:
    invalid = copy.deepcopy(context)
    next(t for t in invalid['targets'] if t['id'] == 'transcript-A1')['target'][member] = value
    context_admit(invalid, True)
invalid = copy.deepcopy(context)
item = next(t for t in invalid['targets'] if t['id'] == 'transcript-A1')
item['evidence'] = copy.deepcopy(next(t for t in invalid['targets'] if t['id'] == 'transcript-A2')['evidence'])
context_admit(invalid, True)

# Matching labels in separate native namespaces do not authorize owner reuse.
spec = json.loads((root / 'inputs/import.json').read_text())
second = copy.deepcopy(spec['assays'][0])
second['id'] = 'another-rna'
second['source_namespace'] = '975400'
spec['assays'].append(second)
spec_path = demo.save(root / 'inputs/two-owners.json', spec)
two = command('import-assays', '--spec', spec_path, '--out', root / 'two-owners')
invalid = copy.deepcopy(context)
invalid['source_digest'] = two['manifest_sha256']
invalid['targets'][0]['target']['assay'] = 'another-rna'
context_admit(invalid, True, root / 'two-owners/workspace.json')
# Explicitly selecting the second owner's actual digest is valid and remains separate.
valid = copy.deepcopy(invalid)
valid['targets'] = valid['targets'][:1]
valid['targets'][0]['target']['source_digest'] = two['assays'][1]['partitions'][0]['source_digest']
second_context = context_admit(valid, owner=root / 'two-owners/workspace.json')
check(second_context['targets'][0]['resolved']['source_namespace'] == '975400')
check(second_context['targets'][0]['binding_sha256'] != catalog['targets'][0]['binding_sha256'])

# Explicit selected views can exceed the former catalog/mention ceilings.
large = copy.deepcopy(context)
large['targets'] = [dict(copy.deepcopy(context['targets'][0]), id='target-' + str(i)) for i in range(301)]
check(len(context_admit(large)['targets']) == 301)
large = copy.deepcopy(context)
large['targets'] = large['targets'][:1]
large['targets'][0]['evidence'] *= 10001
check(len(context_admit(large)['targets'][0]['resolved']['evidence']) == 10001)

resource_view = copy.deepcopy(context)
resource_view['targets'] = [dict(copy.deepcopy(context['targets'][-3]), id='resource-' + str(i)) for i in range(301)]
resource_catalog = context_admit(resource_view)
check(len(resource_catalog['targets']) == 301)
check(all(t['resolved'] == targets['transcript-A1']['resolved'] for t in resource_catalog['targets']))

large_proposal = copy.deepcopy(proposal)
large_proposal['mentions'] = [dict(copy.deepcopy(proposal['mentions'][0]), id='mention-' + str(i)) for i in range(301)]
check(len(annotate(large_proposal)['certificate']['mentions']) == 301)
long_question = copy.deepcopy(context)
long_question['question'] = 'Retain the full source-bound scientific question. ' * 200
check(context_admit(long_question)['question'] == long_question['question'])

response = annotate(proposal)
check(response['status'] == 'prospective')
check(all(m['interpretation_verified'] is False for m in response['certificate']['mentions']))
for key in ['source_digest', 'context_sha256', 'question_sha256']:
    invalid = copy.deepcopy(proposal)
    invalid[key] = 'sha256:' + '0' * 64
    annotate(invalid, True)
for key, value in [('target', 'invented'), ('polarity', 'probably'), ('interpretation', ''), ('operation', 'hodge')]:
    invalid = copy.deepcopy(proposal)
    invalid['mentions'][0][key] = value
    annotate(invalid, True)
invalid = copy.deepcopy(proposal)
invalid['mentions'].append(invalid['mentions'][0])
annotate(invalid, True)
invalid = copy.deepcopy(proposal)
invalid['mentions'] = []
annotate(invalid, True)
clarification = {k: v for k, v in proposal.items() if k != 'mentions'}
clarification.update(action='clarify', reason='The requested clinical severity is not declared in these sources.')
check(annotate(clarification)['status'] == 'clarification')
# Byte admission does not validate interpretation. Retain an intentionally false
# reading with interpretation_verified=false; the numerical source is unchanged.
false_reading = copy.deepcopy(proposal)
false_reading['mentions'][2]['interpretation'] = 'Every RNA coordinate is reported.'
false_response = annotate(false_reading)
check(false_response['status'] == 'prospective')
check(false_response['certificate']['mentions'][2]['interpretation_verified'] is False)
check(context_admit(context) == catalog)
save('false-reading-receipt.json', false_response)

for key in ['annotations', 'requirements']:
    invalid = copy.deepcopy(intent)
    if key == 'annotations':
        invalid['source_context'][key]['question_sha256'] = 'sha256:' + '0' * 64
    else:
        invalid['source_context'][key][0]['mentions'] = ['invented']
    command('recipe-catalog', '--workspace', manifest, '--intent', save('invalid-intent.json', invalid), fail=True)
for field, value in [('id', 'unknown-obligation'), ('mentions', []), ('mentions', ['transcript-A1', 'transcript-A1'])]:
    invalid = copy.deepcopy(intent)
    invalid['source_context']['requirements'][0][field] = value
    command('recipe-catalog', '--workspace', manifest, '--intent', save('invalid-intent.json', invalid), fail=True)
invalid = copy.deepcopy(intent)
invalid['source_context']['requirements'].append(invalid['source_context']['requirements'][0])
command('recipe-catalog', '--workspace', manifest, '--intent', save('invalid-intent.json', invalid), fail=True)
invalid = copy.deepcopy(intent)
invalid.pop('answer_contract')
command('recipe-catalog', '--workspace', manifest, '--intent', save('invalid-intent.json', invalid), fail=True)
result = subprocess.run([str(helper), 'pipeline', str(manifest), str(intent_path), str(selection_path),
                         str(root / 'forged-execution')], capture_output=True, text=True, timeout=180)
check(result.returncode == 0)
checks += json.loads(result.stdout)['checks']

# Tampering with original resource bytes is rejected before publication.
tampered = root / 'tampered-workspace'
shutil.copytree(manifest.parent, tampered)
original = next(r['original']['path'] for r in json.loads(manifest.read_text())['resources']
                if r['declaration']['id'] == 'baseline.fastq')
(tampered / original).write_bytes((tampered / original).read_bytes() + b'corruption')
context_admit(context, True, tampered / 'workspace.json')

report = demo.execute(binary, root, manifest, intent_path, selection_path)
check(report['status'] == 'PASS')
answers = command('record-view','--input',root/'relocated-pipeline/answers.rfcxs','--member','result')
check(len(answers['answers']) == 3)
for answer in answers['answers']:
    check(len(answer['source_mentions']) > 0)
    for artifact in answer['artifacts']:
        for reference in artifact['results']:
            data = (root / 'relocated-pipeline' / reference['path']).read_bytes()
            check(demo.digest(data) == reference['sha256'])
            selection = reference.get('native_selection')
            options = ['--assay', selection['assay']] if selection else []
            certificate = command('record-view','--input',root/'relocated-pipeline'/reference['path'], *options)
            if certificate.get('schema')=='gadget.operation.native.v1':certificate=certificate['result']
            check(certificate['certificate_sha256'] == reference['certificate_sha256'])
check(len(answers['source_context']['annotations']['mentions']) == 10)
report.update(checks=checks, rejected_cases=rejections, native_zero_and_inactive_occurrences='PASS',
              false_interpretation='retained prospectively without semantic verification',
              large_targets=301, large_mentions=301, large_question_bytes=len(long_question['question'].encode()), large_evidence_spans=10001)
save('verification.json', report)
print(json.dumps(report))
