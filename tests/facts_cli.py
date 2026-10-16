"""Literal facts, completed result ownership, falsification, omission and source sized admission."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

binary, root, helper = (Path(p).resolve() for p in sys.argv[1:])
repo = Path(__file__).resolve().parent.parent
shutil.rmtree(root, ignore_errors=True)
spec = importlib.util.spec_from_file_location('fact_example', repo / 'examples/bioinformatics/scientific-claims/replay.py')
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)
manifest, request_path, proposal_path, catalog = demo.prepare(binary, root, repo / 'examples/bioinformatics')
source = demo.load(repo / 'examples/bioinformatics/source-context/replay.py')
request = json.loads(request_path.read_text())
proposal = json.loads(proposal_path.read_text())
checks = rejections = 0


def check(value):
    global checks
    assert value
    checks += 1


def command(*args, fail=False):
    global rejections
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    check(bool(result.returncode) == fail)
    if fail:
        rejections += 1
        return result.stderr
    return json.loads(result.stdout)


def admit(value, fail=False):
    path = source.save(root / 'active-facts.json', value)
    return command('fact-catalog', '--workspace', manifest, '--request', path, '--artifacts', root, fail=fail)


def claim(value, fail=False):
    path = source.save(root / 'active-claims.json', value)
    return command('claim-proposal', '--workspace', manifest, '--facts', request_path,
                   '--request', path, '--artifacts', root, fail=fail)


expected = {'reported-coordinates': 2, 'unobserved-coordinates': 0, 'RNA-units': 'RNA molecules',
            'incidence-0-coefficient': '-1', 'incidence-1-coefficient': '1', 'incidence-2-coefficient': '0',
            'inactive-exists': False, 'transcript-A1-length': 10, 'read-b0-quality-sum': '120',
            'geneA-change-lower': '1', 'geneA-change-upper': '1', 'A1-change-direction': 'unresolved',
            'reported-f1-support': True, 'reported-f1-value': '3', 'missing-f1-support': False}
for fact in catalog['facts']:
    check(fact['value'] == expected[fact['id']])
    check(fact['fact_sha256'].startswith('sha256:'))
check(len(catalog['certificates']) == 4)
check(admit(request) == catalog)
native_facts=root/'selected-facts.rfcxs'
check(command('fact-catalog','--workspace',manifest,'--request',request_path,
    '--artifacts',root,'--out',native_facts)==catalog)
retained=command('record-view','--input',native_facts,'--member','result')
check(retained['catalog']==catalog)

view = command('fact-catalog', '--workspace', manifest, '--request', request_path,
               '--artifacts', root, '--view', 'model')
check(view['schema'] == 'gadget.scientific.facts.model.view.v2')
check(view['task'] == dict(action='report_values', comparison='equal', scope='all_selected_facts'))
check(view['fact_count'] == len(catalog['facts']))
check([f['id'] for f in view['facts']] == [f['id'] for f in catalog['facts']])
owners = {owner['owner_sha256']: owner for owner in view['owners']}
check(len(owners) == len(view['owners']))
check(len(owners) < len(view['facts']))
for shown, fact in zip(view['facts'], catalog['facts']):
    check({k: v for k, v in shown.items() if k != 'owner_sha256'} ==
          {k: v for k, v in fact.items() if k not in ['observed', 'owner']})
    owner = owners[shown['owner_sha256']]
    check(owner['kind'] == fact['from']['kind'])
    if owner['kind'] == 'source':
        check(owner['owner'] == fact['owner'])
    else:
        reference = owner['owner']['reference']
        check('query' not in reference and 'query_sha256' in reference)
        check({k: v for k, v in reference.items() if k != 'query_sha256'} ==
              {k: v for k, v in fact['owner']['reference'].items() if k != 'query'})
check(command('fact-catalog', '--workspace', manifest, '--request', request_path,
              '--artifacts', root, '--view', 'full') == catalog)
command('fact-catalog', '--workspace', manifest, '--request', request_path,
        '--artifacts', root, '--view', 'unknown', fail=True)
response = claim(proposal)
check(response['status'] == 'matched')
check(response['certificate']['matched'] == 15)
check(response['certificate']['interpretation_verified'] is False)
false = json.loads((root / 'contradicted-claims.json').read_text())
response = claim(false)
check(response['status'] == 'contradicted')
check(response['certificate']['contradicted'] == 4)
wrong = {c['claim']['fact']: c['comparison']['observed'] for c in response['certificate']['checks']
         if not c['comparison']['agrees']}
check(wrong == {key: expected[key] for key in ['read-b0-quality-sum', 'A1-change-direction', 'reported-coordinates', 'inactive-exists']})
incomplete = copy.deepcopy(proposal)
incomplete['claims'] = incomplete['claims'][:1]
response = claim(incomplete)
check(response['status'] == 'incomplete')
check(response['certificate']['matched'] == 1)
check(len(response['certificate']['missing']) == 14)
# Repeating a supplied fact does not satisfy a different omitted selection.
incomplete['claims'].append(dict(incomplete['claims'][0], id='second-same-fact'))
check(claim(incomplete)['status'] == 'incomplete')
clarify = {k: v for k, v in proposal.items() if k != 'claims'}
clarify.update(action='clarify', reason='No clinical severity measurement is present in these selected facts')
response = claim(clarify)
check(response['status'] == 'selection_required')
check(len(response['certificate']['missing']) == 15)

for member in ['schema', 'source_digest', 'facts_sha256', 'question_sha256']:
    invalid = copy.deepcopy(proposal)
    invalid[member] = 'wrong'
    claim(invalid, True)
for key, value in [('fact', 'clinical-severity'), ('value', 2.0), ('value', True), ('value', -1),
                   ('comparison', 'approximately_equal'), ('path', ['invented']), ('value_type', 'rational')]:
    invalid = copy.deepcopy(proposal)
    invalid['claims'][0][key] = value
    claim(invalid, True)
for change in ['duplicate', 'empty', 'object', 'extra', 'empty-id', 'text-order', 'boolean-order', 'rational-float', 'rational-object', 'rational-zero-denominator', 'text-null']:
    invalid = copy.deepcopy(proposal)
    if change == 'duplicate': invalid['claims'].append(invalid['claims'][0])
    elif change == 'empty': invalid['claims'] = []
    elif change == 'object': invalid['claims'] = {}
    elif change == 'extra': invalid['new_measurement'] = 'invented'
    elif change == 'empty-id': invalid['claims'][0]['id'] = ''
    elif change in ['text-order', 'boolean-order']:
        target = 'RNA-units' if change == 'text-order' else 'inactive-exists'
        next(c for c in invalid['claims'] if c['fact'] == target)['comparison'] = 'greater'
    elif change == 'rational-float': invalid['claims'][3]['value'] = -1.0
    elif change == 'rational-object': invalid['claims'][3]['value'] = {'rational':'-1','new_measurement':'invented'}
    elif change == 'rational-zero-denominator': invalid['claims'][3]['value'] = '1/0'
    elif change == 'text-null': next(c for c in invalid['claims'] if c['fact'] == 'RNA-units')['value'] = None
    claim(invalid, True)

for member in ['schema', 'source_digest', 'question']:
    invalid = copy.deepcopy(request)
    invalid[member] = 'wrong'
    admit(invalid, True)
for change in ['empty', 'object', 'duplicate', 'missing-context', 'missing-answers', 'invented', 'float-index',
               'negative-index', 'bad-type', 'missing-zero', 'wrong-question', 'absolute-path', 'traversal', 'unsupported-kind']:
    invalid = copy.deepcopy(request)
    if change == 'empty': invalid['facts'] = []
    elif change == 'object': invalid['facts'] = {}
    elif change == 'duplicate': invalid['facts'].append(invalid['facts'][0])
    elif change == 'missing-context': invalid.pop('context')
    elif change == 'missing-answers': invalid.pop('answers')
    elif change == 'invented': invalid['facts'][0]['path'] = ['invented']
    elif change == 'float-index': invalid['facts'][8]['path'][1] = 0.5
    elif change == 'negative-index': invalid['facts'][8]['path'][1] = -1
    elif change == 'bad-type': invalid['facts'][0]['value_type'] = 'boolean'
    elif change == 'missing-zero': invalid['facts'][-1]['path'] = ['observations', 0, 'value']; invalid['facts'][-1]['value_type'] = 'rational'
    elif change == 'wrong-question': invalid['context']['question'] = 'different'
    elif change == 'absolute-path': invalid['answers']['path'] = str(root / 'relocated-pipeline/answers.rfcxs')
    elif change == 'traversal': invalid['answers']['path'] = '../' + invalid['answers']['path']
    elif change == 'unsupported-kind': invalid['facts'][0]['from']['kind'] = 'predicted'
    admit(invalid, True)
for member, value in [('answer', 'wrong'), ('binding', 'wrong'), ('step', 'wrong'), ('result', 1), ('result', 0.5)]:
    invalid = copy.deepcopy(request)
    invalid['facts'][8]['from'][member] = value
    admit(invalid, True)

# Source only views and result only views retain distinct authorities.
source_only = copy.deepcopy(request)
source_only.pop('answers')
source_only['facts'] = source_only['facts'][:8]
check(len(admit(source_only)['facts']) == 8)
check(command('fact-catalog', '--workspace', manifest, '--request', root / 'active-facts.json',
              '--artifacts', root, '--view', 'model')['fact_count'] == 8)
result_only = copy.deepcopy(request)
result_only.pop('context')
result_only['facts'] = result_only['facts'][8:]
check(len(admit(result_only)['facts']) == 7)
check(command('fact-catalog', '--workspace', manifest, '--request', root / 'active-facts.json',
              '--artifacts', root, '--view', 'model')['fact_count'] == 7)
focused = copy.deepcopy(request)
focused['facts'] = [f for f in focused['facts'] if f['id'] in
                    ['geneA-change-lower', 'geneA-change-upper', 'A1-change-direction']]
focus_catalog = admit(focused)
focus_view = command('fact-catalog', '--workspace', manifest, '--request', root / 'active-facts.json',
                     '--artifacts', root, '--view', 'model')
check(focus_view['fact_count'] == 3 and len(focus_view['owners']) == 1)
check(focus_catalog['facts_sha256'] != catalog['facts_sha256'])
check([f['value'] for f in focus_view['facts']] == ['1', '1', 'unresolved'])
focus_proposal = dict(proposal, facts_sha256=focus_catalog['facts_sha256'],
                      claims=[c for c in proposal['claims'] if c['fact'] != 'A1-change-direction'
                              and c['fact'] in [f['id'] for f in focused['facts']]])
focus_response = command('claim-proposal', '--workspace', manifest, '--facts', root / 'active-facts.json',
                         '--request', source.save(root / 'focus-claims.json', focus_proposal), '--artifacts', root)
check(focus_response['status'] == 'incomplete')
check(focus_response['certificate']['missing'] == ['A1-change-direction'])
command('claim-proposal', '--workspace', manifest, '--facts', root / 'active-facts.json',
        '--request', proposal_path, '--artifacts', root, fail=True)
large = copy.deepcopy(source_only)
large['facts'] = [dict(copy.deepcopy(large['facts'][0]), id='fact-' + str(i)) for i in range(301)]
large_catalog = admit(large)
large_view = command('fact-catalog', '--workspace', manifest, '--request', root / 'active-facts.json',
                     '--artifacts', root, '--view', 'model')
check(large_view['fact_count'] == 301 and len(large_view['facts']) == 301)
check(len(large_view['owners']) == 1)
check(len({f['owner_sha256'] for f in large_view['facts']}) == 1)
large_proposal = dict(proposal, facts_sha256=large_catalog['facts_sha256'],
                      claims=[dict(id=f['id'], fact=f['id'], comparison='equal', value=f['value']) for f in large_catalog['facts']])
response = command('claim-proposal', '--workspace', manifest, '--facts', root / 'active-facts.json',
                   '--request', source.save(root / 'large-claims.json', large_proposal), '--artifacts', root)
check(response['certificate']['matched'] == 301)

# Altered outputs fail byte admission. A newly pinned answer file still has to pass its object seal.
answer_path = root / 'relocated-pipeline/answers.rfcxs'
answer_bytes = answer_path.read_bytes()
try:
    changed = command('record-view','--input',answer_path,'--member','result')
    changed['source_digest'] = 'sha256:' + '0' * 64
    forged=source.save(root/'altered-answer-view.json',dict(result=changed))
    temporary=root/'altered-answer.rfcxs'
    subprocess.run([str(binary),'native-record','--input',str(forged),'--out',str(temporary)],check=True)
    answer_path.write_bytes(temporary.read_bytes())
    invalid = copy.deepcopy(request)
    invalid['answers']['sha256'] = 'sha256:' + hashlib.sha256(answer_path.read_bytes()).hexdigest()
    admit(invalid, True)
finally:
    answer_path.write_bytes(answer_bytes)
result_path = root / 'relocated-pipeline/quality/result.rfcxs'
result_bytes = result_path.read_bytes()
try:
    result_path.write_bytes(result_bytes + b'\n')
    admit(request, True)
finally:
    result_path.write_bytes(result_bytes)

p = subprocess.run([str(helper), str(manifest), str(request_path), str(proposal_path), str(root)],
                   capture_output=True, text=True, timeout=180)
if p.returncode:
    raise RuntimeError(p.stderr)
check(p.returncode == 0)
checks += json.loads(p.stdout)['checks']
check(admit(request) == catalog)

# A separate chromatin protocol uses the same result selectors and claim checker.
chromatin = root / 'chromatin'
p = subprocess.run([sys.executable, str(repo / 'examples/bioinformatics/replay-answers.py'), '--gadget', str(binary),
                    '--protocol', str(repo / 'examples/bioinformatics/controlled-coverage'), '--out', str(chromatin)],
                   text=True, capture_output=True, timeout=180)
check(p.returncode == 0)
chrom_answer_path = chromatin / 'relocated-pipeline/answers.rfcxs'
chrom_answers = command('record-view','--input',chrom_answer_path,'--member','result')
chrom_intent = json.loads((chromatin / 'inputs/intent.json').read_text())
first = chrom_answers['answers'][0]
artifact = first['artifacts'][0]
reference = artifact['results'][0]
selection = reference.get('native_selection')
selectors = ['--assay', selection['assay']] if selection else []
certificate = command('record-view','--input',chrom_answer_path.parent/reference['path'], *selectors)
if certificate.get('schema')=='gadget.operation.native.v1':certificate=certificate['result']
# Select a complete declared field/result fact supported by this independent protocol.
def scalar_path(value, path=()):
    if isinstance(value, bool):return list(path), 'boolean', value
    if isinstance(value, int) and value >= 0:return list(path), 'count', value
    if isinstance(value, str):return list(path), 'text', value
    items = value.items() if isinstance(value, dict) else enumerate(value) if isinstance(value, list) else []
    for k, v in items:
        answer = scalar_path(v, (*path, k))
        if answer:return answer
path, type, observed = scalar_path(certificate['result'])
chrom_request = dict(schema=request['schema'], source_digest=chrom_answers['source_digest'],
                     question=chrom_intent['answer_contract']['question'],
                     answers=dict(path='relocated-pipeline/answers.rfcxs', sha256=source.digest(chrom_answer_path.read_bytes())),
                     facts=[dict(id='chromatin-result', meaning=first['observable'], value_type=type, path=path,
                                 **{'from': dict(kind='answer', answer=first['id'], binding=artifact['binding'],
                                                 step=artifact['step'], result=0)})])
chrom_manifest = chromatin / 'relocated-source/workspace.json'
chrom_catalog = command('fact-catalog', '--workspace', chrom_manifest,
                        '--request', source.save(chromatin / 'facts.json', chrom_request), '--artifacts', chromatin)
check(chrom_catalog['facts'][0]['value'] == observed)
chrom_view = command('fact-catalog', '--workspace', chrom_manifest, '--request', chromatin / 'facts.json',
                     '--artifacts', chromatin, '--view', 'model')
check(chrom_view['fact_count'] == 1 and chrom_view['facts'][0]['value'] == observed)
check(chrom_view['owners'][0]['kind'] == 'answer')
chrom_proposal = dict(schema=proposal['schema'], source_digest=chrom_catalog['source_digest'],
                      facts_sha256=chrom_catalog['facts_sha256'], question_sha256=chrom_catalog['question_sha256'],
                      action='claims', claims=[dict(id='reading', fact='chromatin-result', comparison='equal', value=observed)])
check(command('claim-proposal', '--workspace', chrom_manifest, '--facts', chromatin / 'facts.json',
              '--request', source.save(chromatin / 'claims.json', chrom_proposal), '--artifacts', chromatin)['status'] == 'matched')
source.save(root / 'fact-contract-verification.json', dict(status='PASS', checks=checks, rejected=rejections,
             selected_facts=15, captured_certificates=4, large_facts=301, independent_profiles=['RNA', 'chromatin']))
print(json.dumps(dict(status='PASS', checks=checks, rejected=rejections)))
