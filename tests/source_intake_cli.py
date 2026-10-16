"""Source complete admission beyond previous ceilings and key range overlaps."""
import copy
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples/bioinformatics/source-intake'))
import replay
import verify

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True)
replay.experiments(binary, root)
report = verify.verify(root)
request = json.loads((root / 'database-fetch.json').read_text())
refusals = []
for name, edit in [
    ('duplicate-request', lambda q: q['requests'].append(q['requests'][0])),
    ('unsafe-identifier', lambda q: q['requests'][0].update(identifier='../escape')),
    ('foreign-provider-filter', lambda q: q['requests'][1].update(offset=1)),
    ('provider-page-ceiling', lambda q: q['requests'][0].update(limit=1001)),
    ('zero-resource', lambda q: q['transport'].update(max_body_bytes=0)),
    ('unknown-resource', lambda q: q['transport'].update(unknown_bytes=1)),
    ('negative-offset', lambda q: q['requests'][0].update(offset=-1))]:
    invalid = copy.deepcopy(request); edit(invalid)
    replay.run(binary, 'database-request', '--request', replay.save(root / (name + '.json'), invalid), success=False)
    refusals.append(name)
spec = json.loads((root / 'columns/import.json').read_text())
spec['source']['sha256'] = 'sha256:' + '0' * 64
bad = replay.save(root / 'columns/bad-import.json', spec)
replay.run(binary, 'import-table', '--spec', bad, '--out', root / 'bad-study', success=False)
assert not (root / 'bad-study').exists()
replay.run(binary, 'import-table', '--spec', root / 'columns/import.json', '--out', root / 'columns/study', success=False)
assert verify.verify(root) == report
report['refused_requests'] = refusals
# A lexical source identity is retained as a witness identity, not a field.
# Declaring its source column as a covariate must not create a dangling alias.
replay.table(binary, root / 'identity-channel', ['identity', 'group', 'value'],
             [['sample-a', 'a', '0'], ['sample-b', 'b', '1']],
             ['identity', 'group'], ['group'], identity='identity')
identity_spec = json.loads((root / 'identity-channel/import.json').read_text())
identity_spec['columns'][0]['channel'] = 'covariate'
replay.save(root / 'identity-channel/import.json', identity_spec)
replay.run(binary, 'import-table', '--spec', root / 'identity-channel/import.json',
           '--out', root / 'identity-channel/covariate-study')
declaration = json.loads((root / 'identity-channel/covariate-study/study.json').read_text())
assert [c['name'] for c in declaration['design']['covariates']] == ['group']
inspection = replay.run(binary, 'query', '--study', root / 'identity-channel/covariate-study/study.json',
                        '--query', replay.save(root / 'identity-channel/query.json', dict(operation='inspect')))
assert inspection['source_sequence'] == '0'
report['identity_channel_contract'] = 'PASS'
replay.save(root / 'verification.json', report)
print(json.dumps(report))
