import copy
import json
import shutil
import subprocess
import sys
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo/'examples/bioinformatics/linear-observations'))
import replay
import method_answers

binary, producer, session, root = [Path(x).resolve() for x in sys.argv[1:]]
shutil.rmtree(root, ignore_errors=True)
records = replay.examples(binary, producer, root)
base = root/'unfilled-square'
spec = json.loads((base/'spec.json').read_text())
spec['operators'].append(replay.sources.op('full', 3, 'accession', 10, 10, replay.identity(4)))
spec['operators'].append(replay.sources.op('zero_map', 4, 'accession', 10, 20, dict(rows=5,columns=4,entries=[])))
replay.save(base/'spec.json', spec)
shutil.rmtree(base/'study')
subprocess.run([str(producer), str(base/'spec.json'), str(base/'study')], check=True)
q = json.loads((base/'query.json').read_text())
refusals = []
for change in [lambda x: x.update(source_sequence='6'), lambda x: x.update(source_sequence='07'),
               lambda x: x.update(publication='unknown'), lambda x: x.update(fields=['partial']),
               lambda x: x.update(fields=['category']), lambda x: x.update(fields=[]),
               lambda x: x.update(workers=0), lambda x: x.update(accession='R'),
               lambda x: x.update(operator='S'),
               lambda x: x.update(ordering=dict(kind='declared', rows=[0, 1, 2, 2], columns=[0, 1, 2, 3])),
               lambda x: x.update(ordering=dict(kind='declared', rows=[0, 1], columns=[0, 1])),
               lambda x: x.update(observation=dict(kind='coordinates', coordinates=[0, 1, 2, 3, 0]))]:
    changed = copy.deepcopy(q)
    change(changed)
    path = replay.save(root/'refused-query.json', changed)
    p = subprocess.run([str(binary), 'query', '--study', str(base/'study/study.json'),
                        '--query', str(path)], capture_output=True, text=True)
    assert p.returncode, changed
    refusals.append(p.stderr.strip())
p = subprocess.run([str(session), str(base)], capture_output=True, text=True)
assert p.returncode == 0, p.stderr
replay.save(root/'session.json', json.loads(p.stdout))
selected = root/'ordered-coordinates'
spec = json.loads((selected/'spec.json').read_text())
spec['operators'].append(replay.sources.op('full', 3, 'accession', 10, 10, replay.identity(4)))
spec['operators'].append(replay.sources.op('zero_map', 4, 'accession', 10, 20, dict(rows=6, columns=4, entries=[])))
replay.save(selected/'spec.json', spec)
shutil.rmtree(selected/'study')
subprocess.run([str(producer), str(selected/'spec.json'), str(selected/'study')], check=True)
selected_session = subprocess.run([str(session), str(selected)], capture_output=True, text=True)
assert selected_session.returncode == 0, selected_session.stderr
replay.save(root/'selected-session.json', json.loads(selected_session.stdout))
empty = root/'empty-output'
changed = json.loads((empty/'query.json').read_text()); changed['fields'] = ['partial']
refused = subprocess.run([str(binary), 'query', '--study', str(empty/'study/study.json'),
    '--query', str(replay.save(root/'empty-refused-query.json', changed))], capture_output=True, text=True)
assert refused.returncode, 'empty output ignored missing source input'
# Equal local keys do not identify independently admitted publications.
coordinator = json.loads((base/'study/study.json').read_text())
coordinator['native']['path'] = 'unfilled-square/study/'+coordinator['native']['path']
for doc in coordinator['documents']:
    doc['path'] = 'unfilled-square/study/'+doc['path']
coordinator['publications'] = []
for name in ['unfilled-square', 'filled-square']:
    path = root/name/'study/study.json'
    inspected = replay.run(binary, 'query', '--study', path,
                           '--query', replay.save(root/'inspect.json', dict(operation='inspect')))
    coordinator['publications'].append(dict(id=name, study=dict(path=name+'/study/study.json', sha256=replay.sources.digest(path)),
        source_digest=inspected['source_digest'], meaning='Actual independently admitted publication',
        evidence=coordinator['operators'][0]['evidence']))
replay.save(root/'coordinator.json', coordinator)
for name in ['unfilled-square', 'filled-square']:
    query = json.loads((root/name/'query.json').read_text())
    query['publication'] = name
    certificate = replay.run(binary, 'query', '--study', root/'coordinator.json',
                             '--query', replay.save(root/'foreign-query.json', query))
    expected = replay.run(binary, 'query', '--study', root/name/'study/study.json',
                          '--query', replay.save(root/'self-query.json', {**query, 'publication': 'self'}))
    assert certificate['result'] == expected['result']
    assert certificate['field_semantics'] == expected['field_semantics']
# The actual assay-native path, default/custom recipe, DAG and complete answer artifact.
pipeline = root/'pipeline'; pipeline.mkdir()
features = pipeline/'features.tsv'; features.write_text('id\ng1\ng2\ng3\ng4\n')
observations = pipeline/'observations.tsv'; observations.write_text('id\nA\nB\n')
counts = pipeline/'counts.tsv'; counts.write_text('id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\ng4\t1\t2\n')
record = lambda p: dict(path=p.name, sha256=replay.sources.digest(p))
assay = dict(id='counts', kind='expression', format='dense_tsv', source=record(counts), features=record(features),
    observations=record(observations), source_namespace='931100', carrier='features', quantity='read_count', units='counts',
    channel='E', method='Authored counts', arithmetic='exact_decimal', implicit_entries='measured_zero', missing_tokens=[])
spec = dict(schema='gadget.assays.import.v1', id='observation-pipeline', design=dict(question='Observe complete source fields',
    population='authored four-feature source', selection='all', limitations=['Algebra qualification']), assays=[assay], resources=[], links=[])
original = pipeline/'imported'
replay.run(binary, 'import-assays', '--spec', replay.save(pipeline/'import.json', spec), '--out', original)
# Extend the imported actual owner, preserving its fields and original documents.
addition = dict(id='observation-pipeline-declared', namespace='931100', sequence='0',
    carriers=[dict(key='10', coordinate_system='100', size=4)],
    operators=[replay.sources.op('R', 20001, 'resolving', 1, 1, replay.identity(4)),
               replay.sources.op('S', 20002, 'accession', 1, 10, replay.identity(4))], fields=[])
subprocess.run([str(producer), str(replay.save(pipeline/'addition.json', addition)), str(pipeline/'extended'),
                str(original/'assays/0/partitions/0/study.json')], check=True, capture_output=True, text=True)
envelope = json.loads((original/'workspace.json').read_text())
envelope['assays'][0]['partitions'][0]['sha256'] = replay.sources.digest(pipeline/'extended/study.json')
# Workspace artifacts must stay below their declared root.
shutil.copytree(pipeline/'extended', original/'extended')
envelope['assays'][0]['partitions'][0]['path'] = 'extended/study.json'
envelope['assays'][0]['partitions'][0]['source_digest'] = replay.run(binary, 'query',
    '--study', original/'extended/study.json', '--query', replay.save(pipeline/'inspect.json', dict(operation='inspect')))['source_digest']
envelope.pop('manifest_sha256')
# json-c defines the public envelope byte representation; reseal this authored extension.
import ctypes
import hashlib
codec = ctypes.CDLL('libjson-c.so.5')
codec.json_tokener_parse.argtypes = [ctypes.c_char_p]
codec.json_tokener_parse.restype = ctypes.c_void_p
codec.json_object_to_json_string_ext.argtypes = [ctypes.c_void_p, ctypes.c_int]
codec.json_object_to_json_string_ext.restype = ctypes.c_char_p
codec.json_object_put.argtypes = [ctypes.c_void_p]
value = codec.json_tokener_parse(json.dumps(envelope).encode())
encoded = codec.json_object_to_json_string_ext(value, 2 | 16)
envelope['manifest_sha256'] = 'sha256:' + hashlib.sha256(encoded).hexdigest()
codec.json_object_put(value)
replay.save(original/'workspace.json', envelope)
q_assay = dict(operation='linear_response_observation', publication='self', source_sequence='0',
    operator='R', accession='S', fields=['all.A', 'all.B'], ordering=dict(kind='source'),
    observation=dict(kind='linear'), workers=2)
wrapped = dict(schema='gadget.assays.query.v1', operation='native_query', assay='counts', query=q_assay)
result = replay.run(binary, 'assay-query', '--workspace', original/'workspace.json',
                    '--request', replay.save(pipeline/'query.json', wrapped))
method_answers.MEANINGS['linear_response_observation'] = 'linear_response_observation'
answer_checks = method_answers.qualify(binary, pipeline/'answers', original/'workspace.json', {'observations': (wrapped, result)})
replay.save(root/'contracts.json', dict(status='PASS', independent=records, refusals=refusals,
                                     session=json.loads(p.stdout), answer_checks=answer_checks))
print(json.dumps(dict(status="PASS", cases=len(records), refusals=len(refusals), answer_checks=answer_checks)))

# Reopen a result bound to a referenced publication after admission through a coordinator.
foreign_query = json.loads((root/'unfilled-square/query.json').read_text())
foreign_query['publication'] = 'unfilled-square'
replay.save(root/'foreign-native-query.json', foreign_query)
foreign_output = root/'foreign-native.rfcxs'
cli = subprocess.run([str(binary), 'query', '--study', str(root/'coordinator.json'),
                      '--query', str(root/'foreign-native-query.json'), '--out', str(foreign_output)],
                     capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr
cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(foreign_output),
                      '--out-record', str(root/'foreign-replay.rfcxs')], capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr

native_records = []
for case in sorted(root.iterdir()):
    if not (case/'query.json').is_file() or not (case/'study/study.json').is_file():
        continue
    checked = subprocess.run([str(session), str(case), 'native'], capture_output=True, text=True)
    assert checked.returncode == 0, (case.name, checked.stderr)
    native_records.append(dict(case=case.name, **json.loads(checked.stdout)))
    replayed = case/'native/cli-replay.rfcxs'
    cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(case/'native/result.rfcxs'),
                          '--out-record', str(replayed)], capture_output=True, text=True)
    assert cli.returncode == 0, cli.stderr
    changed = case/'native/changed-field.rfcxs'
    if changed.exists():
        refused = case/'native/refused.rfcxs'
        cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(changed),
                              '--out-record', str(refused)], capture_output=True, text=True)
        assert cli.returncode != 0 and not refused.exists(), case.name
replay.save(root/'native-contracts.json', native_records)
print(json.dumps(dict(status='PASS', native_cases=len(native_records), checks=sum(x['checks'] for x in native_records))))

# A single native plan retains different providers over one explicitly declared source.
mixed = root/'mixed'; mixed.mkdir()
spec = json.loads((root/'unfilled-square/spec.json').read_text())
spec['carriers'].append(dict(key='30', grade=0, size=4))
replay.save(mixed/'spec.json', spec)
subprocess.run([str(producer), str(mixed/'spec.json'), str(mixed/'study')], check=True, capture_output=True, text=True)
inspected = replay.run(binary, 'query', '--study', mixed/'study/study.json',
                       '--query', replay.save(mixed/'inspect.json', dict(operation='inspect')))
observation_query = json.loads((root/'unfilled-square/query.json').read_text())
hodge_query = dict(operation='hodge_frame_decomposition', publication='self', source_sequence='7',
                   grade=1, grade_carriers=['30','10'], fields=['x','x2','zero'], metrics=dict(kind='unit'))
plan = dict(schema='gadget.scientific.plan.v1', id='mixed-native-fields', question='Retain sectors and response readings',
            source_digest=inspected['source_digest'], steps=[dict(id='sectors', query=hodge_query),
                                                           dict(id='readings', query=observation_query)])
cli = subprocess.run([str(binary), 'plan', '--study', str(mixed/'study/study.json'),
                      '--request', str(replay.save(mixed/'plan.json', plan)), '--out', str(mixed/'plan.rfcxs')],
                     capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr
cli = subprocess.run([str(binary), 'replay-scientific', '--record', str(mixed/'plan.rfcxs'),
                      '--out-record', str(mixed/'replay.rfcxs')], capture_output=True, text=True)
assert cli.returncode == 0, cli.stderr
