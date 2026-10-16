"""Native only replay after all compatibility views and external inputs disappear."""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

binary, fixture, root = map(Path, sys.argv[1:])
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)

def call(*args, view=True):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, (args, result.stderr)
    return json.loads(result.stdout) if view else None

def native(name, value):
    path = root / (name + '.json')
    path.write_text(json.dumps(value))
    target = root / (name + '.rfcxs')
    call('native-record', '--input', path, '--out', target, view=False)
    path.unlink()
    return target

shutil.copytree(fixture, root / 'input')
call('snapshot', '--input', root/'input', '--measurements', root/'input/measurements.json', '--out', root/'snapshot')
assert {p.name for p in (root/'snapshot').iterdir()} == {'study.rfcxs'}
shutil.rmtree(root/'input')
source = call('audit', '--snapshot', root/'snapshot')
initial = call('store-init', '--snapshot', root/'snapshot', '--out', root/'store')
request = native('publish', dict(schema='gadget.store.publish.v1', transaction='1', expected=initial['root']))
call('store-publish', '--store', root/'store', '--snapshot', root/'snapshot', '--request', request)
assert call('audit', '--store', root/'store')['source_digest'] == source['source_digest']

raw = root/'originals'
raw.mkdir()
def original(name, text):
    (raw/name).write_text(text)
    return dict(path=name, sha256='sha256:'+hashlib.sha256(text.encode()).hexdigest())
features=original('features.tsv', 'id\ng1\ng2\ng3\n')
observations=original('observations.tsv', 'id\nA\nB\n')
matrix=original('counts.tsv', 'id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\n')
reads=original('reads.fastq', '@r1\nACGT\n+\nIIII\n')
spec=dict(schema='gadget.assays.import.v1', id='native-suite',
    design=dict(question='Retain every original count', population='authored cells', selection='all', limitations=['fixture']),
    assays=[dict(id='counts', kind='expression', format='dense_tsv', source=matrix, features=features,
        observations=observations, source_namespace='881100', carrier='features', quantity='read_count', units='counts',
        channel='E', method='authored counts', arithmetic='exact_decimal', implicit_entries='measured_zero', missing_tokens=[])],
    resources=[dict(id='reads', format='fastq', source=reads, method='authored reads')], links=[])
(raw/'import.json').write_text(json.dumps(spec))
call('import-assays', '--spec', raw/'import.json', '--out', root/'assays')
# Only the canonical owning publication is carried to another directory.
standalone=root/'relocated'
standalone.mkdir()
shutil.copyfile(root/'assays/workspace.rfcxs', standalone/'workspace.bin')
shutil.rmtree(raw)
shutil.rmtree(root/'assays')
workspace=standalone/'workspace.bin'
query=native('query', dict(schema='gadget.assays.query.v1', operation='native_query', assay='counts', query=dict(operation='field_summary', field='all.A', support='reported_only')))
result=call('assay-query', '--workspace', workspace, '--request', query)
assert result['result']['answers'][0]['certificate']['result']['mean']['rational']=='3'
physical=native('physical', dict(schema='gadget.assays.query.v1', operation='read_quality', resource='reads', quality_thresholds=[0,20,40,93]))
quality=call('assay-query', '--workspace', workspace, '--request', physical)
assert quality['schema']=='gadget.assay.result.v1'
default=native('defaults', dict(schema='gadget.pipeline.default.v1', profile='rna_seq', question='Check measured support', goals=['coverage'], assay='counts', fields=['all.A']))
index=call('pipeline', '--workspace', workspace, '--request', default, '--out', root/'pipeline')
assert index['status']=='complete'
assert index['results'][0]['result'].endswith('.rfcxs')
assert (root/'pipeline/pipeline.rfcxs').is_file()
print('PASS: native-only source/store, relocated assay, original sequence bytes and native pipeline request/results')

# Preserve the full unsigned JSON integer, rather than int64 saturation.
large=root/'unsigned.json'
large.write_text('18446744073709551615')
call('native-record','--input',large,'--out',root/'unsigned.rfcxs',view=False)
assert call('record-view','--input',root/'unsigned.rfcxs')=='18446744073709551615'
