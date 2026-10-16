"""Execute admitted native pattern requests and archived cold replays through the CLI."""
import json
from pathlib import Path
import subprocess
import sys

binary, consumer, root = map(Path, sys.argv[1:4])
subprocess.run([str(consumer)], check=True)
def run(*args):
    result = subprocess.run([str(binary), *map(str, args)], check=True,
                            capture_output=True, text=True)
    return json.loads(result.stdout)
for name in ('candidates', 'program', 'gluing'):
    view = run('pattern-query', '--study', root/'study.rfcxs', '--model', root/'model.rfcxs',
               '--request', root/(name+'.json'), '--out', root/(name+'-cli.rfcxs'))
    replay = run('replay-pattern', '--record', root/(name+'-cli.rfcxs'),
                 '--out', root/(name+'-replayed.rfcxs'))
    assert view == replay, name
result = run('workflow', '--study', root/'study.rfcxs', '--request', root/'workflow.json',
             '--patterns', root/'roles.json', '--out', root/'workflow-cli.rfcxs')
assert result['status'] == 'complete' and result['result']['matched'] == 1
summary = run('workflow', '--study', root/'study.rfcxs', '--request', root/'workflow.json',
              '--patterns', root/'roles.json', '--out', root/'workflow-summary.rfcxs', '--view', 'receipt')
assert summary['schema'] == 'gadget.scientific.workflow.summary.v1' and summary['status'] == 'complete'
assert 'result' not in summary and 'trace' not in summary and summary['trace_events'] == len(result['trace'])
assert summary['source_digest'] == result['source_digest'] and summary['workflow_sha256'] == result['workflow_sha256']
retained = subprocess.run([str(consumer), '--view-workflow', str(root/'workflow-summary.rfcxs')],
                          check=True, capture_output=True, text=True)
complete = json.loads(retained.stdout)
assert complete['status'] == 'complete' and complete['result']['matched'] == 1
assert len(complete['trace']) == len(result['trace'])
for view, output in [('receipt', False), ('unknown', True)]:
    destination = root/('workflow-view-' + view + '.rfcxs')
    args = [str(binary), 'workflow', '--study', str(root/'study.rfcxs'),
            '--request', str(root/'workflow.json'), '--view', view]
    if output:
        args += ['--out', str(destination)]
    refused = subprocess.run(args, capture_output=True)
    assert refused.returncode != 0 and not destination.exists()
failed_workflow = json.loads((root/'workflow.json').read_text())
failed_workflow['limits']['nodes'] = 1
(root/'failed-workflow.json').write_text(json.dumps(failed_workflow))
failed = subprocess.run([str(binary), 'workflow', '--study', str(root/'study.rfcxs'),
    '--request', str(root/'failed-workflow.json'), '--patterns', str(root/'roles.json'),
    '--out', str(root/'failed-workflow.rfcxs'),
    '--view', 'receipt'], capture_output=True, text=True)
assert failed.returncode == 2 and (root/'failed-workflow.rfcxs').is_file()
assert json.loads(failed.stdout)['status'] == 'failed' and json.loads(failed.stdout)['error']
for name in ('candidates', 'program', 'gluing'):
    query = json.loads((root/(name+'.json')).read_text())
    query['model_sha256'] = 'sha256:foreign'
    path = root/(name+'-foreign.json')
    path.write_text(json.dumps(query)+'\n')
    destination = root/(name+'-refused.rfcxs')
    result = subprocess.run([str(binary), 'pattern-query', '--study', str(root/'study.rfcxs'),
                             '--model', str(root/'model.rfcxs'), '--request', str(path),
                             '--out', str(destination)], capture_output=True)
    assert result.returncode != 0 and not destination.exists(), name
print('Native pattern CLI queries, full replay, dependent workflow and foreign model refusals PASS')
