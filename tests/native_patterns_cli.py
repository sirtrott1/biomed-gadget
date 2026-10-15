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
