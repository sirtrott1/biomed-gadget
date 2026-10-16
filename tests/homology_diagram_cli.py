from pathlib import Path
import json
import shutil
import subprocess
import sys

cli, producer, destination = sys.argv[1:]
root = Path(destination).resolve()
root.mkdir(parents=True, exist_ok=True)
subprocess.run([producer, str(root / 'original')], check=True)
shutil.copytree(root / 'original', root / 'relocated inputs', dirs_exist_ok=True)
for product in ['rational_zigzag', 'integral_modules', 'rational-positive']:
    spec = root / 'relocated inputs' / (product + '.json')
    original = spec.with_suffix('.rfcxs')
    result = root / (product + '-cli.rfcxs')
    replay = root / (product + '-cold.rfcxs')
    subprocess.run([cli, 'homology-diagram', '--spec', str(spec), '--out', str(result)], check=True, cwd='/tmp')
    # The complete archived sources own replay after the input directory moves.
    hidden = root / 'hidden'
    (root / 'relocated inputs').rename(hidden)
    subprocess.run([cli, 'replay-homology', '--record', str(result), '--out', str(replay)], check=True, cwd='/tmp')
    hidden.rename(root / 'relocated inputs')
    bad = json.loads(spec.read_text())
    bad['arrows'][0]['maps'][2]['entries'][0][2] = '1'
    broken = spec.parent / 'invalid-chain.json'
    broken.write_text(json.dumps(bad))
    failed = subprocess.run([cli, 'homology-diagram', '--spec', str(broken), '--out', str(root / 'must-not-exist.rfcxs')], text=True, capture_output=True)
    assert failed.returncode and not (root / 'must-not-exist.rfcxs').exists(), failed
print('relocated CLI source loading, offline cold replay and chain refusal passed')
