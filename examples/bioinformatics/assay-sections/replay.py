#!/usr/bin/env python3
"""Run the complete declared assay-section example using a supplied GADGET binary."""
import argparse
import json
import shutil
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
binary, out = args.binary.resolve(), args.out.resolve()
if out.exists():
    parser.error('output directory already exists')
out.mkdir(parents=True)
inputs = Path(__file__).resolve().parent


def run(name, *arguments):
    result = subprocess.run([str(binary), *map(str, arguments)], check=True,
                            text=True, capture_output=True)
    (out / (name + '.json')).write_text(result.stdout)
    return json.loads(result.stdout)


run('import', 'import-assays', '--spec', inputs/'import.json', '--out', out/'workspace')
manifest = out/'workspace/workspace.json'
run('sections', 'assay-query', '--workspace', manifest, '--request', inputs/'query.json')
run('plan', 'pipeline-prepare', '--workspace', manifest, '--request', inputs/'pipeline.json')
run('pipeline', 'pipeline', '--workspace', manifest, '--request', out/'plan.json', '--out', out/'pipeline')
derived = out/'pipeline/sections/workspace'
shutil.copytree(derived, out/'portable')
# Remove access to the generating workspace. Replay consumes retained lineage alone.
shutil.move(out/'workspace', out/'generating-workspace')
replayed = run('replay', 'replay-derivation', '--workspace', out/'portable/workspace.rfcxs',
               '--assay', 'rna-compatible', '--out', out/'replayed')
if replayed['status'] != 'PASS':
    raise RuntimeError('portable derivation replay failed')
print(json.dumps(dict(status='PASS', output=str(out), sections=5, views=3, portable_replay=True)))
