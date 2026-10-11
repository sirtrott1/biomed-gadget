"""Execute the staged protocol and reproduce its results from portable inputs."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
a=parser.parse_args();binary=a.binary.resolve();out=a.out.resolve();here=Path(__file__).resolve().parent
if out.exists():raise SystemExit('Output directory already exists')
out.mkdir(parents=True)
def run(*args):
 p=subprocess.run([str(binary),*map(str,args)],text=True,capture_output=True)
 if p.returncode:raise RuntimeError(p.stderr)
 return json.loads(p.stdout) if p.stdout.strip() else None
run('import-assays','--spec',here/'import.json','--out',out/'source')
plan=run('pipeline-prepare','--workspace',out/'source/workspace.json','--request',here/'pipeline.json')
(out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
intent=json.loads((here/'intent.json').read_text())
catalog=run('recipe-catalog','--workspace',out/'source/workspace.json','--intent',here/'intent.json')
proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',action='recipes',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],goals=['support'])
(out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
selected=run('recipe-proposal','--workspace',out/'source/workspace.json','--intent',here/'intent.json','--request',out/'proposal.json','--prompt',intent['answer_contract']['question'],'--out',out/'pipeline')
if selected['status']!='admitted' or selected['pipeline']['steps']!=plan['steps']:
 raise RuntimeError('Typed selection does not compile to the authored protocol')
(out/'selection.json').write_text(json.dumps(selected,indent=2)+'\n')
shutil.copytree(out/'pipeline/combined/workspace',out/'bound-assays')
shutil.copytree(out/'pipeline/sections/workspace',out/'compatible-rna')
shutil.move(out/'source',out/'generating-source');shutil.move(out/'pipeline',out/'generating-pipeline')
run('native-workspace','--workspace',out/'bound-assays/workspace.rfcxs','--out',out/'reopened-bound-assays.rfcxs')
if (out/'reopened-bound-assays.rfcxs').read_bytes()!=(out/'bound-assays/workspace.rfcxs').read_bytes():raise RuntimeError('Bound archive bytes differ after reopening')
replays={}
for assay in ['rna','atac']:
 replays[assay]=run('replay-derivation','--workspace',out/'bound-assays/workspace.rfcxs','--assay',assay,'--out',out/('replayed-'+assay))['status']
replays['compatible_rna']=run('replay-derivation','--workspace',out/'compatible-rna/workspace.rfcxs','--assay','rna-compatible','--out',out/'replayed-compatible-rna')['status']
result=run('record-view','--input',out/'compatible-rna/workspace.rfcxs','--assay','rna-compatible')['result']
answers=json.loads((out/'generating-pipeline/answers.json').read_text())
answer_artifacts=[ref for answer in answers['answers'] for artifact in answer['artifacts'] for ref in artifact['results']]
for ref in answer_artifacts:
 data=(out/'generating-pipeline'/ref['path']).read_bytes()
 if ref['sha256']!='sha256:'+hashlib.sha256(data).hexdigest():raise RuntimeError('Answer artifact mismatch')
summary=dict(status='PASS',answer_admission=selected['scientific_admission']['status'],answer_artifacts=len(answer_artifacts),replays=replays,sections=[dict(id=s['id'],residual_quadrance=s['residual_quadrance'],complete=s['complete']) for s in result['sections']])
(out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
