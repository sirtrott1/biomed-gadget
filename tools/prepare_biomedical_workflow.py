#!/usr/bin/env python3
"""Prepare a source-bound C++ binding-evidence workflow; optionally execute it."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('build/presets/cpu/gadget'))
    p.add_argument('--evidence',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--endpoint',choices=['Kd','Ki','IC50','EC50'],default='Ki')
    p.add_argument('--maximum',default='0.5')
    p.add_argument('--model',type=Path)
    p.add_argument('--reading',choices=['affinity_records','coverage','reported_points'],default='affinity_records',
                   help='Complete source query supplied to the observer; all mathematical results remain in the workflow')
    p.add_argument('--context-tokens',type=int)
    p.add_argument('--output-tokens',type=int)
    p.add_argument('--threads',type=int)
    p.add_argument('--run',action='store_true')
    a=p.parse_args(); binary=a.binary.resolve(); evidence=a.evidence.resolve(); study=evidence/'study/study.json'
    catalogue=json.loads(subprocess.check_output([str(binary),'audit','--study',str(study)]))
    store=json.loads((evidence/'store.json').read_text())
    entry={'id':'root','type':'sequence','steps':[
        {'id':'criterion','type':'evidence','query':{'operation':'affinity_evidence','endpoint':a.endpoint,'units':'nM','maximum':a.maximum}},
        {'id':'exact-views','type':'parallel','workers':3,'steps':[
            {'id':'point-support','type':'query','query':{'operation':'support_mask','fields':['all.Kd_nM'],'mode':'all_reported','backend':'cpu_tiles','workers':4,'tile_words':1}},
            {'id':'source-memberships','type':'sequence','steps':[
                {'id':'memberships-cold','type':'query','query':{'operation':'native_boundary','field_set':'category.memberships','backend':'cpu_cached_tiles','workers':4,'tile_rows':1}},
                {'id':'memberships-warm','type':'query','query':{'operation':'native_boundary','field_set':'category.memberships','backend':'cpu_cached_tiles','workers':4,'tile_rows':1}}]},
            {'id':'reported-points','type':'query','query':{'operation':'field_summary','field':'all.Kd_nM','support':'reported_only'}}]},
        {'id':'evidence-followup','type':'branch','from':'criterion','member':'selected','at_least':1,
         'on_true':{'id':'reader','type':'model','role':'reader','prompt':'Explain the supplied per-record threshold certificate. State the selected, contradicted and indeterminate counts and explain why bounds must not be treated as exact point measurements. Describe saved source assertions only; do not infer efficacy, biology or clinical use. Give one concise observation.','context':['criterion']} if a.model else {'id':'retained-coverage','type':'evidence','query':{'operation':'coverage'}},
         'on_false':{'id':'unselected','type':'evidence','query':{'operation':'coverage'}}}]}
    if a.model and a.reading=='coverage':
        entry['steps'][-1]['on_true']={
            'id':'coverage-reading','type':'sequence','steps':[
                {'id':'reported-coverage','type':'evidence','query':{'operation':'coverage'}},
                {'id':'reader','type':'model','role':'reader',
                 'prompt':'Read the complete coverage certificate. State its reported, missing, equality-relation and total binding-record counts. Explain why a reported numerical bound and an equality relation are different, and why equality alone does not establish comparable assays. Give one concise observation about the saved data.',
                 'context':['reported-coverage']}]}
    if a.model and a.reading=='reported_points':
        entry['steps'][-1]['on_true']={
            'id':'checked-point-reading','type':'sequence','steps':[
                {'id':'reader','type':'model','role':'reader',
                 'prompt':'Report exactly four equality comparisons from the supplied field-summary certificate: path ["reported"] as count, ["missing"] as count, ["sum","rational"] as rational, and ["mean","rational"] as rational. Copy the certificate hash and source digest. These describe the saved numerical tokens; make no biological claim.',
                 'context':['reported-points']},
                {'id':'checked-reader','type':'check_claims','from':'reported-points','claims_from':'reader'}]}
    workflow={'schema':'gadget.scientific.workflow.v1','id':'binding-evidence-to-exact-math','source_digest':catalogue['source_digest'],'evidence_digest':store['source_digest'],'limits':{'depth':12,'nodes':64,'model_calls':8 if a.model else 0,'parallel_workers':3},'definitions':[],'entry':entry}
    a.out.mkdir(parents=True,exist_ok=False); request=a.out/'workflow.json'; request.write_text(json.dumps(workflow,indent=2)+'\n')
    command=[str(binary),'workflow','--study',str(study),'--evidence',str(evidence),'--request',str(request.resolve()),'--out',str((a.out/'operation.rfcxs').resolve())]
    if a.model:
        role={'role':'reader','task':'claim_reader' if a.reading=='reported_points' else 'observer','model':str(a.model.resolve())}
        for attribute in ['threads','context_tokens','output_tokens']:
            value=getattr(a,attribute)
            if value is not None:role[attribute]=value
        roles={'schema':'gadget.local.roles.v1','models':[role]}
        path=a.out/'roles.json';path.write_text(json.dumps(roles,indent=2)+'\n');command+=['--models',str(path.resolve())]
    (a.out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    if a.run:
        with (a.out/'receipt.json').open('w') as out,(a.out/'runtime.log').open('w') as err:
            subprocess.run(command,stdout=out,stderr=err,check=True)
    print(json.dumps({'workflow':str(request),'command':command,'executed':a.run}))

if __name__=='__main__':main()
