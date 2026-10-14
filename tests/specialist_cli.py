"""Selected specialist recipes and independent exact assay-map qualification."""
import argparse
import copy
import hashlib
import json
import shutil
import subprocess
from fractions import Fraction as Q
from pathlib import Path


def qualify(binary, root):
    root.mkdir(parents=True, exist_ok=False)
    checks = 0
    calls = 0

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise AssertionError(message)

    def save(name, value):
        path = root / name
        path.write_text(value if isinstance(value, str) else json.dumps(value, indent=2)+'\n')
        return path

    def source(name, value):
        path = save(name, value)
        return dict(path=name, sha256='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest())

    def run(*args, fail=False):
        nonlocal calls
        calls += 1
        p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True, timeout=180)
        check(bool(p.returncode) == fail, str(args)+': '+p.stderr[-3000:])
        return None if fail else json.loads(p.stdout)

    capacity = run('resources')['worker_capacity']
    workers = min(4, capacity)
    design = dict(question='Declared assay maps and retained structure', population='Authored measurements',
                  selection='Every declared source identity', limitations=['Verification fixture'])

    def assay(identity, namespace, features, observations, rows, quantity='read_count'):
        return dict(id=identity, kind='expression', format='dense_tsv',
            source=source(identity+'.tsv', 'id\t'+'\t'.join(observations)+'\n'+''.join(
                f+'\t'+'\t'.join('NA' if v is None else str(v) for v in r)+'\n'
                for f, r in zip(features, rows))),
            features=source(identity+'-features.tsv', 'id\n'+'\n'.join(features)+'\n'),
            observations=source(identity+'-observations.tsv','id\n'+'\n'.join(observations)+'\n'),
            source_namespace=str(namespace), carrier='features', quantity=quantity, units='counts' if quantity=='read_count' else 'units',
            channel='E', method='Authored complete source and explicit absence', arithmetic='exact_rational',
            implicit_entries='unobserved', missing_tokens=['NA'])

    def workspace(name, assays, resources=None):
        spec = dict(schema='gadget.assays.import.v1', id=name, design=design,
                    assays=assays, resources=resources or [], links=[])
        run('import-assays','--spec',save(name+'-import.json',spec),'--out',root/name)
        return root/name/'workspace.json'

    def execute(name, manifest, request):
        out=run('assay-query','--workspace',manifest,'--request',save(name+'-request.json',request))
        save(name+'-result.json',out)
        return out['result']

    def recipe(name, manifest, settings):
        path=save(name+'-default.json',settings)
        prepared=run('pipeline-prepare','--workspace',manifest,'--request',path)
        save(name+'-plan.json',prepared)
        index=run('pipeline','--workspace',manifest,'--request',path,'--out',root/name)
        for item in index['results']:
            p=root/name/item['result']
            check(item['result_sha256']=='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest(), 'Complete result receipt')
        save(name+'-index.json',index)
        return prepared,index

    def output(identity, namespace):
        return dict(id=identity, source_namespace=str(namespace), design=design, channel='E')

    def default(profile, goals, question):
        return dict(schema='gadget.pipeline.default.v1', profile=profile, question=question, goals=goals)

    def leaf(name, step):
        return json.loads((root/name/step/'results/partition-0.json').read_text())['result']

    def generating_result(name, step):
        catalog=json.loads((root/name/step/'result.json').read_text())
        if catalog['schema']=='gadget.assays.workspace.v2':
            return run('record-view','--input',root/name/step/'workspace/workspace.rfcxs',
                       '--assay',catalog['assays'][0]['declaration']['id'])['result']
        record=catalog['assays'][0]['originals']['provenance-0']
        path=root/name/step/'workspace'/record['path']
        check(record['sha256']=='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest(),
              'Producing certificate retained with complete hash')
        return json.loads(path.read_text())['result']

    # The oracle evaluates two declared maps, propagating absence per nonzero
    # dependency at each stage. It never substitutes zero for an unknown input.
    def aggregate_oracle(features, observations, rows, feature_groups, observation_groups):
        fi={f:i for i,f in enumerate(features)}; oi={s:i for i,s in enumerate(observations)}
        first=[]
        for g in feature_groups:
            used=[(fi[m['id']],Q(m['coefficient'])) for m in g['members'] if Q(m['coefficient'])]
            first.append([None if any(rows[f][s] is None for f,c in used)
                          else sum((c*Q(rows[f][s]) for f,c in used),Q(0))
                          for s in range(len(observations))])
        result=[]
        for row in first:
            target=[]
            for g in observation_groups:
                used=[(oi[m['id']],Q(m['coefficient'])) for m in g['members'] if Q(m['coefficient'])]
                target.append('NA' if any(row[s] is None for s,c in used)
                              else str(sum((c*row[s] for s,c in used),Q(0))))
            result.append(target)
        return result

    features=['a','b','c','d','unused']; observations=['x','y','z','unused']
    rows=[[2,5,11,None],[None,7,13,None],[3,None,17,None],[0,0,19,None],[None]*4]
    matrix=assay('cells',891000,features,observations,rows)
    manifest=workspace('maps-workspace',[matrix])
    group=lambda name, members: dict(id=name,members=[dict(id=i,coefficient=str(c)) for i,c in members])
    fg=[group('sum',[('a',1),('b',1)]), group('difference',[('c',1),('d',-1)]),group('zero',[('unused',0)])]
    og=[group('combined',[('x',1),('y',1)]),group('single',[('z',1)]),group('signed',[('x',1),('y',-1)]),group('zero',[('unused',0)])]
    params=dict(assay='cells',feature_groups=fg,observation_groups=og,map_method='Declared linear maps',
                map_evidence='Every coefficient and source identity authored',support_policy='require_all_nonzero_inputs')
    query=dict(schema='gadget.assays.query.v1',operation='aggregate_assay',**params)
    before=execute('aggregate-one',manifest,query)
    after=execute('aggregate-many',manifest,{**query,'workers':workers})
    expected=aggregate_oracle(features,observations,rows,fg,og)
    for actual in [before,after]:
        for a,r in zip(actual['rows'],expected):
            for x,y in zip(a,r): check(x==y,'Independent supported rational aggregation')
        check(actual['execution']['required_source_features']==4,'Zero coefficient does not read missing input')
        check(actual['execution']['required_source_observations']==3,'Unreferenced observations not materialized')
        check(actual['execution']['intermediate_coordinates']==0 and actual['execution']['sparse_map_owner']=='rfc::numeric::PreparedSparseExactMap','Native selected maps retain no dense source intermediate')
        check(actual['execution']['first_map']=='features','Apply the original stored feature carrier map')
    # A smaller target observation axis retains complete-support semantics,
    # including signed cancellation with an absent required input.
    compact={**query,'observation_groups':[group('merged',[('x',1),('y',-1),('z',1)])], 'workers':workers}
    reversed_order=execute('observation-first',manifest,compact)
    check(reversed_order['rows']==aggregate_oracle(features,observations,rows,fg,compact['observation_groups']),
          'Observation-first signed and absent-support oracle')
    check(reversed_order['execution']['first_map']=='features' and
          reversed_order['execution']['intermediate_coordinates']==0,'Contract the smaller observation target without a dense intermediate')
    for invalid in [0,capacity+1]:
        run('assay-query','--workspace',manifest,'--request',save('invalid-workers-'+str(invalid)+'.json',{**query,'workers':invalid}),fail=True)
    batched=execute('aggregate-batch',manifest,dict(schema='gadget.assays.query.v1',operation='query_batch',
        concurrent_jobs=min(2,capacity),queries=[dict(id='aggregate',query={**query,'workers':workers}),
        dict(id='repeat',query={**query,'workers':1})]))
    check(batched['results'][0]['admitted_workers']==workers,'Direct query batch delegates aggregation workers')
    check(batched['results'][0]['answer']['result']['rows']==expected,'Joined batch preserves every coordinate')
    check(batched['results'][1]['id']=='repeat' and batched['results'][1]['answer']['result']['rows']==expected,
          'Independent repeated native publisher retains every original coordinate')
    cancel=[group('annihilate',[('a',1),('b',-1)])]
    cancellation=execute('cancellation',manifest,{**query,'feature_groups':cancel,'observation_groups':[group('x',[('x',1)])]})
    check(cancellation['rows']==[['NA']],'Numerical cancellation cannot certify missing support')
    all_zero=execute('zero-map',manifest,{**query,'feature_groups':[fg[-1]],'observation_groups':[og[-1]],'workers':workers})
    check(all_zero['rows']==[['0']] and all_zero['execution']['feature_stage_completed_tiles']==0,'Empty native map support gives known zero')
    small=default('single_cell',[dict(id='mapped',method='aggregation',parameters={**params,'workers':workers},output=output('mapped',892000)),
        dict(id='support',method='coverage',input='mapped',parameters=dict(fields=['all.combined']))], 'Inspect support after declared feature and observation maps')
    plan,index=recipe('aggregation-pipeline',manifest,small)
    check(generating_result('aggregation-pipeline','mapped')['execution']['requested_workers']==workers,
          'Derived recipe receives requested native CPU allocation')
    if workers > 1:
        invalid=copy.deepcopy(small);invalid['goals'][0]['workers']=1
        run('pipeline','--workspace',manifest,'--request',save('undersized-allocation.json',invalid),
            '--out',root/'undersized-allocation',fail=True)
        check(not (root/'undersized-allocation').exists(),'Insufficient assigned budget has no published result')
    check(leaf('aggregation-pipeline','support')['missing']==2,'Derived source retains unknown coordinates')
    replay=run('replay-derivation','--workspace',root/'aggregation-pipeline/mapped/workspace/workspace.rfcxs','--assay','mapped','--out',root/'aggregation-replay')
    check(replay['status']=='PASS','Parallel aggregation replays from retained source files')

    # More than the former 4096 output groups, complete source ordering and
    # disjoint pseudobulk count semantics, while selecting 3 of 12 observations.
    n=4205; features=['g'+str(i) for i in range(n)]; observations=['cell'+str(i) for i in range(12)]
    rows=[[None if (f==17 and s==1) or s==11 else (7*f+3*s)%97 for s in range(12)] for f in range(n)]
    manifest_large=workspace('large-workspace',[assay('counts',893000,features,observations,rows)])
    fg=[group(f,[(f,1)]) for f in features]
    og=[group('left',[('cell0',1),('cell1',1)]),group('right',[('cell3',1)])]
    large_params=dict(assay='counts',feature_groups=fg,observation_groups=og,map_method='Disjoint specimen sums',
                      map_evidence='Authored cell membership',support_policy='require_all_nonzero_inputs',count_policy='disjoint_count_sum',workers=workers)
    settings=default('single_cell',[dict(id='pseudobulk',method='aggregation',parameters=large_params,output=output('pseudobulk',894000)),
        dict(id='support',method='coverage',input='pseudobulk',parameters=dict(fields=['all.left']),workers=workers)], 'Retain pseudobulk counts and inspect every reported gene')
    plan,index=recipe('large-pseudobulk',manifest_large,settings)
    result=generating_result('large-pseudobulk','pseudobulk')
    expected=aggregate_oracle(features,observations,rows,fg,og)
    check(result['feature_ids']==features,'Every gene identity and original order retained')
    for a,b in zip(result['rows'],expected):
        for x,y in zip(a,b):check(x==y,'Every large pseudobulk source-derived count')
    check(result['quantity']=='read_count' and result['count_contract_verified'],'Count meaning retained through admitted map')
    check(result['execution']['required_source_observations']==3 and result['execution']['intermediate_coordinates']==0 and result['execution']['sparse_map_owner']=='rfc::numeric::PreparedSparseExactMap','Required native count fields use sparse contraction')
    check(len(index['results'])==2,'Unrequested comparison not run')
    support_index=json.loads((root/'large-pseudobulk/support/result.json').read_text())
    source_index=support_index['results'][0]['source_index']
    read=dict(schema='gadget.assays.source.read.v1',source_digest=support_index['source_digest'],
              assay='pseudobulk',partition=0,read=dict(schema='gadget.scientific.source.read.v1',
              source_digest=source_index['source_digest'],source_sequence=source_index['source_sequence'],
              publication='self',stream=dict(kind='catalog',section='cells'),offset=0,limit=len(features)))
    page=run('assay-source-read','--workspace',root/'large-pseudobulk/pseudobulk/workspace/workspace.rfcxs',
             '--request',save('derived-source-cells.json',read))['page']
    check([cell['id'] for cell in page['records']]==features and page['complete'],
          'Native source page covers every output feature in source order')
    check(sum(run('record-view','--input',root/'large-pseudobulk/support/results'/packet['path'],'--member','result')['result']['missing']
              for packet in support_index['results'])==1,'Missing count retained across all actual output storage partitions')

    # Independent carriers and an explicitly declared rectangular map; no
    # label equality, fitting or common embedding supplies the correspondence.
    nl,nr=513,517
    left_ids=['gene_'+str(i) for i in range(nl)];right_ids=['protein_'+str(i) for i in range(nr)]
    left=assay('rna',895000,left_ids,['specimen'],[[2*i+1] for i in range(nl)])
    right=assay('protein',896000,right_ids,['specimen'],[[Q(3*i+1,2)] for i in range(nr)],'measured_abundance')
    multi=workspace('multi-workspace',[left,right])
    mapping=dict(domain_quantity='measured_abundance',domain_units='units',codomain_quantity='read_count',codomain_units='counts',
                 method='Authored quantitative measurement map',evidence='Every cross-assay pair and coefficient declared',
                 entries=[dict(left_id=l,right_id=right_ids[i],coefficient='1/2') for i,l in enumerate(left_ids)])
    transfer=dict(left=dict(assay='rna',partition=0,field='all.specimen'),right=dict(assay='protein',partition=0,field='all.specimen'),metric='identity',map=mapping)
    settings=default('multi_omics',[dict(id='mapped-comparison',method='field_transfer',parameters=transfer)],'Compare independently measured fields through the declared map')
    _,index=recipe('multi-comparison',multi,settings)
    mapped=json.loads((root/'multi-comparison/mapped-comparison/result.json').read_text())['result']
    check(mapped['left_ids']==left_ids and mapped['right_ids']==right_ids,'Independent complete carriers retained')
    for i in range(nl):
        check(Q(mapped['transported_right'][i])==Q(3*i+1,4),'Every mapped coordinate')
        check(Q(mapped['residual'][i])==Q(5*i+3,4),'Every residual coordinate')
    check(mapped['map_rank']==nl and mapped['map_kernel_dimension']==nr-nl,'Exact rectangular rank-nullity')
    check(Q(mapped['residual_quadrance'])==sum(Q(5*i+3,4)**2 for i in range(nl)),'Independent full residual quadrance')
    check(Q(mapped['pairing'])==sum(Q(2*i+1)*Q(3*i+1,4) for i in range(nl)),'Independent full pairing')
    bad=copy.deepcopy(settings);bad['goals'][0]['parameters']['map']['entries'][0]['right_id']=left_ids[0]
    run('pipeline','--workspace',multi,'--request',save('wrong-carrier-default.json',bad),'--out',root/'wrong-carrier',fail=True)
    check(not (root/'wrong-carrier').exists(),'Invalid cross-assay identity has no published result')
    incomplete=copy.deepcopy(right)
    incomplete['id']='incomplete';incomplete['source_namespace']='896100'
    incomplete['source']=source('incomplete.tsv','id\tspecimen\n'+''.join(
        f+'\t'+('NA' if i==nr-1 else str(Q(3*i+1,2)))+'\n' for i,f in enumerate(right_ids)))
    partial=workspace('partial-map-workspace',[left,incomplete])
    invalid=copy.deepcopy(transfer);invalid['right']['assay']='incomplete'
    run('assay-query','--workspace',partial,'--request',save('partial-transfer.json',dict(
        schema='gadget.assays.query.v1',operation='field_transfer',**invalid)),fail=True)

    # A sequence recipe produces a true native grade-one carrier. Dependent
    # recipes retain the measured loop count and original signed occurrences.
    resource=dict(id='reads',format='fasta',source=source('reads.fa','>a\nAACG\n>b\nAACT\n>c\nAAAA\n>d\nAAAA\n'),method='Authored sequences')
    sequence=workspace('sequence-workspace',[],[resource])
    settings=default('sequence',[
        dict(id='assembly',method='sequence_assembly_graph',parameters=dict(resource='reads',k=3,ambiguity='literal_tokens',evidence='All observed k-mers'),output=output('assembly',897000)),
        dict(id='structure',method='chain_certificate',input='assembly',parameters=dict(context='all',grade=1,coefficient_domain='rational',field='all.kmer_occurrence_count')),
        dict(id='response',method='green_action',input='assembly',parameters=dict(context='all',field='all.kmer_occurrence_count',**{'lambda':'1'})),
        dict(id='paths',method='path_comparison',input='assembly',parameters=dict(context='all',grade=1,paths=[dict(id='loop',steps=[dict(index=3,coefficient='1')]),dict(id='empty',steps=[])])),
        dict(id='filtration',method='relational_persistence',input='assembly',parameters=dict(assay='assembly',partition=0,levels=['0','1'],births=[dict(grade=0,index=i,level=0) for i in range(4)]+[dict(grade=1,index=i,level=1) for i in range(4)],filtration_meaning='Declared vertex then observed k-mer stages',evidence='All cell births authored'))
    ],'Inspect the measured loop and complete field on the observed k-mer complex')
    _,index=recipe('sequence-structure',sequence,settings)
    structure=leaf('sequence-structure','structure')
    check(structure['rational_coordinates']==['2','1','1','4'],'Complete observed k-mer counts')
    check(structure['sectors']['harmonic']==['0','0','0','4'],'Loop count remains a complete independent harmonic coordinate')
    check([Q(c['value']['rational']) for c in leaf('sequence-structure','response')]==[Q(1),Q(1,2),Q(1,2),Q(4)],'Full exact Green response')
    paths=leaf('sequence-structure','paths')
    check(paths['pairs'][0]['equal_boundaries'] and not paths['pairs'][0]['equal_fields'],'Equal endpoint view does not erase the observed loop')
    filtration=json.loads((root/'sequence-structure/filtration/result.json').read_text())['result']
    check([[g['betti'] for g in s['grades']] for s in filtration['stages']]==[[4,0],[1,1]],'Declared sequence filtration preserves loop')
    check(filtration['rank_reconstruction_verified'],'Complete interval inclusion-rank identity')
    # The local model may request only structure; admission supplies its exact
    # producing recipe while leaving Green, paths and filtration unselected.
    intent=dict(schema='gadget.pipeline.intent.v1',profile='sequence',recipes=[])
    for original in settings['goals']:
        goal=copy.deepcopy(original);identity=goal.pop('id')
        intent['recipes'].append(dict(id=identity,purpose={'assembly':'Construct the observed k-mer source', 'structure':'Report complete measured counts and gradient, curl, harmonic sectors and cycle witnesses', 'response':'Compute the full Green response', 'paths':'Compare ordered loop and empty histories','filtration':'Compute the declared two-stage persistence'}[identity],goal=goal))
    save('sequence-intent.json',intent)
    catalog=run('recipe-catalog','--workspace',sequence,'--intent',root/'sequence-intent.json')
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',action='recipes',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],goals=['structure'])
    admission=run('recipe-proposal','--workspace',sequence,'--intent',root/'sequence-intent.json','--request',save('sequence-proposal.json',proposal),'--prompt','Report the measured counts and their complete Hodge sectors','--out',root/'selected-structure')
    check(admission['included_goals']==['assembly','structure'],'Source-producing dependency closure only')
    check(not (root/'selected-structure/response').exists(),'Unselected Green operation not executed')
    check(leaf('selected-structure','structure')['sectors']==structure['sectors'],'Selection uses the same native mathematics')
    summary=dict(status='PASS',independent_checks=checks,command_calls=calls,aggregation_output_groups=n,
                 aggregation_workers=workers,cross_assay_carriers=[nl,nr],sequence_grade_one=True,
                 source_digest=index['source_digest'])
    save('verification.json',summary)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary',type=Path)
    parser.add_argument('root',type=Path)
    args=parser.parse_args()
    root=args.root.resolve()
    if root.exists():shutil.rmtree(root)
    print(json.dumps(qualify(args.binary.resolve(),root)))


if __name__=='__main__':main()
