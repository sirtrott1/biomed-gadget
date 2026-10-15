"""Pipeline selection, derived-source dependencies, seals and independent exact fields."""
import copy
import hashlib
import json
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)
checks = 0

def run(*args, fail=False, error=None):
    global checks
    p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True, timeout=180)
    assert bool(p.returncode) == fail, (args, p.stdout, p.stderr)
    if error is not None:
        assert error in p.stderr, (args, p.stderr, error)
    checks += 1
    return None if fail else json.loads(p.stdout)

def save(name, data):
    p = root / name
    p.write_text(data if isinstance(data, str) else json.dumps(data))
    return p

def source(name, data):
    p = save(name, data)
    return dict(path=name, sha256='sha256:' + hashlib.sha256(p.read_bytes()).hexdigest())

features = source('features.tsv', 'id\ng1\ng2\ng3\n')
observations = source('observations.tsv', 'id\nA\nB\n')
matrix = source('counts.tsv', 'id\tA\tB\ng1\t2\t5\ng2\t0\t3\ng3\t7\t11\n')
reads = source('reads.fastq', '@r1\nACGT\n+\nIIII\n@r2\nAGGT\n+\n!III\n')
design = dict(question='Authored source identities and integer measurements', population='three authored features', selection='all', limitations=['nonclinical fixture'])
assay = dict(id='counts', kind='expression', format='dense_tsv', source=matrix, features=features, observations=observations, source_namespace='881100', carrier='features', quantity='read_count', units='counts', channel='E', method='authored counts', arithmetic='exact_decimal', implicit_entries='measured_zero', missing_tokens=[])
spec = dict(schema='gadget.assays.import.v1', id='pipeline-contracts', design=design, assays=[assay], resources=[dict(id='reads', format='fastq', source=reads, method='authored reads')], links=[])
workspace = root / 'workspace'
run('import-assays', '--spec', save('import.json', spec), '--out', workspace)
manifest = workspace / 'workspace.json'
profiles = run('pipelines')
assert {'rna_seq', 'single_cell', 'atac_seq', 'chip_seq', 'spatial', 'multi_omics'} <= {p['profile'] for p in profiles['profiles']}

# Default schemas reference the existing operation arguments. Check that those
# references continue to point to the intended contracts if schemas are reordered.
schema_root=Path(__file__).resolve().parent.parent/'schemas'
recipe_schemas=json.loads((schema_root/'pipeline-default.schema.json').read_text())['properties']['goals']['items']['oneOf'][1:]
recipe_schemas={entry['properties']['method']['const']:entry for entry in recipe_schemas}
assert set(recipe_schemas)=={recipe['method'] for recipe in profiles['recipes']}
for recipe in profiles['recipes']:
    assert recipe['native_publication'], recipe
    assert recipe['native_materialization'] == (recipe['operation'] in ['controlled_coverage','spatial_relations','feature_counts','expression_allocations','aggregate_assay','design_aggregation','variant_genotypes','chromatin_contrast','assay_sections','constrained_allocations','exact_sequence_alignment','sequence_assembly_graph'])
    parameters=recipe_schemas[recipe['method']]['properties']['parameters']
    refs=[value['$ref'] for value in parameters.get('properties',{}).values() if '$ref' in value]
    for ref in refs:
        filename,pointer=ref.split('#',1)
        assert filename==recipe['parameters_schema']
        parts=pointer.strip('/').split('/')
        document=json.loads((schema_root/filename).read_text())
        referenced=document
        for part in parts:
            part=part.replace('~1','/').replace('~0','~')
            referenced=referenced[int(part)] if isinstance(referenced,list) else referenced[part]
        assert isinstance(referenced,dict)
        if recipe.get('workspace_binding'):
            assert document['properties']['schema']['const']=='gadget.assays.bindings.v1'
        else:
            branch=document[parts[0]][int(parts[1])]
            operation=branch['properties']['operation']
            assert recipe['operation'] in operation.get('enum',[operation.get('const')])
assert recipe_schemas['indexed_sequence_alignment']['properties']['parameters']['properties']['execution'] == {
    '$ref': 'assays-query.schema.json#/oneOf/19/properties/execution'}

def request(operation, **kw):
    return dict(schema='gadget.assays.query.v1', operation=operation, **kw)

defaults = dict(schema='gadget.pipeline.default.v1', profile='rna_seq', question='Check measured support only', goals=['coverage'], assay='counts', fields=['all.A'])
plan = run('pipeline-prepare', '--workspace', manifest, '--request', save('default.json', defaults))
assert len(plan['steps']) == 1 and plan['steps'][0]['id'] == 'coverage'
index = run('pipeline', '--workspace', manifest, '--request', save('plan.json', plan), '--out', root/'default-run')
assert len(index['results']) == 1 and index['dispatcher'] == 'rfc::runtime::WorkGraph'
for item in index['results']:
    p = root/'default-run'/item['result']
    assert item['result_sha256'] == 'sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()
    retained = run('record-view', '--input', p, '--member', 'result')
    assert retained['schema'] == 'gadget.assays.run.v1'
    assert run('record-view', '--input', p, '--member', 'request') == plan['steps'][0]['request']

aggregation = request('aggregate_assay', assay='counts',
    feature_groups=[dict(id='total', members=[dict(id=g, coefficient='1') for g in ['g1','g2','g3']])],
    observation_groups=[dict(id=s, members=[dict(id=s, coefficient='1')]) for s in ['A','B']],
    map_method='declared exact feature sum', map_evidence='Authored feature identities', support_policy='require_all_nonzero_inputs')
derive = dict(schema='gadget.assays.materialize.v1', id='summed', source_namespace='883100', design=design, recipe=aggregation, channel='E')
workflow = dict(schema='gadget.scientific.workflow.v1', id='dependent-check', limits=dict(depth=2,nodes=8,model_calls=0,parallel_workers=1), definitions=[], entry=dict(id='summary',type='query',query=dict(operation='field_summary',field='all.A',support='reported_only')))
custom = dict(schema='gadget.pipeline.design.v1', profile='custom', question='Sum retained features, inspect the actual derived source, and verify a native field independently', concurrent_jobs=2, steps=[
    dict(id='sum', type='derive', input='source', request=derive),
    dict(id='after', type='workflow', input='sum', assay='summed', request=workflow),
    dict(id='original', type='query', input='source', request=request('native_query',assay='counts',query=dict(operation='field_summary',field='all.B',support='reported_only'))),
])
index = run('pipeline', '--workspace', manifest, '--request', save('custom.json', custom), '--out', root/'custom-run')
assert [r['id'] for r in index['results']] == ['sum','after','original']
assert index['results'][1]['source_digest'] == index['results'][0]['derived_digest']
workflow_native = run('record-view', '--input', root/'custom-run/after/result.rfcxs', '--member', 'result')
assert workflow_native['schema'] == 'gadget.workflow.native.receipt.v1'
assert workflow_native['native_execution'] is True and workflow_native['status'] == 'complete'
workflow_request = run('record-view', '--input', root/'custom-run/after/result.rfcxs', '--member', 'request')
assert workflow_request == json.loads((root/'custom-run/after/request.json').read_text())
assert workflow_request['entry'] == workflow['entry']
agg = json.loads((root/'custom-run'/'sum'/'workspace'/'lineage.json').read_text()) if (root/'custom-run'/'sum'/'workspace'/'lineage.json').exists() else None
assert (root/'custom-run'/'sum'/'workspace'/'workspace.rfcxs').exists()
native = json.loads((root/'custom-run'/'original'/'results'/'partition-0.json').read_text())
assert Fraction(native['result']['sum']['rational']) == 19
run('replay-derivation','--workspace',root/'custom-run'/'sum'/'workspace'/'workspace.rfcxs','--assay','summed','--out',root/'replayed')

# A model's typed proposal uses the same admission and execution path.
proposal = dict(schema='gadget.pipeline.proposal.v1', action='pipeline', source_digest=plan['source_digest'], profile='rna_seq', question='Inspect measured support', steps=[dict(id='support',type='query',input='source',request_json=json.dumps(request('native_query',assay='counts',query=dict(operation='coverage',field='all.A'))))])
answer = run('pipeline-proposal','--workspace',manifest,'--request',save('proposal.json',proposal),'--out',root/'proposal-run')
assert answer['status'] == 'admitted' and len(answer['execution']['results']) == 1
structured=copy.deepcopy(proposal)
node=structured['steps'][0]
node['request']=json.loads(node.pop('request_json'))
answer=run('pipeline-proposal','--workspace',manifest,'--request',save('structured-proposal.json',structured),'--out',root/'structured-proposal-run')
assert answer['status']=='admitted' and len(answer['execution']['results'])==1
ambiguous=copy.deepcopy(structured);ambiguous['steps'][0]['request_json']=proposal['steps'][0]['request_json']
run('pipeline-proposal','--workspace',manifest,'--request',save('ambiguous-proposal.json',ambiguous),'--out',root/'ambiguous-proposal-run',fail=True)
assert not (root/'ambiguous-proposal-run').exists()
clarify = dict(schema='gadget.pipeline.proposal.v1',action='clarify',source_digest=plan['source_digest'],reason='Supply the actual matched sample identities')
answer = run('pipeline-proposal','--workspace',manifest,'--request',save('clarify.json',clarify),'--out',root/'clarify-run')
assert answer['status'] == 'clarification' and not (root/'clarify-run').exists()

for name, mutate in [
    ('stale', lambda p: p.update(source_digest='sha256:'+'0'*64)),
    ('tampered', lambda p: p.update(question='A changed question under an old seal')),
    ('forward', lambda p: p['steps'][0].update(input='later')),
    ('path', lambda p: p['steps'][0].update(id='../outside')),
]:
    bad=copy.deepcopy(plan); mutate(bad)
    run('pipeline','--workspace',manifest,'--request',save(name+'.json',bad),'--out',root/name,fail=True)
    assert not (root/name).exists()
bad=copy.deepcopy(custom);bad['steps'][2]['request']['query']['field']='all.unknown'
run('pipeline','--workspace',manifest,'--request',save('runtime-failure.json',bad),'--out',root/'runtime-failure',fail=True)
assert not (root/'runtime-failure').exists() and not list(root.glob('.gadget-pipeline-*'))
bad=copy.deepcopy(custom);bad['profile']='rna_seq'
run('pipeline','--workspace',manifest,'--request',save('wrong-profile.json',bad),'--out',root/'wrong-profile',fail=True)
assert not (root/'wrong-profile').exists()

# An unsupported later producer must refuse before an earlier calculation runs.
# The first leaf would fail on its field if execution started before admission.
for name, kind, unsupported, error in [
    ('unknown-query', 'query', request('unregistered_operation'),
     'pipeline query has no native result publisher'),
    ('unknown-derive', 'derive', dict(schema='gadget.assays.materialize.v1',
         id='unregistered', source_namespace='883200', design=design,
         recipe=request('unregistered_operation'), channel='E'),
     'pipeline recipe has no native materialization producer'),
    ('unknown-scientific', 'query', request('native_query', assay='counts',
         query=dict(operation='unregistered_operation')),
     'pipeline scientific query has no native result publisher'),
]:
    bad=dict(schema='gadget.pipeline.v1', source_digest=plan['source_digest'],
        profile='custom', question='Admit all native producers before executing',
        steps=[dict(id='first', type='query', input='source',
            request=request('native_query', assay='counts',
                query=dict(operation='field_summary', field='all.unknown', support='reported_only'))),
            dict(id='unsupported', type=kind, input='source', request=unsupported)])
    run('pipeline', '--workspace', manifest, '--request', save(name+'.json', bad),
        '--out', root/name, fail=True, error=error)
    assert not (root/name).exists() and not list(root.glob('.gadget-pipeline-*'))

# No arbitrary mathematical carrier cap: the source only has isolated measured
# coordinates, so its native exact Green action must return the complete input.
green=request('native_query',assay='counts',query=dict(operation='green_action',context='all',field='all.A',**{'lambda':'3/2'}))
answer=run('assay-query','--workspace',manifest,'--request',save('green.json',green))
certificate=answer['result']['answers'][0]['certificate']
assert certificate['residual_verified'] and certificate['green_iterations']==1
# A requested concurrency is scheduled against actual native workers and the
# selected partitions; it is not rejected by an unrelated four-job policy.
many=request('native_query',assay='counts',query=dict(operation='coverage',field='all.A'),concurrent_jobs=8)
answer=run('assay-query','--workspace',manifest,'--request',save('more-native-jobs.json',many),'--out',root/'more-native-jobs')
assert answer['concurrent_jobs']==1 and answer['query']['concurrent_jobs']==8
# Scientific recipes use actual source declarations and retain resolved requests.
def recipe_run(name, settings, source_manifest=manifest):
    prepared=run('pipeline-prepare','--workspace',source_manifest,'--request',save(name+'.json',settings))
    index=run('pipeline','--workspace',source_manifest,'--request',save(name+'-plan.json',prepared),'--out',root/name)
    for item in index['results']:
        for artifact in ['request','result']:
            data=(root/name/item[artifact]).read_bytes()
            assert item[artifact+'_sha256']=='sha256:'+hashlib.sha256(data).hexdigest()
    return prepared,index

def leaf(name, step):
    return json.loads((root/name/step/'results/partition-0.json').read_text())

paired_parameters=dict(pairs=[dict(id='unit',left='all.A',right='all.B',left_scale='1',right_scale='1')],selection=[],support='joint_reported',normalization=dict(method='Declared identity scales'),output_units='counts',workers=1)
control_parameters=dict(fields=['all.A'],endpoint='all.B',controls=[],selection=[],support='joint_reported',metric='identity',intercept=True,workers=2,gram_fields=['all.A'])
control_recipe=dict(schema='gadget.pipeline.default.v1',profile='rna_seq',question='Pair full residual fields outside the declared intercept',goals=[dict(id='controlled',method='controlled_relations',workers=2,parameters=control_parameters)])
prepared,index=recipe_run('controlled-recipe',control_recipe)
controlled_result=leaf('controlled-recipe','controlled')['result']
assert list(map(Fraction,controlled_result['fields'][0]['residual']))==[-1,-3,4]
assert Fraction(controlled_result['fields'][0]['pairing'])==30
assert controlled_result['residual_gram']['values']==[['26']]
assert controlled_result['execution']['gram_preparations']==1 and controlled_result['execution']['requested_workers']==2
expression=dict(schema='gadget.pipeline.default.v1',profile='rna_seq',question='Retain each matched expression difference',goals=[dict(method='expression_contrast',parameters=paired_parameters,normalization_evidence_field='all.A')])
prepared,index=recipe_run('expression-recipe',expression)
assert [Fraction(f['mean_difference']) for f in leaf('expression-recipe','expression_contrast')['result']['features']]==[-3,-3,-4]
assert len(index['results'])==1
bound=json.loads((root/'expression-recipe/expression_contrast/request.json').read_text())
assert bound['query']['normalization']['evidence'][0]['document']=='assay_declaration'
assert 'evidence' not in prepared['steps'][0]['request']['query']['normalization']

# A declared pseudobulk map preserves supported zero and propagates missing inputs.
cell_axis=source('cells.tsv','id\tx\ty\tcoordinate_frame\tcoordinate_units\nA\t0\t0\tsection\tum\nB\t1\t0\tsection\tum\nC\t2\t0\tsection\tum\nD\t4\t0\tsection\tum\n')
cell_matrix=source('cell-counts.tsv','id\tA\tB\tC\tD\ng1\t1\t3\t8\t4\ng2\t0\t0\t3\t1\ng3\tNA\t2\t5\t5\n')
cell_assay={**assay,'id':'cells','source':cell_matrix,'observations':cell_axis,'source_namespace':'889100','implicit_entries':'unobserved','missing_tokens':['NA']}
cell_spec={**spec,'id':'cell-recipes','assays':[cell_assay]}
cell_workspace=root/'cell-workspace'
run('import-assays','--spec',save('cell-import.json',cell_spec),'--out',cell_workspace)
cell_manifest=cell_workspace/'workspace.json'
groups=request('aggregate_assay',assay='cells',feature_groups=[dict(id=g,members=[dict(id=g,coefficient='1')]) for g in ['g1','g2','g3']],observation_groups=[dict(id=name,members=[dict(id=c,coefficient='1') for c in cells]) for name,cells in [('left',['A','B']),('right',['C','D'])]],map_method='Declared donor sums',map_evidence='The explicitly selected cell identities',support_policy='require_all_nonzero_inputs')
params={k:v for k,v in groups.items() if k not in ['schema','operation']}
params['count_policy']='disjoint_count_sum'
pseudobulk=dict(schema='gadget.pipeline.default.v1',profile='single_cell',question='Compare declared cell groups with complete measurement support',goals=[dict(id='grouped',method='aggregation',parameters=params,output=dict(id='groups',source_namespace='889200',design=design,channel='E')),dict(id='paired',method='expression_contrast',input='grouped',normalization_evidence_field='all.left',parameters={**paired_parameters,'pairs':[dict(id='unit',left='all.left',right='all.right',left_scale='1',right_scale='1')]}),dict(id='support',method='coverage',input='grouped',parameters=dict(fields=['all.left']))])
prepared,index=recipe_run('pseudobulk-recipe',pseudobulk,cell_manifest)
assert prepared['steps'][1]['request']['assay']=='groups'
assert index['results'][1]['source_digest']==index['results'][0]['derived_digest']
result=leaf('pseudobulk-recipe','paired')['result']
assert result['reported']==2 and result['missing']==1
assert [f['entity'] for f in result['features']]==['g1','g2']
assert [Fraction(f['mean_difference']) for f in result['features']]==[-8,-4]
run('replay-derivation','--workspace',root/'pseudobulk-recipe/grouped/workspace/workspace.rfcxs','--assay','groups','--out',root/'pseudobulk-replay')

for name,mutate in [
    ('overlapping-cell-groups',lambda p:p['observation_groups'][1]['members'].append(dict(id='A',coefficient='1'))),
    ('averaged-cell-group',lambda p:p['observation_groups'][0]['members'][0].update(coefficient='1/2')),
    ('merged-features',lambda p:p['feature_groups'][0]['members'].append(dict(id='g2',coefficient='1'))),
]:
    bad=copy.deepcopy(pseudobulk);mutate(bad['goals'][0]['parameters'])
    run('pipeline','--workspace',cell_manifest,'--request',save(name+'.json',bad),'--out',root/name,fail=True)
    assert not (root/name).exists()

# A general weighted map remains available with its derived measurement meaning.
general=copy.deepcopy(pseudobulk)
general['goals']=general['goals'][:1]
general['goals'][0].pop('output')
general['goals'][0]['parameters'].pop('count_policy')
general['goals'][0]['parameters']['observation_groups'][0]['members'][0]['coefficient']='1/2'
prepared,index=recipe_run('general-aggregation',general,cell_manifest)
generic_result=json.loads((root/'general-aggregation/grouped/result.json').read_text())['result']
assert generic_result['quantity']=='declared_linear_aggregate_of_read_count'
assert Fraction(generic_result['rows'][0][0])==Fraction(7,2)
for label,value in [('negative','-1'),('fractional','1/2')]:
    changed_source=source(label+'-counts.tsv','id\tA\tB\tC\tD\ng1\t'+value+'\t3\t8\t4\ng2\t0\t0\t3\t1\ng3\tNA\t2\t5\t5\n')
    changed_assay={**cell_assay,'source':changed_source,'arithmetic':'exact_rational'}
    changed_workspace=root/(label+'-workspace')
    run('import-assays','--spec',save(label+'-import.json',{**cell_spec,'assays':[changed_assay]}),'--out',changed_workspace,fail=True)
    assert not changed_workspace.exists()

# Physical-neighborhood materialization feeds the actual complete Green solve.
spatial=dict(schema='gadget.pipeline.default.v1',profile='spatial',question='Resolve a retained count field on the declared physical neighborhood',goals=[dict(id='neighbors',method='spatial_relations',parameters=dict(assay='cells',coordinates=['x','y'],frame='section',units='um',metric_diagonal=['1','1'],radius_squared='1',relation_type='physical.neighborhood',evidence='Declared coordinate frame and radius'),output=dict(id='physical',source_namespace='889300',design=design,channel='E')),dict(id='response',method='green_action',input='neighbors',parameters=dict(context='all',field='all.g1',**{'lambda':'1'}))])
# Spatial operators require an observation carrier, rather than feature rows.
spatial_assay={**cell_assay,'id':'cells','carrier':'observations','source_namespace':'889400'}
spatial_workspace=root/'spatial-workspace'
run('import-assays','--spec',save('spatial-import.json',{**cell_spec,'assays':[spatial_assay]}),'--out',spatial_workspace)
prepared,index=recipe_run('spatial-recipe',spatial,spatial_workspace/'workspace.json')
response=leaf('spatial-recipe','response')
assert response['residual_verified']
assert [Fraction(c['value']['rational']) for c in response['result']]==[Fraction(19,8),Fraction(15,4),Fraction(47,8),4]

# Exact chromatin interval differences cover every elementary segment.
ref=source('reference.fa','>chrA\nAAAAAAAAAA\n')
frame=dict(assembly='authored',sequence_sha256=ref['sha256'],coordinate_convention='zero_based_half_open',contigs=[dict(id='chrA',length=10)])
chromatin_resources=[dict(id=name,format='bed',source=source(name,body),method='Authored weighted intervals',reference=frame,weight_column=4) for name,body in [('treated.bed','chrA\t0\t4\ta\t2\t+\nchrA\t2\t6\tb\t1\t+\n'),('control.bed','chrA\t1\t5\tc\t1\t+\n')]]
chromatin_workspace=root/'chromatin-workspace'
run('import-assays','--spec',save('chromatin-import.json',{**spec,'resources':chromatin_resources}),'--out',chromatin_workspace)
chromatin=dict(schema='gadget.pipeline.default.v1',profile='atac_seq',question='Retain scaled treated-minus-control interval fields',goals=[dict(id='contrast',method='chromatin_contrast',parameters=dict(samples=[dict(id=n,resource=n,scale='1',scale_provenance='Declared identity scale') for n in ['treated.bed','control.bed']],contrasts=[dict(id='pair',treated='treated.bed',control='control.bed')],strand_policy='separate',thresholds=['0','1','2'],absence_policy='measured_zero_within_union'),output=dict(id='contrast',source_namespace='889500',design=design,channel='observation')),dict(method='coverage',input='contrast',parameters=dict(fields=['all.pair']))])
prepared,index=recipe_run('chromatin-recipe',chromatin,chromatin_workspace/'workspace.json')
result=json.loads((root/'chromatin-recipe/contrast/result.json').read_text())
# Read the retained contrast on its original physical membership complex.
contrast_archive=root/'chromatin-recipe/contrast/workspace/workspace.rfcxs'
assert contrast_archive.exists()
assert not (contrast_archive.parent/'workspace.json').exists()
contrast_result=run('record-view','--input',contrast_archive,'--assay','contrast')['result']
expected=[Fraction(2),Fraction(1),Fraction(2),Fraction(0),Fraction(1)]
assert [Fraction(value) for value in contrast_result['contrasts'][0]['values']]==expected

# Matched-unit counts are bounded by the admitted data, not a 32-pair constant.
pair_count=65
samples=[f'{side}{i}' for i in range(pair_count) for side in ['L','R']]
many_observations=source('many-observations.tsv','id\n'+'\n'.join(samples)+'\n')
many_counts=source('many-counts.tsv','id\t'+'\t'.join(samples)+'\ng1\t'+'\t'.join(str(i+1) if side=='L' else '0' for i in range(pair_count) for side in ['L','R'])+'\n')
many_features=source('many-features.tsv','id\ng1\n')
many_assay={**assay,'id':'many','source':many_counts,'features':many_features,'observations':many_observations,'source_namespace':'889600'}
many_workspace=root/'many-workspace'
run('import-assays','--spec',save('many-import.json',{**spec,'assays':[many_assay]}),'--out',many_workspace)
many={**expression,'goals':[dict(method='expression_contrast',normalization_evidence_field='all.L0',parameters={**paired_parameters,'pairs':[dict(id=str(i),left='all.L'+str(i),right='all.R'+str(i),left_scale='1',right_scale='1') for i in range(pair_count)]})]}
prepared,index=recipe_run('many-pairs-recipe',many,many_workspace/'workspace.json')
f=leaf('many-pairs-recipe','expression_contrast')['result']['features'][0]
assert list(map(Fraction,f['delta']))==list(range(1,pair_count+1))
assert Fraction(f['mean_difference'])==33
assert Fraction(f['heterogeneity_quadrance'])==sum((Fraction(i)-33)**2 for i in range(1,66))

# Selection records are source fields, and need not be fewer than coordinates.
# All 65 indicators intentionally select the same single observed coordinate.
one_observation=source('one-observation.tsv','id\tselected\nunit\tyes\n')
transposed_features=source('transposed-features.tsv','id\n'+'\n'.join(samples)+'\n')
transposed_matrix=source('transposed-counts.tsv','id\tunit\n'+'\n'.join(name+'\t'+(str(i+1) if side=='L' else '0') for i in range(pair_count) for side,name in [('L','L'+str(i)),('R','R'+str(i))])+'\n')
selection_fields=[f'selected_{i}' for i in range(pair_count)]
selected_assay={**assay,'id':'selected','source':transposed_matrix,'features':transposed_features,'observations':one_observation,'carrier':'observations','source_namespace':'889700','observation_fields':[dict(id=name,column='selected',kind='level_indicator',level='yes',quantity='declared_membership',units='dimensionless',method='Authored source membership indicator',missing_tokens=[]) for name in selection_fields]}
selected_workspace=root/'selected-workspace'
run('import-assays','--spec',save('selected-import.json',{**spec,'assays':[selected_assay]}),'--out',selected_workspace)
selected_design=copy.deepcopy(many)
selected_design['goals'][0]['parameters']['selection']=['all.'+name for name in selection_fields]
prepared,index=recipe_run('many-selection-recipe',selected_design,selected_workspace/'workspace.json')
selected_result=leaf('many-selection-recipe','expression_contrast')['result']
assert selected_result['reported']==1 and selected_result['features'][0]['entity']=='unit'
assert Fraction(selected_result['features'][0]['mean_difference'])==33

for name,settings,source_manifest in [
    ('recipe-wrong-profile',{**pseudobulk,'profile':'rna_seq'},cell_manifest),
    ('recipe-future-input',{**expression,'goals':[{**expression['goals'][0],'input':'later'}]},manifest),
    ('recipe-ineligible-output',{**expression,'goals':[{**expression['goals'][0],'output':dict(id='no',source_namespace='99',design=design,channel='E')}]},manifest),
    ('recipe-conflicting-evidence',{**expression,'goals':[{**expression['goals'][0],'parameters':{**paired_parameters,'normalization':dict(method='identity',evidence=[])}}]},manifest),
]:
    run('pipeline-prepare','--workspace',source_manifest,'--request',save(name+'.json',settings),fail=True)
bad=copy.deepcopy(pseudobulk);bad['goals'][1]['normalization_evidence_field']='all.absent'
run('pipeline','--workspace',cell_manifest,'--request',save('recipe-missing-evidence.json',bad),'--out',root/'recipe-missing-evidence',fail=True)
assert not (root/'recipe-missing-evidence').exists() and not list(root.glob('.gadget-pipeline-*'))
# Question-driven selection uses immutable recipe bindings, not model-written parameters.
def intent(settings, purposes=None):
    entries=[]
    for i, original in enumerate(settings['goals']):
        goal=copy.deepcopy(original)
        identity=goal.pop('id', goal['method'])
        entries.append(dict(id=identity, purpose=(purposes or {}).get(identity, identity), goal=goal))
    return dict(schema='gadget.pipeline.intent.v1',profile=settings['profile'],recipes=entries)

def proposal_for(context, source_manifest, goals=None, reason=None):
    context_path=save('intent-active.json',context)
    catalog=run('recipe-catalog','--workspace',source_manifest,'--intent',context_path)
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',action='clarify' if reason else 'recipes',
        source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'])
    proposal.update(dict(reason=reason) if reason else dict(goals=goals))
    return context_path,catalog,proposal

context=intent(pseudobulk)
context_path,catalog,selected=proposal_for(context,cell_manifest,['paired'])
assert len(catalog['recipes'])==3 and catalog['recipes'][1]['dependencies']==['grouped']
answer=run('recipe-proposal','--workspace',cell_manifest,'--intent',context_path,
    '--request',save('recipe-selection.json',selected),'--prompt','Compare the declared matched groups',
    '--out',root/'selected-recipe')
assert answer['included_goals']==['grouped','paired']
assert [node['id'] for node in answer['pipeline']['steps']]==['grouped','paired']
assert not (root/'selected-recipe/support').exists()
assert [Fraction(f['mean_difference']) for f in leaf('selected-recipe','paired')['result']['features']]==[-8,-4]
index=answer['execution']; receipt=root/'selected-recipe'/index['planning']
assert index['planning_sha256']=='sha256:'+hashlib.sha256(receipt.read_bytes()).hexdigest()
record=run("record-view","--input",receipt,"--member","result")
assert record['intent']==context and record['proposal']==selected and record['question']=='Compare the declared matched groups'
assert 'execution' not in record
assert index['results'][1]['source_digest']==index['results'][0]['derived_digest']

spatial_context=intent(spatial)
spatial_path,_,spatial_proposal=proposal_for(spatial_context,spatial_workspace/'workspace.json',['response'])
answer=run('recipe-proposal','--workspace',spatial_workspace/'workspace.json','--intent',spatial_path,
    '--request',save('spatial-selection.json',spatial_proposal),'--prompt','Compute the declared spatial Green response',
    '--out',root/'selected-spatial')
assert answer['included_goals']==['neighbors','response']
assert [Fraction(c['value']['rational']) for c in leaf('selected-spatial','response')['result']]==[Fraction(19,8),Fraction(15,4),Fraction(47,8),4]

# Foreign source, stale design, invented parameters/names and duplicates fail before publication.
context_path,_,selected=proposal_for(context,cell_manifest,['paired'])
for label,change in [('foreign-source',dict(source_digest='sha256:foreign')),
        ('foreign-intent',dict(intent_sha256='sha256:foreign')),('invented-goal',dict(goals=['invented'])),
        ('duplicate-goal',dict(goals=['paired','paired'])),('injected-parameters',dict(parameters={'pairs':[]}))]:
    bad={**selected,**change}
    run('recipe-proposal','--workspace',cell_manifest,'--intent',context_path,
        '--request',save(label+'.json',bad),'--prompt','Compare declared groups','--out',root/label,fail=True)
    assert not (root/label).exists()
stale=copy.deepcopy(context);stale['recipes'][0]['goal']['parameters']['map_evidence']='Changed declaration'
run('recipe-proposal','--workspace',cell_manifest,'--intent',save('stale-intent.json',stale),
    '--request',save('old-selection.json',selected),'--prompt','Compare declared groups','--out',root/'stale-run',fail=True)
assert not (root/'stale-run').exists()
run('recipe-catalog','--workspace',cell_manifest,'--intent',context_path,'--profile','rna_seq',fail=True)

for label,change in [('forward-edge',lambda c:c['recipes'][0]['goal'].update(input='paired')),
        ('cyclic-edge',lambda c:c['recipes'][0]['goal'].update(depends_on=['paired'])),
        ('duplicate-binding',lambda c:c['recipes'][1].update(id='grouped')),
        ('hidden-identity',lambda c:c['recipes'][0]['goal'].update(id='hidden'))]:
    bad=copy.deepcopy(context);change(bad)
    run('recipe-catalog','--workspace',cell_manifest,'--intent',save(label+'.json',bad),fail=True)

missing=dict(schema='gadget.pipeline.intent.v1',profile='rna_seq',recipes=[])
missing_path,_,clarification=proposal_for(missing,manifest,reason='Declare the actual matched sample identities and scales')
answer=run('recipe-proposal','--workspace',manifest,'--intent',missing_path,'--request',save('missing-design-proposal.json',clarification),
    '--prompt','Compare treatment expression','--out',root/'missing-design-run')
assert answer['status']=='clarification' and not (root/'missing-design-run').exists()

# Selection/dependency closure has no former 16-step or 32-field catalog cutoff.
chain=dict(schema='gadget.pipeline.intent.v1',profile='rna_seq',assay='counts',recipes=[
    dict(id='check_'+str(i),purpose='Support checkpoint '+str(i),goal=dict(method='coverage',
         parameters=dict(fields=['all.A']),**(dict(depends_on=['check_'+str(i-1)]) if i else {}))) for i in range(300)])
chain_path,chain_catalog,chain_proposal=proposal_for(chain,manifest,['check_299'])
assert len(chain_catalog['recipes'])==300
answer=run('recipe-proposal','--workspace',manifest,'--intent',chain_path,'--request',save('chain-proposal.json',chain_proposal),
    '--prompt','Run the final declared support checkpoint','--out',root/'chain-run')
assert len(answer['included_goals'])==len(answer['execution']['results'])==300
assert all(leaf('chain-run','check_'+str(i))['result']['reported']==3 for i in range(300))

# Numerical owners still reject invalid source fields; atomic execution leaves no partial publication.
invalid=dict(schema='gadget.pipeline.intent.v1',profile='rna_seq',assay='counts',recipes=[
    dict(id='absent',purpose='Invalid field fixture',goal=dict(method='coverage',parameters=dict(fields=['all.absent'])))])
invalid_path,_,invalid_proposal=proposal_for(invalid,manifest,['absent'])
run('recipe-proposal','--workspace',manifest,'--intent',invalid_path,'--request',save('invalid-field-proposal.json',invalid_proposal),
    '--prompt','Inspect the selected field','--out',root/'invalid-field-run',fail=True)
assert not (root/'invalid-field-run').exists() and not list(root.glob('.gadget-pipeline-*'))

control_bed = source('control.bed', 'chr\t0\t8\tbase\t1\t+\nchr\t8\t12\tlast\t2\t+\n')
treated_bed = source('treated.bed', 'chr\t0\t4\tfirst\t3\t+\nchr\t4\t8\tzero\t0\t+\nchr\t8\t12\tlast\t2\t+\n')
reference = dict(assembly='fixture', sequence_sha256='sha256:'+'0'*64,
    coordinate_convention='zero_based_half_open', contigs=[dict(id='chr', length=12)])
coverage_spec = copy.deepcopy(spec)
coverage_spec['id'] = 'native-producer-fixture'
coverage_spec['resources'] = [dict(id=name, source=data, format='bed', method='Fixture original intervals',
    weight_column=4, reference=reference) for name, data in [('control', control_bed), ('treated', treated_bed)]]
coverage_manifest = root/'coverage-input'/'workspace.rfcxs'
run('import-assays', '--spec', save('coverage-import.json', coverage_spec), '--out', coverage_manifest.parent)
coverage_recipe = request('controlled_coverage',
    samples=[dict(id=name, resource=resource, scale='1', scale_provenance='Fixture')
        for name, resource in [('base', 'control'), ('treated', 'treated')]],
    contrasts=[dict(id='change', terms=[dict(sample='treated', coefficient='1'),
        dict(sample='base', coefficient='-1')], evidence='Fixture')],
    controls=['base'], intercept=True, metric='reference_length', strand_policy='separate',
    thresholds=['1'], absence_policy='measured_zero_within_union', units='count', evidence='Fixture', source_kernel=True)
native_derive = dict(schema='gadget.assays.materialize.v1', id='coverage', source_namespace='889001',
    design=design, recipe=coverage_recipe, channel='observation')
native_binding = dict(schema='gadget.assays.bindings.v1', id='bound', design=design,
    assays=[dict(input='native', assay='coverage', id='alias')], resources=[], links=[],
    method='Retain the selected native owner', evidence='Fixture')
native_program = dict(schema='gadget.pipeline.design.v1', profile='custom',
    question='Retain a native producer across branches and a bound collection', concurrent_jobs=2, steps=[
        dict(id='native', type='derive', input='source', request=native_derive),
        dict(id='bound', type='bind', input='source', request=native_binding),
        dict(id='first', type='query', input='native', request=request('native_query', assay='coverage',
            query=dict(operation='field_summary', field='all.sample:0', support='reported_only'))),
        dict(id='second', type='query', input='bound', request=request('native_query', assay='alias',
            query=dict(operation='field_summary', field='all.sample:0', support='reported_only'))),
        dict(id='ordered', type='query', input='source', depends_on=['native'],
            request=request('native_query', assay='counts', query=dict(operation='coverage', field='all.A'))),
    ])
native_index = run('pipeline', '--workspace', coverage_manifest, '--request',
    save('native-program.json', native_program), '--out', root/'native-program')
assert Fraction(leaf('native-program','first')['result']['sum']['rational']) == 4
assert leaf('native-program','first')['result'] == leaf('native-program','second')['result']
assert native_index['results'][1]['source_digest'] != native_index['results'][0]['derived_digest']
producer_schema = run('record-view', '--input', root/'native-program'/'native/result.rfcxs', '--member', 'schema')
binding_schema = run('record-view', '--input', root/'native-program'/'bound/result.rfcxs', '--member', 'schema')
assert producer_schema == 'gadget.operation.native.v1' and binding_schema == 'gadget.assay.workspace.native.v2'
for item in native_index['results']:
    data = (root/'native-program'/item['result']).read_bytes()
    assert item['result_sha256'] == 'sha256:'+hashlib.sha256(data).hexdigest()
failed_program = copy.deepcopy(native_program)
failed_program['steps'][2]['request']['query']['field'] = 'all.absent'
run('pipeline', '--workspace', coverage_manifest, '--request', save('native-failure.json', failed_program),
    '--out', root/'native-failure', fail=True)
assert not (root/'native-failure').exists() and not list(root.glob('.gadget-pipeline-*'))

sequence_reference = source('indexed-reference.fa', '>reference\nACGTAGGTACGT\n')
sequence_spec = copy.deepcopy(spec)
sequence_spec['id'] = 'indexed-pipeline-contracts'
sequence_spec['resources'].append(dict(id='indexed-reference', format='fasta',
    source=sequence_reference, method='Authored placement source'))
sequence_workspace = root/'sequence-workspace'
run('import-assays', '--spec', save('sequence-import.json', sequence_spec), '--out', sequence_workspace)
sequence_manifest = sequence_workspace/'workspace.json'
parameters = dict(query_resource='reads', reference_resource='indexed-reference', strands=['+'],
    maximum_edits=0, scoring='unit_edit_cost', evidence='Complete authored exact placements')
common = None
for batch, width in [(1, 1), (4, 20)]:
    execution = dict(read_partition_records=batch, reference_block_bases=width)
    settings = dict(schema='gadget.pipeline.default.v1', profile='sequence',
        question='Retain exact placements across execution blocks', goals=[
            dict(id='placements', method='indexed_sequence_alignment', parameters={**parameters,
                'execution': execution})])
    name = f'indexed-pipeline-{batch}-{width}'
    plan = run('pipeline-prepare', '--workspace', sequence_manifest,
        '--request', save(name+'-settings.json', settings))
    assert plan['steps'][0]['request']['execution'] == execution
    run('pipeline', '--workspace', sequence_manifest, '--request', save(name+'-plan.json', plan),
        '--out', root/name)
    result = json.loads((root/name/'placements/result.json').read_text())['result']
    actual = {(row['query_id'], row['reference_id'], row['strand'], int(row['start']),
        int(row['end']), int(row['minimum_edits'])) for row in result['matches']}
    assert actual == {('r1', 'reference', '+', 0, 4, 0), ('r1', 'reference', '+', 8, 12, 0),
        ('r2', 'reference', '+', 4, 8, 0)}
    assert result['work']['peak_partition_reads'] <= batch
    assert result['work']['peak_normalized_reference_block_bases'] <= width
    payload = {key: value for key, value in result.items() if key != 'work'}
    if common is None:
        common = payload
    else:
        assert payload == common
    retained = run('record-view', '--input', root/name/'placements/result.rfcxs', '--member', 'request')
    assert retained['execution'] == execution
bad = copy.deepcopy(settings)
bad['goals'][0]['parameters']['execution']['reference_block_bases'] = 0
bad_plan = run('pipeline-prepare', '--workspace', sequence_manifest,
    '--request', save('invalid-indexed-execution.json', bad))
run('pipeline', '--workspace', sequence_manifest,
    '--request', save('invalid-indexed-plan.json', bad_plan),
    '--out', root/'invalid-indexed-execution', fail=True)
assert not (root/'invalid-indexed-execution').exists()

print(json.dumps(dict(status='PASS',command_checks=checks,derived_source_binding=True,atomic_failure=True,model_proposal_admission=True,assay_recipes=['rna_seq','single_cell','spatial','atac_seq'],matched_units_above_former_cap=pair_count,recipe_planning=True,dependency_chain=300,native_producer_branches=True)))
