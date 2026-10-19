"""Original count fields, conventional models, ownership and native pipelines."""
import copy
from fractions import Fraction as Q
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    source = Path(__file__).resolve().parents[1]/'examples/bioinformatics/relational-ensembles/replay.py'
    spec = importlib.util.spec_from_file_location('count_source', source)
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    ids = [f'sample-{i}' for i in range(6)]
    measurements = dict(gene_a=[1,2,6,4,9,14], gene_b=[2,4,8,6,10,16],
                        zero=[0]*6, missing=[1,2,None,4,5,6], invalid=[1,2,3,4,5,-1],
                        large=[9007199254740993,1,2,3,4,5], separated=[0,0,0,4,9,14])
    demo.stage(root/'inputs', ids, measurements, [], [], units='counts', quantity='declared_count')
    matrix=root/'inputs/measurements.tsv'
    matrix.write_text(matrix.read_text().replace('\tNone','\tNA'))
    declaration=json.loads((root/'inputs/import.json').read_text())
    import hashlib
    declaration['assays'][0]['source']['sha256']='sha256:'+hashlib.sha256(matrix.read_bytes()).hexdigest()
    (root/'inputs/import.json').write_text(json.dumps(declaration))
    def save(name, value):
        path = root/name
        path.write_text(json.dumps(value)+'\n')
        return path
    def run(*args, fail=False):
        result = subprocess.run(list(map(str,args)), capture_output=True, text=True, timeout=120)
        assert bool(result.returncode) == fail, result.stderr
        return json.loads(result.stdout) if result.stdout.strip() and not fail else result
    run(binary, 'import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    workspace = root/'workspace/workspace.rfcxs'
    order = [5,0,4,1,3,2]
    request = dict(schema='gadget.statistics.count.request.v1', operation='count_model_analysis', source_digest='bound by native test',
        samples=[dict(id=ids[i], experimental_unit='unit-'+str(i), coordinate=i) for i in order],
        features=[dict(id=name, field='all.'+name) for name in list(measurements)[:6]], counts=dict(quantity='declared_count', units='counts'),
        design=dict(columns=['intercept','treatment'], rows=[['1',str(int(i>=3))] for i in order], contrast=['0','1']),
        model=dict(family='poisson', dispersion=None), normalization=dict(method='none'), exposure=None,
        inference=dict(method='wald', confidence='0.95', reduced_columns=None), effect_prior=dict(mean='0',variance='1'),
        correction='benjamini_hochberg', outputs=dict(fitted_means=True), workers=1,
        limits=dict(samples=6, features=6, design_columns=2, iterations=200, line_search_steps=40,
                    numerical_work=100000000, declared_buffer_bytes=268435456, request_bytes=1048576))
    def native(name, query):
        return run(consumer, workspace, save(name+'.json',query), root/(name+'.rfcxs'))
    poisson = native('poisson',request)
    rows = poisson['features']
    assert [row['id'] for row in rows] == list(measurements)[:6]
    assert [row['status'] for row in rows] == ['converged','converged','all_zero_counts','missing_counts','invalid_counts','invalid_counts']
    assert rows[3]['missing'] == 1 and rows[3]['p_value'] is None
    assert rows[4]['invalid'] == rows[5]['invalid'] == 1
    assert poisson['tested_family_size'] == 6
    def near(a,b):
        assert math.isclose(float(Q(a)),b,rel_tol=2e-7,abs_tol=2e-8), (a,b)
    near(rows[0]['effect'],math.log(3))
    near(rows[0]['variance'],1/9+1/27)
    assert [sample['coordinate'] for sample in poisson['samples']] == order
    near(rows[0]['effect_posterior']['mean'], math.log(3)/(1+1/9+1/27))
    request['source_digest'] = poisson['source_digest']
    save('cli.json',request)
    cli = run(binary,'statistics','--workspace',workspace,'--assay','measurements','--request',root/'cli.json','--out-record',root/'cli.rfcxs')
    assert cli['features'] == poisson['features']
    for policy in ('none','bonferroni','benjamini_yekutieli'):
        q = copy.deepcopy(request);q['correction']=policy
        output = native(policy,q)
        assert output['features'][3]['adjusted_p_value'] is None
        near(output['features'][0]['p_value'],float(Q(rows[0]['p_value'])))
    q = copy.deepcopy(request)
    q['model'] = dict(family='negative_binomial',dispersion=dict(method='fixed',value='.2'))
    nb = native('nb2',q)
    near(nb['features'][0]['effect'],math.log(3))
    near(nb['features'][0]['variance'],1/9+1/27+2*.2/3)
    q['inference'] = dict(method='likelihood_ratio',confidence='.95',reduced_columns=[0])
    lrt = native('lrt',q)
    assert float(Q(lrt['features'][0]['likelihood_ratio'])) > 0
    q['outputs']['fitted_means'] = False
    output = native('compact',q)
    assert 'fitted_means' not in output['features'][0]
    q=copy.deepcopy(request);q['features']=[dict(id='separated',field='all.separated')]
    output=native('boundary',q)
    assert output['features'][0]['status']=='boundary_likelihood' and output['features'][0]['p_value'] is None
    for method in ('profile','pooled'):
        q = copy.deepcopy(request);q['features']=q['features'][:2]
        d = dict(method=method,minimum='.0001',maximum='3',log_tolerance='.00001',grid_points=9,iterations=80,log_prior=None)
        if method == 'pooled':d['reference_fields']=['all.gene_a','all.gene_b']
        q['model'] = dict(family='negative_binomial',dispersion=d)
        output = native(method,q)
        assert all(row['status']=='converged' for row in output['features']), output
        if method=='pooled':assert output['pooled_dispersion']['converged']
    q=copy.deepcopy(request);q['features']=q['features'][:2]
    q['exposure']=dict(values=['2']*6,meaning='Declared equal observation duration')
    output=native('exposure',q)
    near(output['features'][0]['coefficients'][0],math.log(3)-math.log(2))
    near(output['features'][0]['effect'],math.log(3))
    q['inference']=dict(method='likelihood_ratio',confidence='.95',reduced_columns=[])
    output=native('offset-null',q)
    assert output['likelihood_ratio_degrees']==2
    assert output['features'][0]['reduced_coefficients']==[]
    q=copy.deepcopy(request);q['features']=q['features'][:2]
    q['model']=dict(family='negative_binomial',dispersion=dict(method='profile',minimum='.0001',maximum='3',log_tolerance='.00001',grid_points=9,iterations=80,log_prior=dict(mean='-1',variance='2')))
    output=native('regularized-dispersion',q)
    assert all(row['status']=='converged' for row in output['features']),output
    q['model']['dispersion']['log_prior']=None
    q['inference']=dict(method='likelihood_ratio',confidence='.95',reduced_columns=[0])
    output=native('profile-lrt',q)
    assert all(row['status']=='converged' for row in output['features']),output
    assert all(row['reduced_dispersion_profile']['converged'] for row in output['features'])
    q=copy.deepcopy(request);q['samples']=q['samples'][:4];q['design']['rows']=q['design']['rows'][:4];q['features']=q['features'][:2]
    output=native('selection',q)
    assert len(output['samples'])==4 and len(output['features'][0]['fitted_means'])==4
    q=copy.deepcopy(request);q['limits']['iterations']=1
    output=native('nonconverged',q)
    assert output['features'][0]['status']=='iteration_limit' and output['features'][0]['p_value'] is None
    for method in ('declared_factors','median_ratio'):
        q = copy.deepcopy(request)
        q['normalization'] = dict(method=method,factors=['1','2','1','2','1','2'],center=True) if method=='declared_factors' else dict(method=method,reference_fields=['all.gene_a','all.gene_b'])
        output=native(method,q)
        factors = [float(Q(v)) for v in output['normalization_factors']]
        near(str(math.prod(factors)),1)
        if method=='median_ratio':
            geometric = [math.exp(sum(math.log(v) for v in measurements[name])/6) for name in ('gene_a','gene_b')]
            raw=[sum(measurements[name][i]/geometric[j] for j,name in enumerate(('gene_a','gene_b')))/2 for i in order]
            center=math.exp(sum(math.log(v) for v in raw)/6)
            for a,b in zip(factors,raw):near(str(a),b/center)
    for name, change in [
        ('source',lambda q:q.update(source_digest='sha256:'+'0'*64)),
        ('duplicate-unit',lambda q:q['samples'][1].update(experimental_unit=q['samples'][0]['experimental_unit'])),
        ('coordinate',lambda q:q['samples'][0].update(coordinate=100)),
        ('quantity',lambda q:q['counts'].update(quantity='normalized_signal')),
        ('rank',lambda q:q['design'].update(rows=[['1','1']]*6)),
        ('memory',lambda q:q['limits'].update(declared_buffer_bytes=1024)),
        ('work',lambda q:q['limits'].update(numerical_work=1)),
        ('unknown',lambda q:q.update(unsupported=True))]:
        bad=copy.deepcopy(request);change(bad)
        run(binary,'statistics','--workspace',workspace,'--assay','measurements','--request',save('bad-'+name+'.json',bad),'--out-record',root/('bad-'+name+'.rfcxs'),fail=True)
        assert not (root/('bad-'+name+'.rfcxs')).exists()
    parameters={k:v for k,v in request.items() if k not in ('schema','operation','source_digest')}
    defaults=dict(schema='gadget.pipeline.default.v1',profile='rna_seq',question='Compare declared independent count fields',assay='measurements',
                  goals=[dict(method='count_model_analysis',id='counts',parameters=parameters)])
    plan=run(binary,'pipeline-prepare','--workspace',workspace,'--request',save('default.json',defaults))
    assert plan['steps'][0]['type']=='statistics'
    save('plan.json',plan)
    run(binary,'pipeline','--workspace',workspace,'--request',root/'plan.json','--out',root/'pipeline')
    result=json.loads((root/'pipeline/counts/result.json').read_text())
    assert result['features']==poisson['features']
    def sha(path):return 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()
    gmt='group-a\tDeclared set\tG0\tG1\n'
    mapping='feature_id\tmember\tevidence\n'+''.join(f'{name}\tG{i}\tDeclared correspondence\n' for i,name in enumerate(measurements))
    sets=dict(schema='gadget.pathways.source.v1',id='declared-sets',release='authored.v1',identifier_namespace='declared:gene',format='gmt',data=gmt,sha256='sha256:'+hashlib.sha256(gmt.encode()).hexdigest(),
              mapping=dict(id='declared-map',release='authored.v1',origin='explicit_declared',feature_namespace='declared:count_feature',identifier_namespace='declared:gene',format='tsv',data=mapping,sha256='sha256:'+hashlib.sha256(mapping.encode()).hexdigest()))
    run(binary,'native-record','--input',save('sets.json',sets),'--out',root/'sets.rfcxs')
    query=dict(schema='gadget.pathways.overrepresentation.v1',operation='pathway_overrepresentation',gene_sets_sha256=sha(root/'sets.rfcxs'),feature_namespace='declared:count_feature',
               features_path=['result','view','features'],feature_id_path=['id'],background_indices=list(range(6)),selected_indices=[0],
               mapping_policy=dict(ambiguous='refuse',unmapped='refuse',many_to_one='refuse',duplicate_occurrences='retain_count_once'),set_ids=['group-a'],size_filter=dict(minimum=1,maximum=100,basis='members_in_background'),
               membership_policy=dict(duplicates='retain_count_once',empty='retain_exclude',outside_background='retain_exclude'),correction='none',selection_method='Declared feature index for this control',null_assumptions=['Uniform conditional subset selection'],
               limits=dict(source_bytes=1048576,mapping_rows=1000,membership_cells=1000,sets=10,features=10,field_bytes=4096,kernel_work=1000000,archive_bytes=8388608,declared_buffer_bytes=33554432))
    operation=dict(schema='gadget.result-operation.v1',operation='pathway_overrepresentation',analysis_binding='actual_result',analysis_bytes=8388608,source_record=dict(path=str(root/'sets.rfcxs'),sha256=sha(root/'sets.rfcxs'),max_bytes=8388608),query=query)
    design={k:v for k,v in plan.items() if k not in ('plan_sha256','source_digest')}
    design['schema']='gadget.pipeline.design.v1'
    design['steps'].append(dict(id='sets',type='result',input='source',result_input='counts',request=operation))
    chained=run(binary,'pipeline-prepare','--workspace',workspace,'--request',save('chain.json',design))
    run(binary,'pipeline','--workspace',workspace,'--request',save('chain-plan.json',chained),'--out',root/'chain')
    result=json.loads((root/'chain/sets/result.json').read_text())
    assert Q(result['tests'][0]['p_value'])==Q(1,3),result
    question='Fit the declared independent count model under its conditional assumptions'
    intent=dict(schema='gadget.pipeline.intent.v1',profile='rna_seq',assay='measurements',recipes=[dict(id='count-fit',purpose='Declared statistical comparison',goal=dict(method='count_model_analysis',parameters=parameters))],
                answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',question=question,requirements=[dict(id='fit',observable='count_model_inference',bindings=['count-fit'])]))
    intent_path=save('intent.json',intent)
    catalog=run(binary,'recipe-catalog','--workspace',workspace,'--intent',intent_path)
    assert catalog['scientific_admission']['status']=='satisfied'
    proposal=dict(schema='gadget.pipeline.recipe.proposal.v1',source_digest=catalog['source_digest'],intent_sha256=catalog['intent_sha256'],action='recipes',goals=['count-fit'])
    answer=run(binary,'recipe-proposal','--workspace',workspace,'--intent',intent_path,'--request',save('proposal.json',proposal),'--prompt',question,'--out',root/'answer')
    assert answer['status']=='admitted'
    index=answer['execution']
    receipt=run(binary,'record-view','--input',root/'answer'/index['answers'],'--member','result')
    reference=receipt['answers'][0]['artifacts'][0]['results'][0]
    selection=dict(directory=str(root/'answer'),reference=reference)
    displayed=run(consumer,'--answer',save('answer-view.json',selection))
    near(displayed['result']['features'][0]['effect'],math.log(3))
    bad=copy.deepcopy(selection);bad['reference']['native_selection']['request']['sha256']='sha256:'+'0'*64
    run(consumer,'--answer',save('answer-tamper.json',bad),fail=True)
    bad=copy.deepcopy(selection);bad['reference']['source_digest']='sha256:'+'0'*64
    run(consumer,'--answer',save('answer-source.json',bad),fail=True)
    print(json.dumps(dict(status='PASS',features=len(request['features']),source_coordinates=order,native_restoration=True)))


if __name__=='__main__':
    main(*[Path(v).resolve() for v in sys.argv[1:]])
