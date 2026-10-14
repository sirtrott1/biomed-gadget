"""Source-sized admission and complete native results beyond former profiles."""
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    checks = 0
    cases = []

    def check(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(message)

    def save(name, value):
        path = root / name
        path.write_text(json.dumps(value) + '\n')
        return path

    def file(name, text):
        (root / name).write_text(text)
        return dict(path=name, sha256='sha256:' + hashlib.sha256(text.encode()).hexdigest())

    def run(*args, fail=False):
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
        if fail:
            check(result.returncode != 0, 'Invalid source declaration rejected')
            if isinstance(fail, str):
                check(fail in result.stderr, 'Native publisher refusal names the admission boundary')
            return
        if result.returncode:
            raise RuntimeError(str(args) + ': ' + result.stderr)
        return json.loads(result.stdout)

    design = dict(question='Do complete native results survive larger source declarations?',
                  population='Authored source and operator contracts',
                  selection='All source coordinates and occurrences',
                  limitations=['Computational qualification fixture'])

    def assay(name, features, observations, matrix, carrier='features', **extra):
        return dict(id=name, kind='quantitative', format='matrix_market',
                    source=file(name + '.mtx', matrix),
                    features=file(name + '-features.tsv', 'id\n' + '\n'.join(features) + '\n'),
                    observations=file(name + '-observations.tsv', 'id\n' + '\n'.join(observations) + '\n'),
                    source_namespace='980000', carrier=carrier, channel='observation',
                    quantity='authored_value', units='units', method='Retained fixture coordinates',
                    arithmetic='exact_rational', implicit_entries='unobserved', missing_tokens=[], **extra)

    def admit(name, entries, resources=None, fail=False):
        spec = dict(schema='gadget.assays.import.v1', id=name, design=design,
                    assays=entries, resources=resources or [], links=[])
        run('import-assays', '--spec', save(name + '-import.json', spec), '--out', root / name, fail=fail)
        return root / name / 'workspace.json'

    def native(workspace, name, query):
        request = dict(schema='gadget.assays.query.v1', operation='native_query', assay=name, query=query)
        answer = run('assay-query', '--workspace', workspace, '--request', save('request.json', request))
        return answer['result']['answers'][0]['certificate']['result']

    # Sparse logical shape, complete source axes and more than 10,000 fields.
    n = 10007
    features, observations = [f'g{i}' for i in range(n)], [f'c{i}' for i in range(n)]
    matrix = f'%%MatrixMarket matrix coordinate real general\n{n} {n} 3\n1 1 0\n2 4207 11\n{n} 9999 -7\n'
    large = assay('large', features, observations, matrix, carrier='observations', partition_size=n)
    ws = admit('shape', [large])
    manifest = json.loads(ws.read_text())
    check(len(manifest['assays'][0]['partitions']) == 1, 'Complete requested carrier retained')
    summaries = {}
    for gene, value, squares in [('g0', Q(0), Q(0)), ('g1', Q(11), Q(121)), (f'g{n-1}', Q(-7), Q(49))]:
        summaries[gene] = native(ws, 'large', dict(operation='field_summary', field='all.' + gene, support='reported_only'))
        result = summaries[gene]
        check(result['reported'] == 1 and result['missing'] == n-1, 'Zero and absent support remain distinct')
        check(Q(result['sum']['rational']) == value and Q(result['sum_of_squares']['rational']) == squares,
              'Every retained value has its exact source summary')
    cases.append(dict(case='sparse_shape', features=n, observations=n, logical_entries=n*n,
                      native_fields=n, supported_source_values=3))

    capabilities = run('bioinformatics')
    if capabilities['hdf5_enabled']:
        try:
            import h5py
            import numpy as np
        except ImportError:
            cases.append(dict(case='h5ad_sparse_shape',status='DEPENDENCY_UNAVAILABLE'))
        else:
            path = root/'large.h5ad'
            with h5py.File(path,'w') as image:
                image.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
                for name,ids in [('var',features),('obs',observations)]:
                    group=image.create_group(name)
                    group.attrs['encoding-type']='dataframe'
                    group.attrs['encoding-version']='0.2.0'
                    group.attrs['_index']='id'
                    group.attrs['column-order']=np.array([],dtype=h5py.string_dtype())
                    index=group.create_dataset('id',data=ids,dtype=h5py.string_dtype())
                    index.attrs['encoding-type']='string-array'
                    index.attrs['encoding-version']='0.2.0'
                x=image.create_group('X')
                x.attrs['encoding-type']='csr_matrix'
                x.attrs['encoding-version']='0.1.0'
                x.attrs['shape']=np.array([n,n],dtype=np.int64)
                pointer=np.zeros(n+1,dtype=np.int64)
                pointer[1:]=1;pointer[4207:]=2;pointer[9999:]=3
                x.create_dataset('indptr',data=pointer)
                x.create_dataset('indices',data=np.array([0,1,n-1],dtype=np.int64))
                x.create_dataset('data',data=np.array([0,11,-7],dtype=np.int64))
            h5=copy.deepcopy(large)
            h5.update(format='h5ad',arithmetic='ieee754_exact',nan_policy='reject',
                      source=dict(path=path.name,sha256='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest()))
            h5ws=admit('hdf-shape',[h5])
            for gene,expected_summary in summaries.items():
                check(native(h5ws,'large',dict(operation='field_summary',field='all.'+gene,support='reported_only'))
                      == expected_summary,'Independent H5AD and Matrix Market fields agree completely')
            cases.append(dict(case='h5ad_sparse_shape',status='PASS',logical_entries=n*n,
                              native_fields=n,stored_values=3))

    # Joined native jobs return complete mathematical answers in declaration order.
    capacity = run('resources')['worker_capacity']
    concurrency = min(5, capacity)
    queries = [dict(id=f'job{i}', query=dict(schema='gadget.assays.query.v1', operation='native_query',
                 assay='large', query=dict(operation='field_summary', field='all.' + ['g0','g1',f'g{n-1}'][i % 3],
                                          support='reported_only'))) for i in range(17)]
    batch = dict(schema='gadget.assays.query.v1', operation='query_batch', queries=queries, concurrent_jobs=concurrency)
    native_children = 'native_query' in capabilities['native_publication_methods']
    if native_children:
        answer = run('assay-query', '--workspace', ws, '--request', save('batch.json', batch))['result']
        check(answer['dispatcher'] == 'rfc::runtime::WorkGraph', 'Canonical dispatcher owns jobs')
        for i, entry in enumerate(answer['results']):
            check(entry['id'] == f'job{i}', 'Declared publication order')
            result = entry['answer']['result']['answers'][0]['certificate']['result']
            check(result == summaries[['g0','g1',f'g{n-1}'][i % 3]], 'Complete batch result equals direct result')
    else:
        run('assay-query', '--workspace', ws, '--request', save('batch.json', batch),
            fail='source query child has no native publisher')
    invalid = copy.deepcopy(batch)
    invalid['queries'][1]['id'] = 'job0'
    run('assay-query', '--workspace', ws, '--request', save('invalid-batch.json', invalid), fail=True)
    cases.append(dict(case='joined_source_queries', jobs=17, requested_concurrency=concurrency,
                      capacity=capacity, native_children_admitted=native_children))

    # Partition selection follows the actual source, including researcher order.
    partition_cells = 129
    partitioned = assay('partitioned', [f'v{i}' for i in range(partition_cells)], ['x'],
                        f'%%MatrixMarket matrix coordinate integer general\n{partition_cells} 1 0\n',
                        partition_size=1)
    partition_ws = admit('selected-partitions', [partitioned])
    selection = list(reversed(range(partition_cells)))
    selected = run('assay-query', '--workspace', partition_ws, '--request', save('partition-query.json',
                   dict(schema='gadget.assays.query.v1', operation='native_query', assay='partitioned',
                        partitions=selection,
                        query=dict(operation='field_summary', field='all.x', support='reported_only'))))
    check([answer['partition'] for answer in selected['result']['answers']] == selection,
          'Every selected partition retained in declared order')
    cases.append(dict(case='selected_partition_admission', partitions=partition_cells))

    # A primary membership relation exceeds the former whole-carrier/arity limits.
    m = 4111
    values = [Q(i % 13 - 6) for i in range(m)]
    relations = dict(relations=[dict(id='library', meaning='Declared source membership',
                     participants=[dict(cell=f'v{i}', exists=True, negative=False, magnitude='1') for i in range(m)])])
    members = assay('membership', [f'v{i}' for i in range(m)], ['x'],
                    f'%%MatrixMarket matrix coordinate real general\n{m} 1 {m}\n' +
                    ''.join(f'{i+1} 1 {x}\n' for i, x in enumerate(values)),
                    relations=file('membership-relations.json', json.dumps(relations)))
    membership_ws = admit('membership-source', [members])
    packet = membership_ws.parent/'assays/0/partitions/0'
    declaration = json.loads((packet/'study.json').read_text())
    document = (packet/'relations.json').read_bytes()
    check(json.loads(document)==relations, 'Canonical source document retains every declared value and order')
    for i, occurrence in enumerate(declaration['occurrences']):
        evidence = occurrence['evidence'][0]
        quote = document[evidence['start']:evidence['end']].decode()
        check(quote==evidence['text'] and json.loads(quote)==relations['relations'][0]['participants'][i],
              'Precise original occurrence evidence in retained order')
    response = native(membership_ws, 'membership', dict(operation='green_action', context='all', field='all.x', **{'lambda':'1'}))
    total = sum(values)
    expected = [x-total/Q(m+1) for x in values]
    response_total = sum(expected)
    check([row['entity'] for row in response] == [f'v{i}' for i in range(m)], 'Complete membership carrier order')
    for row, x, y in zip(response, values, expected):
        check(Q(row['value']['rational']) == y and y + response_total == x, 'Exact full Green equation')
    cases.append(dict(case='whole_relation_green', source_cells=m, primary_arity=m, response_coordinates=m))
    invalid = copy.deepcopy(members)
    invalid['partition_size'] = m-1
    admit('split-relation', [invalid], fail=True)

    # More than 256 source cells; inactive and cancelling occurrences remain.
    p = 257
    cancellation = dict(relations=[dict(id='loop', meaning='Retained cancelled relation', participants=[
        dict(cell='v0',exists=True,negative=True,magnitude='1'),
        dict(cell='v0',exists=True,negative=False,magnitude='1'),
        dict(cell='v0',exists=False,negative=False,magnitude='7')])])
    cells = assay('chain', [f'v{i}' for i in range(p)], ['x'],
                  f'%%MatrixMarket matrix coordinate real general\n{p} 1 0\n',
                  relations=file('chain-relations.json',json.dumps(cancellation)))
    chain_ws = admit('chain-source', [cells])
    certificate = native(chain_ws, 'chain', dict(operation='chain_certificate',context='all',grade=1,
                         coefficient_domain='rational',values=['1']))
    save('large-chain-result.json',certificate)
    check(list(map(Q,certificate['boundary'])) == [Q(0)]*p, 'All boundary coordinates cancel exactly')
    check(certificate['facts']['cycle'] and not certificate['facts']['rational_boundary'], 'Retained unfilled cycle')
    check(certificate['sectors']['harmonic']==['1'] and certificate['integral_homology']['free_rank']==1,
          'Exact harmonic field and integral free class')
    check(len(certificate['occurrences'])==3, 'Every inactive and cancelling occurrence retained')
    check(list(map(Q,[v['contribution'] for v in certificate['occurrences']]))==[Q(-1),Q(1),Q(0)],
          'Independent occurrence contributions')

    paths = [dict(id=f'p{i}', steps=[dict(index=0, coefficient='1')] * (i+1)) for i in range(33)]
    compared = native(chain_ws, 'chain', dict(operation='realized_path_comparison', context='all',
                                             grade=1, paths=paths))
    check(len(compared['paths']) == len(paths), 'Complete path family beyond former count ceiling')
    for i, path in enumerate(compared['paths']):
        check(path['steps'] == paths[i]['steps'] and list(map(Q, path['field'])) == [Q(i+1)],
              'Ordered repeated occurrences retain their complete field realization')
    gram = compared['shared_field_form']
    check(gram['rows'] == len(paths) and gram['columns'] == len(paths) and
          {(i,j):Q(x) for i,j,x in gram['entries']} ==
          {(i,j):Q((i+1)*(j+1)) for i in range(len(paths)) for j in range(len(paths))},
          'Independent complete field Gram for every declared path')
    check(compared['realization_rank'] == 1 and compared['endpoint_rank'] == 0 and
          len(compared['pairs']) == len(paths)*(len(paths)-1)//2,
          'All closed path differences remain available beyond endpoint view')
    cases.append(dict(case='complete_realized_paths', paths=len(paths), steps=sum(len(p['steps']) for p in paths)))

    # Grade is not arity: adjacent cancelling boundaries define an eight-grade tower.
    higher = []
    for grade in range(1, 9):
        lower = 'v' if grade == 1 else f'r{grade-1}'
        higher.append(dict(id=f'r{grade}', grade=grade, meaning='Authored adjacent zero boundary',
                      participants=[dict(cell=lower, exists=True, negative=sign, magnitude='1')
                                    for sign in [True,False]]))
    tower = assay('tower', ['v'], ['x'], '%%MatrixMarket matrix coordinate integer general\n1 1 0\n',
                  relations=file('tower-relations.json', json.dumps(dict(relations=higher, fields=[
                      dict(id='top', grade=8, channel='observation', quantity='authored_top_value',
                           units='units', method='Measured source operand',
                           values=[dict(relation='r8',value='2')])]))))
    tower_ws = admit('higher-source', [tower])
    top = native(tower_ws, 'tower', dict(operation='chain_certificate', context='all', grade=8,
                                       coefficient_domain='rational', values=['1']))
    check(top['boundary'] == ['0'] and top['facts']['cycle'] and not top['facts']['rational_boundary'],
          'Higher grade admits the exact authored adjacent boundary')
    check(top['integral_homology']['free_rank'] == 1 and top['sectors']['harmonic'] == ['1'],
          'Native rational and integral higher-grade class is retained')
    check(len(top['occurrences']) == 2, 'Adjacent higher-grade cancelling occurrences stay distinct')
    measured = native(tower_ws, 'tower', dict(operation='chain_certificate', context='all', grade=8,
                       coefficient_domain='rational', field='all.top'))
    check(measured['rational_coordinates'] == ['2'] and measured['sectors']['harmonic'] == ['2'],
          'Higher-grade source field and carrier descend through the complete native declaration')
    tower_source = json.loads((tower_ws.parent/'assays/0/partitions/0/study.json').read_text())
    glued = native(tower_ws, 'tower', dict(operation='local_sections', context='all',
                   construction='cellular_sheaf', metric='identity',
                   stalks=[dict(grade=c['grade'],key=c['key'],axes=['a']) for c in tower_source['cells']],
                   maps=[dict(grade=o['grade'],key=o['key'],matrix=[['1']]) for o in tower_source['occurrences']],
                   sections=[dict(id=f's{g}',grade=g,values=['1'],quantity='authored_value',units='units')
                             for g in range(9)]))
    check([(g['grade'],g['dimension']) for g in glued['cohomology']] == [(g,1) for g in range(9)],
          'Independent scalar sheaf cohomology through every authored grade')
    check(all(s['compatible'] and s['glued']==['1'] and s['residual']==['0'] for s in glued['sections']),
          'Native higher-grade restriction and section gluing preserve all supplied fields')
    cases.append(dict(case='authored_higher_grade', maximum_grade=8, retained_occurrences=16))
    for levels in [[Q(0)], list(map(Q,range(33)))]:
        selected_ws, selected_id = (chain_ws,'chain') if len(levels)==1 else (membership_ws,'membership')
        # The 33-stage case uses a small complex to isolate the stage contract.
        if len(levels)>1:
            small = assay('stages',['a','b'],['x'],'%%MatrixMarket matrix coordinate integer general\n2 1 0\n')
            selected_ws,selected_id = admit('stage-source',[small]),'stages'
        births = [dict(grade=0,index=i,level=0) for i in range(p if len(levels)==1 else 2)]
        if len(levels)==1:births.append(dict(grade=1,index=0,level=0))
        request = dict(schema='gadget.assays.query.v1',operation='relational_persistence',assay=selected_id,
                       partition=0,levels=list(map(str,levels)),births=births,
                       filtration_meaning='Every retained coordinate present at level zero',evidence='Authored complete births')
        persistent = run('assay-query','--workspace',selected_ws,'--request',save('persistence.json',request))['result']
        save(f'persistence-{len(levels)}-result.json',persistent)
        check(len(persistent['coordinates']) == (p+1 if len(levels)==1 else 2), 'Complete persistence coordinates')
        check(len(persistent['stages']) == len(levels) and persistent['rank_reconstruction_verified'],
              'Every supplied stage and complete rank invariant retained')
        check([(b['grade'],b['birth_index'],b['multiplicity'],b['survives_final_level']) for b in persistent['barcode']]
              == ([(0,0,p,True),(1,0,1,True)] if len(levels)==1 else [(0,0,2,True)]),
              'Independent complete rational interval multiplicities')
    cases.append(dict(case='chain_and_persistence', source_cells=p+1, retained_occurrences=3, levels=33))

    # Source-backed requests use the codec representation, not a fixed CLI MiB cap.
    padded = root / 'complete-request.json'
    padded.write_text(' '*(8*1024*1024) + json.dumps(request) + '\n')
    transported = run('assay-query', '--workspace', selected_ws, '--request', padded)['result']
    check(transported == persistent, 'Complete larger request yields the identical persistence certificate')
    cases.append(dict(case='complete_query_transport', bytes=padded.stat().st_size))

    # A complete retained matrix exceeds the former 64 MiB source ceiling.
    large_file = assay('bytes',['a','b'],['x'],'%%MatrixMarket matrix coordinate integer general\n%'+
                       'c'*(64*1024*1024)+'\n2 1 1\n1 1 7\n')
    byte_ws = admit('byte-source',[large_file])
    retained = byte_ws.parent/'assays/0/source-matrix'
    check(hashlib.sha256(retained.read_bytes()).hexdigest() == large_file['source']['sha256'][7:], 'Every original byte retained')
    check(Q(native(byte_ws,'bytes',dict(operation='field_summary',field='all.x',support='reported_only'))['sum']['rational'])==7,
          'Large retained source gives exact measured result')
    cases.append(dict(case='retained_source_bytes',bytes=retained.stat().st_size))

    # HTS records exceed the previous 100,000-record intake profile.
    if capabilities['htslib_enabled']:
        count = 100001
        reference = dict(assembly='authored',sequence_sha256='sha256:'+hashlib.sha256(b'>chr\nA\n').hexdigest(),
                         coordinate_convention='zero_based_half_open',contigs=[dict(id='chr',length=1)])
        sam = '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:1\n'+''.join(f'r{i}\t0\tchr\t1\t60\t1M\t*\t0\t0\tA\tI\n' for i in range(count))
        resource = dict(id='reads',format='sam',source=file('reads.sam',sam),reference=reference,method='All authored SAM records')
        hts_ws = admit('hts-source',[],[resource])
        manifest = json.loads(hts_ws.read_text())
        check(manifest['resources'][0]['inspection']['record_count']==count, 'Every alignment record decoded')
        check('records' not in manifest['resources'][0]['inspection'], 'Import catalog omits redundant record expansion')
        inspection = run('inspect-file', '--spec', save('hts-inspection.json',
                         dict(schema='gadget.bioinformatics.file.v1', format='sam',
                              source=resource['source'], reference=reference)))
        check([r['source_record_index'] for r in inspection['records']]==list(range(count)),
              'Complete original record order')
        check(hashlib.sha256((hts_ws.parent/manifest['resources'][0]['original']['path']).read_bytes()).hexdigest()==resource['source']['sha256'][7:],
              'All original SAM records retained')
        cases.append(dict(case='alignment_records',records=count))
    for name in ['features','observations','logical_matrix_entries','whole_relational_carrier_cells',
                 'transpose_feature_fields','interval_records','structural_certificate_cells','persistence_cells','persistence_levels']:
        check(name not in capabilities['limits'], 'Removed ceiling absent from capability contract: '+name)
    receipt = dict(status='PASS',independent_checks=checks,cases=cases)
    save('verification.json',receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
