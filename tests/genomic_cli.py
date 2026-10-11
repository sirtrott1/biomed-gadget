"""Independent finite oracles for genomic, spatial and derived-assay contracts."""
import copy
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys


def main():
    binary, root = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    codec = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else None
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    checks = 0

    def check(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(message)

    def run(*args, fail=False):
        p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
        if fail:
            check(p.returncode != 0, f'expected rejection: {args}')
            return p.stderr
        if p.returncode:
            raise RuntimeError(f'{args}: {p.stderr}')
        return json.loads(p.stdout)

    def save(name, value):
        p = root / name
        p.write_text(json.dumps(value, indent=2) + '\n')
        return p

    def source(name, data):
        p = root / name
        p.write_bytes(data.encode() if isinstance(data, str) else data)
        return dict(path=name, sha256='sha256:' + hashlib.sha256(p.read_bytes()).hexdigest())

    reference = source('reference.fa', '>chrA\n' + 'A' * 100 + '\n')
    frame = dict(assembly='authored-reference', sequence_sha256=reference['sha256'],
                 coordinate_convention='zero_based_half_open', contigs=[dict(id='chrA', length=100)])
    resources = []

    def resource(name, fmt, data, **kwargs):
        item = dict(id=name, format=fmt, source=source(name, data), method='Authored independent contract fixture', **kwargs)
        resources.append(item)
        return item

    fastq = resource('reads.fastq', 'fastq', '@read1\nACGN\n+\n!5I~\n@read2\nSSWA\n+\n5555\n')
    resource('query.fa', 'fasta', '>q\nACGT\n>qtie\nA\n')
    resource('target.fa', 'fasta', '>r\nAGT\n>rtie\nC\n')
    gtf = ('chrA\tfixture\texon\t1\t4\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
           'chrA\tfixture\texon\t8\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
           'chrA\tfixture\texon\t3\t6\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n')
    resource('annotation.gtf', 'gtf', gtf, reference=frame)
    resource('annotation.gff3', 'gff3', '##gff-version 3\nchrA\tfixture\texon\t1\t4\t.\t+\t.\tID=exon%3A1;Parent=T1,T2\n', reference=frame)
    chip = [('chip1.bed', 'chrA\t0\t4\ta\t2\t+\nchrA\t2\t6\tb\t1\t+\n'),
            ('input1.bed', 'chrA\t1\t5\tc\t1\t+\n'),
            ('chip2.bed', 'chrA\t0\t3\td\t3\t+\nchrA\t3\t6\te\t2\t+\n'),
            ('input2.bed', 'chrA\t1\t5\tf\t2\t+\n')]
    for name, data in chip:
        resource(name, 'bed', data, reference=frame, weight_column=4)
    header = '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chrA\tLN:100\n@RG\tID:lane1\tSM:sample\n@RG\tID:lane2\tSM:sample\n'
    def samrow(q, flag, pos, cigar, seq, group='lane1', mapq=60):
        return f'{q}\t{flag}\tchrA\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*len(seq))+f'\tRG:Z:{group}\n'
    sam = header + ''.join([
        samrow('unique', 0, 1, '2M', 'AA'),
        samrow('ambiguous', 0, 3, '2M', 'AA'),
        samrow('multiple', 0, 1, '2M', 'AA'),
        samrow('multiple', 256, 20, '2M', 'AA'),
        samrow('splice', 0, 1, '2M5N3M', 'AAAAA'),
        samrow('deletion', 0, 1, '2M3D2M', 'AAAA'),
        samrow('paired', 65, 1, '2M', 'AA'),
        samrow('paired', 129, 3, '2M', 'AA'),
        samrow('unique', 0, 3, '2M', 'AA', group='lane2'),
        samrow('unknownq', 0, 1, '2M', 'AA', mapq=255),
        samrow('duplicate', 1024, 1, '2M', 'AA')])
    if codec:
        resource('alignments.sam', 'sam', sam, reference=frame)
        for fmt in ('bam', 'cram'):
            p = root / f'alignments.{fmt}'
            subprocess.run([str(codec), fmt, str(root/'alignments.sam'), str(p), str(root/'reference.fa')], check=True, capture_output=True)
            resource(p.name, fmt, p.read_bytes(), reference=frame, **({'reference_source': reference} if fmt == 'cram' else {}))
        vcf = ('##fileformat=VCFv4.3\n##contig=<ID=chrA,length=100>\n'
               '##FILTER=<ID=q10,Description="Low quality">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n'
               'chrA\t3\tv1\tA\tC,G\t.\tPASS\t.\tGT\t0/1\t1|2\t./1\t2\n'
               'chrA\t8\tv2\tA\tT\t.\tq10\t.\tGT\t0/0\t1/1\t.\t1\n'
               'chrA\t10\tv3\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/0\t0\n')
        resource('calls.vcf', 'vcf', vcf, reference=frame)
        for kind, name, fmt in [('bcf', 'calls.bcf', 'bcf'), ('vcfgz', 'calls.vcf.gz', 'vcf')]:
            subprocess.run([str(codec), kind, str(root/'calls.vcf'), str(root/name)], check=True, capture_output=True)
            resource(name, fmt, (root/name).read_bytes(), reference=frame)
    features = source('features.tsv', 'id\ng1\ng2\n')
    observations = source('observations.tsv', 'id\tx\ty\tcoordinate_frame\tcoordinate_units\nA\t0\t0\tsection1\tum\nB\t1\t0\tsection1\tum\nC\t0\t1\tsection1\tum\nD\t4\t4\tsection1\tum\n')
    matrix = source('matrix.tsv', 'id\tA\tB\tC\tD\ng1\t2\t4\t6\t8\ng2\t1\tNA\t0\t3\n')
    rels = source('relations.json', json.dumps(dict(relations=[dict(id='authored-contact', grade=1, type='authored.contact', meaning='Authored comparison between A and D', participants=[dict(cell='A', exists=True, negative=True, magnitude='1'),dict(cell='D', exists=True, negative=False, magnitude='1')])])) + '\n')
    assay = dict(id='spatial', kind='spatial_expression', format='dense_tsv', source=matrix, features=features, observations=observations,
                 source_namespace='801000', carrier='observations', quantity='read_count', units='counts', channel='E',
                 method='Authored finite counts', arithmetic='exact_rational', implicit_entries='unobserved', missing_tokens=['NA'], relations=rels)
    design = dict(question='Inspect source-supported scientific objects and native invariants', population='Authored contract fixture',
                  selection='All retained records with declared filters', limitations=['Finite verification data; no biological finding'])
    spec = dict(schema='gadget.assays.import.v1', id='suite-contracts', design=design, assays=[assay], resources=resources, links=[])
    save('import.json', spec)
    run('import-assays', '--spec', root/'import.json', '--out', root/'workspace')
    workspace = root/'workspace/workspace.json'
    requests = {}
    results = {}
    def query(name, op, **kwargs):
        request = dict(schema='gadget.assays.query.v1', operation=op, **kwargs)
        requests[name] = request
        save(name+'-request.json', request)
        response = run('assay-query', '--workspace', workspace, '--request', root/(name+'-request.json'))
        save(name+'-result.json', response)
        results[name] = response['result']
        return response['result']

    qc = query('quality', 'read_quality', resource='reads.fastq', quality_thresholds=[20,40])
    check(qc['total_bases'] == 8, 'QC bases')
    check([r['gc_lower_bases'] for r in qc['reads']] == ['2','2'], 'GC lower')
    check([r['gc_upper_bases'] for r in qc['reads']] == ['3','2'], 'GC upper')
    check(qc['reads'][0]['quality_score_sum'] == '153', 'stored quality sum')
    check([c['bases_below'] for c in qc['reads'][0]['quality_thresholds']] == ['1','2'], 'quality thresholds')
    gff = query('gff', 'annotation_features', resource='annotation.gff3', feature_type='exon', group_attribute='Parent')
    check([f['id'] for f in gff['features']] == ['T1','T2'], 'multi-parent grouping')
    check(all(f['union_intervals'] == [[0,4]] for f in gff['features']), 'GFF coordinate conversion')
    spatial = query('spatial', 'spatial_relations', assay='spatial', coordinates=['x','y'], frame='section1', units='um', metric_diagonal=['1','1'],
                    radius_squared='1', relation_type='physical.neighborhood', evidence='Retained coordinate columns and the declared section frame')
    check(len(spatial['relations']['relations']) == 2, 'spatial geometry carrier')
    check([d['squared_distance'] for d in spatial['distances']] == ['1','1'], 'spatial exact quadrance')
    agg = query('aggregate', 'aggregate_assay', assay='spatial',
                feature_groups=[dict(id='g1', members=[dict(id='g1',coefficient='1')]), dict(id='both',members=[dict(id='g1',coefficient='1'),dict(id='g2',coefficient='1')])],
                observation_groups=[dict(id='AB',members=[dict(id='A',coefficient='1'),dict(id='B',coefficient='1')]), dict(id='CDmean',members=[dict(id='C',coefficient='1/2'),dict(id='D',coefficient='1/2')])],
                map_method='Declared group sum and explicit arithmetic mean', map_evidence='Authored identity dictionary', support_policy='require_all_nonzero_inputs')
    check(agg['rows'] == [['6','7'],['NA','17/2']], 'native aggregate with missing dependencies')
    samples = [dict(id=name,resource=name,scale='1',scale_provenance='Declared identity scale') for name, _ in chip]
    contrast = query('chromatin', 'chromatin_contrast', samples=samples,
                     contrasts=[dict(id='pair1',treated='chip1.bed',control='input1.bed'),dict(id='pair2',treated='chip2.bed',control='input2.bed')],
                     strand_policy='separate', thresholds=['0','1','2'], absence_policy='measured_zero_within_union')
    perbase = {}
    for name, data in chip:
        values = [Q(0)]*100
        for line in data.splitlines():
            _, begin, end, _, weight, _ = line.split('\t')
            for base in range(int(begin),int(end)):
                values[base] += Q(weight)
        perbase[name] = values
    for c, (treated, control) in zip(contrast['contrasts'], [('chip1.bed','input1.bed'),('chip2.bed','input2.bed')]):
        expected = [a-b for a,b in zip(perbase[treated],perbase[control])]
        check(Q(c['integral']) == sum(expected), 'coverage integral oracle')
        for i, segment in enumerate(contrast['segments']):
            for base in range(segment['start'], segment['end']):
                check(Q(c['values'][i]) == expected[base], 'full coverage coordinate oracle')
        for threshold in c['threshold_components']:
            selected = {b for s in contrast['segments'] for b in range(s['start'],s['end']) if expected[b] >= Q(threshold['threshold'])}
            reported = {b for s in threshold['regions'] for b in range(s['start'],s['end'])}
            check(selected == reported and Q(threshold['length']) == len(selected), 'threshold physical components')
    alignment = query('alignment', 'exact_sequence_alignment', query_resource='query.fa', query_id='q', reference_resource='target.fa', reference_id='r',
                      scoring=dict(match='0',mismatch='1',query_gap='1',reference_gap='1',units='declared_edit_cost',evidence='Authored cost contract'), ambiguity='definite',mode='global_linear_gap')
    check(alignment['minimum_cost'] == '1' and alignment['primal_dual_verified'], 'native alignment certificate')
    check(alignment['aligned_query'] == 'ACGT' and alignment['aligned_reference'] == 'A-GT', 'native alignment path')
    tied = copy.deepcopy(requests['alignment']); tied.update(query_id='qtie',reference_id='rtie'); tied['scoring']['mismatch']='2'
    tie = query('alignment-ties','exact_sequence_alignment',**{k:v for k,v in tied.items() if k not in ['schema','operation']})
    check(tie['minimum_cost'] == '2' and len(tie['all_optimal_relations']) == 5,'three optimal edit paths retain their complete relation support')
    # Exhaustively enumerate every tiny global edit word, independently of RexGraph.
    rng = random.Random(1729)
    def paths(a,b,i=0,j=0):
        if i == len(a) and j == len(b):
            yield Q(0), []
        if i < len(a):
            for value,p in paths(a,b,i+1,j): yield Q(2,3)+value,[(i,j,'I')]+p
        if j < len(b):
            for value,p in paths(a,b,i,j+1): yield Q(3,4)+value,[(i,j,'D')]+p
        if i < len(a) and j < len(b):
            op = '=' if a[i] == b[j] else 'X'
            for value,p in paths(a,b,i+1,j+1): yield (Q(0) if op == '=' else Q(5,6))+value,[(i,j,op)]+p
    for trial in range(12):
        a = ''.join(rng.choice('AC') for _ in range(rng.randrange(1,5)))
        b = ''.join(rng.choice('AC') for _ in range(rng.randrange(1,5)))
        # Import new independently hashed sequence inputs through the same methods.
        extra = copy.deepcopy(spec)
        extra['resources'] = [dict(id='a',format='fasta',source=source(f'a-{trial}.fa','>q\n'+a+'\n'),method='Finite alignment oracle'),dict(id='b',format='fasta',source=source(f'b-{trial}.fa','>r\n'+b+'\n'),method='Finite alignment oracle')]
        extra['assays'] = []
        save('oracle-import.json',extra)
        trial_ws = root/f'oracle-{trial}'
        run('import-assays','--spec',root/'oracle-import.json','--out',trial_ws)
        request = copy.deepcopy(requests['alignment']); request.update(query_resource='a',reference_resource='b')
        request['scoring'].update(query_gap='3/4',reference_gap='2/3',mismatch='5/6')
        save('oracle-request.json',request)
        out = run('assay-query','--workspace',trial_ws/'workspace.json','--request',root/'oracle-request.json')['result']
        exhaustive = list(paths(a,b)); minimum = min(x[0] for x in exhaustive)
        check(Q(out['minimum_cost']) == minimum,'exhaustive alignment minimum')
        optimal = {step for cost,path in exhaustive if cost == minimum for step in path}
        found = {(p['query_prefix'],p['reference_prefix'],p['operation']) for p in out['all_optimal_relations']}
        check(found == optimal,'all optimal alignment relation support')
    if codec:
        policy = dict(exclude_flags=1024,minimum_mapping_quality=0,unknown_mapping_quality='exclude')
        search = dict(method='Retained authored SAM mappings and annotation overlap',completeness='Candidate universe is exactly the retained mappings; no whole-genome completeness asserted',reference_sha256=frame['sequence_sha256'])
        counts = []
        for fmt in ('sam','bam','cram'):
            evidence = query('evidence-'+fmt,'alignment_evidence',resource='alignments.'+fmt,filter=policy)
            check(len(evidence['selected_record_indices']) == 9 and len(evidence['excluded_record_indices']) == 2,'alignment filters')
            check([(j['intron_start'],j['intron_end'],j['alignment_record_count']) for j in evidence['junctions']] == [(2,7,'1')],'splice skips')
            check(evidence['records'][5]['aligned_blocks'] == [[0,2],[5,7]],'deletion excludes sequenced coverage')
            counts.append(query('counts-'+fmt,'feature_counts',resource='alignments.'+fmt,annotation='annotation.gtf',feature_type='exon',group_attribute='gene_id',unit='sequenced_segment',overlap='any_aligned_base',strand_policy='ignore',filter=policy,candidate_search=search))
        expected = {'G1':(Q(3),Q(8)),'G2':(Q(0),Q(4))}
        for out in counts:
            bounds = {b['id']:(Q(b['lower']),Q(b['upper'])) for b in out['allocation_certificate']['result']['samples'][0]['bounds']}
            check(bounds == expected,'segment assignment bounds')
            check(len(out['segments']) == 8,'RG/query/mate retained identity')
            check(out['allocation_request']['units'] == 'assigned_segment_count','segment quantity preserved')
        for fmt in ('calls.vcf','calls.bcf','calls.vcf.gz'):
            calls = query('genotypes-'+fmt,'variant_genotypes',resource=fmt,filter_policy='pass')
            check(len(calls['allele_fields']) == 2,'PASS differs from unfiltered dot')
            check([v.get('copy_number') for v in calls['allele_fields'][0]['values']] == ['1','1',None,'0'],'called C allele copies')
            check([v.get('copy_number') for v in calls['allele_fields'][1]['values']] == ['0','1',None,'1'],'haploid/partial GT support')
            check(calls['records'][0]['genotypes'][1]['phased_separators'] == [True],'phase retained')
        bad = copy.deepcopy(resources[-1]); bad.update(schema='gadget.bioinformatics.file.v1'); bad.pop('id'); bad.pop('method')
        bad['format'] = 'bam'; save('bad.json',bad); run('inspect-file','--spec',root/'bad.json',fail=True)
        bad = dict(schema='gadget.bioinformatics.file.v1',format='cram',source=resources[[r['id'] for r in resources].index('alignments.cram')]['source'],reference=frame)
        save('bad.json',bad); run('inspect-file','--spec',root/'bad.json',fail=True)

    for name in ['spatial','aggregate','chromatin','alignment','alignment-ties'] + (['counts-sam','genotypes-calls.vcf'] if codec else []):
        request = dict(schema='gadget.assays.materialize.v1', id='derived-'+name, source_namespace=str(900000+1000*len(name)),design=design,recipe=requests[name],channel='E' if name in ['aggregate','spatial','counts-sam'] else 'observation')
        save(name+'-derive.json',request)
        destination = root/('derived-'+name)
        admitted = run('derive-assays','--workspace',workspace,'--request',root/(name+'-derive.json'),'--out',destination)
        if name not in ('spatial', 'counts-sam'):
            check(len(admitted['assays'][0]['declaration']['provenance']) >= 3,'derived source lineage')
        if name == 'spatial':
            check(admitted['schema'] == 'gadget.assays.workspace.v2','native spatial publication')
            reference = run('record-view','--input',destination/'workspace.rfcxs','--assay','derived-spatial')['result']
            check(reference['neighbor_begin'] == 1 and reference['neighbor_count'] == 2,'authored topology and appended physical relations retained')
            check(reference['source_index']['catalog_sections']['cells'] == 7 and
                  reference['source_index']['catalog_sections']['occurrences'] == 6,'complete native spatial cells and occurrences')
            q = dict(schema='gadget.assays.query.v1',operation='native_query',assay='derived-spatial',query=dict(operation='hodge',context='all',field='all.g1'))
            save('derived-query.json',q)
            native = run('assay-query','--workspace',destination/'workspace.rfcxs','--request',root/'derived-query.json')
            check('answers' in native['result'],'derived spatial native field pipeline')
            values = native['result']['answers'][0]['certificate']['result']
            check([v['value']['rational'] for v in values] == ['-12','2','4','6'],'native Hodge uses authored and physical relations together')
            save('derived-spatial-hodge-result.json',native)
            save('derived-spatial-hodge-request.json',q)
            workflow = dict(schema='gadget.scientific.workflow.v1',id='derived-native-graph',
                            source_digest=admitted['assays'][0]['partitions'][0]['source_digest'],
                            limits=dict(depth=8,nodes=16,model_calls=0,parallel_workers=2), definitions=[],
                            entry=dict(id='analysis',type='graph',workers=2,steps=[
                                dict(id='hodge',after=[],node=dict(id='hodge',type='query',query=q['query'])),
                                dict(id='summary',after=[],node=dict(id='summary',type='query',query=dict(operation='field_summary',field='all.g1',support='reported_only')))]))
            save('derived-spatial-workflow.json',workflow)
            graph = run('workflow','--workspace',destination/'workspace.rfcxs','--assay','derived-spatial','--request',root/'derived-spatial-workflow.json')
            save('derived-spatial-workflow-result.json',graph)
            check(graph['status'] == 'complete','derived assay joined parallel workflow')
        if name == 'alignment-ties':
            q = dict(schema='gadget.assays.query.v1',operation='native_query',assay='derived-alignment-ties',query=dict(operation='green_response',context='all',field='all.chosen_path',**{'lambda':'1'},interior=[0,1],coordinates=dict(id='8800',cells=[str(1000000+i) for i in range(5)])))
            save('tied-green-request.json',q)
            green = run('assay-query','--workspace',destination/'workspace.json','--request',root/'tied-green-request.json')
            save('tied-green-result.json',green)
            g = green['result']['answers'][0]['certificate']['result']
            check([Q(x['value']['rational']) for x in g['field']] == [Q(-1,5),Q(-1,5),Q(3,5),Q(-1,5),Q(-1,5)],'native Green acts on the full optimal-alignment relation carrier')
            sections=dict(operation='local_sections',context='all',construction='cellular_sheaf',
                          stalks=[dict(grade=0,key=str(i),axes=['all.minimum_prefix_cost']) for i in range(1,5)]+[dict(grade=1,key=str(1000000+i),axes=['edit-coordinate']) for i in range(5)],
                          maps=[dict(grade=1,key=str(i),matrix=[['1']]) for i in range(1,11)],
                          sections=[dict(id='prefix-cost',grade=0,fields=['all.minimum_prefix_cost'],quantity='declared_minimum_edit_cost',units='declared_edit_cost')],metric='identity')
            q=dict(schema='gadget.assays.query.v1',operation='native_query',assay='derived-alignment-ties',query=sections)
            save('tied-sections-request.json',q)
            cohom=run('assay-query','--workspace',destination/'workspace.json','--request',root/'tied-sections-request.json')
            save('tied-sections-result.json',cohom)
            check([c['dimension'] for c in cohom['result']['answers'][0]['certificate']['result']['cohomology']] == [1,2],'full retained edit-grid sheaf cohomology')
        moved = root/('relocated '+name)
        shutil.copytree(destination,moved)
        reimported = root/('reimported-'+name)
        if name in ('spatial', 'counts-sam'):
            check((moved/'workspace.rfcxs').read_bytes() == (destination/'workspace.rfcxs').read_bytes(),
                  'complete portable native spatial publication')
        else:
            run('import-assays','--spec',moved/'reimport.json','--out',reimported)
            check((reimported/'workspace.json').exists(),'portable derived reimport')
        original_location = root/'workspace unavailable during generation replay'
        (root/'workspace').rename(original_location)
        try:
            replay = run('replay-derivation','--workspace',moved/('workspace.rfcxs' if name in ('spatial', 'counts-sam') else 'workspace.json'),'--assay','derived-'+name,'--out',root/('regenerated-'+name))
            check(replay['status'] == 'PASS','original generating method re-executed from relocated lineage alone')
            save(name+'-generation-replay.json',replay)
        finally:
            original_location.rename(root/'workspace')
        run('derive-assays','--workspace',workspace,'--request',root/(name+'-derive.json'),'--out',destination,fail=True)
    wrong = copy.deepcopy(requests['spatial']); wrong['frame'] = 'unrelated-section'; save('bad.json',wrong)
    run('assay-query','--workspace',workspace,'--request',root/'bad.json',fail=True)
    wrong = copy.deepcopy(requests['chromatin']); wrong['absence_policy'] = 'infer'; save('bad.json',wrong)
    run('assay-query','--workspace',workspace,'--request',root/'bad.json',fail=True)
    check(hashlib.sha256((root/'workspace/resources/0/source').read_bytes()).hexdigest() == fastq['source']['sha256'][7:],'original FASTQ retained')
    from method_answers import qualify
    checks += qualify(binary, root/'typed-answers', workspace,
                      {name: (request, results[name]) for name, request in requests.items()})
    receipt = dict(schema='gadget.genomic.contracts.v1',status='PASS',checks=checks,htslib=bool(codec),methods=list(requests))
    save('verification.json',receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
