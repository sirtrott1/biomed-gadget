"""Independent finite-allocation, physical-coverage, axis/support and native-dispatch contracts."""
import copy
import hashlib
import itertools
import json
import random
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def main():
    binary, root = str(Path(sys.argv[1]).resolve()), Path(sys.argv[2]).resolve()
    shutil.rmtree(root, ignore_errors=True); root.mkdir(parents=True)
    checks = 0
    def command(*args, success=True):
        nonlocal checks
        p = subprocess.run([binary, *map(str, args)], text=True, capture_output=True, timeout=120)
        assert (p.returncode == 0) == success, (args, p.stderr)
        checks += 1
        return json.loads(p.stdout) if success else None
    def save(name, content):
        data = content.encode() if isinstance(content, str) else content
        (root/name).write_bytes(data)
        return dict(path=name, sha256='sha256:'+hashlib.sha256(data).hexdigest())
    def request(name, value, cmd, *args, success=True):
        (root/name).write_text(json.dumps(value)); return command(cmd,*args,'--request',root/name,success=success)
    feature=save('features.tsv','id\tkind\ng1\tgene\ng2\tgene\ng3\tgene\n')
    obs=save('observations.tsv','id\tspecimen\nL1\tperson-A\nR1\tperson-A\nL2\tperson-B\nR2\tperson-B\n')
    counts=[[12,4,7,3],[0,0,3,3],[1,5,9,1]]
    dense=save('counts.tsv','id\tL1\tR1\tL2\tR2\n'+'\n'.join('g'+str(i+1)+'\t'+'\t'.join(map(str,row)) for i,row in enumerate(counts))+'\n')
    entries=[(r+1,c+1,v) for r,row in enumerate(counts) for c,v in enumerate(row) if v]
    sparse=save('counts.mtx','%%MatrixMarket matrix coordinate integer general\n% retained authored counts\n3 4 '+str(len(entries))+'\n'+'\n'.join(' '.join(map(str,e)) for e in entries)+'\n')
    def assay(id,format,source,ns,carrier='features'):
        return dict(id=id,kind='expression',format=format,source=source,features=feature,observations=obs,source_namespace=str(ns),carrier=carrier,quantity='read_count',units='counts',channel='E',method='reported authored integer counts',arithmetic='exact_decimal',implicit_entries='measured_zero',missing_tokens=[])
    spec=dict(schema='gadget.assays.import.v1',id='unrelated-assay-contracts',design=dict(question='Verify explicitly authored assay and physical-source contracts',population='authored nonclinical specimens',selection='all retained records',limitations=['authored fixture']),assays=[assay('dense','dense_tsv',dense,901000),assay('sparse','matrix_market',sparse,902000),assay('cells','dense_tsv',dense,903000,'observations')],resources=[],links=[dict(id='L1',left_assay='dense',left_observation='L1',right_assay='cells',right_observation='L1',relation='same_specimen',evidence='Authored stable specimen IDs; no numerical map implied')])
    reference=dict(assembly='authored-reference',sequence_sha256='sha256:'+'0'*64,coordinate_convention='zero_based_half_open',contigs=[dict(id='chrA',length=100)])
    bed=save('coverage.bed','chrA\t0\t4\ta\t2\t+\nchrA\t2\t6\tb\t3\t+\nchrA\t2\t3\tc\t1\t+\nchrA\t2\t4\td\t5\t-\n')
    fasta=save('sequence.fa','>seq1 retained ambiguous base\nACNTA\n>seq2\nGGG\n')
    eq=save('classes.txt','3\n3\nA1\nA2\nB1\n2\t0\t1\t6\n2\t1\t2\t8\n1\t2\t3\n')
    relations=save('contacts.json',json.dumps(dict(relations=[dict(id='pair-A',meaning='authored physical contact',participants=[dict(cell='L1',exists=True,negative=True,magnitude='1'),dict(cell='R1',exists=True,negative=False,magnitude='1')]),dict(id='pair-B',meaning='authored physical contact',participants=[dict(cell='L2',exists=True,negative=True,magnitude='1'),dict(cell='R2',exists=True,negative=False,magnitude='1')])])) )
    spec['assays'][2]['relations']=relations
    spec['assays'][2]['observation_fields']=[dict(id='person-A',column='specimen',kind='level_indicator',level='person-A',quantity='declared_specimen_membership',units='dimensionless',method='authored specimen membership indicator',missing_tokens=[])]
    cycle=[dict(id='e01',meaning='authored oriented comparison',grade=1,participants=[dict(cell='L1',exists=True,negative=True,magnitude='1'),dict(cell='R1',exists=True,negative=False,magnitude='1')]),dict(id='e12',meaning='authored oriented comparison',grade=1,participants=[dict(cell='R1',exists=True,negative=True,magnitude='1'),dict(cell='L2',exists=True,negative=False,magnitude='1')]),dict(id='e20',meaning='authored oriented comparison',grade=1,participants=[dict(cell='L2',exists=True,negative=True,magnitude='1'),dict(cell='L1',exists=True,negative=False,magnitude='1')])]
    flux=[dict(id='flux',grade=1,channel='observation',quantity='authored_relation_signal',units='dimensionless',method='independently authored exact relation field',values=[dict(relation=e['id'],value='1') for e in cycle])]
    open_rel=save('open.json',json.dumps(dict(relations=cycle,fields=flux)))
    face=dict(id='closure',meaning='explicitly authored higher closure identity',grade=2,participants=[dict(cell=e['id'],exists=True,negative=False,magnitude='1') for e in cycle])
    filled_rel=save('filled.json',json.dumps(dict(relations=cycle+[face],fields=flux)))
    for id,rel,ns in [('open',open_rel,904000),('filled',filled_rel,905000)]:
        a=assay(id,'dense_tsv',dense,ns,'observations');a['relations']=rel;spec['assays'].append(a)
    spec['resources']=[dict(id='coverage',format='bed',source=bed,method='authored weighted intervals',reference=reference,weight_column=4),dict(id='sequence',format='fasta',source=fasta,method='authored IUPAC sequence'),dict(id='classes',format='salmon_eq',source=eq,method='authored unweighted Salmon export')]
    (root/'spec.json').write_text(json.dumps(spec));catalog=command('import-assays','--spec',root/'spec.json','--out',root/'workspace')
    workspace=root/'workspace/workspace.json'
    answer_cases = {}
    def query(op,success=True,**kwargs):
        declaration=dict(schema='gadget.assays.query.v1',operation=op,**kwargs)
        result=request('query.json',declaration,'assay-query','--workspace',workspace,success=success)
        if success and op in {'field_transfer','interval_coverage','sequence_motifs','allocation_bounds'} and op not in answer_cases:
            answer_cases[op]=(declaration,result)
        return result
    def native(a,q):return query('native_query',assay=a,query=q)['result']['answers'][0]['certificate']
    for name in ['dense','sparse']:
        r=native(name,dict(operation='field_summary',field='all.L1',support='reported_only'))['result']
        assert Q(r['sum']['rational'])==13 and Q(r['sum_of_squares']['rational'])==145 and r['reported']==3
    matrix_spec=json.dumps(spec['assays'][0],separators=(',',':')) # Evidence must use actual packet bytes.
    packet=root/'workspace/assays/0/partitions/0'
    evbytes=(packet/'assay.json').read_text()
    paired=dict(operation='paired_field_contrasts',pairs=[dict(id='A',left='all.L1',right='all.R1',left_scale='1',right_scale='1'),dict(id='B',left='all.L2',right='all.R2',left_scale='1',right_scale='1')],selection=[],support='joint_reported',normalization=dict(method='identity',evidence=[dict(document='assay_declaration',start=0,end=len(evbytes.encode()),text=evbytes)]),output_units='counts',workers=1)
    p=native('dense',paired)['result']
    for row,input_values in zip(p['features'],counts):
        d=[Q(input_values[0]-input_values[1]),Q(input_values[2]-input_values[3])];m=sum(d)/2
        assert list(map(Q,row['delta']))==d and Q(row['mean_difference'])==m and list(map(Q,row['residual']))==[v-m for v in d]
        checks+=4
    # Native spatial sheaf gluing retains all gene coordinates and exact departures.
    sheaf=dict(operation='local_sections',context='all',construction='cellular_sheaf',stalks=[dict(grade=0,key=str(i),axes=['all.g1']) for i in range(1,5)]+[dict(grade=1,key=str(1000000+i),axes=['contact']) for i in range(2)],maps=[dict(grade=1,key=str(i),matrix=[['1']]) for i in range(1,5)],sections=[dict(id='g1',grade=0,fields=['all.g1'],quantity='read_count',units='counts')],metric='identity')
    glued=native('cells',sheaf)['result'];section=glued['sections'][0]
    assert list(map(Q,section['glued']))==list(map(Q,[8,8,5,5])) and list(map(Q,section['residual']))==list(map(Q,[4,-4,2,-2])) and Q(section['residual_quadrance'])==40
    assert [x['dimension'] for x in glued['cohomology']]==[2,0]
    green=dict(operation='green_response',context='all',field='all.g1',**{'lambda':'1'},interior=[0,2],coordinates=dict(id='7000',cells=['1','2','3','4']))
    green_answer=native('cells',green)['result']
    assert [Q(x['value']['rational']) for x in green_answer['field']]==[Q(28,3),Q(20,3),Q(17,3),Q(13,3)]
    checks+=8
    controls=dict(operation='controlled_pairings',fields=['all.g1','all.g2'],endpoint='all.g3',controls=['all.person-A'],selection=[],support='joint_reported',metric='identity',intercept=True)
    controlled=native('cells',controls)['result'];assert controlled['n']==4 and controlled['control_rank']==2
    checks+=2
    # Multi-core native field work and independent partition results are equal.
    runq=dict(schema='gadget.assays.query.v1',operation='native_query',assay='dense',query={**paired,'workers':2},concurrent_jobs=2)
    run=request('run.json',runq,'assay-query','--workspace',workspace,'--out',root/'run')
    assert len(run['results'])==1 and run['results'][0]['workers']==2
    leaf=command('record-view','--input',root/'run'/run['results'][0]['path'],'--member','result');assert leaf['result']['features']==p['features']
    # The same B1 and relation field have different Green response under an
    # explicitly authored higher attachment. No face is inferred from overlap.
    higher_green=dict(operation='green_response',context='all',field='all.flux',**{'lambda':'1'},interior=[0],coordinates=dict(id='7100',cells=['1000000','1000001','1000002']))
    for a,expected in [('open',Q(1)),('filled',Q(1,4))]:
        response=native(a,higher_green)['result'];assert [Q(v['value']['rational']) for v in response['field']]==[expected]*3;checks+=3
    # Cross-assay transport uses declared identity maps; no equal-shape shortcut.
    mapping=dict(domain_quantity='read_count',domain_units='counts',codomain_quantity='read_count',codomain_units='counts',method='authored identity on the same feature identities',evidence='Declared fixture correspondence',entries=[dict(left_id=f'g{i}',right_id=f'g{i}',coefficient='1') for i in range(1,4)])
    transfer=query('field_transfer',left=dict(assay='dense',partition=0,field='all.L1'),right=dict(assay='sparse',partition=0,field='all.L1'),map=mapping,metric='identity')['result']
    assert Q(transfer['residual_quadrance'])==0 and transfer['map_rank']==3 and transfer['map_kernel_dimension']==0
    bad=copy.deepcopy(mapping);bad['domain_units']='TPM'
    query('field_transfer',False,left=dict(assay='dense',partition=0,field='all.L1'),right=dict(assay='sparse',partition=0,field='all.L1'),map=bad,metric='identity')
    relation_map=dict(domain_quantity='authored_relation_signal',domain_units='dimensionless',codomain_quantity='authored_relation_signal',codomain_units='dimensionless',method='Authored correspondence between the retained lower relation identities',evidence='The open and filled declarations share the same ordered grade-one relations',entries=[dict(left_id='relation:'+e['id'],right_id='relation:'+e['id'],coefficient='1') for e in cycle])
    transported=query('field_transfer',left=dict(assay='open',partition=0,field='all.flux'),right=dict(assay='filled',partition=0,field='all.flux'),map=relation_map,metric='identity')['result']
    assert transported['left_carrier_grade']==transported['right_carrier_grade']==1
    assert transported['left_ids']==transported['right_ids']==['relation:'+e['id'] for e in cycle]
    assert list(map(Q,transported['residual']))==[Q(0)]*3 and transported['map_rank']==3 and transported['map_kernel_dimension']==0
    checks+=8
    wrong=copy.deepcopy(relation_map)
    wrong['entries']=[dict(left_id=f'g{i}',right_id=f'g{i}',coefficient='1') for i in range(1,4)]
    query('field_transfer',False,left=dict(assay='open',partition=0,field='all.flux'),right=dict(assay='filled',partition=0,field='all.flux'),map=wrong,metric='identity')
    # Physical coverage is checked independently at each base on each strand.
    intervals=[(0,4,Q(2),'+'),(2,6,Q(3),'+'),(2,3,Q(1),'+'),(2,4,Q(5),'-')]
    r=query('interval_coverage',resource='coverage',thresholds=['0','3','5','6'],strand_policy='separate')['result']
    integral=Q(0)
    for atom in r['segments']:
        expected=[sum(w for a,b,w,s in intervals if a<=x<b and s==atom['strand']) for x in range(atom['start'],atom['end'])]
        assert all(v==Q(atom['coverage']) for v in expected);integral+=sum(expected);checks+=len(expected)
    assert integral==Q(r['coverage_integral'])==31
    for family in r['threshold_components']:
        threshold=Q(family['threshold']);length=sum(1 for s in '+-' for x in range(6) if any(a<=x<b and strand==s for a,b,w,strand in intervals) and sum(w for a,b,w,strand in intervals if a<=x<b and strand==s)>=threshold)
        assert Q(family['retained_length'])==length;checks+=1
    # IUPAC ambiguity is preserved rather than treating N as a reported A/C/G/T.
    motifs=query('sequence_motifs',resource='sequence',motifs=[dict(id='AC',pattern='AC'),dict(id='CT',pattern='CT')],strands='forward')['result']
    s={ (v['sequence'],v['motif']):(v['definite_occurrences'],v['possible_occurrences']) for v in motifs['summaries']}
    assert s['seq1','AC']==(1,1) and s['seq1','CT']==(0,2)
    search=dict(method='authored finite candidate list',completeness='only the listed candidate sets',reference_sha256='sha256:'+'0'*64)
    a=query('allocation_bounds',resource='classes',targets=[dict(id='G',candidates=['A1','A2']),dict(id='H',candidates=['B1'])],candidate_search=search,scale='1',scale_provenance='raw class counts')['result']['samples'][0]['bounds']
    assert [(Q(x['lower']),Q(x['upper'])) for x in a]==[(Q(6),Q(14)),(Q(3),Q(11))]
    # Independent enumeration across many varied candidate sets and targets.
    rng=random.Random(173)
    for case in range(20):
        names=['a','b','c'];classes=[dict(id=str(i),candidates=rng.sample(names,rng.randint(1,3))) for i in range(3)];class_counts=[rng.randint(0,3) for _ in classes]
        targets=[dict(id=t,candidates=[t]) for t in names]+[dict(id='ab',candidates=['a','b'])]
        alloc=dict(schema='gadget.assay.allocations.v1',operation='allocation_bounds',candidates=names,classes=classes,samples=[dict(id='sample',counts=list(map(str,class_counts)),scale='1',scale_provenance='raw')],targets=targets,contrasts=[],class_partition='disjoint_observations',candidate_search=search,evidence=['Authored independent finite allocations'],units='assigned_fragment_count')
        def compositions(n,k):
            if k==1:return [(n,)]
            return [(i,*rest) for i in range(n+1) for rest in compositions(n-i,k-1)]
        possibilities=[]
        for joint in itertools.product(*(compositions(n,len(c['candidates'])) for c,n in zip(classes,class_counts))):
            values={t:0 for t in names}
            for c,x in zip(classes,joint):
                for t,n in zip(c['candidates'],x):values[t]+=n
            possibilities.append([sum(values[t] for t in target['candidates']) for target in targets])
        answer=request('allocation.json',alloc,'allocation-bounds')['result']['samples'][0]['bounds']
        for i,bound in enumerate(answer):
            assert Q(bound['lower'])==min(x[i] for x in possibilities) and Q(bound['upper'])==max(x[i] for x in possibilities);checks+=2
    # Genomic convention conversion, FASTQ validation, coordinate and integrity failures.
    gtf=save('gene.gtf','chrA\tannotation\texon\t1\t4\t.\t+\t0\tgene_id "g1";\n')
    physical=dict(schema='gadget.bioinformatics.file.v1',format='gtf',source=gtf,reference=reference)
    (root/'physical.json').write_text(json.dumps(physical));inspection=command('inspect-file','--spec',root/'physical.json');assert inspection['records'][0]['start']==0 and inspection['records'][0]['end']==4
    fastq=save('reads.fq','@read1\nACGT\n+\n!!!!\n');physical=dict(schema='gadget.bioinformatics.file.v1',format='fastq',source=fastq)
    (root/'physical.json').write_text(json.dumps(physical));assert command('inspect-file','--spec',root/'physical.json')['record_count']==1
    invalid=save('bad.fq','@read1\nACGT\n+\n!!!\n');physical['source']=invalid;(root/'physical.json').write_text(json.dumps(physical));command('inspect-file','--spec',root/'physical.json',success=False)
    duplicated=save('duplicate.mtx','%%MatrixMarket matrix coordinate integer general\n3 4 2\n1 1 4\n1 1 8\n')
    bad=copy.deepcopy(spec);bad['assays'][1]['source']=duplicated;(root/'bad-spec.json').write_text(json.dumps(bad));command('import-assays','--spec',root/'bad-spec.json','--out',root/'failed',success=False);assert not (root/'failed').exists()
    # Explicit zeros on sparse support differ from unreported coordinates.
    partial=save('partial.mtx','%%MatrixMarket matrix coordinate integer general\n3 4 2\n1 1 0\n3 1 2\n')
    bad=copy.deepcopy(spec);bad['assays']=[assay('partial','matrix_market',partial,910000)];bad['assays'][0]['implicit_entries']='unobserved';bad['links']=[]
    (root/'partial.json').write_text(json.dumps(bad));command('import-assays','--spec',root/'partial.json','--out',root/'partial')
    q=dict(schema='gadget.assays.query.v1',operation='native_query',assay='partial',query=dict(operation='field_summary',field='all.L1',support='reported_only'))
    summary=request('query.json',q,'assay-query','--workspace',root/'partial/workspace.json')['result']['answers'][0]['certificate']['result'];assert summary['reported']==2 and summary['missing']==1
    # Replay from a relocated workspace requires no original paths.
    relocated=root/'relocated space/workspace';shutil.copytree(root/'workspace',relocated)
    command('import-assays','--spec',relocated/'reimport.json','--out',root/'replayed')
    q=dict(schema='gadget.assays.query.v1',operation='native_query',assay='dense',query=dict(operation='field_summary',field='all.L1',support='reported_only'))
    replayed=request('query.json',q,'assay-query','--workspace',root/'replayed/workspace.json')['result']['answers'][0]['certificate']['result'];assert Q(replayed['sum']['rational'])==13
    source_digest=native('dense',dict(operation='field_summary',field='all.L1',support='reported_only'))['source_digest']
    node=dict(id='summary',type='query',query=dict(operation='field_summary',field='all.L1',support='reported_only'))
    branch=dict(id='branch',type='branch',**{'from':'summary'},condition=dict(path=['reported'],value_type='count',comparison='equal',value=3),on_true=dict(id='followup',type='query',query=dict(operation='coverage',field='all.L1')),on_false=dict(id='other',type='query',query=dict(operation='inspect')))
    workflow=dict(schema='gadget.scientific.workflow.v1',id='workspace-branch',source_digest=source_digest,limits=dict(depth=8,nodes=16,model_calls=0,parallel_workers=2),definitions=[],entry=dict(id='jobs',type='graph',workers=2,steps=[dict(id='summary',after=[],node=node),dict(id='branch',after=['summary'],node=branch)]))
    w=request('workflow.json',workflow,'workflow','--workspace',workspace,'--assay','dense');assert w['status']=='complete'
    checks+=3
    # Mutation of an original source invalidates the complete workspace.
    original=root/'workspace/assays/0/source-matrix';old=original.read_bytes();original.write_bytes(old+b'changed')
    query('inspect',False);original.write_bytes(old)
    caps=command('bioinformatics')
    assert 'stored_matrix_values' not in caps['limits']
    # UMI observations are integer counts with their own source units. Their
    # channel does not determine whether disjoint count summation is valid.
    for encoding,source in [('dense_tsv',dense),('matrix_market',sparse)]:
        uspec=copy.deepcopy(spec);u=assay('umi',encoding,source,930000)
        u.update(quantity='umi_count',units='antibody_assigned_umi',channel='observation')
        uspec['assays']=[u];uspec['links']=[]
        (root/'umi.json').write_text(json.dumps(uspec));destination=root/('umi-'+encoding)
        command('import-assays','--spec',root/'umi.json','--out',destination)
        aggregation=dict(schema='gadget.assays.query.v1',operation='aggregate_assay',assay='umi',
            feature_groups=[dict(id='g1',members=[dict(id='g1',coefficient='1')])],
            observation_groups=[dict(id='combined',members=[dict(id=v,coefficient='1') for v in ['L1','L2']])],
            count_policy='disjoint_count_sum',support_policy='require_all_nonzero_inputs',
            map_method='Declared disjoint observation counts',map_evidence='Authored retained source IDs')
        a=request('umi-query.json',aggregation,'assay-query','--workspace',destination/'workspace.json')['result']
        assert a['rows']==[['19']] and a['quantity']=='umi_count' and a['units']=='antibody_assigned_umi'
        assert a['count_contract_verified']
        for token in ['1/2','-1']:
            invalid=copy.deepcopy(uspec)
            content=('id\tL1\tR1\tL2\tR2\ng1\t'+token+'\t4\t7\t3\ng2\t0\t0\t3\t3\ng3\t1\t5\t9\t1\n'
                if encoding=='dense_tsv' else '%%MatrixMarket matrix coordinate real general\n3 4 1\n1 1 '+token+'\n')
            bad=save('invalid-umi-'+encoding,content)
            invalid['assays'][0].update(source=bad,arithmetic='exact_rational')
            (root/'invalid-umi.json').write_text(json.dumps(invalid))
            command('import-assays','--spec',root/'invalid-umi.json','--out',root/('invalid-'+encoding+token.replace('/','_')),success=False)
        checks+=4
    for token in ['1/2','-1']:
        invalid=copy.deepcopy(spec);a=invalid['assays'][3]
        values=json.loads((root/'open.json').read_text());values['fields'][0].update(quantity='umi_count',units='relation_assigned_umi')
        values['fields'][0]['values'][0]['value']=token
        a['relations']=save('invalid-relation-umi.json',json.dumps(values));invalid['assays']=[a];invalid['links']=[]
        (root/'invalid-relation-import.json').write_text(json.dumps(invalid))
        command('import-assays','--spec',root/'invalid-relation-import.json','--out',root/('invalid-relation-'+token.replace('/','_')),success=False)
    if caps['hdf5_enabled']:
        import h5py
        import numpy as np
        for encoding in ['dense','csr_matrix','csc_matrix']:
            path=root/(encoding+'.h5ad')
            with h5py.File(path,'w') as f:
                f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
                for name,axis_ids in [('obs',['L1','R1','L2','R2']),('var',['g1','g2','g3'])]:
                    group=f.create_group(name);group.attrs.update({'encoding-type':'dataframe','encoding-version':'0.2.0','_index':'_index','column-order':np.array([],dtype=h5py.string_dtype())})
                    d=group.create_dataset('_index',data=np.array(axis_ids,dtype=h5py.string_dtype()));d.attrs.update({'encoding-type':'string-array','encoding-version':'0.2.0'})
                data=np.array(counts,dtype=np.int64).T
                if encoding=='dense':
                    d=f.create_dataset('X',data=data);d.attrs.update({'encoding-type':'array','encoding-version':'0.2.0'})
                else:
                    g=f.create_group('X');g.attrs.update({'encoding-type':encoding,'encoding-version':'0.1.0','shape':data.shape});v=[];ix=[];ptr=[0]
                    major=data if encoding=='csr_matrix' else data.T
                    for row in major:
                        for j,x in enumerate(row):
                            if x:v.append(x);ix.append(j)
                        ptr.append(len(v))
                    g.create_dataset('data',data=np.array(v,dtype=np.int64));g.create_dataset('indices',data=np.array(ix,dtype=np.int32));g.create_dataset('indptr',data=np.array(ptr,dtype=np.int64))
            hspec=copy.deepcopy(spec);h=assay(encoding,'h5ad',save(path.name,path.read_bytes()),920000);h.update(arithmetic='ieee754_exact',nan_policy='reject',quantity='umi_count',units='gene_assigned_umi');hspec['assays']=[h];hspec['links']=[]
            (root/'h5spec.json').write_text(json.dumps(hspec));command('import-assays','--spec',root/'h5spec.json','--out',root/('h5-'+encoding))
            q=dict(schema='gadget.assays.query.v1',operation='native_query',assay=encoding,query=dict(operation='field_summary',field='all.L1',support='reported_only'))
            answer=request('query.json',q,'assay-query','--workspace',root/('h5-'+encoding)/'workspace.json')['result']['answers'][0]['certificate']['result'];assert Q(answer['sum']['rational'])==13 and answer['reported']==3
        for dtype in [np.float32,np.float64]:
            # Check exact stored IEEE values and NaN support independently.
            path=root/'dense.h5ad'
            with h5py.File(path,'a') as f:
                del f['X'];data=np.array(counts,dtype=dtype).T;data[0,0]=dtype('0.1');data[0,1]=np.nan
                d=f.create_dataset('X',data=data);d.attrs.update({'encoding-type':'array','encoding-version':'0.2.0'})
            hs=copy.deepcopy(hspec);hs['assays'][0].update(id='float',source=save(path.name,path.read_bytes()),quantity='recorded_intensity',units='assay_intensity',method='Authored IEEE-754 intensity values with a missing observation',nan_policy='unobserved')
            (root/'float.json').write_text(json.dumps(hs));destination=root/('float-'+dtype.__name__)
            command('import-assays','--spec',root/'float.json','--out',destination)
            q=dict(schema='gadget.assays.query.v1',operation='native_query',assay='float',query=dict(operation='field_summary',field='all.L1',support='reported_only'))
            answer=request('query.json',q,'assay-query','--workspace',destination/'workspace.json')['result']['answers'][0]['certificate']['result']
            assert answer['reported']==2 and answer['missing']==1 and Q(answer['sum']['rational'])==Q(float(dtype('0.1')))+1
            checks+=3
            invalid=copy.deepcopy(hs);invalid['assays'][0]['quantity']='umi_count'
            (root/'float-umi.json').write_text(json.dumps(invalid))
            command('import-assays','--spec',root/'float-umi.json','--out',root/('fractional-umi-'+dtype.__name__),success=False)
    from method_answers import qualify
    checks += qualify(binary, root/'typed-answers', workspace, answer_cases)
    print(json.dumps(dict(status='PASS',independent_checks=checks,hdf5=caps['hdf5_enabled'])))

if __name__=='__main__': main()
