"""Independent finite law, paired design and weighted contraction qualifications."""
import copy
import importlib.util
import itertools
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path


def main(binary, root):
    shutil.rmtree(root, ignore_errors=True)
    spec = importlib.util.spec_from_file_location('finite_demo', Path(__file__).resolve().parents[1]/
                                                'examples/bioinformatics/finite-calibration/replay.py')
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    observations = [dict(id=f'specimen-{i}', block=f'pair-{i//2}', treated=str(i % 2)) for i in range(8)]
    pairs = [[row['id'] for row in observations[i:i+2]] for i in range(0, 8, 2)]
    fields = dict(A=['10','12','20','24','30','36','40','48'],
                  B=['7','10','15','12','20','21','30','29'],
                  flat=['10','10','20','20','30','30','40','40'],
                  partial=['NA','12','20','24','30','36','40','48'])
    query = demo.stage(root/'inputs', observations, fields, pairs)
    relations=json.loads((root/'inputs/relations.json').read_text())
    for name,values in [('rel_A',[1,2,3,4]),('rel_B',[1,2,3,4]),('rel_C',[4,3,2,1])]:
        relations['fields'].append(dict(id=name,grade=1,channel='E',quantity='declared_measurement',
            units='counts',method='Authored scalar field on the original primary relations',
            values=[dict(relation=f'pair-{i}',value=str(v)) for i,v in enumerate(values)]))
    demo.save(root/'inputs/relations.json',relations)
    imported=json.loads((root/'inputs/import.json').read_text())
    imported['assays'][0]['relations']['sha256']=demo.digest(root/'inputs/relations.json')
    demo.save(root/'inputs/import.json',imported)
    query['model']['responses'] = ['all.A','all.B']
    query['hypotheses'] = query['hypotheses'][:2]
    for state in query['law']['states']:
        state['responses'] = state['responses'][:2]
    checks = refused = 0
    def check(condition):
        nonlocal checks
        assert condition
        checks += 1
    def command(*args, fail=False):
        nonlocal refused
        p = subprocess.run([str(binary), *map(str,args)], text=True, capture_output=True, timeout=180)
        check(bool(p.returncode) == fail)
        if fail:
            refused += 1
            return None
        return json.loads(p.stdout)
    command('import-assays','--spec',root/'inputs/import.json','--out',root/'workspace')
    demo.bind_evidence(query,root/'workspace')
    def execute(value, fail=False):
        demo.save(root/'request.json',dict(schema='gadget.assays.query.v1',operation='native_query',
                                         assay='measurements',query=value))
        answer = command('assay-query','--workspace',root/'workspace/workspace.json',
                         '--request',root/'request.json',fail=fail)
        if answer:
            check(answer['result']['answers'][0]['certificate']['query']==value)
            return answer['result']['answers'][0]['certificate']['result']

    measurement_rows = (root/'inputs/measurements.tsv').read_text().splitlines()[1:]
    source = {row.split('\t')[0]:[None if v=='NA' else Q(v) for v in row.split('\t')[1:]] for row in measurement_rows}
    def verify(result, request):
        check(result['total_mass']=='1' and result['law']==request['law'])
        states=request['law']['states'];hypotheses=request['hypotheses']
        weights=[Q(state['mass']) for state in states]
        expected=[]
        for state,reading in zip(states,result['states']):
            check(reading['id']==state['id'] and Q(reading['mass'])==Q(state['mass']))
            values=[]
            for hypothesis in hypotheses:
                j=request['model']['responses'].index(hypothesis['response'])
                values_for_state=source[state['responses'][j].removeprefix('all.')]
                chosen=list(range(8)) if hypothesis['stratum']=='all' else [] if hypothesis['stratum']=='absent' else [0,1]
                differences=[values_for_state[i+1]-values_for_state[i] for i in chosen[::2]]
                denominator=len(differences)*sum((d*d for d in differences),Q(0))
                values.append(sum(differences,Q(0))**2/denominator if denominator else None)
            check([Q(v) if v is not None else None for v in reading['statistics']]==values)
            expected.append(values)
        observed=next(i for i,state in enumerate(states) if state['id']==request['law']['observed_state'])
        complete=all(v is not None for values in expected for v in values)
        check((result['status']=='completed')==complete)
        maxima=[max(values) for values in expected] if complete else None
        for h,reading in enumerate(result['hypotheses']):
            available=all(values[h] is not None for values in expected)
            check((reading['status']=='completed')==available)
            check(reading['observed']==(str(expected[observed][h]) if expected[observed][h] is not None else None))
            if not available:
                check(reading['reference_tail_mass'] is None and reading['max_statistic_tail_mass'] is None)
                continue
            tail=sum((w for w,values in zip(weights,expected) if values[h]>=expected[observed][h]),Q(0))
            ties=sum((w for w,values in zip(weights,expected) if values[h]==expected[observed][h]),Q(0))
            check(Q(reading['reference_tail_mass'])==tail and Q(reading['tie_mass'])==ties)
            check(Q(reading['strict_tail_mass'])==tail-ties and tail>=weights[observed])
            distribution={v:sum((w for w,values in zip(weights,expected) if values[h]==v),Q(0)) for v in set(values[h] for values in expected)}
            check({Q(r['value']):Q(r['mass']) for r in reading['distribution']}==distribution)
            # Exhaustively check superuniformity for every possible observation.
            pvalues=[sum((w for w,values in zip(weights,expected) if values[h]>=trial[h]),Q(0)) for trial in expected]
            for alpha in set(pvalues):
                check(sum((w for w,p in zip(weights,pvalues) if p<=alpha),Q(0))<=alpha)
            if complete:
                adjusted=sum((w for w,v in zip(weights,maxima) if v>=expected[observed][h]),Q(0))
                check(Q(reading['max_statistic_tail_mass'])==adjusted and adjusted>=tail)
            else:
                check(reading['max_statistic_tail_mass'] is None)
        if not complete:
            check(result['joint_family']['global_tail_mass'] is None)
            return
        check(Q(result['joint_family']['global_tail_mass'])==sum((w for w,v in zip(weights,maxima) if v>=maxima[observed]),Q(0)))
        # Independent coordinate covariance and sample Gram contractions agree.
        mean=[sum((w*row[h] for w,row in zip(weights,expected)),Q(0)) for h in range(len(hypotheses))]
        centered=[[row[h]-mean[h] for h in range(len(hypotheses))] for row in expected]
        covariance=[[sum((w*row[i]*row[j] for w,row in zip(weights,centered)),Q(0)) for j in range(len(hypotheses))] for i in range(len(hypotheses))]
        gram=[[sum((a*b for a,b in zip(left,right)),Q(0)) for right in centered] for left in centered]
        variance=sum((covariance[h][h] for h in range(len(hypotheses))),Q(0))
        second=sum((value**2 for row in covariance for value in row),Q(0))
        check(variance==sum((w*gram[s][s] for s,w in enumerate(weights)),Q(0)))
        check(second==sum((weights[s]*weights[t]*gram[s][t]**2 for s in range(len(weights)) for t in range(len(weights))),Q(0)))
        ensemble=result['joint_family']['ensemble']
        check(list(map(Q,ensemble['mean']))==mean)
        check([list(map(Q,row)) for row in ensemble['centered_statistics']]==centered)
        check(Q(ensemble['variance'])==variance and Q(ensemble['second_contraction'])==second)
        check(ensemble['contraction_space']==('reference_states' if len(states)<len(hypotheses) else 'hypotheses'))
        if variance:
            check(Q(ensemble['effective_dimension_2'])==variance**2/second)
            check(Q(ensemble['collision_probability'])==second/variance**2)
        else:
            check(ensemble['status']=='constant_statistics' and ensemble['effective_dimension_2'] is None)
        models=result['model_comparison']['strata'][0]['models']
        for model in models:
            check(model['execution']['gram_preparations']==1)
            check(model['execution']['applications']==len(result['expanded_model_query']['responses']))
            check(model['execution']['completed_tiles']==model['execution']['applications'])
    baseline=execute(query);verify(baseline,query)
    serial=copy.deepcopy(query);serial['workers']=1
    serial_result=execute(serial);verify(serial_result,serial)
    def science(result):
        result=copy.deepcopy(result);result['expanded_model_query'].pop('workers')
        for stratum in result['model_comparison']['strata']:
            for model in stratum.get('models',[]):model.pop('execution')
        return result
    check(science(baseline)==science(serial_result))
    weighted=copy.deepcopy(query)
    for i,state in enumerate(weighted['law']['states']):state['mass']=str(Q(i+1,136))
    weighted['law']['derivation']='Authored unequal joint finite law; not the uniform paired-swap law'
    verify(execute(weighted),weighted)
    small=copy.deepcopy(query);small['law']['states']=small['law']['states'][:2]
    for state in small['law']['states']:state['mass']='1/2'
    small['law']['derivation']='Explicit authored two-state law'
    small['model']['strata']=[dict(id='first',selection=['all.block_0'])]
    for h in query['hypotheses']:
        small['hypotheses'].append(dict(h,id=h['id']+'-first',stratum='first'))
    verify(execute(small),small)
    constant=copy.deepcopy(small)
    constant['law']['states'][1]['responses']=constant['law']['states'][0]['responses']
    verify(execute(constant),constant)
    undefined=copy.deepcopy(query)
    undefined['model']['responses'].append('all.flat')
    undefined['hypotheses'].append(dict(query['hypotheses'][0],id='flat',response='all.flat'))
    for i,state in enumerate(undefined['law']['states']):state['responses'].append('all.flat' if i==0 else f'all.ref_{i}_flat')
    verify(execute(undefined),undefined)
    in_reference=copy.deepcopy(query);in_reference['law']['states'][1]['responses'][0]='all.flat'
    verify(execute(in_reference),in_reference)
    empty=copy.deepcopy(query)
    empty['model']['strata']=[dict(id='absent',selection=['all.block_0','all.block_1'])]
    empty['hypotheses'].append(dict(query['hypotheses'][0],id='absent',stratum='absent'))
    verify(execute(empty),empty)
    relational=copy.deepcopy(query)
    relational['model'].update(responses=['all.rel_A'],terms=[dict(id='x',factors=['all.rel_B'])],controls=[],ablations=[])
    relational['law']['states']=[dict(id='forward',mass='1/2',responses=['all.rel_A']),dict(id='reverse',mass='1/2',responses=['all.rel_C'])]
    relational['law']['observed_state']='forward'
    relational['law']['derivation']='Authored two-state law on the original grade-one primary relation field'
    relational['hypotheses']=[dict(id='relation-response',stratum='all',response='all.rel_A',reduced='controls',statistic='fraction_of_reduced_residual')]
    relation_result=execute(relational)
    check(relation_result['model_comparison']['grade']==1 and relation_result['hypotheses'][0]['reference_tail_mass']=='1')
    check(relation_result['joint_family']['ensemble']['status']=='constant_statistics')
    # Refuse incomplete/ill typed declarations instead of changing the law.
    failures=[]
    def bad(mutator):
        value=copy.deepcopy(query);mutator(value);failures.append(value)
    bad(lambda q:q['law']['states'][0].update(mass='0'))
    bad(lambda q:q['law']['states'][0].update(mass='-1/16'))
    bad(lambda q:q['law']['states'][0].update(mass='1/15'))
    bad(lambda q:q['law']['states'][0].update(mass='0.0625'))
    bad(lambda q:q['law']['states'][1].update(id='state-0'))
    bad(lambda q:q['law'].update(observed_state='unknown'))
    bad(lambda q:q['law'].update(observed_state='state-1'))
    bad(lambda q:q['law']['states'][1].update(responses=['all.A']))
    bad(lambda q:q['law']['states'][1].update(responses=['all.A','all.A']))
    bad(lambda q:q['law']['states'][1]['responses'].__setitem__(0,'all.partial'))
    bad(lambda q:q['law']['states'][1]['responses'].__setitem__(0,'all.treated'))
    bad(lambda q:q['law']['states'][1]['responses'].__setitem__(0,'all.rel_A'))
    bad(lambda q:q['law']['states'][1]['responses'].__setitem__(0,'all.unknown'))
    bad(lambda q:q['law'].update(assumptions=[]))
    bad(lambda q:q['law'].update(evidence=[]))
    bad(lambda q:q['law']['evidence'][0].update(text='invented quote'))
    bad(lambda q:q['law'].update(kind='sampled_permutations'))
    bad(lambda q:q['hypotheses'][0].update(stratum='unknown'))
    bad(lambda q:q['hypotheses'][0].update(reduced='full'))
    bad(lambda q:q['hypotheses'][0].update(statistic='predicted_probability'))
    bad(lambda q:q['hypotheses'][0].update(response='all.B'))
    bad(lambda q:q['hypotheses'].append(q['hypotheses'][0]))
    bad(lambda q:q.update(alpha='0'))
    bad(lambda q:q.update(alpha='1'))
    bad(lambda q:q.update(p_value='0.001'))
    bad(lambda q:q['model'].update(operation='controlled_pairings'))
    bad(lambda q:q['model'].update(workers=1))
    for value in failures:execute(value,fail=True)
    demo.save(root/'verification.json',dict(status='PASS',checks=checks,refusals=refused,
                                            states=16,contraction_routes=['hypotheses','reference_states']))
    print(json.dumps(dict(status='PASS',checks=checks,refusals=refused)))


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
