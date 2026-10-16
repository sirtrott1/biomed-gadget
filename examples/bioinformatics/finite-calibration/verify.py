#!/usr/bin/env python3
"""Independent paired law verification and optional conventional benchmarks."""
import argparse
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

from replay import digest, save


def verify(root, benchmark=False):
    certificate = json.loads((root/'certificate.json').read_text())
    query, result = certificate['query'], certificate['result']
    construction = json.loads((root/'inputs/reference-construction.json').read_text())
    observations = construction['observations']
    ids = [r['id'] for r in observations]
    position = {name:i for i,name in enumerate(ids)}
    pairs = construction['pairs']
    source = {row.split('\t')[0]:list(map(Q,row.split('\t')[1:]))
              for row in (root/'inputs/measurements.tsv').read_text().splitlines()[1:]}
    checks = 0
    def check(condition):
        nonlocal checks
        assert condition
        checks += 1
    check(len(query['law']['states'])==2**len(pairs))
    check([s['swaps'] for s in construction['states']]==[list(s) for s in itertools.product((0,1),repeat=len(pairs))])
    check(query['law']['observed_state']==construction['observed_state'])
    masses = [Q(state['mass']) for state in query['law']['states']]
    check(all(mass==Q(1,len(masses)) for mass in masses))
    fields = list(construction['original_measurements'])
    expected, differences = [], {}
    for declaration, map_record, state in zip(query['law']['states'],construction['states'],result['states']):
        check(declaration['id']==map_record['id']==state['id'])
        statistics=[]
        for name,reference in zip(fields,declaration['responses']):
            original = list(map(Q,construction['original_measurements'][name]))
            values = source[reference.removeprefix('all.')]
            check(values==[original[position[origin]] for origin in map_record['destination_to_source']])
            delta=[values[position[b]]-values[position[a]] for a,b in pairs]
            numerator=sum(delta,Q(0))**2
            denominator=len(pairs)*sum((d*d for d in delta),Q(0))
            check(denominator>0)
            statistics.append(numerator/denominator)
            if declaration['id']==query['law']['observed_state']:
                differences[name]=[str(d) for d in delta]
        expected.append(statistics)
        check(list(map(Q,state['statistics']))==statistics)
    observed=next(i for i,s in enumerate(query['law']['states']) if s['id']==query['law']['observed_state'])
    maxima=[max(state) for state in expected]
    for h,reading in enumerate(result['hypotheses']):
        threshold=expected[observed][h]
        tail=sum((mass for mass,row in zip(masses,expected) if row[h]>=threshold),Q(0))
        adjusted=sum((mass for mass,value in zip(masses,maxima) if value>=threshold),Q(0))
        check(Q(reading['observed'])==threshold)
        check(Q(reading['reference_tail_mass'])==tail and Q(reading['max_statistic_tail_mass'])==adjusted)
        check(Q(reading['tie_mass'])==sum((mass for mass,row in zip(masses,expected) if row[h]==threshold),Q(0)))
        for candidate in {row[h] for row in expected}:
            p=sum((mass for mass,row in zip(masses,expected) if row[h]>=candidate),Q(0))
            rejection=sum((mass for mass,row in zip(masses,expected)
                           if sum((w for w,other in zip(masses,expected) if other[h]>=row[h]),Q(0))<=p),Q(0))
            check(rejection<=p)
    check(Q(result['joint_family']['global_tail_mass'])==sum((mass for mass,v in zip(masses,maxima) if v>=maxima[observed]),Q(0)))
    mean=[sum((w*row[h] for w,row in zip(masses,expected)),Q(0)) for h in range(len(fields))]
    centered=[[row[h]-mean[h] for h in range(len(fields))] for row in expected]
    covariance=[[sum((w*row[i]*row[j] for w,row in zip(masses,centered)),Q(0)) for j in range(len(fields))] for i in range(len(fields))]
    variance=sum((covariance[i][i] for i in range(len(fields))),Q(0))
    second=sum((entry**2 for row in covariance for entry in row),Q(0))
    ensemble=result['joint_family']['ensemble']
    check(list(map(Q,ensemble['mean']))==mean)
    check(Q(ensemble['variance'])==variance and Q(ensemble['second_contraction'])==second)
    check(Q(ensemble['effective_dimension_2'])==variance**2/second)
    for model in result['model_comparison']['strata'][0]['models']:
        check(model['execution']['gram_preparations']==1)
        check(model['execution']['applications']==len(result['expanded_model_query']['responses']))
        for response in model['responses']:
            check(all(Q(v)==0 for v in response['normal_equations']))
            check([Q(a)+Q(b) for a,b in zip(response['represented'],response['residual'])]==list(map(Q,response['input'])))
            if model['id']=='full':
                values=source[response['field'].removeprefix('all.')]
                expected_coefficient=sum((values[position[b]]-values[position[a]] for a,b in pairs),Q(0))/len(pairs)
                check(Q(response['particular_coefficients'][model['columns'].index('treatment')])==expected_coefficient)
    verification=dict(status='PASS',checks=checks,source_digest=certificate['source_digest'],
                      certificate_sha256=digest(root/'certificate.json'),differences=differences,
                      states=len(expected),hypotheses=len(fields),reference_law='Complete joint within-pair swaps')
    save(root/'independent.json',verification)
    if benchmark:
        conventional(root,result,fields,source,pairs,position)
    return verification


def conventional(root,result,fields,source,pairs,position):
    import numpy as np
    import scipy
    from scipy import stats
    fits=[]
    for model in result['model_comparison']['strata'][0]['models']:
        design=model['design'];X=np.zeros((design['rows'],design['columns']))
        for entry in design['entries']:X[entry['row'],entry['column']]=float(Q(entry['value']))
        Y=np.array([[float(Q(v)) for v in response['input']] for response in model['responses']]).T
        coefficients,_,rank,singular=np.linalg.lstsq(X,Y,rcond=None)
        represented=X@coefficients
        error=max((abs(represented[i,j]-float(Q(response['represented'][i])))
                   for j,response in enumerate(model['responses']) for i in range(len(represented))),default=0)
        fits.append(dict(model=model['id'],response_fields=Y.shape[1],rank=int(rank),native_rank=model['rank'],
                         max_fitted_field_difference=float(error),singular_values=singular.tolist(),rcond='None'))
    tests=[]
    for name,reading in zip(fields,result['hypotheses']):
        untreated=np.array([float(source[name][position[a]]) for a,b in pairs])
        treated=np.array([float(source[name][position[b]]) for a,b in pairs])
        def statistic(a,b,axis):return np.mean(a-b,axis=axis)
        permutation=stats.permutation_test((treated,untreated),statistic,permutation_type='samples',
                                          vectorized=True,n_resamples=np.inf,alternative='two-sided')
        parametric=stats.ttest_rel(treated,untreated,alternative='two-sided')
        native=Q(reading['reference_tail_mass'])
        if float(native)!=float(permutation.pvalue):
            raise ValueError('Conventional complete paired swap tail disagrees with exact native tail')
        tests.append(dict(field=name,scipy_complete_swap_p=float(permutation.pvalue),
                          native_tail=str(native),null_means=permutation.null_distribution.tolist(),
                          paired_t_statistic=float(parametric.statistic),paired_t_p=float(parametric.pvalue),
                          paired_t_df=float(parametric.df),
                          interpretation='Paired t benchmark uses a Student reference for independent normal differences; this is a different null model from the declared finite exchangeability law'))
    save(root/'benchmark.json',dict(numpy=np.__version__,scipy=scipy.__version__,fits=fits,tests=tests,
         references=['https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html',
                     'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html'],
         scope='Same count tokens, paired blocks and designs. Complete permutation benchmark; separate parametric model. No benchmark-fitted dispersion or expression normalization is imported.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--benchmark',action='store_true')
    a=p.parse_args();print(json.dumps(verify(a.run,a.benchmark)))
