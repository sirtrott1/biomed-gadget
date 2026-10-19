"""Independent likelihood, reference tail and correction checks."""
import argparse
from fractions import Fraction
import json
import math
import random
import subprocess


def run(binary, cases):
    result = subprocess.run([str(binary), '--evaluate'], input=''.join(json.dumps(q)+'\n' for q in cases),
                            text=True, capture_output=True, timeout=120)
    if result.returncode:
        raise AssertionError(result.stderr)
    output = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(output) == len(cases)
    return output


def near(a, b, tolerance=2e-8):
    assert math.isfinite(a) and abs(a-b) <= tolerance*(1+abs(b)), (a, b)


def main(binary, scipy=False):
    cases, references = [], []
    for y in (0, 1, 3, 17, 101):
        for mu in (.01, 1, 19, 200):
            for alpha in (0, .01, .25, 2):
                cases.append(dict(operation='mass', y=str(y), mu=str(mu), alpha=str(alpha)))
                if alpha:
                    reference = sum(math.log1p(alpha*k) for k in range(y)) + y*math.log(mu) - math.lgamma(y+1) - (y+1/alpha)*math.log1p(alpha*mu)
                else:
                    reference = y*math.log(mu)-mu-math.lgamma(y+1)
                references.append(reference)
    for value, reference in zip(run(binary, cases), references):
        near(value['value'], reference, 1e-10)
    checks = len(cases)
    corrections = ('none', 'bonferroni', 'benjamini_hochberg', 'benjamini_yekutieli')
    rng = random.Random(9401)
    for n in (1, 2, 7, 33):
        probabilities = [Fraction(rng.randrange(101), 100) for _ in range(n)]
        positive = [max(p, Fraction(1, 10000)) for p in probabilities]
        cases = [dict(operation='adjust', values=[str(math.log(float(p))) for p in positive], policy=policy) for policy in corrections]
        for policy, output in zip(corrections, run(binary, cases)):
            multiplier = n*sum((Fraction(1, k) for k in range(1, n+1)), Fraction()) if policy == 'benjamini_yekutieli' else n
            expected = list(positive)
            if policy == 'bonferroni':
                expected = [min(Fraction(1), p*n) for p in positive]
            elif policy != 'none':
                order = sorted(range(n), key=positive.__getitem__)
                following = Fraction(1)
                for j in range(n, 0, -1):
                    following = min(following, positive[order[j-1]]*multiplier/j)
                    expected[order[j-1]] = following
            for value, truth in zip(output['values'], expected):
                near(value, math.log(float(truth)), 1e-12)
                checks += 1
    if scipy:
        import mpmath as mp
        import numpy as np
        from scipy.optimize import minimize, minimize_scalar
        from scipy.special import gammaln
        from scipy.stats import chi2, norm
        mp.mp.dps = 90
        cases, references = [], []
        for y in (0, 1, 7, 8, 32, 1000, 1000000, 9007199254740992):
            for alpha in (0, 1e-18, 1e-10, .001, 1, 100):
                mu = y+.75
                cases.append(dict(operation='mass', y=str(y), mu=str(mu), alpha=str(alpha)))
                m, k = mp.mpf(mu), mp.mpf(y)
                if alpha:
                    r = 1/mp.mpf(alpha)
                    reference = mp.loggamma(k+r)-mp.loggamma(r)-mp.loggamma(k+1)-r*mp.log1p(m/r)+k*(mp.log(m/r)-mp.log1p(m/r))
                else:
                    reference = k*mp.log(m)-m-mp.loggamma(k+1)
                references.append(float(reference))
        for value, reference in zip(run(binary, cases), references):
            near(value['value'], reference, 2e-9)
            checks += 1
        cases = [dict(operation='tails', z=str(z), statistic=str(z*z), degrees=df) for z in (.01, .8, 2, 5, 12, 40) for df in (1, 2, 3, 8, 17)]
        for q, value in zip(cases, run(binary, cases)):
            z, statistic = float(q['z']), float(q['statistic'])
            near(value['normal'], math.log(2)+norm.logsf(z), 1e-12)
            shape = mp.mpf(q['degrees'])/2
            near(value['chi_square'], float(mp.log(mp.gammainc(shape, statistic/2, mp.inf)/mp.gamma(shape))), 2e-12)
            checks += 2
        design = np.array([[1, i >= 6, (i % 3)-1] for i in range(12)], float)
        y = np.array([1, 11, 3, 7, 2, 14, 4, 19, 6, 23, 3, 28], float)
        offset = np.log(np.linspace(.6, 1.5, 12))
        contrast = np.array([0, 1, 0], float)
        def fit(alpha):
            def objective(beta):
                eta = design@beta+offset
                mu = np.exp(eta)
                if alpha:
                    r = 1/alpha
                    ll = gammaln(y+r)-gammaln(r)-gammaln(y+1) - r*np.log1p(mu/r) + y*(eta-math.log(r)-np.log1p(mu/r))
                else:
                    ll = y*eta-mu-gammaln(y+1)
                return -sum(ll), -design.T@((y-mu)/(1+alpha*mu))
            result = minimize(objective, np.array([2., .5, .1]), jac=True, method='BFGS', options=dict(gtol=1e-8))
            assert np.linalg.norm(result.jac) < 1e-5, result
            return result
        cases = [dict(operation='fit', design=[[str(int(v)) for v in row] for row in design], y=list(map(str,y)), offset=list(map(str,offset)), contrast=list(map(str,contrast)), alpha=str(alpha)) for alpha in (0, .1, .7, 2)]
        for q, output in zip(cases, run(binary, cases)):
            assert output['converged'], output
            alpha = float(q['alpha'])
            reference = fit(alpha)
            for a, b in zip(output['coefficients'], reference.x):
                near(a, b, 2e-7)
            mu = np.exp(design@reference.x+offset)
            info = design.T@((mu/(1+alpha*mu))[:, None]*design)
            near(output['variance'], contrast@np.linalg.solve(info,contrast), 2e-7)
            near(output['log_likelihood'], -reference.fun, 2e-8)
            checks += 5
        q = cases[-1] | dict(operation='profile')
        output = run(binary, [q])[0]
        assert output['profile_converged'] and output['converged'], output
        reference = minimize_scalar(lambda log_alpha: fit(math.exp(log_alpha)).fun, bounds=(math.log(1e-5), math.log(10)), method='bounded', options=dict(xatol=2e-6))
        near(output['alpha'], math.exp(reference.x), 2e-5)
        near(output['log_likelihood'], -reference.fun, 2e-8)
        checks += 2
        from scipy.optimize import linprog
        candidates=[];truth=[]
        for _ in range(90):
            x=np.array([[rng.randrange(-2,3) for _ in range(3)] for _ in range(6)],float)
            if np.linalg.matrix_rank(x)!=3:continue
            y=np.array([rng.choice((0.,0.,1.,3.)) for _ in range(6)])
            positive=x[y>0];zero=x[y==0]
            if not len(zero):
                boundary=False
            else:
                result=linprog(np.zeros(3),A_ub=np.vstack([zero,zero.sum(axis=0)]),b_ub=np.r_[np.zeros(len(zero)),-1.],A_eq=positive if len(positive) else None,b_eq=np.zeros(len(positive)) if len(positive) else None,bounds=[(None,None)]*3,method='highs',options={'threads':1})
                assert result.status in (0,2),result
                boundary=result.status==0
            candidates.append(dict(operation='support',design=[[str(int(v)) for v in row] for row in x],y=list(map(str,y))))
            truth.append('all_zero_counts' if boundary and not any(y) else 'boundary_likelihood' if boundary else 'finite_coefficients')
        for output,reference in zip(run(binary,candidates),truth):
            assert output['status']==reference,(output,reference)
            checks+=1
    print(json.dumps(dict(status='PASS', comparisons=checks, independent_scipy=scipy)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary')
    parser.add_argument('--scipy', action='store_true')
    args = parser.parse_args()
    main(args.binary, args.scipy)
