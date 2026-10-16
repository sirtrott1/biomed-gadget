#!/usr/bin/env python3
"""Compare reported observation bounds with SciPy's Kaplan-Meier estimate."""
import csv,json,sys
from pathlib import Path
import scipy
from scipy.stats import ecdf,CensoredData
from fractions import Fraction as Q

root=Path(sys.argv[1]);rows=list(csv.DictReader((root/'source/lung.csv').open()));certificate=json.loads((root/'certificate.json').read_text());r=certificate['result'];times=[float(Q(x['time'])) for x in r['readings']];results=[]
for g,group in enumerate(r['groups']):
    selected=rows if g==0 else [row for row in rows if int(Q(row['sex']))==g]
    uncensored=[float(Q(row['time'])) for row in selected if int(Q(row['status']))==2];right=[float(Q(row['time'])) for row in selected if int(Q(row['status']))==1]
    estimate=ecdf(CensoredData(uncensored=uncensored,right=right)).cdf
    values=estimate.evaluate(times).tolist()
    bounds=[dict(time=x['time'],lower=x['cumulative_event_fraction_lower'][g],upper=x['cumulative_event_fraction_upper'][g],kaplan_meier=value) for x,value in zip(r['readings'],values)]
    for x in bounds:
        assert float(Q(x['lower']))-1e-12<=x['kaplan_meier']<=float(Q(x['upper']))+1e-12
    results.append(dict(group=group['id'],records=len(selected),deaths=len(uncensored),censored=len(right),readings=bounds))
out=dict(comparator='scipy.stats.ecdf(CensoredData)',scipy_version=scipy.__version__,documentation='https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ecdf.html',estimand='Kaplan-Meier event distribution under its censoring and population assumptions',native_estimand='Sharp reported-observation bounds for the original finite population',groups=results)
(root/'kaplan-meier.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='PASS',groups=len(results),comparisons=len(results)*len(times))))
