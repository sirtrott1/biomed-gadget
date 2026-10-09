#!/usr/bin/env python3
"""Independent rational checks of shared event assignments and temporal fields."""
import argparse,json
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

def load_history():
    path=Path(__file__).resolve().parents[1]/'versioned-history/verify.py';spec=importlib.util.spec_from_file_location('history_verify',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
checks=0

def check(value):
    global checks
    assert value;checks+=1

def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def vector(v):return list(map(Q,v))
def sparse(m):return {(e['row'],e['column']):Q(e['value']) for e in m['entries']}
def apply(m,x):
    out=[Q(0)]*m['rows']
    for (i,j),v in sparse(m).items():out[i]+=v*x[j]
    return out

def coupled(root):
    cert=json.loads((root/'joint.json').read_text());out=cert['result'];q=cert['query'];joint=out['joint_family'];coords=out['coordinate_carrier'];entities={c['id']:i for i,c in enumerate(coords)};n=len(coords);times=vector(q['times']);states=q['joint_family']['states'];weights=vector([s['weight'] for s in states]);check(sum(weights)==1);check(joint['declaration']==q['joint_family']);check(joint['native_constraints_checked']==1)
    # This fixture's native comparison relations explicitly retain alternating
    # group memberships. Their exact grade-zero curl projector is the mean on
    # each disjoint four-record membership relation; all other original
    # directions remain in the harmonic field.
    groups=[list(range(0,n,2)),list(range(1,n,2))]
    def sectors(x):
        curl=[Q(0)]*n
        for group in groups:
            mean=sum((x[i] for i in group),Q(0))/len(group)
            for i in group:curl[i]=mean
        return [x,[Q(0)]*n,curl,[a-b for a,b in zip(x,curl)]]
    fields={'event':[],'restricted_time':[]}
    for state,actual in zip(states,joint['states']):
        event_times=[None]*n
        for entry in state['values']:event_times[entities[entry['entity']]]=None if entry['event_time'] is None else Q(entry['event_time'])
        check(actual['id']==state['id']);check(Q(actual['weight'])==Q(state['weight']));events=[];durations=[]
        for t,reading in zip(times,actual['readings']):
            event=[Q(int(v is not None and v<=t)) for v in event_times];duration=[min(v,t) if v is not None else t for v in event_times]
            check(Q(reading['time'])==t)
            for kind,x,destination in [('event',event,events),('restricted_time',duration,durations)]:
                parts=sectors(x);destination.extend(parts)
                for name,part in zip(['full','gradient','curl','harmonic'],parts):check(vector(reading[kind+'_sectors'][name])==part)
                key='event_fraction' if kind=='event' else 'restricted_time';check(vector(reading[key])==apply(out['readout'],x));check(vector(reading[key+'_contrast'])==apply(out['contrast_readout'],x))
        fields['event'].append(events);fields['restricted_time'].append(durations)
    for kind,samples in fields.items():
        result=joint[kind+'_moments'];channels=len(samples[0]);means=[[sum((w*s[c][i] for w,s in zip(weights,samples)),Q(0)) for i in range(n)] for c in range(channels)]
        check(len(result['channels'])==channels)
        for c,mean in enumerate(means):check(vector(result['channels'][c]['mean_field'])==mean)
        actual={key:sparse(result[key]) for key in ['expected_gram','mean_gram','centered_gram']}
        for a in range(channels):
            for b in range(channels):
                second=sum((w*dot(s[a],s[b]) for w,s in zip(weights,samples)),Q(0));mean=dot(means[a],means[b]);centered=sum((w*dot([x-y for x,y in zip(s[a],means[a])],[x-y for x,y in zip(s[b],means[b])]) for w,s in zip(weights,samples)),Q(0))
                for key,value in [('expected_gram',second),('mean_gram',mean),('centered_gram',centered)]:check(actual[key].get((a,b),Q(0))==value)
                check(centered==second-mean)
        check(any(a//4!=b//4 and v!=0 for (a,b),v in actual['centered_gram'].items()))
        centered=[[[x-y for x,y in zip(state[c],means[c])] for c in range(channels)] for state in samples]
        variances=[sum((w*dot(s[c],s[c]) for w,s in zip(weights,centered)),Q(0)) for c in range(channels)]
        for c,family in enumerate(result['field_families']):
            check(Q(family['variance'])==variances[c]);m=sparse(family['centered_fields'])
            for state in range(len(samples)):
                for i in range(n):check(m.get((state,i),Q(0))==centered[state][c][i])
        positive_orthogonal=False
        for pair in result['ensemble_coupling']:
            a,b=pair['left'],pair['right'];defined=variances[a]!=0 and variances[b]!=0;check(pair['defined']==defined)
            if not defined:check(pair['value'] is None);continue
            # Full original-coordinate cross covariance, computed directly from
            # the stated law rather than using native sample-Gram contractions.
            covariance=[[sum((w*s[a][i]*s[b][j] for w,s in zip(weights,centered)),Q(0)) for j in range(n)] for i in range(n)]
            gamma=sum((x*x for row in covariance for x in row),Q(0))/(variances[a]*variances[b]);check(Q(pair['value'])==gamma);check(0<=gamma<=1)
            if {a%4,b%4}=={2,3} and gamma>0:
                check(actual['centered_gram'].get((a,b),Q(0))==0);positive_orthogonal=True
        check(positive_orthogonal)

    # All original marginal bounds remain exactly the uncoupled observation
    # certificate; the finite law neither overwrites nor narrows that source.
    baseline=json.loads((root/'certificate.json').read_text())['result'];without={k:v for k,v in out.items() if k!='joint_family'};check(without==baseline)
    return dict(status='PASS',checks=checks)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args();v=load_history();history=v.authored(a.directory/'native');joint=coupled(a.directory/'coupled');print(json.dumps(dict(status='PASS',native=history,joint=joint)))
