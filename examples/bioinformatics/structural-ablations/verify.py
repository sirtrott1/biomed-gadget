#!/usr/bin/env python3
"""Independent exact reference from staged inputs, using SymPy rational algebra."""
import argparse
import csv
import json
from pathlib import Path
import sympy as S


def sparse(record):
    m = S.zeros(record['rows'], record['columns'])
    for entry in record['entries']:
        i, j, q = ((entry['row'], entry['column'], entry['value'])
                   if isinstance(entry, dict) else entry)
        m[i,j] = S.Rational(q)
    return m


def column(values): return S.Matrix(len(values),1,list(map(S.Rational,values)))
def frame(m):
    basis=m.columnspace()
    return S.Matrix.hstack(*basis) if basis else S.zeros(m.rows,0)
def kernel(m):
    basis=m.nullspace()
    return S.Matrix.hstack(*basis) if basis else S.zeros(m.cols,0)
def projector(f, metric):
    return f*(f.T*metric*f).inv()*f.T*metric if f.cols else S.zeros(metric.rows)


def verify(root, certificate_name='certificate.json', report_name='verification.json'):
    certificate=json.loads((root/certificate_name).read_text());r=certificate['result'];q=certificate['query']
    study=json.loads(next((root/'workspace').rglob('study.json')).read_text())
    rows=list(csv.reader((root/'inputs/measurements.tsv').open(),delimiter='\t'))
    declared=json.loads((root/'inputs/relations.json').read_text())
    spec=json.loads((root/'inputs/import.json').read_text())['assays'][0]
    cells=r['source']['cells'];ids=[[c['id'] for c in grade] for grade in cells];n=[len(g) for g in ids]
    index=[{name:i for i,name in enumerate(g)} for g in ids]
    measured={}
    if spec['carrier']=='observations':
        for row in rows[1:]:measured['all.'+row[0]]=column(row[1:])
        for name in q['roles']:
            if name['kind']=='field_span' and 'all.treatment_contrast' in name['fields']:
                samples={s['id']:s for s in csv.DictReader((root/'inputs/observations.tsv').open(),delimiter='\t')}
                measured['all.treatment_contrast']=column([samples[name]['treatment_contrast'] for name in ids[0]])
    else:
        for j,name in enumerate(rows[0][1:]):measured['all.'+name]=column([row[j+1] for row in rows[1:]])
    for field in declared['fields']:
        v=S.zeros(n[field['grade']],1)
        for value in field['values']:v[index[field['grade']]['relation:'+value['relation']],0]=S.Rational(value['value'])
        measured['all.'+field['id']]=v
    boundaries=[S.zeros(0,n[0])]+[S.zeros(n[k-1],n[k]) for k in range(1,len(n))]+[S.zeros(n[-1],0)]
    flags=[]
    # Retained occurrence order is explicitly authored, including zeros,
    # repeated participants, cancelling contributions and inactive flags.
    for relation in declared['relations']:
        k=relation['grade'];upper=index[k]['relation:'+relation['id']]
        for p in relation['participants']:
            lower=index[k-1][('relation:' if k>1 else '')+p['cell']]
            coefficient=S.Rational(p['magnitude'])*(-1 if p['negative'] else 1) if p['exists'] else S.Rational(0)
            boundaries[k][lower,upper]+=coefficient;flags.append((k,lower,upper,coefficient,p['exists']))
    metrics=[]
    for k,m in enumerate(study['contexts'][0]['metrics']):
        metrics.append(S.eye(n[k]) if m['kind']=='identity' else S.diag(*map(S.Rational,m['values'])))
    selected_grade=r['grade'];checks=0
    def check(value):
        nonlocal checks
        assert value
        checks+=1
    def equal(actual,expected):check(actual==expected)
    def operators(bs,ms,k):
        size=bs[k].cols
        gradient=projector(frame(ms[k].inv()*bs[k].T),ms[k]) if k else S.zeros(size)
        curl=projector(frame(bs[k+1]),ms[k])
        harmonic=S.eye(size)-gradient-curl
        down=ms[k].inv()*bs[k].T*ms[k-1]*bs[k] if k else S.zeros(size)
        up=bs[k+1]*ms[k+1].inv()*bs[k+1].T*ms[k] if k+1<len(ms) else S.zeros(size)
        return [S.eye(size),gradient,curl,harmonic],down+up
    original,lbase=operators(boundaries,metrics,selected_grade);lam=S.Rational(q['lambda'])
    gbase=(S.eye(n[selected_grade])+lam*lbase).inv()
    names=('full','gradient','curl','harmonic')
    for k,g in enumerate(r['original_grades']):equal(g['betti'],n[k]-boundaries[k].rank()-boundaries[k+1].rank())
    for field in r['original_fields']:
        x=measured[field['field']]
        for name,op in zip(names,original):
            equal(column(field[name]),op*x);equal(S.Rational(field['quadrances'][name]),(x.T*op.T*metrics[selected_grade]*op*x)[0])
        equal(column(field['green']),gbase*x)
    roles={role['id']:role for role in q['roles']};variant_refs={};field_refs={}
    for ablation in r['ablations']:
        removed=set();span=set()
        for name in ablation['remove_roles']:
            role=roles[name];removed.update(role.get('cells',[]));span.update(role.get('fields',[]))
        keep=[[i for i,name in enumerate(grade) if name not in removed] for grade in ids]
        inclusions=[S.eye(n[k])[:,keep[k]] for k in range(len(n))]
        bs=[S.zeros(0,len(keep[0]))]+[inclusions[k-1].T*boundaries[k]*inclusions[k] for k in range(1,len(n))]+[S.zeros(len(keep[-1]),0)]
        ms=[j.T*m*j for j,m in zip(inclusions,metrics)]
        ps,lnew=operators(bs,ms,selected_grade);j=inclusions[selected_grade];restriction=j.T
        span_frame=S.Matrix.hstack(*(measured[name] for name in sorted(span))) if span else S.zeros(n[selected_grade],0)
        removed_projection=projector(frame(span_frame),metrics[selected_grade]);equal(ablation['removed_span_rank'],span_frame.rank())
        green=(S.eye(j.cols)+lam*lnew).inv()
        for k,g in enumerate(ablation['grades']):
            equal([c['id'] for c in ablation['coordinates'][k]],[ids[k][i] for i in keep[k]])
            equal(sparse(g['inclusion']),inclusions[k]);equal(sparse(g['cochain_restriction']),inclusions[k].T)
            equal(sparse(g['boundary']),bs[k]);equal(sparse(g['pullback_metric']),ms[k])
            equal(g['betti'],bs[k].cols-bs[k].rank()-bs[k+1].rank())
            cycles=kernel(bs[k]);mapped=inclusions[k]*cycles
            image=S.Matrix.hstack(boundaries[k+1],mapped).rank()-boundaries[k+1].rank()
            equal(g['homology_image_rank'],image);equal(g['killed_by_inclusion'],g['betti']-image)
            equal(g['ambient_classes_outside_image'],r['original_grades'][k]['betti']-image)
            excluded=[[i for i in range(n[t]) if i not in keep[t]] for t in range(len(n))]
            relative=S.zeros(0,len(excluded[0])) if k==0 else boundaries[k].extract(excluded[k-1],excluded[k])
            upper_relative=S.zeros(len(excluded[k]),0) if k+1==len(n) else boundaries[k+1].extract(excluded[k],excluded[k+1])
            equal(sparse(g['relative_boundary']),relative);equal(g['relative_betti'],relative.cols-relative.rank()-upper_relative.rank())
            if k:
                equal(boundaries[k]*inclusions[k],inclusions[k-1]*bs[k]);equal(bs[k-1]*bs[k],S.zeros(bs[k-1].rows,bs[k].cols))
            if 'harmonic_representatives' in g:
                h=sparse(g['harmonic_representatives']);equal(h.cols,g['betti']);equal(ps[3]*h,h)
                included=sparse(g['included_representatives']);equal(included,j*h)
                images=sparse(g['ambient_harmonic_images']);equal(images,original[3]*included);equal(images.rank(),image)
                killed=sparse(g['killed_class_coefficients']);equal(images*killed,S.zeros(images.rows,killed.cols));equal(killed.rank(),g['killed_by_inclusion'])
                equal(sparse(g['killed_representatives']),included*killed)
                equal(boundaries[k+1]*sparse(g['ambient_fillings']),included*killed)
        for field in ablation['fields']:
            x=measured[field['field']];removed_field=removed_projection*x;y=x-removed_field;z=restriction*y;innovation=restriction*(y-x)
            equal(column(field['removed_span_field']),removed_field);equal(column(field['retained_span_field']),y)
            equal(column(field['restricted_input']),z);equal(column(field['discarded_source_field']),y-j*z);equal(column(field['field_innovation']),innovation)
            for name,p,old in zip(names,ps,original):
                equal(column(field['sectors'][name]),p*z);equal(S.Rational(field['quadrances'][name]),(z.T*p.T*ms[selected_grade]*p*z)[0])
                if name!='full':equal(column(field['sector_change'][name]),p*z-restriction*old*x)
            for t in field['sector_transport']:equal(column(t['value']),ps[names.index(t['to'])]*restriction*original[names.index(t['from'])]*x)
            gr=field['green'];equal(column(gr['value']),green*z);equal(column(gr['restricted_original']),restriction*gbase*x)
            equal(column(gr['structure_term']),lam*green*(restriction*lbase-lnew*restriction)*gbase*x)
            equal(column(gr['field_term']),green*innovation);equal(column(gr['change']),green*z-restriction*gbase*x)
            field_refs[ablation['id'],field['field']]={name:p*z for name,p in zip(names,ps)}
        variant_refs[ablation['id']]=(keep,inclusions,bs,ms,ps)
    if 'sheaf' in q:
        sheaf=q['sheaf'];dims=[[0]*size for size in n]
        for stalk in sheaf['stalks']:dims[stalk['grade']][next(i for i,c in enumerate(cells[stalk['grade']]) if c['key']==stalk['key'])]=len(stalk['axes'])
        offsets=[]
        for grade_dims in dims:
            off=[0]
            for d in grade_dims:off.append(off[-1]+d)
            offsets.append(off)
        ds=[S.zeros(offsets[k+1][-1],offsets[k][-1]) for k in range(len(n)-1)]+[S.zeros(0,offsets[-1][-1])]
        map_records={(m['grade'],m['key']):S.Matrix([[S.Rational(q) for q in row] for row in m['matrix']]) for m in sheaf['maps']}
        for flag,witness in zip(flags,r['source']['occurrences']):
            k,lower,upper,coef,_=flag;rho=map_records[k,witness['key']]
            for i in range(rho.rows):
                for a in range(rho.cols):ds[k-1][offsets[k][upper]+i,offsets[k-1][lower]+a]+=coef*rho[i,a]
        for k in range(len(n)-1):equal(ds[k+1]*ds[k],S.zeros(ds[k+1].rows,ds[k].cols))
        for k,native in enumerate(r['original_sheaf']['cochain_maps']):equal(S.Matrix(native['values']) if native['rows'] else S.zeros(0,native['columns']),ds[k])
        def reference_sections(record,keep):
            stalk_maps=[S.eye(offsets[k][-1])[:,[a for i in keep[k] for a in range(offsets[k][i],offsets[k][i+1])]] for k in range(len(n))]
            local_delta=[stalk_maps[k+1].T*ds[k]*stalk_maps[k] for k in range(len(n)-1)]+[S.zeros(0,stalk_maps[-1].cols)]
            result={}
            for s,section in enumerate(sheaf['sections']):
                k=section['grade']
                value=column([str(measured[name][i]) for i in range(n[k]) for name in section['fields']]) if 'fields' in section else column(section['values'])
                z=stalk_maps[k].T*value;p=projector(kernel(local_delta[k]),S.eye(z.rows));fitted=p*z
                actual=record['sections'][s];equal(column(actual['glued']),fitted);equal(column(actual['residual']),z-fitted);equal(column(actual['coboundary_defect']),local_delta[k]*z)
                equal(S.Rational(actual['residual_quadrance']),((z-fitted).T*(z-fitted))[0]);equal(actual['compatible'],local_delta[k]*z==S.zeros(local_delta[k].rows,1))
                result[k]=stalk_maps[k]*fitted
            for k,g in enumerate(record['cohomology']):equal(g['dimension'],local_delta[k].cols-local_delta[k].rank()-(local_delta[k-1].rank() if k else 0))
            return result,stalk_maps
        base_sections,_=reference_sections(r['original_sheaf'],[list(range(size)) for size in n]);section_refs={}
        for ablation in r['ablations']:
            keep=variant_refs[ablation['id']][0];record=ablation['sheaf'];values,smaps=reference_sections(record,keep);section_refs[ablation['id']]=values
            for section in record['sections']:
                k=section['grade'];observed=smaps[k].T*base_sections[k]
                equal(column(section['restricted_original_compatible']),observed);equal(column(section['projection_order_defect']),smaps[k].T*values[k]-observed)
            for k,c in enumerate(record['cochain_restriction_certificates']):equal(sparse(c['defect']),S.zeros(smaps[k+1].cols,smaps[k].rows))
        for comparison in r['comparisons']:
            left,right=comparison['left'],comparison['right'];a,b=variant_refs[left][0],variant_refs[right][0]
            equal(comparison['sheaf_overlap_compatible'],all(section_refs[left][k][offsets[k][i]:offsets[k][i+1],:]==section_refs[right][k][offsets[k][i]:offsets[k][i+1],:] for k in range(len(n)) for i in set(a[k])&set(b[k])))
            discrepancies=[]
            for k in range(len(n)):
                for i in sorted(set(a[k])&set(b[k])):
                    av=section_refs[left][k][offsets[k][i]:offsets[k][i+1],:];bv=section_refs[right][k][offsets[k][i]:offsets[k][i+1],:]
                    if av!=bv:discrepancies.append((ids[k][i],av,bv))
            equal(len(comparison['sheaf_overlap_disagreements']),len(discrepancies))
            for actual,(name,av,bv) in zip(comparison['sheaf_overlap_disagreements'],discrepancies):
                equal(actual['cell']['id'],name);equal(column(actual['left']),av);equal(column(actual['right']),bv);equal(column(actual['difference']),bv-av)
            if comparison['sheaf_overlap_compatible']:
                expected={}
                for k in range(len(n)):
                    for i in sorted(set(a[k])|set(b[k])):
                        owner=left if i in a[k] else right;expected[ids[k][i]]=section_refs[owner][k][offsets[k][i]:offsets[k][i+1],:]
                equal(len(comparison['glued_sections']),len(expected))
                for actual in comparison['glued_sections']:equal(column(actual['value']),expected[actual['cell']['id']])
    for comparison in r['comparisons']:
        a,b=variant_refs[comparison['left']][0],variant_refs[comparison['right']][0]
        for k in range(len(n)):equal([c['id'] for c in comparison['common_cells'][k]],[ids[k][i] for i in sorted(set(a[k])&set(b[k]))])
        equal(comparison['union_covers_ambient'],all(set(a[k])|set(b[k])==set(range(n[k])) for k in range(len(n))))
        common=sorted(set(a[selected_grade])&set(b[selected_grade]));m=metrics[selected_grade].extract(common,common)
        for field in comparison['fields']:
            for name in names:
                av=field_refs[comparison['left'],field['field']][name];bv=field_refs[comparison['right'],field['field']][name]
                delta=bv.extract([b[selected_grade].index(i) for i in common],[0])-av.extract([a[selected_grade].index(i) for i in common],[0])
                equal(column(field['difference_on_common_carrier'][name]),delta);equal(S.Rational(field['quadrances'][name]),(delta.T*m*delta)[0])
    report=dict(status='PASS',checks=checks,reference='SymPy exact rational algebra from staged measurements, retained occurrence declarations, role bindings, metrics and sheaf maps')
    (root/report_name).write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);args=p.parse_args()
    results=[verify(args.run.resolve())]
    if (args.run/'compatible-cover-certificate.json').exists():
        results.append(verify(args.run.resolve(),'compatible-cover-certificate.json','compatible-cover-verification.json'))
    print(json.dumps(dict(status='PASS',checks=sum(r['checks'] for r in results),references=results)))
