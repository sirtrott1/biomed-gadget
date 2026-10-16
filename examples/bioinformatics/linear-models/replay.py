#!/usr/bin/env python3
"""Reproduce declared linear comparisons; NumPy is an optional benchmark only."""
import argparse
import hashlib
import itertools
import json
import shutil
import subprocess
import time
from fractions import Fraction as Q
from pathlib import Path


def digest(path):
    return 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value if isinstance(value, str) else json.dumps(value, indent=2) + '\n')
    return path


def run(binary, *args):
    start = time.perf_counter()
    process = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True,
                             timeout=240)
    if process.returncode:
        raise RuntimeError(process.stderr)
    return json.loads(process.stdout), time.perf_counter() - start


def query(terms, responses=None, **options):
    return dict(operation='linear_model_comparison', responses=responses or ['all.response'],
                terms=[dict(id=name, factors=['all.' + f for f in factors])
                       for name, factors in terms], controls=[], selection=[],
                support='joint_reported', metric='identity', intercept=True, **options)


def stage(root):
    """Authored nonclinical factorial observations with a declared rational error field."""
    root.mkdir(parents=True)
    epsilon = Q(1, 10**20)
    rows, measurements = [], []
    for i, (x, z, arm, subtype) in enumerate(itertools.product((-1, 0, 1), (-1, 1), ('A', 'B'), ('L', 'M'))):
        b, m = int(arm == 'B'), int(subtype == 'M')
        response = 5 + 2*x - 3*z + 4*b + 5*b*x + 2*m + Q(x*x) - Q(2, 3)
        rows.append([f'unit-{i}', str(x), str(z), str(x), str(Q(x) + epsilon*z),
                     arm, subtype, 'negative' if x < 0 else 'positive' if x > 0 else 'zero',
                     str(1 + (i % 3)), 'NA' if i == 0 else str(x)])
        measurements.append(str(response))
    observations = save(root/'observations.tsv', '\t'.join(['id', 'x', 'z', 'x_copy', 'near_x', 'arm', 'subtype', 'category', 'weight', 'partial_x']) + '\n' +
                        ''.join('\t'.join(row) + '\n' for row in rows))
    features = save(root/'features.tsv', 'id\nresponse\n')
    source = save(root/'response.tsv', 'id\t' + '\t'.join(row[0] for row in rows) + '\nresponse\t' + '\t'.join(measurements) + '\n')
    def ref(path):
        return dict(path=path.name, sha256=digest(path))
    covariates = [dict(id=name, column=name, kind='exact_value', quantity=name,
                       units='dimensionless', method='Authored rational observation', missing_tokens=['NA'])
                  for name in ('x', 'z', 'x_copy', 'near_x', 'weight', 'partial_x')]
    for name, column, level in [('armA', 'arm', 'A'), ('armB', 'arm', 'B'),
                                ('subtypeM', 'subtype', 'M'), ('negative', 'category', 'negative'),
                                ('zero', 'category', 'zero'), ('positive', 'category', 'positive'),
                                ('absent_category', 'category', 'absent')]:
        covariates.append(dict(id=name, column=column, kind='level_indicator', level=level,
                               quantity=f'{column} equals {level}', units='dimensionless',
                               method='Declared lexical level membership', missing_tokens=['NA']))
    relations = [dict(id=f'relation-{i}', grade=1, meaning='Authored relation carrying independently declared scalar fields',
                      participants=[dict(cell=f'unit-{2*i}', exists=True, negative=True, magnitude='1'),
                                    dict(cell=f'unit-{2*i+1}', exists=True, negative=False, magnitude='1'),
                                    dict(cell=f'unit-{2*i+2}', exists=i != 2, negative=False, magnitude='0')])
                 for i in range(3)]
    relation_fields = [dict(id=name, grade=1, channel='observation', quantity=name, units='dimensionless',
                           method='Authored relation field', values=[dict(relation=f'relation-{i}',value=str(v)) for i,v in enumerate(values)])
                       for name,values in [('rel_x',[0,1,2]),('rel_y',[3,5,7])]]
    relation_source = save(root/'relations.json', dict(relations=relations,fields=relation_fields))
    specification = dict(schema='gadget.assays.import.v1', id='declared-linear-comparisons',
        design=dict(question='Which declared field directions represent these observations?',
                    population='24 authored nonclinical observations', selection='all',
                    limitations=['Factorial construction for numerical and interpretability checks']),
        assays=[dict(id='measurements', kind='expression', format='dense_tsv', source=ref(source),
                     features=ref(features), observations=ref(observations), source_namespace='991701',
                     carrier='observations', quantity='authored_response', units='dimensionless',
                     channel='E', method='Explicit factorial equation', arithmetic='exact_rational',
                     implicit_entries='unobserved', missing_tokens=['NA'], partition_size=len(rows),
                     observation_fields=covariates, relations=ref(relation_source))], resources=[], links=[])
    save(root/'import.json', specification)
    terms = [('x', ['x']), ('z', ['z']), ('armB', ['armB']), ('x_armB', ['x', 'armB'])]
    continuous = query(terms)
    continuous.update(controls=['all.subtypeM'], workers=2,
        strata=[dict(id='armA', selection=['all.armA']), dict(id='armB', selection=['all.armB']),
                dict(id='absent', selection=['all.absent_category'])],
        ablations=[dict(id='without_interaction', remove=['x_armB']),
                   dict(id='without_arm', remove=['armB', 'x_armB']),
                   dict(id='without_z', remove=['z'])],
        coefficient_contrasts=[dict(id='slope_A', weights=[dict(column='x', value='1')]),
                               dict(id='slope_B', weights=[dict(column='x', value='1'), dict(column='x_armB', value='1')])])
    categories = query([('x', ['x']), ('z', ['z'])], ['all.negative', 'all.zero', 'all.positive'], workers=2,
                       ablations=[dict(id='without_x', remove=['x']), dict(id='without_z', remove=['z'])])
    quadratic = {**categories, 'terms': categories['terms'] + [dict(id='x_squared', factors=['all.x', 'all.x'])],
                 'ablations': [dict(id='without_quadratic', remove=['x_squared'])]}
    redundant = query([('x', ['x']), ('x_copy', ['x_copy']), ('z', ['z']), ('armA', ['armA']), ('armB', ['armB']), ('x_armB', ['x', 'armB'])])
    redundant.update(controls=['all.subtypeM'], workers=2,
        coefficient_contrasts=[dict(id='total_x', weights=[dict(column='x', value='1'), dict(column='x_copy', value='1')]),
                               dict(id='x_alone', weights=[dict(column='x', value='1')]),
                               dict(id='arm_B_minus_A', weights=[dict(column='armB', value='1'), dict(column='armA', value='-1')])])
    weighted = {**continuous, 'metric': dict(kind='diagonal_field', field='all.weight', meaning='Authored positive metric; these are not estimated measurement precisions')}
    near = query([('x', ['x']), ('near_x', ['near_x'])], ['all.z'], workers=2,
                 ablations=[dict(id='without_near_x', remove=['near_x'])])
    queries = dict(continuous=continuous, categories=categories, category_quadratic=quadratic,
                   redundant=redundant, weighted=weighted, near_collinear=near,
                   relation_fields=query([('rel_x',['rel_x'])],['all.rel_y']))
    for name, declaration in queries.items():
        save(root/f'{name}.query.json', declaration)
    return queries


def dense(record):
    out = [[Q(0)] * record['columns'] for _ in range(record['rows'])]
    for entry in record['entries']:
        out[entry['row']][entry['column']] = Q(entry['value'])
    return out


def benchmark(receipts):
    import numpy as np
    out = dict(backend='numpy.linalg.lstsq', numpy_version=np.__version__, dtype='float64',
               rcond=None, rank_policy='NumPy default: machine precision times max(matrix.shape)',
               reference='https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html',
               cases=[])
    for case, receipt in receipts.items():
        for stratum in receipt['result']['strata']:
            if stratum['status'] != 'completed':
                continue
            weights = np.array([float(Q(v)) for v in stratum['metric_diagonal']])
            for model in stratum['models']:
                exact = dense(model['design'])
                x = np.array([[float(v) for v in row] for row in exact], dtype='float64')
                x = x.reshape(model['design']['rows'], model['design']['columns'])
                wx = np.sqrt(weights)[:, None] * x
                input_errors = [abs(Q(float(v)) - v) for row in exact for v in row]
                for response in model['responses']:
                    yq = [Q(v) for v in response['input']]
                    y = np.array([float(v) for v in yq])
                    start = time.perf_counter()
                    coefficients, _, rank, singular = np.linalg.lstsq(wx, np.sqrt(weights)*y, rcond=None)
                    elapsed = time.perf_counter() - start
                    fitted = x @ coefficients
                    residual = y-fitted
                    expected = np.array([float(Q(v)) for v in response['represented']])
                    out['cases'].append(dict(case=case, stratum=stratum['id'], model=model['id'],
                        response=response['field'], exact_rank=model['rank'], float_rank=int(rank),
                        coefficient_nullity=model['coefficient_nullity'],
                        coefficients=coefficients.tolist(), represented=fitted.tolist(), residual=residual.tolist(),
                        float_residual_quadrance=float(np.dot(weights*residual, residual)),
                        exact_residual_quadrance=response['residual_quadrance'],
                        max_represented_error=float(np.max(np.abs(expected-fitted))),
                        max_design_conversion_error=str(max(input_errors, default=Q(0))),
                        max_response_conversion_error=str(max((abs(Q(float(v))-v) for v in yq), default=Q(0))),
                        singular_values=singular.tolist(), solve_seconds=elapsed,
                        comparison_basis='Fitted fields on the declared carrier; coefficients remain elements of their declared fibres'))
    return out


def replay(binary, root, with_benchmark=False):
    relocated = root.with_name(root.name + '-relocated')
    if root.exists() or relocated.exists():
        raise ValueError('Choose a new output directory; retained evidence is never overwritten')
    declarations = stage(root/'inputs')
    imported, seconds = run(binary, 'import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    save(root/'import-receipt.json', imported)
    receipts, timings = {}, dict(import_seconds=seconds)
    for name in declarations:
        receipts[name], timings[name] = run(binary, 'assay-query', '--workspace', root/'workspace/workspace.json',
                                            '--request', save(root/f'{name}.assays-query.json',
                                              dict(schema='gadget.assays.query.v1', operation='native_query',
                                                   assay='measurements', query=declarations[name])))
        # The owning assay certificate contains the source bound scientific leaf.
        receipts[name] = receipts[name]['result']['answers'][0]['certificate']
        save(root/f'{name}.certificate.json', receipts[name])
    question = 'How do the two continuous fields and the declared category terms represent the observed response within and across arms?'
    intent = dict(schema='gadget.pipeline.intent.v1', profile='custom', assay='measurements',
                  recipes=[dict(id='continuous-comparison', purpose=question,
                    goal=dict(method='linear_model_comparison', parameters={k:v for k,v in declarations['continuous'].items() if k!='operation'}))],
                  answer_contract=dict(schema='gadget.pipeline.answer.contract.v1', question=question,
                    requirements=[dict(id='joint-fields', observable='linear_model_comparison', bindings=['continuous-comparison'])]))
    save(root/'intent.json', intent)
    catalog, _ = run(binary, 'recipe-catalog', '--workspace', root/'workspace/workspace.json', '--intent', root/'intent.json')
    save(root/'recipe-catalog.json', catalog)
    proposal = dict(schema='gadget.pipeline.recipe.proposal.v1', source_digest=catalog['source_digest'],
                    intent_sha256=catalog['intent_sha256'], action='recipes', goals=['continuous-comparison'])
    answer, _ = run(binary, 'recipe-proposal', '--workspace', root/'workspace/workspace.json', '--intent', root/'intent.json',
                     '--request', save(root/'proposal.json', proposal), '--prompt', question, '--out', root/'pipeline')
    save(root/'pipeline-receipt.json', answer)
    answers_path = root/'pipeline/answers.json'
    answers = json.loads(answers_path.read_text())
    artifact = answers['answers'][0]['artifacts'][0]
    selections = [
        ('full-rank', ['strata',0,'models',1,'rank'], 'count', 'Rank of the six declared full-model columns'),
        ('retained-observations', ['strata',0,'n'], 'count', 'Number of jointly observed coordinates in the full population'),
        ('residual-quadrance', ['strata',0,'models',1,'responses',0,'residual_quadrance'], 'rational', 'Full-model residual quadrance in the declared identity metric'),
        ('slope-A', ['strata',0,'models',1,'responses',0,'coefficient_contrasts',0,'value'], 'rational', 'Identified x slope for arm A'),
        ('slope-B', ['strata',0,'models',1,'responses',0,'coefficient_contrasts',1,'value'], 'rational', 'Identified x slope for arm B'),
        ('interaction-quadrance', ['strata',0,'comparisons',1,'responses',0,'added_quadrance'], 'rational', 'Added representation quadrance from the x-by-arm interaction'),
        ('arm-block-quadrance', ['strata',0,'comparisons',2,'responses',0,'added_quadrance'], 'rational', 'Added representation quadrance from both declared arm terms'),
        ('z-quadrance', ['strata',0,'comparisons',3,'responses',0,'added_quadrance'], 'rational', 'Added representation quadrance from the declared z term')]
    facts = dict(schema='gadget.scientific.facts.v1', source_digest=catalog['source_digest'], question=question,
                 answers=dict(path='pipeline/answers.json', sha256=digest(answers_path)),
                 facts=[dict(id=id, path=path, value_type=kind, meaning=meaning,
                             **{'from':dict(kind='answer',answer='joint-fields',binding=artifact['binding'],
                                            step=artifact['step'],result=0)}) for id,path,kind,meaning in selections])
    save(root/'facts.json', facts)
    fact_catalog, _ = run(binary,'fact-catalog','--workspace',root/'workspace/workspace.json',
                          '--request',root/'facts.json','--artifacts',root)
    save(root/'facts-catalog.json',fact_catalog)
    reading_view, _ = run(binary,'fact-catalog','--workspace',root/'workspace/workspace.json',
                          '--request',root/'facts.json','--artifacts',root,'--view','model')
    save(root/'facts-model-view.json',reading_view)
    claims = dict(schema='gadget.scientific.claims.proposal.v1',action='claims',
                  source_digest=fact_catalog['source_digest'],facts_sha256=fact_catalog['facts_sha256'],
                  question_sha256=fact_catalog['question_sha256'],
                  claims=[dict(id=f['id'],fact=f['id'],comparison='equal',value=f['value']) for f in fact_catalog['facts']])
    check, _ = run(binary,'claim-proposal','--workspace',root/'workspace/workspace.json','--facts',root/'facts.json',
                   '--request',save(root/'claims.json',claims),'--artifacts',root)
    if check['status'] != 'matched' or check['certificate']['missing']:
        raise ValueError('Literal reading did not cover every selected value')
    save(root/'claims-receipt.json',check)
    save(root/'timings.json', timings)
    if with_benchmark:
        save(root/'benchmark.json', benchmark(receipts))
    root.rename(relocated)
    after, _ = run(binary,'fact-catalog','--workspace',relocated/'workspace/workspace.json',
                   '--request',relocated/'facts.json','--artifacts',relocated)
    view_after, _ = run(binary,'fact-catalog','--workspace',relocated/'workspace/workspace.json',
                        '--request',relocated/'facts.json','--artifacts',relocated,'--view','model')
    if after != fact_catalog or view_after != reading_view:
        raise ValueError('Relocation changed the admitted values or their complete owners')
    return receipts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--benchmark', action='store_true', help='Also run the optional NumPy float64 benchmark')
    args = parser.parse_args()
    result = replay(args.binary.resolve(), args.out.resolve(), args.benchmark)
    print(json.dumps(dict(status='completed', cases=list(result), selected_facts=8, relocation='PASS',
                          out=str(args.out.resolve().with_name(args.out.name + '-relocated')))))
