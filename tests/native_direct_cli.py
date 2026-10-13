"""Native direct actions and summaries with complete original source oracles."""
import importlib.util
import json
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
import sympy as S


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    path = Path(__file__).resolve().parents[1] / 'examples/bioinformatics/relational-ensembles/replay.py'
    spec = importlib.util.spec_from_file_location('direct_source', path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    source.authored(root / 'inputs', field_values=dict(A=[6, 2, 4, 8], B=[2, 4, 8, 6],
        C=[3, 7, 5, 1], zero=[0, 0, 0, 0], missing=[None]*4, partial=[None, 0, 8, None]))
    subprocess.run([str(binary), 'import-assays', '--spec', str(root / 'inputs/import.json'),
                    '--out', str(root / 'workspace')], check=True, capture_output=True, text=True)
    queries = {operation: dict(operation=operation, context='all', field='all.object_A' if operation == 'adjoint' else 'all.A')
               for operation in ('boundary', 'adjoint', 'hodge', 'green_action')}
    queries['green_action']['lambda'] = '3/2'
    queries['boundary_top'] = dict(operation='boundary', context='all', field='all.top_B')
    queries['hodge_top'] = dict(operation='hodge', context='all', field='all.top_B')
    for name in ('A', 'zero', 'missing', 'partial'):
        queries['summary_' + name] = dict(operation='field_summary', field='all.' + name, support='reported_only')
    query_path = root / 'queries.json'
    query_path.write_text(json.dumps(queries))
    study_path = next((root / 'workspace').rglob('study.json'))
    weighted = root / 'weighted'
    shutil.copytree(study_path.parent, weighted)
    declaration = json.loads(study_path.read_text())
    declaration['contexts'][0]['metrics'] = [dict(kind='diagonal', values=v)
        for v in (['1', '2', '3', '4'], ['2', '3', '5', '7'], ['11'])]
    (weighted / 'study.json').write_text(json.dumps(declaration))
    checks = 0
    for name, study, metrics in [('identity', study_path, [S.eye(4), S.eye(4), S.eye(1)]),
                               ('weighted', weighted / 'study.json', [S.diag(1, 2, 3, 4), S.diag(2, 3, 5, 7), S.diag(11)])]:
        boundary = S.Matrix([[-1]*4, [S.Rational(1, 3)]*4, [S.Rational(1, 3)]*4, [S.Rational(1, 3)]*4])
        upper = S.Matrix([1, 1, -1, -1])
        laplacian = metrics[1].inv()*boundary.T*metrics[0]*boundary + upper*metrics[2].inv()*upper.T*metrics[1]
        observed = S.Matrix([6, 2, 4, 8])
        expected = dict(boundary=boundary*observed, adjoint=metrics[1].inv()*boundary.T*metrics[0]*S.Matrix([0, 1, 2, 3]),
                        hodge=laplacian*observed, green_action=(S.eye(4)+S.Rational(3, 2)*laplacian).inv()*observed,
                        boundary_top=upper, hodge_top=metrics[2].inv()*upper.T*metrics[1]*upper)
        for operation, reference in expected.items():
            control = root / (name + '-' + operation + '.json')
            control.write_text(json.dumps(queries[operation]))
            process = subprocess.run([str(binary), 'query', '--study', str(study), '--query', str(control)],
                                     capture_output=True, text=True, check=True)
            actual = json.loads(process.stdout)['result']
            assert [S.Rational(row['value']['rational']) for row in actual] == list(reference)
            checks += len(reference)
        for field, expected_values in [('A', [6, 2, 4, 8]), ('zero', [0]*4), ('missing', []), ('partial', [0, 8])]:
            control = root / (name + '-summary-' + field + '.json')
            control.write_text(json.dumps(queries['summary_' + field]))
            process = subprocess.run([str(binary), 'query', '--study', str(study), '--query', str(control)],
                                     capture_output=True, text=True, check=True)
            actual = json.loads(process.stdout)['result']
            assert actual['reported'] == len(expected_values) and actual['missing'] == 4 - len(expected_values)
            assert Fraction(actual['sum']['rational']) == sum(expected_values)
            assert Fraction(actual['sum_of_squares']['rational']) == sum(x*x for x in expected_values)
            assert actual['mean'] is None if not expected_values else Fraction(actual['mean']['rational']) == Fraction(sum(expected_values), len(expected_values))
            checks += 5
    for name, source_path in [('identity', root / 'workspace/workspace.rfcxs'), ('weighted', weighted / 'study.json')]:
        process = subprocess.run([str(consumer), str(source_path), str(query_path), str(root / (name + '-results')), 'direct'],
                                 capture_output=True, text=True, timeout=600)
        assert process.returncode == 0, process.stderr
        print(process.stdout.strip())
    print(json.dumps(dict(status='PASS', independent_checks=checks, reference='Exact original declared source matrices and metric forms; original raw reported counts and missing support')))


if __name__ == '__main__':
    main(*[Path(value).resolve() for value in sys.argv[1:]])
