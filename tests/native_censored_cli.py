"""Native interval families, exact joint moments and original field owners."""
import copy
import csv
import importlib.util
import json
import shutil
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

repository = Path(__file__).resolve().parents[1]
binary, consumer, output = (Path(value).resolve() for value in sys.argv[1:])
shutil.rmtree(output, ignore_errors=True)
output.mkdir(parents=True)

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

intervals = load('native_intervals', repository / 'examples/bioinformatics/longitudinal-endpoints/replay.py')
interval_oracle = load('native_interval_oracle', repository / 'examples/bioinformatics/longitudinal-endpoints/verify.py')
joint = load('native_joint', repository / 'examples/bioinformatics/temporal-execution/replay.py')
joint_oracle = load('native_joint_oracle', repository / 'examples/bioinformatics/temporal-execution/verify.py')
receipts = []
for name in ['authored', 'weighted', 'joint']:
    case = output / name
    if name == 'joint':
        query, certificate = joint.joint(binary, case)
        independent = joint_oracle.coupled(case)
    else:
        query, certificate = intervals.endpoints(binary, case)
        independent = interval_oracle.verify_endpoints(case)
    if name == 'weighted':
        path = case / 'study/study.json'
        study = json.loads(path.read_text())
        weights = list(map(str, range(1, 9)))
        study['contexts'][0]['metrics'][0] = dict(kind='diagonal', values=weights)
        intervals.save(path, study)
        certificate, _ = intervals.run(binary, 'query', '--study', path, '--query', case / 'inputs/query.json')
        intervals.save(case / 'weighted-certificate.json', certificate)
        with (case / 'inputs/observations.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        groups = {label: [i for i, row in enumerate(rows) if row['group'] == label] for label in ['A', 'B']}
        weighted_checks = 0
        for reading in certificate['result']['readings']:
            family = reading['admissible_event_family']
            offset = list(map(Q, reading['event_field_lower']))
            expected = [Q(0)] * len(rows)
            for members in groups.values():
                mean = sum(Q(weights[i]) * offset[i] for i in members) / sum(Q(weights[i]) for i in members)
                for i in members:
                    expected[i] = mean
            assert list(map(Q, family['sectors']['curl']['offset'])) == expected
            assert list(map(Q, family['sectors']['harmonic']['offset'])) == [a - b for a, b in zip(offset, expected)]
            weighted_checks += 2 * len(rows)
        independent['checks'] += weighted_checks
    queries = {'original': query}
    if name == 'authored':
        empty = copy.deepcopy(query)
        empty['selection'] = ['all.groupA']
        empty['strata'] = ['all.groupB']
        empty.pop('contrasts')
        queries['empty-stratum'] = empty
        shared = copy.deepcopy(query)
        shared['contrasts'] = [dict(id='same', left='all.groupA', right='all.groupA', meaning='The same original readout')]
        queries['same-readout'] = shared
        zero = copy.deepcopy(query)
        entities = [value['id'] for value in certificate['result']['coordinate_carrier']]
        assignments = [['2', None, '0', '4', '7', None, '21/2', '12'],
                       ['2', '5', '5', '8', '7', '10', '21/2', '12']]
        zero['joint_family'] = dict(
            meaning='Two authored complete event assignments including supported zero time and no event',
            evidence=query['evidence'],
            states=[dict(id='assignment-' + str(i), weight=weight,
                         values=[dict(entity=entity, event_time=value) for entity, value in zip(entities, values)])
                    for i, (values, weight) in enumerate(zip(assignments, ['1/3', '2/3']))])
        queries['zero-event-time'] = zero
    queries_path = intervals.save(case / 'queries.json', queries)
    process = subprocess.run([str(consumer), str(case / 'study/study.json'), str(queries_path), str(case / 'native'), 'censored'],
                             capture_output=True, text=True, timeout=900)
    assert process.returncode == 0, process.stdout + process.stderr
    receipts.append(dict(case=name, independent=independent, native=json.loads(process.stdout)))
print(json.dumps(dict(status='PASS', cases=receipts)))
