#!/usr/bin/env python3
"""Execute complete declared paired reference laws through GADGET's native model engine."""
import argparse
import csv
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
    process = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True)
    if process.returncode:
        raise RuntimeError(process.stderr)
    return json.loads(process.stdout), time.perf_counter() - start


def stage(root, observations, measurements, pairs, *, units='counts', origin=None):
    """Stage every joint within-pair swap. No numerical model is calculated here."""
    if root.exists():
        raise ValueError('Choose a new input directory')
    root.mkdir(parents=True)
    ids = [row['id'] for row in observations]
    positions = {name: i for i, name in enumerate(ids)}
    if len(positions) != len(ids) or sorted(i for pair in pairs for i in pair) != sorted(ids):
        raise ValueError('Pairs must partition the unique observation identities')
    if any(len(pair) != 2 for pair in pairs):
        raise ValueError('This example authors paired swaps; production admits other finite laws')
    if any(len(values) != len(ids) for values in measurements.values()):
        raise ValueError('Measurement axes must retain all original observations')
    if any(observations[positions[a]]['block'] != observations[positions[b]]['block'] or
           {observations[positions[a]]['treated'], observations[positions[b]]['treated']} != {'0', '1'}
           for a, b in pairs):
        raise ValueError('Each pair must retain one untreated and one treated observation in its block')
    fields = dict(measurements)
    states, maps = [], []
    probability = Q(1, 2**len(pairs))
    for state_number, swaps in enumerate(itertools.product((0, 1), repeat=len(pairs))):
        permutation = list(range(len(ids)))
        for swapped, (left, right) in zip(swaps, pairs):
            if swapped:
                a, b = positions[left], positions[right]
                permutation[a], permutation[b] = permutation[b], permutation[a]
        axis = []
        for name, values in measurements.items():
            ref = name if state_number == 0 else f'ref_{state_number}_{name}'
            if state_number:
                fields[ref] = [values[i] for i in permutation]
            axis.append('all.' + ref)
        identity = f'state-{state_number}'
        states.append(dict(id=identity, mass=str(probability), responses=axis))
        maps.append(dict(id=identity, swaps=list(swaps), destination_to_source=[ids[i] for i in permutation]))
    witnesses = dict(pairs=pairs, observed_state='state-0', states=maps,
                     construction='All joint within-pair swaps; equal positive mass; original field tokens copied without rounding',
                     original_measurements=measurements, observations=observations, origin=origin)
    save(root/'reference-construction.json', witnesses)
    source = save(root/'measurements.tsv', 'id\t' + '\t'.join(ids) + '\n' +
                  ''.join(name + '\t' + '\t'.join(map(str, values)) + '\n' for name, values in fields.items()))
    features = save(root/'features.tsv', 'id\n' + ''.join(name + '\n' for name in fields))
    metadata = save(root/'observations.tsv', 'id\tblock\ttreated\n' +
                    ''.join('\t'.join(row[k] for k in ('id', 'block', 'treated')) + '\n' for row in observations))
    def ref(path):
        return dict(path=path.name, sha256=digest(path))
    blocks = list(dict.fromkeys(row['block'] for row in observations))
    controls = [f'block_{i}' for i in range(len(blocks)-1)]
    covariates = [dict(id='treated', column='treated', kind='level_indicator', level='1',
                      quantity='recorded_treatment_membership', units='dimensionless',
                      method='Source condition membership', missing_tokens=['NA'])]
    for name, block in zip(controls, blocks[:-1]):
        covariates.append(dict(id=name, column='block', kind='level_indicator', level=block,
                               quantity='recorded_block_membership', units='dimensionless',
                               method='Source block membership', missing_tokens=['NA']))
    relations = [dict(id=f'pair-{i}', grade=1, meaning='Declared comparison of the recorded paired specimens',
                      participants=[dict(cell=left, exists=True, negative=True, magnitude='1'),
                                    dict(cell=right, exists=True, negative=False, magnitude='1')])
                 for i, (left, right) in enumerate(pairs)]
    relation_file = save(root/'relations.json', dict(relations=relations, fields=[]))
    spec = dict(schema='gadget.assays.import.v1', id='paired-reference-law',
                design=dict(question='How extreme is the declared treatment representation under complete joint paired exchangeability?',
                            population='Recorded paired specimens', selection='All declared pairs and preselected response fields',
                            limitations=['Exchangeability is an explicit reference assumption; pair matching does not prove it']),
                assays=[dict(id='measurements', kind='expression', format='dense_tsv', source=ref(source),
                             features=ref(features), observations=ref(metadata), source_namespace='992101',
                             carrier='observations', quantity='declared_measurement', units=units, channel='E',
                             method='Original measurements and explicitly prospective joint reference rearrangements',
                             arithmetic='exact_rational', implicit_entries='unobserved', missing_tokens=['NA'],
                             partition_size=max(len(ids), len(fields)), observation_fields=covariates,
                             relations=ref(relation_file))], resources=[], links=[])
    save(root/'import.json', spec)
    query = dict(operation='finite_reference_calibration', workers=2,
                 model=dict(operation='linear_model_comparison', responses=['all.'+name for name in measurements],
                            terms=[dict(id='treatment', factors=['all.treated'])],
                            controls=['all.'+name for name in controls], selection=[], support='joint_reported',
                            metric='identity', intercept=True,
                            ablations=[dict(id='without_treatment', remove=['treatment'])]),
                 law=dict(kind='declared_joint_finite', observed_state='state-0', states=states,
                          null_hypothesis='The joint response vectors are invariant in law to the declared within-pair condition exchanges',
                          population='These retained paired specimens, conditional on their joint response orbit',
                          assumptions=['Every declared within-pair swap is equally likely under the joint null',
                                       'The same swaps act on every response; blocks remain fixed',
                                       'The supplied complete paired measurements and preprocessing are the conditional observations'],
                          derivation='Complete Cartesian product of the two within-pair assignments; no sampled permutations',
                          evidence=[]),
                 hypotheses=[dict(id=name, stratum='all', response='all.'+name, reduced='controls',
                                  statistic='fraction_of_reduced_residual') for name in measurements], alpha='1/20')
    # Execution binds the law to the original retained assay-declaration bytes.
    save(root/'query.json', query)
    return query


def authored(root):
    observations = [dict(id=f'specimen-{i}', block=f'pair-{i//2}', treated=str(i % 2)) for i in range(8)]
    pairs = [[row['id'] for row in observations[i:i+2]] for i in range(0, len(observations), 2)]
    measurements = dict(response_A=['10', '12', '20', '24', '30', '36', '40', '48'],
                        response_B=['7', '10', '15', '12', '20', '21', '30', '29'])
    return stage(root, observations, measurements, pairs,
                 origin=dict(scope='Authored nonclinical experiment; all measurement tokens are declared'))


def expression(root, source, feature_file):
    """Retain the complete count packet and select only researcher-declared gene IDs."""
    if root.exists():
        raise ValueError('Choose a new output directory')
    root.mkdir(parents=True)
    source = source.resolve()
    retained = root/'source'
    retained.mkdir()
    receipt = json.loads((source/'source.json').read_text())
    required = list(dict.fromkeys(['source.json', 'counts.tsv', 'samples.tsv'] + list(receipt['artifacts'])))
    for name in required:
        path = (source/name).resolve()
        if not path.is_relative_to(source) or not path.is_file():
            raise ValueError('Source artifact escapes the original packet')
        if name in receipt['artifacts'] and digest(path).removeprefix('sha256:') != receipt['artifacts'][name]:
            raise ValueError('Original source packet hash changed')
        shutil.copyfile(path, retained/name)
    metadata = list(csv.DictReader((retained/'samples.tsv').open(), delimiter='\t'))
    selected = json.loads(feature_file.read_text())
    feature_ids = [item['id'] for item in selected]
    if not feature_ids or len(set(feature_ids)) != len(feature_ids):
        raise ValueError('Declare unique feature identities before execution')
    with (retained/'counts.tsv').open() as handle:
        rows = csv.reader(handle, delimiter='\t')
        header = next(rows)
        measurements = {row[0]: row[1:] for row in rows if row[0] in feature_ids}
    if set(measurements) != set(feature_ids):
        raise ValueError('Selected features must exist in the complete original count matrix')
    by_id = {row['sample']: row for row in metadata}
    observations = [dict(id=name, block=by_id[name]['cell'], treated=str(int(by_id[name]['dex'] == 'trt')))
                    for name in header[1:]]
    pairs = []
    for block in dict.fromkeys(row['block'] for row in observations):
        group = [row for row in observations if row['block'] == block]
        group.sort(key=lambda row: row['treated'])
        pairs.append([row['id'] for row in group])
    origin = dict(provider='Bioconductor airway; GSE52778', representation='Original assigned gene-count tokens',
                  source_files={name: digest(retained/name) for name in required},
                  feature_selection=selected, normalization='none',
                  interpretation='Conditional exchangeability analysis; documented pairing alone does not establish random assignment')
    save(root/'selected-features.json', selected)
    query = stage(root/'inputs', observations, {name: measurements[name] for name in feature_ids}, pairs,
                  origin=origin)
    return query


def bind_evidence(query, workspace):
    paths = list(workspace.rglob('assay.json'))
    if len(paths) != 1:
        raise ValueError('The example requires one complete observation carrier')
    original = paths[0].read_text()
    query['law']['evidence'] = [dict(document='assay_declaration', start=0, end=len(original.encode()), text=original)]
    return query


def execute(binary, root, query):
    method = root/'method'
    method.mkdir()
    for name in ['replay.py', 'verify.py', 'pipeline.py', 'README.md']:
        shutil.copyfile(Path(__file__).resolve().parent/name, method/name)
    imported, import_seconds = run(binary, 'import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    save(root/'import-receipt.json', imported)
    bind_evidence(query, root/'workspace')
    save(root/'inputs/query.json', query)
    certificate, seconds = run(binary, 'assay-query', '--workspace', root/'workspace/workspace.json',
        '--request', save(root/'request.json', dict(schema='gadget.assays.query.v1', operation='native_query',
                                                 assay='measurements', query=query)))
    certificate = certificate['result']['answers'][0]['certificate']
    save(root/'certificate.json', certificate)
    save(root/'timings.json', dict(import_seconds=import_seconds, query_seconds=seconds,
                                 scope='Complete CLI source admission, model preparation, exact actions, checks and serialization'))
    return certificate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--expression-source', type=Path)
    p.add_argument('--features', type=Path)
    a = p.parse_args()
    if a.out.exists():
        raise ValueError('Choose a fresh reproduction directory')
    if bool(a.expression_source) != bool(a.features):
        raise ValueError('Expression input and declared features are both required')
    if a.expression_source:
        query = expression(a.out, a.expression_source, a.features)
    else:
        query = authored(a.out/'inputs')
    certificate = execute(a.binary.resolve(), a.out, query)
    print(json.dumps(dict(status=certificate['result']['status'], states=len(query['law']['states']),
                          hypotheses=len(query['hypotheses']), out=str(a.out.resolve()))))


if __name__ == '__main__':
    main()
