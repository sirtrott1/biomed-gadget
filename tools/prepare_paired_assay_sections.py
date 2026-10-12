#!/usr/bin/env python3
"""Compile a complete captured RNA/antibody library-balance workflow."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return 'sha256:' + value.hexdigest()


def source_contract(folder):
    folder = folder.resolve(strict=True)
    specification = folder / 'paired-import.json'
    spec = json.loads(specification.read_text())
    if spec['schema'] != 'gadget.assays.import.v1' or len(spec['assays']) != 2:
        raise ValueError('paired capture requires two declared assay sources')
    assays = {a['id']: a for a in spec['assays']}
    if set(assays) != {'rna', 'antibody'}:
        raise ValueError('paired capture requires declared rna and antibody sources')
    hashes = {}

    def artifact(record):
        path = (folder / record['path']).resolve(strict=True)
        if folder not in path.parents or not path.is_file():
            raise ValueError('capture artifact is outside its source directory')
        actual = digest(path)
        if actual != record['sha256']:
            raise ValueError('capture artifact hash differs: ' + record['path'])
        hashes[record['path']] = actual
        return path

    def axis(record):
        with artifact(record).open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        ids = [row['id'] for row in rows]
        if not ids or any(not value for value in ids) or len(set(ids)) != len(ids):
            raise ValueError('capture axis must retain distinct nonempty identities')
        return ids

    features, observations, provenance = {}, {}, {}
    for name, assay in assays.items():
        if (assay['format'] != 'matrix_market' or assay['carrier'] != 'features' or
                assay['quantity'] != 'umi_count' or assay['implicit_entries'] != 'measured_zero'):
            raise ValueError('paired workflow requires the declared feature UMI capture profile')
        artifact(assay['source'])
        features[name] = axis(assay['features'])
        observations[name] = axis(assay['observations'])
        provenance[name] = {p['path']: p for p in assay['provenance']}
        for record in assay['provenance']:
            artifact(record)
    if assays['rna']['source_namespace'] == assays['antibody']['source_namespace']:
        raise ValueError('paired assays require separate native source identities')
    if assays['rna']['units'] == assays['antibody']['units']:
        raise ValueError('paired capture must retain separate modality unit declarations')
    for name in ('barcodes.tsv', 'features.tsv', 'matrix.mtx'):
        left, right = provenance['rna'].get(name), provenance['antibody'].get(name)
        if left is None or left != right:
            raise ValueError('matching requires one shared original capture artifact: ' + name)
    original_barcodes = artifact(provenance['rna']['barcodes.tsv']).read_text().splitlines()
    if observations['rna'] != original_barcodes or observations['antibody'] != original_barcodes:
        raise ValueError('modality axes must preserve the same original capture column order')
    with artifact(provenance['rna']['features.tsv']).open(newline='') as handle:
        original_features = list(csv.reader(handle, delimiter='\t'))
    by_type = {'Gene Expression': [], 'Antibody Capture': []}
    for row in original_features:
        if len(row) != 3 or row[2] not in by_type:
            raise ValueError('capture feature type is outside the selected paired profile')
        by_type[row[2]].append(row[0])
    if features['rna'] != by_type['Gene Expression'] or features['antibody'] != by_type['Antibody Capture']:
        raise ValueError('modality features must retain the complete original type partition')
    correspondence = dict(schema='gadget.paired.capture.correspondence.v1',
        original_barcode_source=provenance['rna']['barcodes.tsv'],
        original_matrix_source=provenance['rna']['matrix.mtx'],
        method='Each modality preserves the same original captured matrix column; labels alone do not identify independent sources',
        pairs=[dict(column=i, section=barcode,
            rna=dict(source_namespace=assays['rna']['source_namespace'], field='all.'+barcode),
            antibody=dict(source_namespace=assays['antibody']['source_namespace'], field='all.'+barcode))
            for i, barcode in enumerate(original_barcodes)])
    return dict(specification=specification, design=spec['design'], assays=assays,
        features=features, observations=original_barcodes, hashes=hashes,
        correspondence=correspondence)


def pipeline_design(contract, namespace, workers, correspondence_sha):
    assays, features, observations = contract['assays'], contract['features'], contract['observations']
    design = dict(question='Describe the declared RNA library-total versus antibody assignment balance for every captured barcode',
        population=contract['design']['population'], selection='All captured barcodes, every RNA feature and every antibody feature',
        limitations=['Library-assignment contrast with explicit unit domains and declared positive costs; no biochemical calibration or cell-type inference',
            'One captured donor and upstream barcode/UMI selection; no independence, treatment or clinical endpoint assertion'])
    aggregation = dict(schema='gadget.assays.materialize.v1', id='rna-library',
        source_namespace=namespace, design=design, channel=assays['rna']['channel'], recipe=dict(
            schema='gadget.assays.query.v1', operation='aggregate_assay', assay='rna',
            feature_groups=[dict(id='RNA_library_total', members=[dict(id=f, coefficient='1') for f in features['rna']])],
            observation_groups=[dict(id=s, members=[dict(id=s, coefficient='1')]) for s in observations],
            map_method='Sum every original RNA feature, preserving each original barcode',
            map_evidence='Complete original captured feature partition and barcode column order; no feature or cell filtering',
            support_policy='require_all_nonzero_inputs', workers=workers))
    binding = dict(schema='gadget.assays.bindings.v1', id='matched-capture', design=design,
        assays=[dict(input='rna-library', assay='rna-library', id='rna-library'),
                dict(input='source', assay='antibody', id='antibody')], resources=[], links=[],
        method='Retain independent native owners on the explicitly shared original capture columns',
        evidence='Original barcode and matrix hashes plus complete correspondence ' + correspondence_sha)
    penalty = dict(quantity='declared_library_balance_cost', units='declared squared count',
        method='Declared positive identity cost in each separate measurement domain',
        evidence='Comparison choice; not fitted sampling variance or modality calibration')
    interface = dict(id='library-balance', axes=['assigned_count_balance'],
        quantity='declared_library_assignment_balance', units='declared count', channel='observation',
        method='Signed RNA library-total minus the sum of antibody assignments',
        evidence='Declared library-count contrast, without identifying the assay units', participants=[])
    views = []
    for i, name in enumerate(('rna', 'antibody')):
        assay = 'rna-library' if name == 'rna' else 'antibody'
        ids = ['RNA_library_total'] if name == 'rna' else features['antibody']
        quantity = 'declared_linear_aggregate_of_umi_count' if name == 'rna' else assays[name]['quantity']
        meaning = dict(quantity=quantity, units=assays[name]['units'], channel=assays[name]['channel'])
        metric = dict(diagonal=['1']*len(ids), domain_quantity=quantity, domain_units=meaning['units'],
            quadrance_quantity=penalty['quantity'], quadrance_units=penalty['units'],
            method=penalty['method'], evidence=penalty['evidence'])
        views.append(dict(id=name, assay=assay, partition=0,
            fields=[dict(section=s, field='all.'+s) for s in observations], metric=metric))
        transport = dict(**{'domain_'+k: v for k, v in meaning.items()},
            **{'codomain_'+k: interface[k] for k in ('quantity', 'units', 'channel')},
            method=interface['method'], evidence=interface['evidence'],
            entries=[dict(observable_id='assigned_count_balance', source_id=f, coefficient='1') for f in ids])
        interface['participants'].append(dict(key=str(i+1), view=name, exists=True,
            orientation='positive' if name == 'rna' else 'negative', magnitude='1', map=transport))
    sections = dict(schema='gadget.assays.query.v1', operation='assay_sections', views=views,
        sections=observations, interfaces=[interface], metric=penalty,
        support_policy='require_all_nonzero_occurrences', workers=workers, complete_families=True)
    return dict(schema='gadget.pipeline.design.v1', profile='multi_omics', question=design['question'],
        steps=[dict(id='rna-library', type='derive', input='source', request=aggregation),
               dict(id='matched', type='bind', input='source', request=binding),
               dict(id='sections', type='query', input='matched', request=sections)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source-namespace', default='976300')
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output directory already exists')
    if args.workers < 1:
        parser.error('worker count must be positive')
    contract = source_contract(args.source)
    args.out.mkdir(parents=True)
    correspondence = args.out/'correspondence.json'
    correspondence.write_text(json.dumps(contract['correspondence'], indent=2)+'\n')
    design = pipeline_design(contract, args.source_namespace, args.workers, digest(correspondence))
    (args.out/'design.json').write_text(json.dumps(design, indent=2)+'\n')
    binary = args.binary.resolve(strict=True)

    def run(name, *arguments):
        result = subprocess.run([str(binary), *map(str, arguments)], capture_output=True, text=True)
        (args.out/(name+'.log')).write_text(result.stderr)
        if result.returncode:
            raise RuntimeError(name+' refused: '+result.stderr.strip())
        (args.out/(name+'.json')).write_text(result.stdout)

    run('import', 'import-assays', '--spec', contract['specification'], '--out', args.out/'source')
    workspace = args.out/'source/workspace.rfcxs'
    run('plan', 'pipeline-prepare', '--workspace', workspace, '--request', args.out/'design.json')
    run('execution', 'pipeline', '--workspace', workspace, '--request', args.out/'plan.json', '--out', args.out/'analysis')
    receipt = dict(schema='gadget.paired.capture.workflow.v1', source_hashes=contract['hashes'],
        sections=len(contract['observations']), features={k: len(v) for k,v in contract['features'].items()},
        correspondence_sha256=digest(correspondence), native_result='analysis/sections/result.rfcxs',
        scope='Complete captured inputs and declared library-assignment balance; no clinical or biochemical calibration claim')
    (args.out/'source-receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
