#!/usr/bin/env python3
"""Prepare complete published spatial counts, coordinates and native requests."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import urllib.request

BASE = ('https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/'
        'V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_')
ARCHIVES = (
    ('spatial-matrix.tar.gz', 'filtered_feature_bc_matrix.tar.gz',
     'e619857b8505ae0090063b2d63211d6315cfa424d5534a396bc617aaa2cb1a2a'),
    ('spatial.tar.gz', 'spatial.tar.gz',
     '415530a9d31172f63fd358d5d1cefb0798698ef9c0e614f52825894cd45c61bf'),
)


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def extract(path, target, names):
    found = set()
    with tarfile.open(path) as archive:
        for member in archive:
            name = Path(member.name).name
            plain = name[:-3] if name.endswith('.gz') else name
            if not member.isfile() or plain not in names:
                continue
            if plain in found:
                raise ValueError(f'{path.name}: duplicate archive member {plain}')
            with archive.extractfile(member) as source:
                data = source.read()
            (target / plain).write_bytes(gzip.decompress(data) if name.endswith('.gz') else data)
            found.add(plain)
    if found != set(names):
        raise ValueError(f'{path.name}: original members missing {set(names) - found}')


def prepare(args):
    data = args.data.resolve()
    originals, inputs = data / 'originals', data / 'input'
    for name in ('import', 'spatial', 'materialize', 'sections'):
        if (data / (name + '.rfcxs')).exists():
            raise FileExistsError(data / (name + '.rfcxs'))
    originals.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(exist_ok=True)
    for name, remote, expected in ARCHIVES:
        path = originals / name
        if not path.exists():
            if not args.download:
                raise FileNotFoundError(f'{path}; supply the original archive or use --download')
            temporary = path.with_name(path.name + '.download')
            try:
                urllib.request.urlretrieve(BASE + remote + '?download=1', temporary)
                if digest(temporary) != expected:
                    raise ValueError(f'{name}: original byte hash differs')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        if digest(path) != expected:
            raise ValueError(f'{name}: original byte hash differs')
    extract(originals / 'spatial-matrix.tar.gz', inputs, ['matrix.mtx', 'barcodes.tsv', 'features.tsv'])
    extract(originals / 'spatial.tar.gz', inputs,
            ['tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_lowres_image.png'])
    features = list(csv.reader((inputs / 'features.tsv').open(), delimiter='\t'))
    barcodes = (inputs / 'barcodes.tsv').read_text().splitlines()
    positions = list(csv.reader((inputs / 'tissue_positions_list.csv').open()))
    by_barcode = {row[0]: row[1:] for row in positions}
    if len(by_barcode) != len(positions) or any(by_barcode[b][0] != '1' for b in barcodes):
        raise ValueError('Original barcode membership or tissue flag differs')
    (inputs / 'axis-features.tsv').write_text('id\tsymbol\ttype\n' +
                                             ''.join('\t'.join(row) + '\n' for row in features))
    (inputs / 'axis-observations.tsv').write_text(
        'id\tin_tissue\tarray_row\tarray_col\ty\tx\tcoordinate_frame\tcoordinate_units\n' +
        ''.join(b + '\t' + '\t'.join(by_barcode[b]) + '\tV19L29-097-B1\tfull_resolution_pixel\n'
                for b in barcodes))

    def artifact(path):
        return dict(path=str(path.relative_to(data)), sha256='sha256:' + digest(path))

    def native(name, value):
        authoring = data / (name + '.json')
        try:
            authoring.write_text(json.dumps(value) + '\n')
            subprocess.run([str(args.gadget.resolve()), 'native-record', '--input', str(authoring),
                            '--out', str(data / (name + '.rfcxs'))], check=True)
        finally:
            authoring.unlink(missing_ok=True)

    design = dict(question='Measure spatial expression variation on the published coordinate carrier',
                  population='10x Genomics breast cancer block A section 1',
                  selection='Every released filtered gene and tissue barcode',
                  limitations=['Publisher filtering and gene assignment are retained'])
    assay = dict(id='measurements', kind='spatial_expression', format='matrix_market',
                 source=artifact(inputs / 'matrix.mtx'), features=artifact(inputs / 'axis-features.tsv'),
                 observations=artifact(inputs / 'axis-observations.tsv'),
                 source_namespace=str(args.source_namespace), carrier='observations', quantity='umi_count',
                 units='gene_assigned_umi', channel='E', method='Published gene assigned UMI counts',
                 arithmetic='exact_decimal', implicit_entries='measured_zero', missing_tokens=[],
                 provenance=[artifact(originals / name) for name, _, _ in ARCHIVES] +
                            [artifact(inputs / name) for name in ['features.tsv', 'barcodes.tsv',
                             'tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_lowres_image.png']])
    native('import', dict(schema='gadget.assays.import.v1', id='spatial-expression', design=design,
                          assays=[assay], resources=[], links=[]))
    recipe = dict(schema='gadget.assays.query.v1', operation='spatial_relations', assay='measurements',
                  coordinates=['x', 'y'], frame='V19L29-097-B1', units='full_resolution_pixel',
                  metric_diagonal=['1', '1'], radius_squared=args.radius_squared,
                  relation_type='physical.spot.neighborhood',
                  evidence='Declared Euclidean pixel radius over every original tissue coordinate',
                  execution_certificate=True, charts='grade_one_stars')
    native('spatial', recipe)
    native('materialize', dict(schema='gadget.assays.materialize.v1', id='spatial-fields',
                               source_namespace=str(args.source_namespace + 1), design=design,
                               recipe=recipe, channel='E'))
    selected = []
    for symbol in args.gene:
        matches = [row[0] for row in features if row[1] == symbol]
        if len(matches) != 1:
            raise ValueError(f'{symbol}: expected one unambiguous original feature id, found {matches}')
        selected.append('all.' + matches[0])
    if len(set(selected)) != len(selected):
        raise ValueError('Repeated requested gene')
    native('sections', dict(schema='gadget.scientific.sequence.v1', assay='spatial-fields',
                            context_template=dict(fields=selected, methods=['local_sections', 'green_action', 'hodge'],
                                                  restriction='scalar_identity', metric='identity', green_lambda='1')))
    print(f'Prepared all {len(features)} genes and {len(barcodes)} tissue spots')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--gadget', required=True, type=Path)
    parser.add_argument('--source-namespace', required=True, type=int)
    parser.add_argument('--radius-squared', required=True)
    parser.add_argument('--gene', required=True, action='append')
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    if not 0 < args.source_namespace < 2 ** 64 - 1:
        parser.error('--source-namespace must leave room for the derived source namespace')
    prepare(args)
