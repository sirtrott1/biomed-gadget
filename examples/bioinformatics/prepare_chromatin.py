#!/usr/bin/env python3
"""Prepare complete published chr2 and chr3 fragments and native declarations."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request


SOURCES = {
    'fragments.tsv.gz': (
        'https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz',
        '47cf3697d4a1a8afb68debf346c4e7aeb302add102d9478c36a236ace77713ef'),
    'chr2.fa.gz': (
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr2.fa.gz',
        'f831824c5bdcc650c62ae2db0812d04bc2933c20c7cce43e43b25c3167d84d55'),
    'chr3.fa.gz': (
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz',
        '1edd1dc9c69454208b246f9cacfb427b4f79ac18c776186a58dbc9816fb1fad4'),
    'web-summary.html': (
        'https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_web_summary.html',
        'bbc432565ea4fea6d1e4209e4fb53938b0ffc9d9f03edf4c8c4bd9f47021338a'),
}
CONTIGS = ('chr2', 'chr3')


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def prepare(data, gadget, download, fragments):
    data = data.resolve()
    originals = data/'originals'
    originals.mkdir(parents=True, exist_ok=True)
    for name in ('import', 'coverage', 'materialize'):
        destination = data/(name+'.rfcxs')
        if destination.exists():
            raise FileExistsError(destination)
    if fragments and not (originals/'fragments.tsv.gz').exists():
        try:
            os.link(fragments, originals/'fragments.tsv.gz')
        except OSError:
            shutil.copyfile(fragments, originals/'fragments.tsv.gz')
    for name, (url, expected) in SOURCES.items():
        path = originals/name
        if not path.exists():
            if not download:
                raise FileNotFoundError(f'{path}; supply the published file or use --download')
            temporary = path.with_name(path.name+'.download')
            try:
                urllib.request.urlretrieve(url, temporary)
                if digest(temporary) != expected:
                    raise ValueError(f'{name}: original byte hash differs')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        if digest(path) != expected:
            raise ValueError(f'{name}: original byte hash differs')
    inputs = data/'input'
    inputs.mkdir(exist_ok=True)
    bed = inputs/'fragments.bed'
    counts = {contig: [0, 0, 0] for contig in CONTIGS}
    with gzip.open(originals/'fragments.tsv.gz', 'rb') as source, bed.open('wb') as output:
        for line in source:
            fields = line.rstrip(b'\n').split(b'\t')
            if len(fields) != 5:
                raise ValueError('Unexpected publisher fragment row')
            contig = fields[0].decode()
            if contig in counts:
                start, end, support = map(int, (fields[1], fields[2], fields[4]))
                if not 0 <= start < end or support <= 0:
                    raise ValueError('Invalid published fragment')
                output.write(line)
                counts[contig][0] += 1
                counts[contig][1] += support
                counts[contig][2] += (end-start)*support
    fasta = inputs/'reference.fa'
    lengths = {}
    with fasta.open('wb') as output:
        for contig in CONTIGS:
            length = 0
            with gzip.open(originals/(contig+'.fa.gz'), 'rb') as source:
                first = source.readline()
                if first.strip() != ('>'+contig).encode():
                    raise ValueError('Reference contig header differs')
                output.write(first)
                for line in source:
                    if line.startswith(b'>'):
                        raise ValueError('Unexpected reference sequence')
                    length += len(line.strip())
                    output.write(line)
            lengths[contig] = length

    def artifact(path):
        return dict(path=str(path.relative_to(data)), sha256='sha256:'+digest(path))

    def native(name, value):
        temporary = data/(name+'.compat.json')
        destination = data/(name+'.rfcxs')
        if destination.exists():
            raise FileExistsError(destination)
        try:
            temporary.write_text(json.dumps(value, indent=2)+'\n')
            subprocess.run([str(gadget), 'native-record', '--input', str(temporary),
                            '--out', str(destination)], check=True)
        finally:
            temporary.unlink(missing_ok=True)

    frame = dict(assembly='GRCh37 / UCSC hg19; complete chr2 and chr3',
                 sequence_sha256=artifact(fasta)['sha256'],
                 coordinate_convention='zero_based_half_open',
                 contigs=[dict(id=contig, length=lengths[contig]) for contig in CONTIGS])
    design = dict(question='Account for every reported fragment and physical incidence on chr2 and chr3',
                  population='Published Next GEM v1.1 PBMC library',
                  selection='Every chr2 and chr3 row; all barcodes and reported multiplicities',
                  limitations=['The question is conditional on supplied mapped fragment coordinates'])
    resource = dict(id='fragments', format='bed', source=artifact(bed), reference=frame,
                    reference_source=artifact(fasta), weight_column=4,
                    method='Original interval and barcode; weight is reported read pair multiplicity',
                    provenance=[artifact(originals/name) for name in SOURCES])
    native('import', dict(schema='gadget.assays.import.v1', id='chromatin-fragments',
                          design=design, assays=[], resources=[resource], links=[]))
    recipe = dict(schema='gadget.assays.query.v1', operation='controlled_coverage',
                  samples=[dict(id='published', resource='fragments', scale='1',
                                scale_provenance='Identity divisor preserving reported read support')],
                  contrasts=[dict(id='reported_support', terms=[dict(sample='published', coefficient='1')],
                                  evidence='Every source interval with its original integer multiplicity')],
                  controls=[], intercept=False, metric='reference_length', strand_policy='ignore',
                  thresholds=['1', '2', '10'], absence_policy='measured_zero_within_union',
                  units='reported read pair support',
                  evidence='Complete chr2 and chr3; outside the observed interval union is unobserved',
                  source_kernel=False)
    native('coverage', recipe)
    native('materialize', dict(schema='gadget.assays.materialize.v1', id='fragment-coverage',
                               source_namespace='908104', channel='E', design=design, recipe=recipe))
    print('contig\trecords\treported_support\tweighted_base_integral\treference_length')
    for contig, values in counts.items():
        print(contig, *values, lengths[contig], sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--gadget', required=True, type=Path)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--fragments', type=Path,
                        help='Reuse the complete publisher archive, checked against its pinned hash')
    args = parser.parse_args()
    prepare(args.data, args.gadget.resolve(), args.download, args.fragments)
