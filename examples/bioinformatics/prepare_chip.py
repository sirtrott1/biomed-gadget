#!/usr/bin/env python3
"""Prepare complete published EGR1 and input control intervals on chromosome 21."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

BASE = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/'
LIBRARIES = (
    ('egr1_1', 'wgEncodeHudsonalphaChipSeqAlignmentsRep1K562Egr1.tagAlign.gz',
     'cd795b8bde13da1f649be9fe54ec9fa9d62120f2b729cd4c0932077f1ec15b7b'),
    ('egr1_2', 'wgEncodeHudsonalphaChipSeqAlignmentsRep2K562Egr1.tagAlign.gz',
     'e99d23965229041b88c3f1192c82d709dc1febb4f9d03bec1ca38502ffbfd052'),
    ('input_1', 'wgEncodeHudsonalphaChipSeqAlignmentsRep1K562ControlPcr2x.tagAlign.gz',
     'fb2f6ce3e15be8acb07abcf5da1fedd0603063b63498b62a71b35d0fdd751a53'),
    ('input_2', 'wgEncodeHudsonalphaChipSeqAlignmentsRep2K562ControlPcr2x.tagAlign.gz',
     '4b9a84bb619bbeadcaede8bd3c184741307092cde8de2e5a374c85dcad136b06'),
)
REFERENCE_SHA256 = '142c895916d4975bd8bc787c091f3af12c3fc04cb0a6d4ec08e777a4449c6b25'


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def fetch(path, url, expected, download):
    if not path.exists():
        if not download:
            raise FileNotFoundError(f'{path}; supply the original file or use --download')
        temporary = path.with_name(path.name + '.download')
        try:
            urllib.request.urlretrieve(url, temporary)
            if expected and digest(temporary) != expected:
                raise ValueError(f'{path.name}: original byte hash differs')
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
    if expected and digest(path) != expected:
        raise ValueError(f'{path.name}: original byte hash differs')


def prepare(data, gadget, source_namespace, download):
    data = data.resolve()
    for name in ('import', 'coverage', 'materialize'):
        if (data / (name + '.rfcxs')).exists():
            raise FileExistsError(data / (name + '.rfcxs'))
    originals, inputs = data / 'originals', data / 'input'
    originals.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(exist_ok=True)
    for name, url in [('index.html', BASE), ('files.txt', BASE + 'files.txt'),
                      ('md5sum.txt', BASE + 'md5sum.txt')]:
        fetch(originals / name, url, None, download)
    published = {line.split()[1]: line.split()[0]
                 for line in (originals / 'md5sum.txt').read_text().splitlines()}
    metadata = (originals / 'index.html').read_text()
    for sample, name, expected in LIBRARIES:
        fetch(originals / name, BASE + name, expected, download)
        with (originals / name).open('rb') as source:
            if hashlib.file_digest(source, 'md5').hexdigest() != published[name]:
                raise ValueError(f'{name}: publisher checksum differs')
        row = next(line for line in metadata.splitlines() if f'HREF="{name}"' in line)
        for fact in ('cell=K562;', 'protocol=PCR2x;', 'treatment=None;'):
            if fact not in row:
                raise ValueError(f'{name}: publisher design differs')
        if sample.startswith('egr1'):
            if 'antibody=Egr-1;' not in row or 'controlId=K562/None/Input/PCR2x;' not in row:
                raise ValueError(f'{name}: target or control binding differs')
        elif 'antibody=Input;' not in row:
            raise ValueError(f'{name}: input control declaration differs')
    fetch(originals / 'chr21.fa.gz',
          'https://hgdownload.soe.ucsc.edu/goldenPath/hg18/chromosomes/chr21.fa.gz',
          REFERENCE_SHA256, download)
    fasta = inputs / 'reference.fa'
    length = 0
    with gzip.open(originals / 'chr21.fa.gz', 'rb') as source, fasta.open('wb') as output:
        header = source.readline()
        if header.strip() != b'>chr21':
            raise ValueError('Reference contig header differs')
        output.write(header)
        for line in source:
            if line.startswith(b'>'):
                raise ValueError('Unexpected reference contig')
            length += len(line.strip())
            output.write(line)
    counts = []
    for sample, name, _ in LIBRARIES:
        complete = selected = integral = 0
        with gzip.open(originals / name, 'rb') as source, (inputs / (sample + '.bed')).open('wb') as output:
            for line in source:
                columns = line.rstrip(b'\r\n').split(b'\t')
                if len(columns) != 6:
                    raise ValueError('Publisher tagAlign must contain six BED columns')
                complete += 1
                if columns[0] == b'chr21':
                    start, end = int(columns[1]), int(columns[2])
                    if not 0 <= start < end <= length or columns[5] not in (b'+', b'-'):
                        raise ValueError('Invalid published interval')
                    output.write(line)
                    selected += 1
                    integral += end - start
        counts.append((sample, complete, selected, integral))

    def artifact(path):
        return dict(path=str(path.relative_to(data)), sha256='sha256:' + digest(path))

    def native(name, value):
        temporary = data / (name + '.compat.json')
        try:
            temporary.write_text(json.dumps(value, indent=2) + '\n')
            subprocess.run([str(gadget), 'native-record', '--input', str(temporary),
                            '--out', str(data / (name + '.rfcxs'))], check=True)
        finally:
            temporary.unlink(missing_ok=True)

    frame = dict(assembly='NCBI36 / UCSC hg18; complete chromosome 21',
                 sequence_sha256=artifact(fasta)['sha256'],
                 coordinate_convention='zero_based_half_open', contigs=[dict(id='chr21', length=length)])
    design = dict(question='Which observed EGR1 read coverage lies outside the declared input control span?',
                  population='K562 EGR1 ChIP replicates 1 and 2, with PCR2x input controls 1 and 2',
                  selection='Every chromosome 21 interval from all four complete publisher archives',
                  limitations=['Mapped read intervals are the measured quantity; no fragment extension is inferred',
                               'Publisher replicate labels are retained without inferring additional pairing',
                               'Coverage residuals describe these libraries and their declared reference domain'])
    resources = []
    for sample, name, _ in LIBRARIES:
        resources.append(dict(id=sample, format='bed', source=artifact(inputs / (sample + '.bed')),
                              reference=frame, reference_source=artifact(fasta),
                              method='One mapped read per original BED row; BED score is not a count',
                              provenance=[artifact(originals / name), artifact(originals / 'index.html'),
                                          artifact(originals / 'files.txt'), artifact(originals / 'md5sum.txt'),
                                          artifact(originals / 'chr21.fa.gz')]))
    native('import', dict(schema='gadget.assays.import.v1', id='chip-coverage', design=design,
                          assays=[], resources=resources, links=[]))
    contrasts = [dict(id=sample, terms=[dict(sample=sample, coefficient='1')],
                      evidence='Original EGR1 library and its measured read intervals')
                 for sample in ('egr1_1', 'egr1_2')]
    contrasts.append(dict(id='replicate_difference',
                          terms=[dict(sample='egr1_1', coefficient='1'), dict(sample='egr1_2', coefficient='-1')],
                          evidence='Difference between publisher labeled replicates on the same physical domain'))
    recipe = dict(schema='gadget.assays.query.v1', operation='controlled_coverage',
                  samples=[dict(id=sample, resource=sample, scale='1',
                                scale_provenance='Identity divisor preserving original read count units')
                           for sample, _, _ in LIBRARIES],
                  contrasts=contrasts, controls=['input_1', 'input_2'], intercept=False,
                  metric='reference_length', strand_policy='ignore', thresholds=['0'],
                  absence_policy='measured_zero_within_union', units='mapped read count',
                  evidence='Publisher matched K562 PCR2x input family; separate original control columns',
                  source_kernel=False)
    native('coverage', recipe)
    native('materialize', dict(schema='gadget.assays.materialize.v1', id='chip-coverage',
                               source_namespace=str(source_namespace), channel='E', design=design, recipe=recipe))
    print('sample\tcomplete_archive_records\tchr21_records\tweighted_base_integral')
    for row in counts:
        print(*row, sep='\t')
    print('reference_length\t' + str(length))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('--gadget', required=True, type=Path)
    parser.add_argument('--source-namespace', required=True, type=int)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    if not 0 < args.source_namespace < 2**64:
        parser.error('--source-namespace must be a nonzero unsigned 64 bit identity')
    prepare(args.data, args.gadget.resolve(), args.source_namespace, args.download)
