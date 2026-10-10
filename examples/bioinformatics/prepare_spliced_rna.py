#!/usr/bin/env python3
"""Prepare complete published chromosome 4 alignment inputs and native requests."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request

SOURCES = {
    'pasillaBamSubset_0.48.0.tar.gz': (
        'https://bioconductor.posit.co/packages/3.22/data/experiment/src/contrib/pasillaBamSubset_0.48.0.tar.gz',
        '3ea580497f5d65ffe2f73370ce69d1ab91952559ebca9c468992f2fdab90fa41'),
    'dm3.ensGene.gtf.gz': (
        'https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/genes/dm3.ensGene.gtf.gz',
        '5da9156a192528ccecfaee08633d0a3c1c8a7265148d5b08393c3a4b025c62d1'),
    'dm3.fa.gz': (
        'https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz',
        'd623d985911191982f26d38eb04c681fe0fc6b0c1808404c62d0450cb4e47065'),
}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def prepare(data, gadget, download):
    data = data.resolve()
    originals = data/'originals'
    originals.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in SOURCES.items():
        path = originals/name
        if not path.exists():
            if not download:
                raise FileNotFoundError(f'{path}; use --download or supply the published file')
            temporary = path.with_name(path.name+'.download')
            try:
                urllib.request.urlretrieve(url, temporary)
                if digest(temporary) != expected:
                    raise ValueError(f'{name}: published byte hash differs')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        if digest(path) != expected:
            raise ValueError(f'{name}: published byte hash differs')
    with tarfile.open(originals/'pasillaBamSubset_0.48.0.tar.gz') as archive:
        archive.extractall(originals, filter='data')
    fasta = originals/'dm3.fa'
    if not fasta.exists():
        with gzip.open(originals/'dm3.fa.gz', 'rb') as source, fasta.open('wb') as output:
            shutil.copyfileobj(source, output)
    if digest(fasta) != 'c87fc4fae9c8576396bddada408b92dc384756721fd84171ddf25ca658784736':
        raise ValueError('Complete reference FASTA differs')
    inputs = data/'input'
    inputs.mkdir(exist_ok=True)
    annotation = inputs/'chr4.ensGene.gtf'
    with gzip.open(originals/'dm3.ensGene.gtf.gz', 'rt') as source, annotation.open('w') as output:
        for line in source:
            if line.startswith('chr4\t'):
                output.write(line)
    contigs, chromosome = {}, []
    for line in fasta.open():
        if line.startswith('>'):
            key = line[1:].strip().split()[0]
            contigs[key] = 0
        else:
            contigs[key] += len(line.strip())
            if key == 'chr4':
                chromosome.append(line.strip().upper())
    bundled = ''.join(line.strip().upper() for line in
                      (originals/'pasillaBamSubset/inst/extdata/dm3_chr4.fa').open()
                      if not line.startswith('>'))
    if ''.join(chromosome) != bundled:
        raise ValueError('Published chromosome 4 sequence and dm3 reference differ')

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

    frame = dict(assembly='BDGP Release 5 / UCSC dm3', sequence_sha256=artifact(fasta)['sha256'],
                 coordinate_convention='zero_based_half_open',
                 contigs=[dict(id=k, length=v) for k, v in contigs.items()])
    resources = [dict(id='reference', format='fasta', source=artifact(fasta),
                      method='Complete UCSC dm3 reference; chr4 matches the published BAM package',
                      provenance=[artifact(originals/'dm3.fa.gz'),
                                  artifact(originals/'pasillaBamSubset_0.48.0.tar.gz')]),
                 dict(id='annotation', format='gtf', source=artifact(annotation), reference=frame,
                      method='All chr4 records from the complete UCSC dm3 ensGene GTF',
                      provenance=[artifact(originals/'dm3.ensGene.gtf.gz')])]
    for name in ('untreated1', 'untreated3'):
        path = originals/('pasillaBamSubset/inst/extdata/'+name+'_chr4.bam')
        resources.append(dict(id=name, format='bam', source=artifact(path), reference=frame,
                              method='Complete published chromosome 4 BAM; original TopHat 1.2.0 header retained'))
    design = dict(question='Account for every supplied alignment, junction and assignment candidate',
                  population='Published Pasilla untreated1 and untreated3 chromosome 4 files',
                  selection='All alignments and all chromosome 4 annotation records',
                  limitations=['Published chromosome subsets, conditional on supplied mappings and annotation',
                               'Paired ends are separate sequenced segments'])
    native('import', dict(schema='gadget.assays.import.v1', id='spliced-rna', design=design,
                          assays=[], resources=resources, links=[]))
    policy = dict(exclude_flags=3588, minimum_mapping_quality=0, unknown_mapping_quality='include')
    for name in ('untreated1', 'untreated3'):
        native(name+'-alignment', dict(schema='gadget.assays.query.v1', operation='alignment_evidence',
                                       resource=name, filter=policy))
        for mode in ('gene', 'transcript'):
            request = dict(schema='gadget.assays.query.v1', operation='feature_counts', resource=name,
                           annotation='annotation', feature_type='exon',
                           group_attribute='gene_id' if mode == 'gene' else 'transcript_id',
                           unit='sequenced_segment', overlap='any_aligned_base' if mode == 'gene' else 'all_aligned_bases',
                           strand_policy='ignore', filter=policy,
                           candidate_search=dict(method='Original TopHat 1.2.0 alignments',
                                                 completeness='All mappings in the complete published chr4 BAM',
                                                 reference_sha256=frame['sequence_sha256']))
            if mode == 'transcript':
                request['junction_policy'] = 'declared_exon_chain'
            native(name+'-'+mode, request)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--gadget', type=Path, required=True)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    prepare(args.data, args.gadget.resolve(), args.download)


if __name__ == '__main__':
    main()
