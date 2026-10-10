#!/usr/bin/env python3
"""Compare native alignment tables with every original BAM record."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
from pathlib import Path
import struct


def exact_read(stream, size):
    value = stream.read(size)
    if len(value) != size:
        raise ValueError('Truncated BAM record')
    return value


def read_group(data):
    position = 0
    group = ''
    widths = {'A': 1, 'c': 1, 'C': 1, 's': 2, 'S': 2,
              'i': 4, 'I': 4, 'f': 4, 'd': 8}
    while position < len(data):
        tag = data[position:position+2]
        kind = chr(data[position+2])
        position += 3
        if kind in widths:
            position += widths[kind]
        elif kind in ('Z', 'H'):
            end = data.index(0, position)
            if tag == b'RG':
                if kind != 'Z':
                    raise ValueError('Read group is not text')
                group = data[position:end].decode()
            position = end+1
        elif kind == 'B':
            subtype = chr(data[position])
            count = struct.unpack_from('<I', data, position+1)[0]
            position += 5+count*widths[subtype]
        else:
            raise ValueError(f'Unsupported BAM auxiliary type {kind}')
        if position > len(data):
            raise ValueError('Truncated BAM auxiliary field')
    return group


def bam_records(path):
    with gzip.open(path, 'rb') as stream:
        if exact_read(stream, 4) != b'BAM\1':
            raise ValueError('BAM magic')
        size = struct.unpack('<i', exact_read(stream, 4))[0]
        exact_read(stream, size)
        contigs = []
        count = struct.unpack('<i', exact_read(stream, 4))[0]
        for _ in range(count):
            size = struct.unpack('<i', exact_read(stream, 4))[0]
            name = exact_read(stream, size)
            length = struct.unpack('<i', exact_read(stream, 4))[0]
            contigs.append((name[:-1].decode(), length))
        index = 0
        while size_bytes := stream.read(4):
            if len(size_bytes) != 4:
                raise ValueError('Truncated BAM record length')
            size = struct.unpack('<i', size_bytes)[0]
            data = exact_read(stream, size)
            reference, start, bmq, fc, length, _, _, _ = struct.unpack_from('<iiIIiiii', data)
            flags, mapq, name_size, cigar_size = fc >> 16, (bmq >> 8) & 255, bmq & 255, fc & 65535
            query = data[32:32+name_size-1].decode()
            offset = 32+name_size
            cigar = struct.unpack_from('<'+'I'*cigar_size, data, offset)
            auxiliary = offset+4*cigar_size+(length+1)//2+length
            group = read_group(data[auxiliary:])
            blocks, junctions = [], []
            position = start
            if flags & 4:
                contig, strand, start, position = '*', '.', None, None
            else:
                contig = contigs[reference][0]
                strand = '-' if flags & 16 else '+'
                for word in cigar:
                    span, operation = word >> 4, word & 15
                    if operation not in range(9):
                        raise ValueError('Unsupported CIGAR operation')
                    if operation in (0, 2, 3, 7, 8):
                        if operation in (0, 7, 8):
                            blocks.append((position, position+span))
                        elif operation == 3:
                            junctions.append((position, position+span))
                        position += span
            yield dict(index=index, query=query, group=group, flags=flags, mapq=mapq,
                       contig=contig, strand=strand, start=start, end=position,
                       query_length=length, blocks=blocks, junctions=junctions)
            index += 1


def verify(bam, prefix, exclude_flags, minimum_mapq, include_unknown):
    junctions = defaultdict(list)
    blocks = []
    selected = []
    flags = Counter()
    qualities = Counter()
    spliced = 0
    with Path(str(prefix)+'-records.tsv').open() as output:
        count = 0
        for record in bam_records(bam):
            keep = not record['flags'] & exclude_flags and (
                include_unknown if record['mapq'] == 255 else record['mapq'] >= minimum_mapq)
            selected.append(int(keep))
            expected = [record[name] for name in ('index', 'query', 'group', 'flags', 'mapq',
                                                  'contig', 'strand', 'start', 'end', 'query_length')]
            expected.append(int(keep))
            expected = [str(value) if value is not None else 'NA' for value in expected]
            if output.readline().rstrip('\n').split('\t') != expected:
                raise AssertionError(f"Original BAM record {record['index']} differs")
            for start, end in record['blocks']:
                blocks.append((record['contig'], record['strand'], start, end, record['index']))
            for start, end in record['junctions']:
                junctions[(record['contig'], record['strand'], start, end)].append(record['index'])
            flags[record['flags']] += 1
            qualities[record['mapq']] += 1
            spliced += bool(record['junctions'])
            count += 1
        if output.read():
            raise AssertionError('Extra native alignment records')
    with Path(str(prefix)+'-relations.tsv').open() as output:
        index = 0
        for contig, strand, start, end, record in blocks:
            expected = list(map(str, [index, 0, contig, strand, start, end, selected[record], record]))
            if output.readline().rstrip('\n').split('\t') != expected:
                raise AssertionError(f'Aligned block {index} differs')
            index += 1
        for (contig, strand, start, end), records in sorted(junctions.items()):
            expected = list(map(str, [index, 1, contig, strand, start, end,
                                      sum(selected[record] for record in records), *records]))
            if output.readline().rstrip('\n').split('\t') != expected:
                raise AssertionError(f'Splice junction {index} differs')
            index += 1
        if output.read():
            raise AssertionError('Extra native relations')
    with bam.open('rb') as original:
        source_hash = hashlib.file_digest(original, 'sha256').hexdigest()
    return dict(records=count, aligned_blocks=len(blocks), spliced_records=spliced,
                junction_occurrences=sum(map(len, junctions.values())),
                distinct_stranded_junctions=len(junctions), selected_records=sum(selected),
                flags=dict(flags), mapping_quality=dict(qualities),
                original_sha256=source_hash)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--exclude-flags', type=int, required=True)
    parser.add_argument('--minimum-mapq', type=int, required=True)
    parser.add_argument('--unknown-mapq', choices=('include', 'exclude'), required=True)
    args = parser.parse_args()
    result = verify(args.bam, args.prefix, args.exclude_flags, args.minimum_mapq,
                    args.unknown_mapq == 'include')
    for name, value in result.items():
        print(f'{name}\t{value}')


if __name__ == '__main__':
    main()
