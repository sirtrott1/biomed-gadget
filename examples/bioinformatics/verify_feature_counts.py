#!/usr/bin/env python3
"""Verify every native assignment class, overlap witness and marginal bound."""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path
import re

from verify_alignment_fields import bam_records


def annotation_features(path, attribute):
    features = {}
    for line in path.open():
        if not line.strip() or line.startswith('#'):
            continue
        values = line.rstrip('\n').split('\t')
        if values[2] != 'exon':
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', values[8]))
        key = attrs[attribute]
        if key not in features:
            features[key] = dict(contig=values[0], strand=values[6], intervals=[])
        feature = features[key]
        if (feature['contig'], feature['strand']) != (values[0], values[6]):
            raise ValueError('Feature spans different contigs or strands')
        feature['intervals'].append((int(values[3])-1, int(values[4])))
    for feature in features.values():
        union = []
        for start, end in sorted(feature['intervals']):
            if union and start <= union[-1][1]:
                union[-1] = (union[-1][0], max(end, union[-1][1]))
            else:
                union.append((start, end))
        feature['intervals'] = union
    return features


def verify(bam, annotation, prefix, attribute, all_bases, exon_chain,
           strand_policy, exclude_flags, minimum_mapq, include_unknown):
    features = annotation_features(annotation, attribute)
    ordered = list(features)
    feature_index = {name: i for i, name in enumerate(ordered)}
    intervals = defaultdict(list)
    for name, feature in features.items():
        for start, end in feature['intervals']:
            intervals[feature['contig']].append((start, end, name))
    indices = {}
    for contig, entries in intervals.items():
        entries.sort()
        prefix_ends, end = [], 0
        for start, stop, name in entries:
            end = max(end, stop)
            prefix_ends.append(end)
        indices[contig] = (entries, [e[0] for e in entries], prefix_ends)
    candidates = defaultdict(set)
    source_records = defaultdict(list)
    overlaps = 0
    selected = 0
    with Path(str(prefix)+'-overlaps.tsv').open() as reported:
        for record in bam_records(bam):
            if record['flags'] & exclude_flags or (
                not include_unknown if record['mapq'] == 255 else record['mapq'] < minimum_mapq):
                continue
            selected += 1
            mate = 1 if record['flags'] & 64 else 2 if record['flags'] & 128 else 0
            identity = (record['group'], record['query'], mate)
            source_records[identity].append(record['index'])
            possible = set()
            if record['contig'] in indices:
                entries, starts, ends = indices[record['contig']]
                for start, end in record['blocks']:
                    begin = bisect_right(ends, start)
                    stop = bisect_left(starts, end)
                    for entry in entries[begin:stop]:
                        if entry[1] > start:
                            possible.add(entry[2])
            assigned = False
            for name in sorted(possible, key=feature_index.get):
                feature = features[name]
                if strand_policy != 'ignore' and (
                    feature['strand'] == '.' or
                    ((strand_policy == 'same') != (feature['strand'] == record['strand']))):
                    continue
                common = sum(max(0, min(end, right)-max(start, left))
                             for start, end in record['blocks']
                             for left, right in feature['intervals'])
                aligned = sum(end-start for start, end in record['blocks'])
                if not common:
                    continue
                compatible = True
                if exon_chain:
                    union = feature['intervals']
                    gaps = {(union[i-1][1], union[i][0]) for i in range(1, len(union))}
                    observed = set(record['junctions'])
                    compatible = all(gap in gaps for gap in observed)
                    compatible &= all(gap in observed for gap in gaps
                                      if record['start'] < gap[1] and record['end'] > gap[0])
                    compatible &= any(a <= record['start'] < b for a, b in union)
                    compatible &= any(a <= record['end']-1 < b for a, b in union)
                eligible = (not all_bases or common == aligned) and compatible
                expected = [record['index'], name, common, aligned, int(eligible)]
                if exon_chain:
                    expected.append(int(compatible))
                if reported.readline().rstrip('\n').split('\t') != list(map(str, expected)):
                    raise AssertionError(f"Overlap witness for record {record['index']}, {name} differs")
                overlaps += 1
                if eligible:
                    candidates[identity].add(name)
                    assigned = True
            if not assigned:
                candidates[identity].add('__unassigned__')
        if reported.read():
            raise AssertionError('Extra overlap witnesses')
    lower, upper = Counter(), Counter()
    class_counts = Counter()
    for identity in source_records:
        choices = candidates[identity]
        class_counts[frozenset(choices)] += 1
        for name in choices:
            if name != '__unassigned__':
                upper[name] += 1
                if choices == {name}:
                    lower[name] += 1
    observed_bounds = {}
    with Path(str(prefix)+'-bounds.tsv').open() as reported:
        for line in reported:
            name, lo, hi = line.rstrip('\n').split('\t')
            if name in observed_bounds:
                raise AssertionError('Duplicate feature bound')
            observed_bounds[name] = (lo, hi)
    expected_bounds = {name: (str(lower[name]), str(upper[name])) for name in ordered}
    expected_bounds['__unassigned__'] = ('NA', 'NA')
    if observed_bounds != expected_bounds:
        raise AssertionError('Complete marginal bound fields differ')
    classes = {}
    with Path(str(prefix)+'-classes.tsv').open() as reported:
        for line in reported:
            index, count, *choices = line.rstrip('\n').split('\t')
            key = frozenset(choices)
            if len(key) != len(choices) or key not in class_counts or int(count) != class_counts[key]:
                raise AssertionError('Native class boundary or count differs')
            classes[int(index)] = key
    if len(classes) != len(class_counts):
        raise AssertionError('Incomplete native classes')
    with Path(str(prefix)+'-segments.tsv').open() as reported:
        for i, identity in enumerate(sorted(source_records)):
            values = reported.readline().rstrip('\n').split('\t')
            _, query, group, mate, class_index, *records = values
            if int(values[0]) != i or (group, query, int(mate)) != identity:
                raise AssertionError('Segment identity differs')
            if classes[int(class_index)] != candidates[identity] or list(map(int, records)) != source_records[identity]:
                raise AssertionError('Segment class or original records differ')
        if reported.read():
            raise AssertionError('Extra sequenced segments')
    return dict(features=len(features), selected_alignment_records=selected,
                sequenced_segments=len(source_records), overlap_witnesses=overlaps,
                assignment_classes=len(classes),
                ambiguous_segments=sum(count for key, count in class_counts.items() if len(key) > 1),
                unassigned_only_segments=class_counts.get(frozenset({'__unassigned__'}), 0),
                identifiable_features=sum(lower[name] == upper[name] for name in ordered),
                features_with_positive_lower_bound=sum(lower[name] > 0 for name in ordered))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, required=True)
    parser.add_argument('--annotation', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--attribute', required=True)
    parser.add_argument('--overlap', choices=('any', 'all'), required=True)
    parser.add_argument('--exon-chain', action='store_true')
    parser.add_argument('--strand-policy', choices=('ignore', 'same', 'opposite'), required=True)
    parser.add_argument('--exclude-flags', type=int, required=True)
    parser.add_argument('--minimum-mapq', type=int, required=True)
    parser.add_argument('--unknown-mapq', choices=('include', 'exclude'), required=True)
    args = parser.parse_args()
    result = verify(args.bam, args.annotation, args.prefix, args.attribute,
                    args.overlap == 'all', args.exon_chain, args.strand_policy,
                    args.exclude_flags, args.minimum_mapq, args.unknown_mapq == 'include')
    for name, value in result.items():
        print(f'{name}\t{value}')


if __name__ == '__main__':
    main()
