#!/usr/bin/env python3
"""Reproduce physically constrained count, ordinal and measured dose models."""
import argparse
import importlib.util
from pathlib import Path

p = Path(__file__).resolve().parents[1] / 'observed-responses/replay.py'
spec = importlib.util.spec_from_file_location('observed_source', p)
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
save, digest, run = source.save, source.digest, source.run


def constrained(query, dose_responses=()):
    query['operation'] = 'constrained_response_model'
    query['physical_shape'] = dict(
        meaning='Nonnegative counts, original trial bounds and jointly ordered category fields; selected measured dose directions',
        evidence=[], dose_constraints=[dict(response=r, comparison='measured-dose') for r in dose_responses])
    return query


def execute(binary, root, query):
    source.run(binary, 'import-assays', '--spec', root/'inputs/import.json', '--out', root/'workspace')
    source.source.bind_evidence(query, root/'workspace')
    evidence = query['interpretation']['evidence']
    for response in query['responses']:
        response['evidence'] = evidence
    for dose in query.get('dose_comparisons', []):
        dose['evidence'] = evidence
    if 'reference_law' in query:
        query['reference_law']['evidence'] = evidence
    query['physical_shape']['evidence'] = evidence
    save(root/'inputs/query.json', query)
    request = dict(schema='gadget.assays.query.v1', operation='native_query', assay='measurements', query=query)
    out, elapsed = run(binary, 'assay-query', '--workspace', root/'workspace/workspace.json',
                       '--request', save(root/'request.json', request))
    save(root/'result.json', out)
    save(root/'certificate.json', out['result']['answers'][0]['certificate'])
    save(root/'timing.json', dict(query_seconds=elapsed))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--panel', type=Path)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--dose-response', action='append', default=[])
    args = parser.parse_args()
    root = args.out.resolve()
    if args.source:
        if not args.panel:
            parser.error('A declared expression panel is required')
        query = source.expression(root, args.source.resolve(), args.panel.resolve(), args.workers)
    else:
        query = source.authored(root, args.workers)
    execute(args.binary.resolve(), root, constrained(query, args.dose_response))
