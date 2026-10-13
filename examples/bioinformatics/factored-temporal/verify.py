#!/usr/bin/env python3
"""Expand every factor and compare all original map entries and scientific evidence."""
import argparse,json
from pathlib import Path
from maps import check_source_bindings,expand_result
import sympy as S


def verify_pair(full, factored):
    a, b = full['result'], factored['result']
    assert full['source_digest'] == factored['source_digest']
    expanded = expand_result(b)
    def numerical_maps(value):
        if isinstance(value, dict):
            if all(key in value for key in ('rows', 'columns', 'entries')):
                return {key: value[key] for key in ('rows', 'columns', 'entries')}
            return {key: numerical_maps(item) for key, item in value.items()}
        if isinstance(value, list):
            return [numerical_maps(item) for item in value]
        return value
    assert numerical_maps(expanded) == numerical_maps(a), 'expanded maps or original evidence differ'
    for step in expanded['steps']:
        matrices=[factor['contribution'] for factor in step['word_factors']]
        from maps import sparse
        direct=sparse(step['word_defect'])
        assert sum((sparse(value) for value in matrices),S.zeros(direct.rows,direct.cols))==direct
    count = check_source_bindings(b['factored_maps'], b['endpoints'])
    for t, step in enumerate(b['steps']):
        count += check_source_bindings(step['structure']['factored_maps'], b['endpoints'][t:t+2])
    def entries(value):
        if isinstance(value, dict):
            return (value['rows'] * value['columns'] if 'entries' in value else
                    sum(entries(v) for v in value.values()))
        return sum(entries(v) for v in value) if isinstance(value, list) else 0
    return dict(status='PASS', complete_map_coordinates=entries(a), source_binding_coordinates=count)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('enumerated', type=Path)
    p.add_argument('factored', type=Path)
    a = p.parse_args()
    print(json.dumps(verify_pair(json.loads(a.enumerated.read_text()),json.loads(a.factored.read_text()))))
