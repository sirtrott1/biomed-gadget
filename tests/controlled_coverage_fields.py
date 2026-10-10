#!/usr/bin/env python3
"""Verify complete control spans, rational contrasts and original byte bindings."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


gadget, verifier, root = map(Path, sys.argv[1:])
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)


def run(*command, success=True):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True)
    if (result.returncode == 0) != success:
        raise AssertionError(f'{command}: {result.stdout}\n{result.stderr}')
    return result


def native(name, value):
    compatibility, output = root / (name + '.json'), root / (name + '.rfcxs')
    compatibility.write_text(json.dumps(value) + '\n')
    run(gadget, 'native-record', '--input', compatibility, '--out', output)
    compatibility.unlink()
    return output


def artifact(path):
    return dict(path=path.name, sha256='sha256:' + hashlib.sha256(path.read_bytes()).hexdigest())


contents = {
    'base': b'b\t1\t5\tsame\t3\t-\n'
            b'a\t0\t4\tfirst\t2\t+\n'
            b'a\t0\t4\tduplicate\t2\t+\n'
            b'b\t10\t12\tisolated_zero\t0\t+\n',
    'treated': b'a\t2\t6\tfirst\t1\t+\n'
               b'b\t0\t4\tsecond\t2\t-\n'
               b'a\t1\t3\topposite\t4\t-\n',
    'replicate': b'a\t0\t3\tfirst\t2\t+\n'
                 b'b\t2\t6\tsecond\t3\t-\n'
                 b'a\t1\t3\topposite\t0\t-\n',
}
for name, data in contents.items():
    (root / (name + '.bed')).write_bytes(data)
fasta = root / 'reference.fa'
fasta.write_bytes(b'>a\nACGTACGTACGT\n>b\nACGTACGTACGT\n')
design = dict(question='Compare original measured fields with a declared control span',
              population='Authored fixture', selection='Every original row',
              limitations=['Supplied physical intervals'])
reference = dict(assembly='fixture', sequence_sha256=artifact(fasta)['sha256'],
                 coordinate_convention='zero_based_half_open',
                 contigs=[dict(id=name, length=12) for name in ('a', 'b')])
specification = native('import', dict(
    schema='gadget.assays.import.v1', id='controlled', design=design, assays=[], links=[],
    resources=[dict(id=name, source=artifact(root / (name + '.bed')), format='bed',
                    reference=reference, reference_source=artifact(fasta), weight_column=4,
                    method='Original measured interval weights') for name in contents]))
run(gadget, 'import-assays', '--spec', specification, '--out', root / 'imported')
workspace = root / 'imported/workspace.rfcxs'
samples = [dict(id=name, resource=name, scale=scale, scale_provenance='Declared exact divisor')
           for name, scale in [('base', '3/2'), ('treated', '2'), ('replicate', '1')]]
samples.append(dict(id='base_copy', resource='base', scale='3/2',
                    scale_provenance='Same retained resource, separately declared control'))
originals = [root / (sample['resource'] + '.bed') for sample in samples]
contrasts = [dict(id='change', terms=[dict(sample='treated', coefficient='1'),
                                    dict(sample='base', coefficient='-1')], evidence='Declared difference'),
             dict(id='combination', terms=[dict(sample='replicate', coefficient='-1/2'),
                                          dict(sample='treated', coefficient='2/3')], evidence='Declared contrast')]
recipe = dict(schema='gadget.assays.query.v1', operation='controlled_coverage',
              samples=samples, contrasts=contrasts, metric='reference_length',
              thresholds=['-1', '0', '2/3'], absence_policy='measured_zero_within_union',
              units='declared interval weight per divisor', evidence='Complete original union',
              source_kernel=False)
for case, controls, intercept, strand, rank in [
        ('dependent', ['base', 'base_copy'], True, 'separate', 2),
        ('empty', [], False, 'ignore', 0),
        ('intercept', [], True, 'ignore', 1)]:
    value = dict(recipe, controls=controls, intercept=intercept, strand_policy=strand)
    request = native(case + '-request', value)
    publication = root / (case + '-result.rfcxs')
    run(gadget, 'assay-query', '--workspace', workspace, '--request', request,
        '--out-record', publication)
    verified = run(verifier, publication, workspace, request, *originals)
    assert 'status\tPASS' in verified.stdout
    assert 'original_records\t14' in verified.stdout
    assert f'control_rank\t{rank}' in verified.stdout
    wrong = native(case + '-wrong', dict(value, units='Different unit declaration'))
    assert 'retained original bytes differ' in run(
        verifier, publication, workspace, wrong, *originals, success=False).stderr
    modified = root / 'modified.bed'
    modified.write_bytes(contents['base'].replace(b'same', b'else'))
    assert 'retained original bytes differ' in run(
        verifier, publication, workspace, request, modified, *originals[1:], success=False).stderr
    materialize = native(case + '-materialize', dict(
        schema='gadget.assays.materialize.v1', id='coverage', source_namespace=str(7190 + rank),
        channel='E', design=design, recipe=value))
    derived, replayed = root / (case + '-derived'), root / (case + '-replayed')
    run(gadget, 'derive-assays', '--workspace', workspace, '--request', materialize, '--out', derived)
    run(verifier, '--assay', 'coverage', derived / 'workspace.rfcxs', workspace, materialize, *originals)
    run(gadget, 'replay-derivation', '--workspace', derived / 'workspace.rfcxs',
        '--assay', 'coverage', '--out', replayed)
    run(verifier, '--assay', 'coverage', replayed / 'workspace.rfcxs', workspace, materialize, *originals)
print('Complete declared control spans, exact contrasts, source occurrences and restart PASS')
