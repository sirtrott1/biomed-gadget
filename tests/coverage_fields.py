#!/usr/bin/env python3
"""Check original interval accounting across contigs, zero support and restart."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


gadget, verifier, root = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
shutil.rmtree(root, ignore_errors=True)
root.mkdir(parents=True)


def run(*command, success=True):
    result = subprocess.run([str(part) for part in command], capture_output=True, text=True)
    if (result.returncode == 0) != success:
        raise AssertionError(f'{command}: {result.stdout}\n{result.stderr}')
    return result


def native(name, value):
    compat, output = root/(name+'.json'), root/(name+'.rfcxs')
    compat.write_text(json.dumps(value)+'\n')
    run(gadget, 'native-record', '--input', compat, '--out', output)
    compat.unlink()
    return output


def artifact(path):
    return dict(path=path.name, sha256='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest())


bed = root/'original.bed'
bed.write_bytes(b'b\t2\t6\tsame\t0\n'
                b'a\t0\t4\tfirst\t2\n'
                b'b\t0\t4\tother\t3\n'
                b'a\t0\t4\tduplicate\t2\n'
                b'a\t3\t8\tlast\t1\n'
                b'b\t10\t12\tisolated_zero\t0\n')
fasta = root/'reference.fa'
fasta.write_bytes(b'>a\nACGTACGTACGT\n>b\nACGTACGTACGT\n')
design = dict(question='Account for every original interval', population='Authored fixture',
              selection='All original rows', limitations=['Supplied coordinate intervals'])
frame = dict(assembly='fixture', sequence_sha256=artifact(fasta)['sha256'],
             coordinate_convention='zero_based_half_open',
             contigs=[dict(id=c, length=12) for c in ('a', 'b')])
spec = native('import', dict(schema='gadget.assays.import.v1', id='fixture', design=design,
                             assays=[], links=[], resources=[dict(id='source', format='bed',
                             source=artifact(bed), reference=frame, reference_source=artifact(fasta),
                             weight_column=4, method='Original intervals')]))
run(gadget, 'import-assays', '--spec', spec, '--out', root/'imported')
workspace = root/'imported/workspace.rfcxs'
recipe = dict(schema='gadget.assays.query.v1', operation='controlled_coverage',
              samples=[dict(id='observed', resource='source', scale='1', scale_provenance='Identity')],
              contrasts=[dict(id='identity', terms=[dict(sample='observed', coefficient='1')],
                              evidence='Original rows')], controls=[], intercept=False,
              metric='reference_length', strand_policy='ignore', thresholds=['0', '2', '7'],
              absence_policy='measured_zero_within_union', units='count',
              evidence='Full interval union', source_kernel=False)
request = native('coverage', recipe)
publication = root/'coverage-result.rfcxs'
run(gadget, 'assay-query', '--workspace', workspace, '--request', request,
    '--out-record', publication)
result = run(verifier, publication, bed, workspace, request)
assert 'status\tPASS' in result.stdout
assert 'original_records\t6' in result.stdout
assert 'physical_segments\t7' in result.stdout
assert 'reported_support\t8' in result.stdout
assert 'weighted_base_integral\t33' in result.stdout
assert 'quadrance\t113' in result.stdout
modified = root/'modified.bed'
modified.write_bytes(bed.read_bytes().replace(b'other', b'wrong'))
failure = run(verifier, publication, modified, success=False)
assert 'original byte owner changed' in failure.stderr
wrong = native('wrong-request', {**recipe, 'units': 'different count declaration'})
failure = run(verifier, publication, bed, workspace, wrong, success=False)
assert 'input request publication bytes changed' in failure.stderr
materialize = native('materialize', dict(schema='gadget.assays.materialize.v1',
                     id='coverage', source_namespace='2187', channel='E', design=design, recipe=recipe))
run(gadget, 'derive-assays', '--workspace', workspace, '--request', materialize,
    '--out', root/'derived')
derived_workspace = root/'derived/workspace.rfcxs'
run(verifier, '--assay', 'coverage', derived_workspace, bed, workspace, materialize)
replayed = root/'replayed'
run(gadget, 'replay-derivation', '--workspace', derived_workspace, '--assay', 'coverage',
    '--out', replayed)
run(verifier, '--assay', 'coverage', replayed/'workspace.rfcxs', bed, workspace, materialize)
print('Complete contig separation, source multiplicity, supported zeros and provenance PASS')
