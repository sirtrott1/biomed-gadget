"""Complete native section and reference ownership with portable replay."""
import copy
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(binary, consumer, root):
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    demo = load('native_structural', Path(__file__).resolve().parents[1] /
                'examples/bioinformatics/structural-calibration/replay.py')
    structural = root / 'structural'
    query = demo.authored(structural, workers=2)
    demo.execute(binary, structural, query)
    queries = {'complete': query}
    compact = copy.deepcopy(query)
    compact['reference_output'] = 'statistics'
    queries['compact'] = compact
    undefined = copy.deepcopy(query)
    undefined.pop('uncertainty')
    undefined['fields'] = ['all.zero']
    undefined['hypotheses'] = [dict(id='zero', field='all.zero', statistic='sector_fraction', sector='harmonic')]
    queries['undefined'] = undefined
    unstable = copy.deepcopy(query)
    unstable.pop('uncertainty')
    unstable['fields'] = ['all.probe']
    unstable['hypotheses'] = [dict(id='curl-probe', field='all.probe', statistic='probe_fraction',
                                  sector='curl', probes=['all.A'])]
    queries['partial'] = unstable
    zero = copy.deepcopy(query)
    zero['hypotheses'] = [dict(id='zero-probe', field='all.A', statistic='probe_fraction',
                              sector='harmonic', probes=['all.zero'])]
    queries['supported-zero'] = zero
    planned = copy.deepcopy(query)
    planned['execution'] = 'plan'
    queries['planned'] = planned
    receipts = []

    def verify(case, inputs, mode):
        path = root / (mode + '-queries.json')
        path.write_text(json.dumps(inputs))
        process = subprocess.run([str(consumer), str(case / 'workspace/workspace.rfcxs'),
                                  str(path), str(root / (mode + '-results')), mode],
                                 capture_output=True, text=True, timeout=300)
        if process.returncode:
            raise AssertionError(mode + ': ' + process.stderr)
        receipts.append(dict(case=mode, **json.loads(process.stdout)))
        results = root / (mode + '-results')
        workflow = subprocess.run([str(binary), 'workflow', '--workspace',
                                   str(case / 'workspace.moved/workspace.rfcxs'), '--assay',
                                   'measurements', '--request', str(results / 'workflow-request.rfcxs'),
                                   '--out', str(results / 'cli-workflow.rfcxs')],
                                  capture_output=True, text=True, timeout=300)
        if workflow.returncode:
            raise AssertionError(mode + ' workflow CLI: ' + workflow.stderr)
        view = json.loads(workflow.stdout)
        expected = json.loads((results / 'workflow-display.json').read_text())
        assert view['schema'] == 'gadget.scientific.workflow.receipt.v1' and view['status'] == 'complete'
        assert view['result'] == expected['result']
        assert len(view['result']) == 2 and (results / 'cli-workflow.rfcxs').is_file()
        display = subprocess.run(workflow.args[:-2], capture_output=True, text=True, timeout=300)
        assert display.returncode == 0, display.stderr
        assert json.loads(display.stdout)['result'] == expected['result']
        legacy = copy.deepcopy(expected['workflow'])
        legacy['entry'] = dict(id='legacy-summary', type='query',
                               query=dict(operation='coverage', field='all.A'))
        legacy_path = results / 'legacy-workflow.json'
        legacy_path.write_text(json.dumps(legacy))
        command = [str(binary), 'workflow', '--workspace',
                   str(case / 'workspace.moved/workspace.rfcxs'), '--assay', 'measurements',
                   '--request', str(legacy_path)]
        display = subprocess.run(command, capture_output=True, text=True, timeout=300)
        assert display.returncode == 0, display.stderr
        assert json.loads(display.stdout)['result']['result']['reported'] > 0
        target = results / 'unconverted-workflow.rfcxs'
        refused = subprocess.run([*command, '--out', str(target)], capture_output=True, text=True)
        assert refused.returncode != 0 and 'native result publisher' in refused.stderr
        assert not target.exists()

    verify(structural, queries, 'structural')
    sections = root / 'sections'
    query = demo.authored(sections)
    demo.run(binary, 'import-assays', '--spec', sections / 'inputs/import.json', '--out', sections / 'workspace')
    paths = list((sections / 'workspace').rglob('study.json'))
    assert len(paths) == 1
    declaration = json.loads(paths[0].read_text())
    cells, occurrences = declaration['cells'], declaration['occurrences']
    # Independent stalk axes have equal dimensions but distinct cell owners.
    stalks = [dict(grade=c['grade'], key=c['key'], axes=['all.A'] if c['grade'] == 1 else ['scalar'])
              for c in cells]
    maps = [dict(grade=o['grade'], key=o['key'], matrix=[['1']]) for o in occurrences]
    # The complete retained chain supplies all two-step cancellations.
    local = dict(operation='local_sections', context='all', construction='cellular_sheaf',
                 stalks=stalks, maps=maps, metric='identity',
                 sections=[dict(id='source', grade=1, fields=['all.A'], quantity='declared_amplitude',
                                units='arbitrary_units')])
    source_meaning = next(f for f in declaration['fields'] if f['id'] == 'all.A')
    local['sections'][0].update(quantity=source_meaning['quantity'], units=source_meaning['units'])
    declared = copy.deepcopy(local)
    declared['sections'] = [dict(id='declared', grade=1, values=['0', '2', '4', '0'],
                                quantity='authored_test', units='test')]
    fibers = copy.deepcopy(declared)
    fibers['construction'] = 'boundary_fibers'
    heterogeneous = copy.deepcopy(local)
    heterogeneous['stalks'] = [dict(grade=c['grade'], key=c['key'],
                                  axes=['scalar'] if c['grade'] == 0 else
                                  ['all.A', 'all.B'] if c['grade'] == 1 else
                                  ['terminal', 'unused', 'diagnostic']) for c in cells]
    heterogeneous['maps'] = [dict(grade=o['grade'], key=o['key'],
                                 matrix=[['1'], ['0']] if o['grade'] == 1 else
                                 [['0', '0'], ['0', '0'], ['0', '0']]) for o in occurrences]
    heterogeneous['sections'][0]['fields'] = ['all.A', 'all.B']
    downward = copy.deepcopy(heterogeneous)
    downward['construction'] = 'boundary_fibers'
    downward['maps'] = [dict(grade=o['grade'], key=o['key'],
                            matrix=[['1', '0']] if o['grade'] == 1 else
                            [['0', '1', '0'], ['0', '0', '1']]) for o in occurrences]
    downward['sections'] = [dict(id='upper', grade=2, values=['2', '3', '4'],
                                quantity='authored_test', units='test')]
    verify(sections, {'source': local, 'declared': declared, 'fibers': fibers,
                      'heterogeneous': heterogeneous, 'downward': downward}, 'sections')
    print(json.dumps(dict(status='PASS', cases=receipts)))


if __name__ == '__main__':
    main(*[Path(v).resolve() for v in sys.argv[1:]])
