"""Original formula tokens, supported zero, complete native fields and portable replay."""
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def digest(path):
    return 'sha256:' + hashlib.sha256(path.read_bytes()).hexdigest()


def main(binary, collector, consumer, parent):
    parent.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='formula inputs # ', dir=parent))
    tool = Path(__file__).resolve().parents[1] / 'tools/prepare_chemical_formula_sections.py'
    module_spec = importlib.util.spec_from_file_location('formula_adapter', tool)
    adapter = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(adapter)
    checks = 0

    def save(path, value):
        path.write_text(json.dumps(value, indent=2) + '\n')
        return path

    def fixture(path, left='CH4', right='CH4O', key='AUTHORED-IDENTIFIER'):
        path.mkdir()
        bodies = [dict(molecule_chembl_id='CHEMBLTEST', molecule_structures=dict(standard_inchi_key=key),
                       molecule_properties=dict(full_molformula=left)),
                  dict(PropertyTable=dict(Properties=[dict(CID=123, InChIKey=key, MolecularFormula=right)]))]
        rows = []
        for provider, resource, identifier, body in zip(['chembl', 'pubchem'], ['molecule', 'compound'],
                                                      ['CHEMBLTEST', '123'], bodies):
            raw = save(path / (provider + '.json'), body)
            headers = path / (provider + '.headers')
            headers.write_bytes(b'Offline authored contract input\r\nX-Retained: complete\r\n\r\n')
            rows.append(dict(request=dict(id=provider, provider=provider, resource=resource, identifier=identifier),
                path=raw.name, sha256=digest(raw), headers_path=headers.name, headers_sha256=digest(headers),
                url='offline:authored-formula', data_terms='Authored contract input'))
        return save(path/'snapshot.json', dict(schema='gadget.database.snapshot.v1', records=rows))

    def execute(snapshot, destination, fail=False):
        result = subprocess.run([sys.executable, '-B', str(tool), '--snapshot', str(snapshot),
            '--chembl-response', 'chembl', '--pubchem-response', 'pubchem', '--binary', str(binary),
            '--collection-tool', str(collector), '--out', str(destination)], text=True, capture_output=True,
            cwd=root, timeout=180)
        assert bool(result.returncode) == fail, (result.stdout, result.stderr)
        if fail:
            assert not destination.exists(), result.stderr
        return result

    manifest_path = fixture(root/'original')
    snapshot = manifest_path.parent
    original_hashes = {p.name:digest(p) for p in snapshot.iterdir()}
    output = root/'native result # full paths'
    execute(snapshot, output)
    declaration = json.loads((output/'inputs/collection.json').read_text())
    assert declaration['assays'] == ['chembl', 'pubchem']
    for original in json.loads(manifest_path.read_text())['records']:
        for path_key in ['path', 'headers_path']:
            retained = output/'inputs'/(original['request']['provider']+'-'+path_key+'.original')
            assert retained.read_bytes() == (snapshot/original[path_key]).read_bytes()
            checks += 1
    expected = save(root/'expected.json', dict(axes=['C','H','O'], left=['1','4','0'], right=['1','4','1'],
        sources=declaration['sources']))
    view = json.loads((output/'result.json').read_text())['result']
    assert view['sections'][0]['complete'] and Q(view['sections'][0]['residual_quadrance']) == Q(1,2)
    for row in view['sections'][0]['fields']:
        assert row['ids'] == ['C','H','O']
        assert list(map(Q,row['compatible'])) == [Q(1),Q(4),Q(1,2)]
        assert list(map(Q,row['residual'])) == [Q(0),Q(0), Q(-1,2) if row['view']=='chembl' else Q(1,2)]
        checks += 3
    assert view['prepared_groups'][0]['compatible_dimension'] == 3
    for cold in [False, True]:
        if cold:
            snapshot.rename(root/'original.offline')
            (output/'inputs').rename(output/'inputs.offline')
            (output/'imported').rename(output/'imported.offline')
        receipt = subprocess.run([str(consumer),str(output),str(expected)], capture_output=True,text=True,timeout=180)
        assert receipt.returncode == 0, (receipt.stdout,receipt.stderr)
        checks += json.loads(receipt.stdout)['native_checks']
    snapshot = root/'original.offline'
    assert {p.name:digest(p) for p in snapshot.iterdir()} == original_hashes
    checks += len(original_hashes)
    for formula in ['CH4','H2CH2','Og2','NaCl','C1000000H2']:
        counts,tokens = adapter.formula_counts(formula)
        assert ''.join(t['text'] for t in tokens) == formula and all(n>0 for n in counts.values())
        checks += 1
    assert adapter.formula_counts('H2CH2')[0] == {'H':4,'C':1}
    for formula in ['',None,'Xx2','C0','C01','[13C]H4','CH4+','CH4.O','CH4-2','2H','C1.2']:
        try:
            adapter.formula_counts(formula)
        except (ValueError,TypeError):
            checks += 1
        else:
            raise AssertionError('unsupported original formula admitted: '+str(formula))
    for index, kind in enumerate(['missing','key','hash','escape','duplicate-id','duplicate-key','grammar']):
        broken = root/('refusal-'+kind)
        shutil.copytree(snapshot,broken)
        manifest = json.loads((broken/'snapshot.json').read_text())
        if kind in ['missing','key','grammar','duplicate-key']:
            raw = json.loads((broken/'chembl.json').read_text())
            if kind == 'missing': raw['molecule_properties']['full_molformula'] = None
            if kind == 'key': raw['molecule_structures']['standard_inchi_key'] = 'DIFFERENT'
            if kind == 'grammar': raw['molecule_properties']['full_molformula'] = 'CH4+'
            save(broken/'chembl.json',raw)
            if kind == 'duplicate-key':
                (broken/'chembl.json').write_text('{"molecule_chembl_id":"CHEMBLTEST","molecule_chembl_id":"CHEMBLTEST"}')
            manifest['records'][0]['sha256'] = digest(broken/'chembl.json')
        elif kind == 'hash': manifest['records'][0]['sha256'] = 'sha256:'+'0'*64
        elif kind == 'escape':
            manifest['records'][0]['path'] = '../expected.json'
            manifest['records'][0]['sha256'] = digest(expected)
        else: manifest['records'].append(copy.deepcopy(manifest['records'][0]))
        save(broken/'snapshot.json',manifest)
        execute(broken,root/('refused-'+str(index)),fail=True)
        checks += 1
    prior = {str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}
    existing = subprocess.run([sys.executable,'-B',str(tool),'--snapshot',str(snapshot),
        '--chembl-response','chembl','--pubchem-response','pubchem','--binary',str(binary),
        '--collection-tool',str(collector),'--out',str(output)],capture_output=True,text=True)
    assert existing.returncode and prior == {str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}
    checks += 1
    print(json.dumps(dict(status='PASS',checks=checks,scope='Authored two source original formula counts, supported zero, exact defects, complete families, full paths, cold native replay and source refusals')))


if __name__ == '__main__':
    main(*(Path(p).resolve() for p in sys.argv[1:]))
