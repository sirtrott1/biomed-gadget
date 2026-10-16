#!/usr/bin/env python3
"""Compare complete producer reported molecular formula counts on a declared interface."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


ELEMENTS = set('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split())


def digest(path):
    with path.open('rb') as stream:
        return 'sha256:' + hashlib.file_digest(stream, 'sha256').hexdigest()


def pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('duplicate original JSON key: ' + key)
        result[key] = value
    return result


def formula_counts(value):
    if not isinstance(value, str) or not value:
        raise ValueError('reported molecular formula is missing')
    counts, tokens, position = Counter(), [], 0
    for match in re.finditer(r'([A-Z][a-z]?)([0-9]*)', value):
        symbol, number = match.groups()
        if match.start() != position or symbol not in ELEMENTS:
            raise ValueError('formula profile requires ordinary element symbols and integer counts')
        if number and (number.startswith('0') or int(number) < 1):
            raise ValueError('formula profile requires positive element counts')
        count = int(number) if number else 1
        counts[symbol] += count
        tokens.append(dict(element=symbol, count=str(count), start=match.start(), end=match.end(), text=match.group()))
        position = match.end()
    if position != len(value) or not tokens:
        raise ValueError('formula profile does not admit charges, isotopes, brackets or mixture notation')
    return dict(counts), tokens


def prepare(snapshot, chembl_response, pubchem_response):
    snapshot = snapshot.resolve(strict=True)
    manifest = json.loads((snapshot/'snapshot.json').read_text(), object_pairs_hook=pairs)
    if manifest.get('schema') != 'gadget.database.snapshot.v1':
        raise ValueError('a frozen database snapshot is required')
    sources = []
    values = []
    selected = {}
    for name, provider, resource in [(chembl_response,'chembl','molecule'), (pubchem_response,'pubchem','compound')]:
        candidates = [r for r in manifest['records'] if r['request']['id'] == name]
        if len(candidates) != 1:
            raise ValueError('response identifier is missing or ambiguous: ' + name)
        row = candidates[0]
        if (row['request']['provider'],row['request']['resource']) != (provider,resource):
            raise ValueError('response provider or resource differs from selected formula profile')
        for path_key,hash_key in [('path','sha256'),('headers_path','headers_sha256')]:
            path = (snapshot/row[path_key]).resolve(strict=True)
            if snapshot not in path.parents or not path.is_file() or digest(path) != row[hash_key]:
                raise ValueError('original response or header owner differs from frozen declaration')
            sources.append(dict(id=provider+'-'+path_key,original=path,sha256=row[hash_key],
                                source_url=row['url'],license=row.get('data_terms','Terms not supplied by this source declaration')))
        raw = json.loads((snapshot/row['path']).read_bytes(),object_pairs_hook=pairs)
        if provider == 'chembl':
            identity = 'chembl:' + raw['molecule_chembl_id']
            if raw['molecule_chembl_id'] != row['request']['identifier']:
                raise ValueError('original molecule identity differs from selected request')
            key = raw['molecule_structures']['standard_inchi_key']
            formula = raw['molecule_properties']['full_molformula']
            pointer = '/molecule_properties/full_molformula'
        else:
            properties = raw['PropertyTable']['Properties']
            if len(properties) != 1 or str(properties[0]['CID']) != row['request']['identifier']:
                raise ValueError('formula profile requires the complete single selected PubChem compound')
            identity = 'pubchem:cid:' + str(properties[0]['CID'])
            key, formula = properties[0]['InChIKey'], properties[0]['MolecularFormula']
            pointer = '/PropertyTable/Properties/0/MolecularFormula'
        if not isinstance(key,str) or not key:
            raise ValueError('reported InChIKey is missing')
        counts,tokens = formula_counts(formula)
        selected[provider] = dict(identity=identity,reported_inchikey=key,formula=formula,
                                  formula_pointer=pointer,counts=counts,tokens=tokens,receipt=row)
        values.append(key)
    if values[0] != values[1]:
        raise ValueError('producer reported InChIKeys differ; no common chemical interface admitted')
    return selected,sources


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,required=True)
    parser.add_argument('--chembl-response',required=True)
    parser.add_argument('--pubchem-response',required=True)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--collection-tool',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    selected,sources = prepare(args.snapshot,args.chembl_response,args.pubchem_response)
    binary,collector = args.binary.resolve(strict=True),args.collection_tool.resolve(strict=True)
    output = args.out.resolve()
    output.mkdir(parents=True,exist_ok=False)
    inputs = output/'inputs'
    inputs.mkdir()
    provenance = []
    for source in sources:
        path = inputs/(source['id']+'.original')
        shutil.copyfile(source['original'],path)
        provenance.append({k:v for k,v in source.items() if k not in ['original'] } | dict(path=path.name))
    shutil.copyfile(args.snapshot.resolve()/'snapshot.json',inputs/'snapshot.original.json')
    provenance.append(dict(id='source-snapshot',path='snapshot.original.json',sha256=digest(inputs/'snapshot.original.json'),
                           source_url='Declared frozen database snapshot metadata',license='Original provider terms retained in snapshot records'))
    save(inputs/'correspondence.json',selected)
    provenance.append(dict(id='formula-correspondence',path='correspondence.json',sha256=digest(inputs/'correspondence.json'),
                           source_url='Declared correspondence to retained original database formula fields',
                           license='GADGET application declaration'))
    axes = sorted(set().union(*(set(v['counts'])for v in selected.values())))
    (inputs/'elements.tsv').write_text('id\n'+'\n'.join(axes)+'\n')
    (inputs/'observations.tsv').write_text('id\nreported_formula\n')
    design = dict(question='Do the two frozen producer reported formulas agree on every declared element count?',
        population='The selected producer compound records',selection='Every element token in both complete formulas',
        limitations=['Common interface is conditional on equal producer reported InChIKey; no independent structure or stereochemistry validation',
                     'Unbracketed ordinary element counts only; mixtures, charges and isotope notation refuse',
                     'Identity costs measure source count discrepancy; no assay calibration, clinical or causal inference'])
    spec = dict(schema='gadget.assays.import.v1',id='reported-chemical-formula-counts',design=design,assays=[],resources=[],links=[])
    views = []
    for index,(provider,record)in enumerate(selected.items()):
        table = inputs/(provider+'.tsv')
        table.write_text('id\treported_formula\n'+''.join(element+'\t'+str(record['counts'].get(element,0))+'\n'for element in axes))
        file_record = lambda path:dict(path=path.name,sha256=digest(path))
        spec['assays'].append(dict(id=provider,kind='paired_measurement',format='dense_tsv',source=file_record(table),
            features=file_record(inputs/'elements.tsv'),observations=file_record(inputs/'observations.tsv'),
            source_namespace=str(977300+index),carrier='features',quantity='chemical_element_count',units='atoms_per_molecule',
            channel='observation',method='Exact counts from every original formula token; missing formulas refuse',
            arithmetic='exact_rational',implicit_entries='unobserved',missing_tokens=['NA']))
        views.append(dict(id=provider,assay=provider,partition=0,fields=[dict(section='reported_formula',field='all.reported_formula')],
            metric=dict(diagonal=['1']*len(axes),domain_quantity='chemical_element_count',domain_units='atoms_per_molecule',
                quadrance_quantity='declared_count_discrepancy',quadrance_units='atoms_per_molecule_squared',
                method='Identity cost on all original declared element counts',evidence='Declared comparison cost, not experimental error variances')))
    save(inputs/'import.json',spec)
    subprocess.run([str(binary),'import-assays','--spec',str(inputs/'import.json'),'--out',str(output/'imported')],check=True,stdout=(output/'import.log').open('w'))
    control = dict(schema='gadget.biomedical.collection.request.v1',id=spec['id'],design=design,
                   assays=list(selected),method='Full original provider bytes, lexical formula tokens and declared element count correspondence',sources=provenance)
    save(inputs/'collection.json',control)
    subprocess.run([str(collector),str(output/'imported/workspace.json'),str(inputs/'collection.json'),str(output/'workspace.rfcxs')],
                   cwd=inputs,check=True,stdout=(output/'collection.log').open('w'))
    participants = []
    for index,provider in enumerate(selected):
        participants.append(dict(key=str(index+1),view=provider,exists=True,orientation='positive'if index==0 else'negative',
            magnitude='1',map=dict(domain_quantity='chemical_element_count',domain_units='atoms_per_molecule',domain_channel='observation',
                codomain_quantity='chemical_element_count',codomain_units='atoms_per_molecule',codomain_channel='observation',
                method='Declared element identity under the retained producer chemical identifier assertions',evidence='See retained original formula and InChIKey byte owners',
                entries=[dict(observable_id=e,source_id=e,coefficient='1')for e in axes])))
    query = dict(schema='gadget.assays.query.v1',operation='assay_sections',views=views,sections=['reported_formula'],
        interfaces=[dict(id='reported-chemical-identity',axes=axes,quantity='chemical_element_count',units='atoms_per_molecule',
            channel='observation',method='Compare every element count across independently admitted source studies',
            evidence='Producer identity assertions and full original formula correspondence retained',participants=participants)],
        metric=dict(quantity='declared_count_discrepancy',units='atoms_per_molecule_squared',
                    method='Direct sum of complete identity forms on reported element counts',evidence='Declared source comparison costs, not precision estimates'),
        support_policy='require_all_nonzero_occurrences',workers=1,complete_families=True)
    save(output/'query.json',query)
    with (output/'execution.log').open('w') as stream:
        subprocess.run([str(binary),'assay-query','--workspace',str(output/'workspace.rfcxs'),'--request',str(output/'query.json'),
                        '--out-record',str(output/'result.rfcxs')],stdout=stream,check=True)
    with (output/'result.json').open('w') as stream:
        subprocess.run([str(binary),'record-view','--input',str(output/'result.rfcxs'),'--member','result'],stdout=stream,check=True)
    print(json.dumps(dict(workspace=str(output/'workspace.rfcxs'),native_result=str(output/'result.rfcxs'),elements=axes)))


if __name__ == '__main__':
    main()
