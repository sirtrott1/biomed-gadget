"""Connect independently verified producers to complete typed answer artifacts."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess


MEANINGS = {
    'quality': 'read_quality',
    'alignment_evidence': 'alignment_records',
    'annotation_features': 'annotation_features',
    'feature_counts': 'feature_allocations',
    'allocation_bounds': 'allocation_ranges',
    'expression_allocations': 'expression_allocation_ranges',
    'interval_coverage': 'interval_coverage',
    'chromatin_contrast': 'chromatin_contrast',
    'controlled_coverage': 'controlled_coverage',
    'spatial_relations': 'spatial_relations',
    'field_transfer': 'field_transfer',
    'genotypes': 'genotype_observations',
    'motifs': 'motif_matches',
    'exact_sequence_alignment': 'global_edit_routes',
    'reference_sequence_matches': 'reference_placements',
    'indexed_sequence_alignment': 'indexed_edit_placements',
    'sequence_splice_discovery': 'splice_chains',
    'sequence_assembly_graph': 'kmer_structure',
    'read_overlap_reconstruction': 'read_cover_reconstructions',
    'constrained_allocations': 'constrained_allocations',
    'path_comparison': 'realized_path_comparison',
}


def qualify(binary, root, workspace, cases):
    """cases maps a fixture ID to (actual request, independently checked result)."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    checks = 0

    def check(value):
        nonlocal checks
        assert value
        checks += 1

    def run(*args):
        p = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True, timeout=180)
        if p.returncode:
            raise AssertionError((args, p.stderr))
        check(True)
        return json.loads(p.stdout)

    def save(name, value):
        path = root / name
        path.write_text(json.dumps(value, indent=2) + '\n')
        return path

    definitions = run('pipelines')['recipes']
    by_operation = {r['operation']: r for r in definitions}
    recipes, requirements, expected, actual_requests = [], [], {}, {}
    for name, (request, result) in cases.items():
        native = request['operation'] == 'native_query'
        operation = request['query']['operation'] if native else request['operation']
        if operation not in by_operation or by_operation[operation]['method'] not in MEANINGS:
            continue
        label, name = name, 'case-' + str(len(recipes))
        method = by_operation[operation]['method']
        params = {k: v for k, v in (request['query'] if native else request).items()
                  if k not in {'schema', 'operation'}}
        goal = dict(method=method, parameters=params)
        if native:
            goal['assay'] = request['assay']
            if 'partitions' in request:
                goal['partitions'] = request['partitions']
            if 'concurrent_jobs' in request:
                goal['partition_jobs'] = request['concurrent_jobs']
        recipes.append(dict(id=name, purpose='Retained output for independently verified ' + label, goal=goal))
        requirements.append(dict(id=name, observable=MEANINGS[method], bindings=[name]))
        expected[name] = result
        actual_requests[name] = request
    check(bool(recipes))
    question = 'Execute the declared methods and retain their complete scientific results.'
    intent = dict(schema='gadget.pipeline.intent.v1', profile='custom', recipes=recipes,
                  answer_contract=dict(schema='gadget.pipeline.answer.contract.v1',
                                       question=question, requirements=requirements))
    intent_path = save('intent.json', intent)
    catalog = run('recipe-catalog', '--workspace', workspace, '--intent', intent_path)
    check(catalog['scientific_admission']['status'] == 'satisfied')
    check(catalog['scientific_admission']['native']['assignment'] == [1] * len(recipes))
    for r in catalog['scientific_admission']['requirements']:
        check(r['eligible_bindings'] == [r['id']])
    proposal = dict(schema='gadget.pipeline.recipe.proposal.v1', source_digest=catalog['source_digest'],
                    intent_sha256=catalog['intent_sha256'], action='recipes', goals=list(expected))

    def admit(proposed, declaration=intent_path, out=None):
        args = ['recipe-proposal', '--workspace', workspace, '--intent', declaration,
                '--request', save('proposal.json', proposed), '--prompt', question]
        if out:
            args += ['--out', root / out]
        return run(*args)

    admitted = admit(proposal, out='complete')
    index = admitted['execution']
    directory = root / 'complete'
    check(index['answers_sha256'] == 'sha256:' + hashlib.sha256((directory / index['answers']).read_bytes()).hexdigest())
    answers = run('record-view','--input',directory/index['answers'],'--member','result')
    check(len(answers['answers']) == len(requirements))
    checked_files = 0
    for answer in answers['answers']:
        name = answer['id']
        check(answer['observable'] == next(r['observable'] for r in requirements if r['id'] == name))
        check(len(answer['artifacts']) == 1 and answer['artifacts'][0]['binding'] == name)
        raw = expected[name]
        # Direct assay queries wrap native partition certificates in answers.
        if isinstance(raw, dict) and raw.get('schema') == 'gadget.assay.result.v1':
            raw = raw['result']
        expected_results = [r['certificate']['result'] for r in raw['answers']] if actual_requests[name]['operation'] == 'native_query' else [raw]
        refs = answer['artifacts'][0]['results']
        check(len(refs) == len(expected_results))
        for ref, expected_result in zip(refs, expected_results):
            check(not Path(ref['path']).is_absolute() and '..' not in Path(ref['path']).parts)
            data = (directory / ref['path']).read_bytes()
            check(ref['sha256'] == 'sha256:' + hashlib.sha256(data).hexdigest())
            certificate = run('record-view','--input',directory/ref['path'])
            if certificate.get('schema')=='gadget.operation.native.v1':certificate=certificate['result']
            check(ref['pointer'] == '/result')
            check(ref['source_digest'] == certificate['source_digest'])
            check(ref['query'] == certificate['query'])
            check(ref['certificate_sha256'] == certificate['certificate_sha256'])
            check(certificate['query'] == (actual_requests[name]['query'] if actual_requests[name]['operation'] == 'native_query' else actual_requests[name]))
            check(certificate['result'] == expected_result)
            checked_files += 1
    if len(recipes) > 1:
        missing = copy.deepcopy(proposal)
        missing['goals'].pop(0)
        rejected = admit(missing, out='missing')
        check(rejected['status'] == 'incompatible' and not (root / 'missing').exists())
        check(rejected['scientific_admission']['native']['unsat_core']['subset_minimal'])
    # Purpose text cannot turn enumerated ranges, counts or routes into p-values.
    invalid = copy.deepcopy(intent)
    invalid['recipes'][0]['purpose'] = 'Calibrated population p-values and statistical confidence'
    invalid['answer_contract']['requirements'][0]['observable'] = 'population_p_value'
    invalid_path = save('mislabelled-intent.json', invalid)
    unavailable = run('recipe-catalog', '--workspace', workspace, '--intent', invalid_path)
    changed = {**proposal, 'intent_sha256': unavailable['intent_sha256']}
    rejected = admit(changed, declaration=invalid_path, out='mislabelled')
    check(rejected['status'] == 'unsupported' and not (root / 'mislabelled').exists())
    check(rejected['scientific_admission']['native']['unsat_core']['subset_minimal'])
    methods = sorted({r['goal']['method'] for r in recipes})
    save('verification.json', dict(status='PASS', checks=checks, methods=methods,
                                  independently_verified_cases=len(recipes), complete_certificates=checked_files))
    return checks
