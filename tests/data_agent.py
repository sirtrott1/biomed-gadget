"""Optional native GADGET providers under the existing Agent/Hive boundaries."""
import os
import subprocess

import pytest
from agent.mcp_tools import Context

from agent import RelationalSession
from gadget_data.agent import GadgetSession
from rcdb import MemoryStore, VersionConflictError


@pytest.fixture
def native(tmp_path):
    adapter = pytest.importorskip("gadget_rcdb")
    binary = os.environ.get("GADGET_TEST_BINARY")
    if not binary:
        pytest.skip("GADGET_TEST_BINARY selects the optional compiled application")
    from rcdb.header import DEFAULT_RECORD_CODECS

    from rcdb import RecordCodec
    budget = RecordCodec.__dataclass_fields__["max_bytes"].default
    codec = adapter.install_codec(gadget=binary, max_bytes=budget)
    store = MemoryStore(record_codecs=(*DEFAULT_RECORD_CODECS, codec))
    context = Context(workspace="gadget-tests", identity="alice", auth_enabled=True)
    session = RelationalSession(store, context=context)
    service = GadgetSession(session, gadget=binary, max_bytes=budget)
    packets = []
    for index in range(2):
        source = tmp_path / (str(index) + ".json")
        target = source.with_suffix(".rfcxs")
        source.write_text('{"original_observation":' + str(index) + '}')
        subprocess.run([binary, "native-record", "--input", str(source), "--out", str(target)],
                       check=True, capture_output=True, timeout=30)
        packets.append((target, target.read_bytes()))
    return service, store, packets, binary, budget


def test_optional_setup_rejects_invalid_capabilities_and_requires_a_session():
    session = RelationalSession(MemoryStore())
    with pytest.raises(TypeError, match="RelationalSession"):
        GadgetSession(object(), gadget="unused", max_bytes=1)
    with pytest.raises(ValueError, match="capabilities"):
        GadgetSession(session, gadget="unused", max_bytes=1, capabilities=("shell",))
    with pytest.raises(ValueError, match="byte budget"):
        GadgetSession(session, gadget="unused", max_bytes=True)
    with pytest.raises(ValueError, match="selector"):
        GadgetSession(session, gadget="unused", max_bytes=1, capabilities=("select",))


def test_operation_inspection_keeps_native_execution_ungranted(native):
    service, _, _, _, _ = native
    report = service.operation_capabilities(("GADGET_DATA_FIELD", "DEMAND_RUN", "C_EVALUATE"))
    assert report["declared_permissions"] == ["identity", "read", "records"]
    assert all(row["platform_qualification"] == "unverified" for row in report["operations"])
    assert all(row["biomedical_qualification"] == "unverified" for row in report["operations"])
    assert "gadget_execute" not in report["declared_permissions"]


def test_gadget_optional_dependency_is_loaded_only_when_requested(monkeypatch):
    import builtins
    previous = builtins.__import__
    def missing(name, *args, **kwargs):
        if name == "gadget_rcdb":
            raise ImportError("optional application absent")
        return previous(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", missing)
    session = RelationalSession(MemoryStore())
    with pytest.raises(ImportError, match="biomed-gadget-data\\[agent\\]"):
        GadgetSession(session, gadget="unused", max_bytes=1)


def test_gadget_hive_worker_preserves_original_bytes_versions_and_scope(native, tmp_path, monkeypatch):
    from agent.hive import Hive
    from agent.hive_workflow import WorkflowStep

    from rcql import Result
    service, store, packets, _, _ = native
    monkeypatch.setenv("REXGRAPH_CONFIG_DIR", str(tmp_path / "config"))
    service.publish_file("publication", packets[0][0], expected_version=0)
    service.publish("publication", packets[1][1], expected_version=1)
    hive = Hive("gadget-tests")
    try:
        hive.add_worker("gadget", service, capability="analyze")
        flow = hive.run_workflow([WorkflowStep("historical", "gadget", {
            "query": 'FROM $db RETURN GADGET_PUBLICATION_BYTES("publication",1)',
            "save_as": "result", "expected_version": 0})])
        assert flow.ok, flow.steps
        result = flow.outputs["historical"]
        assert result.values == (packets[0][1],)
        assert Result.from_bytes(result.to_bytes()).values == result.values
        assert service.session.read("result", 1).value.to_bytes() == result.to_bytes()
        assert store.get_record("publication").meta["stored_by"] == "alice"
        assert store.get_record("result").meta["workspace"] == "gadget-tests"
        assert result.execution[0]["methods"][0]["record_version"] == 1
    finally:
        hive.stop_all()
    foreign = RelationalSession(store, context=Context(workspace="other", identity="bob", auth_enabled=True))
    with pytest.raises(KeyError):
        GadgetSession(foreign, gadget=service.gadget, max_bytes=service.max_bytes).query(
            'FROM $db RETURN GADGET_PUBLICATION_BYTES("publication",1)')


def test_gadget_execution_requires_explicit_capability_before_dispatch(native, monkeypatch):
    from rcql import call, query, source
    service, _, _, binary, budget = native
    enabled = GadgetSession(service.session, gadget=binary, max_bytes=budget, capabilities=("execute",))
    def must_not_execute(*args, **kwargs):
        raise AssertionError("ungranted native computation executed")
    monkeypatch.setattr(enabled._adapter, "execute_publication", must_not_execute)
    declaration = query(source("db"), call("GADGET_ASSAY_QUERY", "workspace", 1, "request", 1))
    with pytest.raises(PermissionError):
        service.query(declaration)
    with pytest.raises(PermissionError):
        service.session.query(declaration)


def test_stored_read_program_keeps_historical_publication_and_complete_steps(native):
    from rcql import Program, ProgramStep, parse
    service, _, packets, _, _ = native
    service.publish("publication", packets[0][1], expected_version=0)
    service.publish("publication", packets[1][1], expected_version=1)
    program = Program("gadget_history", (
        ProgramStep("original", parse('FROM $db RETURN GADGET_PUBLICATION_BYTES("publication",1)')),
        ProgramStep("current", parse('FROM $db RETURN GADGET_PUBLICATION_BYTES("publication",2)'))))
    service.session.publish_program("reader", program, expected_version=0)
    result = service.run_program("reader", version=1)
    assert result.values == (packets[1][1],)
    assert result.native_plan["step_results"][0]["values"] == (packets[0][1],)
    assert len(result.native_plan["step_results"]) == 2


@pytest.fixture
def native_data(native, tmp_path):
    import numpy as np
    from gadget_data import DataStore, DataWriter
    service, _, _, _, _ = native
    with DataWriter(tmp_path / "counts", name="counts") as writer:
        writer.axis("rows", ("a", "b"), offset=0)
        writer.axis("columns", ("x", "y"), offset=0)
        writer.field("counts", axes=("rows", "columns"), shape=(2, 2), quantity="count", units="count")
        writer.dense("counts", np.array(((0, 2), (3, 4)), dtype=np.int32),
                     np.ones((2, 2), dtype=bool), offset=(0, 0))
        writer.finish()
    with DataStore(tmp_path / "counts") as store:
        reference = next(iter(store.manifest["fields"]["counts"]["blocks"]))
        owner = service.bind_native("counts", store, reference)
        yield service, store, reference, owner


def test_native_dataset_query_and_stored_program_use_the_same_actual_owner(native_data):
    from rcql import Program, ProgramInput, ProgramStep, Result, parse
    from unittest.mock import patch
    service, store, _, owner = native_data
    with patch.object(store, "read", side_effect=AssertionError("native source decoded again")):
        field, = service.query("FROM $counts RETURN GADGET_DATA_FIELD([0])").values
        assert field.source.source is owner.value
        program = Program("read_counts", (ProgramStep("selected",
            parse("FROM $counts RETURN GADGET_DATA_FIELD($columns)")),), (ProgramInput("columns"),))
        service.session.publish_program("reader", program, expected_version=0)
        service.session.save_result("columns", Result(((1, 0),)), expected_version=0)
        output = service.run_program("reader", version=1,
            record_parameters={"columns": {"record_id": "columns", "version": 1, "output": 0}},
            save_as="readout")
        assert output.values[0].source.source is owner.value
        assert output.values[0].values.tolist() == [[2, 0], [4, 3]]
        assert output.native_plan["record_parameters"]["columns"]["version"] == 1
        restored = service.session.read("readout", 1).value
        assert restored.values[0].coefficient_digest == output.values[0].coefficient_digest


def test_native_program_capture_replays_with_fresh_explicit_dataset_bindings(native_data, tmp_path):
    from gadget_data import DataStore
    from rcql import Program, ProgramStep, parse
    service, store, reference, owner = native_data
    program = Program("read_counts", (ProgramStep("selected",
        parse("FROM $counts RETURN GADGET_DATA_FIELD([0,1])")),))
    service.session.publish_program("reader", program, expected_version=0)
    service.capture_program("reader", tmp_path / "execution", version=1)
    store.close()
    fresh = GadgetSession(service.session, gadget=service.gadget, max_bytes=service.max_bytes)
    with DataStore(tmp_path / "counts") as reopened:
        rebound = fresh.bind_native("counts", reopened, reference)
        result = fresh.replay_program(tmp_path / "execution")
        assert result.values[0].source.source is rebound.value
        assert result.values[0].source.source is not owner.value
        assert result.values[0].values.tolist() == [[0, 2], [3, 4]]
    with pytest.raises(ValueError, match="store is closed"):
        fresh.replay_program(tmp_path / "execution")


def test_native_binding_revocation_during_query_prevents_saved_results(native_data):
    from rcql import Executor
    from unittest.mock import patch
    service, store, _, _ = native_data
    execute = Executor.execute
    def close_after_query(executor, query):
        result = execute(executor, query)
        store.close()
        return result
    with patch.object(Executor, "execute", close_after_query):
        with pytest.raises(ValueError, match="store is closed"):
            service.query("FROM $counts RETURN GADGET_DATA_FIELD([0])", save_as="invalid")
    with pytest.raises(KeyError):
        service.session.read("invalid", 1)


def test_worker_arguments_cannot_bind_datasets_or_raise_source_permissions(native_data):
    service, store, reference, _ = native_data
    for key in ("bindings", "sources", "store"):
        with pytest.raises(ValueError, match="worker"):
            service({"query": "FROM $counts RETURN GADGET_DATA_FIELD([0])", key: {"counts": reference}})
    with pytest.raises(ValueError, match="read and identity"):
        service.bind_native("other", store, reference, permissions=("train",))
    service.unbind_native("counts")
    with pytest.raises((KeyError, ValueError)):
        service.query("FROM $counts RETURN GADGET_DATA_FIELD([0])")


def test_gadget_corrupt_publication_and_stale_writes_preserve_prior_state(native):
    service, _, packets, _, _ = native
    service.publish("publication", packets[0][1], expected_version=0)
    with pytest.raises(ValueError):
        service.publish("publication", packets[0][1][:-1], expected_version=1)
    with pytest.raises(VersionConflictError):
        service.publish("publication", packets[1][1], expected_version=0)
    result = service.query('FROM $db RETURN GADGET_PUBLICATION_BYTES("publication",1)')
    assert result.values == (packets[0][1],)
    assert len(service.session.store.history("publication")) == 1
    with pytest.raises(ValueError, match="worker"):
        service({"query": 'FROM $db RETURN 1', "gadget": "another executable"})


def test_native_filesystem_actions_stay_outside_pure_programs(native):
    from rcql import Program, ProgramStep, call, query, source
    service, _, _, binary, budget = native
    GadgetSession(service.session, gadget=binary, max_bytes=budget, capabilities=("execute",))
    with pytest.raises(ValueError, match="read calculations"):
        Program("native_action", (ProgramStep("compute", query(source("db"),
            call("GADGET_ASSAY_QUERY", "workspace", 1, "request", 1))),))


def test_missing_agent_extra_is_reported_on_construction(monkeypatch):
    import builtins
    previous = builtins.__import__
    def missing(name, *args, **kwargs):
        if name == 'agent.relational':
            raise ImportError('optional Agent absent')
        return previous(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', missing)
    with pytest.raises(ImportError, match=r'biomed-gadget-data\[agent\]'):
        GadgetSession(object(), gadget='unused', max_bytes=1)


@pytest.fixture
def recipe(native, tmp_path):
    """Produce portable native inputs using the existing full-owner contract."""
    import importlib.util
    import json
    from pathlib import Path
    service, store, _, binary, budget = native
    consumer = Path(binary).with_name('gadget_native_biomedical_fields_tests')
    if not consumer.exists():
        pytest.fail('the native biomedical contract consumer is required')
    path = Path(__file__).with_name('native_biomedical_fields.py')
    spec = importlib.util.spec_from_file_location('gadget_recipe_inputs', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    parent = tmp_path / 'recipe-fixture'
    module.main(Path(binary), consumer, parent)
    root = next(parent.glob('biomedical-fields-*'))
    enabled = GadgetSession(service.session, gadget=binary, max_bytes=budget, capabilities=('execute',))
    enabled.publish_file('workspace', root / 'moved-workspace/workspace.rfcxs', expected_version=0)
    intent = json.loads((root / 'intent.json').read_text())
    for name in ('typed', 'untyped'):
        value = dict(intent)
        if name == 'untyped':
            del value['answer_contract']
        path = tmp_path / (name + '.json')
        path.write_text(json.dumps(value))
        subprocess.run([binary, 'native-record', '--input', str(path), '--out', str(path.with_suffix('.rfcxs'))],
                       check=True, capture_output=True)
        enabled.publish_file(name, path.with_suffix('.rfcxs'), expected_version=0)
    return enabled, store, root, intent['answer_contract']['question']


def recipe_catalog(service, intent):
    import json
    result = service.query('FROM $db RETURN GADGET_RECIPE_CATALOG("workspace",1,"' + intent + '",1)')
    return json.loads(result.values[0])


def recipe_proposal(catalog, *, action='recipes'):
    return dict(schema='gadget.pipeline.recipe.proposal.v1', source_digest=catalog['source_digest'],
                intent_sha256=catalog['intent_sha256'], action=action,
                **({'goals': [catalog['recipes'][0]['id']]} if action == 'recipes' else {'reason': 'More input is needed'}))


def test_native_recipe_receipt_keeps_whole_fields_and_pinned_sources(recipe):
    import json
    from rcql import Result, call, query, source
    service, _, _, question = recipe
    catalog = recipe_catalog(service, 'typed')
    result = service.query(query(source('db'), call('GADGET_RECIPE_PROPOSAL', 'workspace', 1, 'typed', 1,
        json.dumps(recipe_proposal(catalog)), question, True)), save_as='recipe-result')
    value = result.values[0]
    assert value['observation']['executed'] is True
    assert value['observation']['scientific_status'] == 'satisfied'
    assert 'admission.rfcxs' in value['publications']
    assert 'execution/glue/result.rfcxs' in value['publications']
    assert all(type(raw) is bytes for raw in value['publications'].values())
    assert Result.from_bytes(result.to_bytes()).to_bytes() == result.to_bytes()
    assert service.session.read('recipe-result', 1).value.to_bytes() == result.to_bytes()
    method = result.execution[0]['methods'][0]
    assert method['intent_version'] == method['workspace_version'] == 1
    assert method['execution_requested'] is True


def test_valid_but_untyped_model_choice_cannot_execute(recipe, tmp_path, monkeypatch):
    import json
    from agent.hive import Hive, ChatResult
    from agent.hive_workflow import WorkflowStep
    service, _, _, _ = recipe
    catalog = recipe_catalog(service, 'untyped')
    proposal = recipe_proposal(catalog)
    monkeypatch.setenv('REXGRAPH_CONFIG_DIR', str(tmp_path / 'hive-config'))
    hive = Hive('unchecked-model-choice')
    # A valid, deliberately unrelated model choice must still retain an untyped
    # admission family, never execute merely because its goal name is available.
    monkeypatch.setattr(hive, 'ask_full', lambda *a, **kw: ChatResult(
        content=json.dumps(proposal), finish_reason='stop', raw={'proposal': proposal}))
    worker = service.recipe_worker(hive, 'model', workspace={'record_id': 'workspace', 'version': 1},
        intent={'record_id': 'untyped', 'version': 1})
    hive.add_worker('recipe', worker, capability='analyze')
    try:
        flow = hive.run_workflow([WorkflowStep('proposal', 'recipe', {
            'question': 'Book an airline ticket for tomorrow.', 'save_as': 'prospective'})])
        assert flow.ok
        result = flow.outputs['proposal']
        value = result.values[0]
        assert value['observation']['status'] == 'admitted'
        assert value['observation']['scientific_status'] == 'untyped'
        assert value['observation']['executed'] is False
        assert set(value['publications']) == {'admission.rfcxs'}
        assert result.native_plan['gadget_recipe_model']['proposal'] == proposal
        assert service.session.read('prospective', 1).value.to_bytes() == result.to_bytes()
    finally:
        hive.stop_all()


def test_wrong_question_or_source_is_refused_without_saved_result(recipe):
    import json
    from rcql import call, query, source
    service, _, _, question = recipe
    proposal = recipe_proposal(recipe_catalog(service, 'typed'))
    for control, text in [(proposal, question + ' changed'),
                          (dict(proposal, source_digest='sha256:foreign'), question)]:
        with pytest.raises(ValueError, match='admission failed'):
            service.query(query(source('db'), call('GADGET_RECIPE_PROPOSAL', 'workspace', 1, 'typed', 1,
                json.dumps(control), text, True)), save_as='must-not-exist')
    with pytest.raises(KeyError):
        service.session.read('must-not-exist', 1)


def test_recipe_proposals_require_capability_and_stay_outside_pure_programs(recipe):
    import json
    from rcql import Program, ProgramStep, call, query, source
    service, _, _, question = recipe
    declaration = query(source('db'), call('GADGET_RECIPE_PROPOSAL', 'workspace', 1, 'typed', 1,
        json.dumps(recipe_proposal(recipe_catalog(service, 'typed'))), question, True))
    plain = GadgetSession(service.session, gadget=service.gadget, max_bytes=service.max_bytes)
    with pytest.raises(PermissionError):
        plain.query(declaration)
    with pytest.raises(ValueError, match='read calculations'):
        Program('recipe_action', (ProgramStep('execute', declaration),))


def test_recipe_worker_rejects_unpinned_sources_and_model_execution_overrides(native):
    service, _, _, _, _ = native
    with pytest.raises(ValueError, match='pins'):
        service.recipe_worker(object(), 'model', workspace={'record_id': 'workspace'},
                              intent={'record_id': 'intent', 'version': 1})
    worker = service.recipe_worker(object(), 'model', workspace={'record_id': 'workspace', 'version': 1},
        intent={'record_id': 'intent', 'version': 1})
    with pytest.raises(ValueError, match='question'):
        worker({'question': 'Inspect the source', 'execute': True})


def test_incomplete_model_output_remains_inspectable_without_native_dispatch(recipe, monkeypatch):
    from agent.hive import ChatResult
    from rcql import Result
    service, _, _, question = recipe
    class Model:
        def ask_full(self, *args, **kwargs):
            return ChatResult(content='{"action":"clarify","reason":"', finish_reason='length',
                              raw={'unfinished_original_text': '{"action":"clarify","reason":"'})
    worker = service.recipe_worker(Model(), 'model', workspace={'record_id': 'workspace', 'version': 1},
        intent={'record_id': 'typed', 'version': 1})
    def no_native_execution(*args, **kwargs):
        raise AssertionError('incomplete model output reached native execution')
    monkeypatch.setattr(service._adapter, 'recipe_publications', no_native_execution)
    result = worker({'question': question, 'save_as': 'unfinished'})
    assert result.values[0]['observation'] == {'status': 'incomplete_model_output', 'executed': False}
    assert result.values[0]['publications'] == {}
    evidence = result.native_plan['gadget_recipe_model']
    assert evidence['proposal'] is None and evidence['raw_reply']['unfinished_original_text']
    assert Result.from_bytes(evidence['catalog_packet']).values
    assert service.session.read('unfinished', 1).value.to_bytes() == result.to_bytes()


def test_valid_wrong_goal_is_blocked_by_actual_typed_native_obligations(recipe, tmp_path):
    import json
    from rcql import call, query, source
    service, _, root, _ = recipe
    intent = json.loads((root / 'intent.json').read_text())
    control = json.loads((root / 'requests.json').read_text())['transfer']
    intent['recipes'].append({'id': 'transfer', 'purpose': 'Apply the declared map',
        'goal': {'method': 'field_transfer', 'parameters': {k: v for k, v in control.items()
            if k not in ('schema', 'operation')}}})
    question = 'Apply the original declared map'
    intent['answer_contract'] = {'schema': 'gadget.pipeline.answer.contract.v1', 'question': question,
        'requirements': [{'id': 'map', 'observable': 'field_transfer', 'bindings': ['transfer']}]}
    path = tmp_path / 'map-intent.json'
    path.write_text(json.dumps(intent))
    subprocess.run([service.gadget, 'native-record', '--input', str(path), '--out', str(path.with_suffix('.rfcxs'))],
                   capture_output=True, check=True)
    service.publish_file('map-intent', path.with_suffix('.rfcxs'), expected_version=0)
    catalog = recipe_catalog(service, 'map-intent')
    proposal = recipe_proposal(catalog)  # existing glue goal, valid name, wrong observable
    result = service.query(query(source('db'), call('GADGET_RECIPE_PROPOSAL', 'workspace', 1,
        'map-intent', 1, json.dumps(proposal), question, True)))
    value = result.values[0]
    assert value['observation']['status'] == 'incompatible'
    assert value['observation']['scientific_status'] == 'unsatisfied'
    assert not value['observation']['executed']
    assert set(value['publications']) == {'admission.rfcxs'}
    from pathlib import Path
    archive = tmp_path / 'incompatible-admission.rfcxs'
    archive.write_bytes(value['publications']['admission.rfcxs'])
    proof = subprocess.run([str(Path(service.gadget).with_name('gadget_native_biomedical_fields_tests')),
        '--admission', str(archive)], capture_output=True, text=True, check=True)
    assert json.loads(proof.stdout)['status'] == 'PASS'


def test_model_goal_budget_matches_native_unique_selection(recipe):
    import json
    from agent.hive import ChatResult
    from rcql import call, query, source
    service, _, _, question = recipe
    catalog = recipe_catalog(service, 'typed')
    names = [entry['id'] for entry in catalog['recipes']]
    assert names
    proposal = dict(recipe_proposal(catalog), goals=names)
    class Model:
        def ask_full(self, *args, **kwargs):
            branches = kwargs['response_format']['json_schema']['schema']['oneOf']
            assert all(next(iter(row['properties'])) == 'action' for row in branches)
            choices = next(row['properties']['goals'] for row in branches if 'goals' in row['properties'])
            # The full legitimate selection must fit, while the repeated-goal
            # sequence that previously consumed the output budget cannot fit.
            assert choices['maxItems'] == len(names) == len(choices['items']['enum'])
            return ChatResult(content=json.dumps(proposal), finish_reason='stop', raw={'proposal': proposal})
    worker = service.recipe_worker(Model(), 'model', workspace={'record_id': 'workspace', 'version': 1},
        intent={'record_id': 'typed', 'version': 1}, execute=False)
    result = worker({'question': question})
    assert result.values[0]['observation']['status'] == 'admitted'
    assert not result.values[0]['observation']['executed']
    duplicate = dict(proposal, goals=[names[0], names[0]])
    with pytest.raises(ValueError, match='admission failed'):
        service.query(query(source('db'), call('GADGET_RECIPE_PROPOSAL', 'workspace', 1, 'typed', 1,
            json.dumps(duplicate), question, False)))
