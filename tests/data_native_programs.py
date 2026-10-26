"""Retain exact transport families, function programs and parameter domains."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/rcdb"))
from gadget_data import (DataStore, DataWriter, NativeBindings, capture_program, record_link,
                         replay_program, native_owner)
from gadget_data._records import graph


class NativePrograms(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.addCleanup(self.directory.cleanup)

    def test_whole_function_derivative_and_structural_revision_replay(self):
        from rcql import (TemporalModel, Program, ProgramInput, ProgramStep, OutputRef,
                          Executor, parse, BoundSource, SourcePolicy)
        model = TemporalModel.create(81).declare("coordinate",
            {"name": "time", "clock": "event", "unit": "seconds", "origin": "zero"})
        coordinate = model.selected
        model = model.declare("function", {"coordinate": coordinate, "pieces": ({"numerator": (1, 2, 3)},)})
        function = model.selected
        model = model.declare("program", {"opcode": "input"})
        input_program = model.selected
        model = model.declare("program", {"opcode": "derivative", "operands": (input_program,)})
        derivative = model.selected
        with DataWriter(self.root / "owners", name="owners") as writer:
            ref = writer.native("context", graph("context", ("original",)))
            writer.finish()
        program = Program("whole_function", (ProgramStep("function",
            parse("FROM $source RETURN FUNCTION_APPLY($model,$operation,$inputs)")),),
            tuple(ProgramInput(name) for name in ("model", "operation", "inputs")),
            (("function", OutputRef("function")),))
        parameters = {"model": model, "operation": derivative, "inputs": {1: function}}
        with DataStore(self.root / "owners") as owners:
            sources = {"source": BoundSource(owners.native(ref), SourcePolicy.allow("read", "identity"))}
            capture_program(program, self.root / "execution", sources=sources, parameters=parameters,
                            contributors=((owners, ref),))
            result, = replay_program(self.root / "execution", sources=sources, parameters=parameters,
                                     datasets=(owners,)).values
            self.assertEqual(result.value(Q(1, 2)), 5)
            self.assertEqual(result.value(2), 14)
            revised = model.revise_program(derivative, {"opcode": "multiply",
                "operands": (input_program, input_program)})
            parameters = dict(parameters, model=revised, operation=revised.selected)
            actual, = Executor(sources=sources, params=parameters).execute_program(program).values
            self.assertEqual(actual.value(2), 289)
            self.assertEqual(result.value(2), 14)

    def test_parameter_rank_and_original_denominator_hole_survive_retention(self):
        from rcql import CWorkflow
        from rexgraph_cpp import _c_workflow as backend
        payload = backend.Payload
        def row(kind, key, value, dependencies=()):
            return {"kind": kind, "key": key, "revision": 1,
                    "dependencies": dependencies, "payload": value.take()}
        records = (row("Coordinate", 1, payload().string("parameter").u64(0)),
            row("Domain", 2, payload().rational(-2).rational(2).flag(True).flag(True).vectors(()), (0,)),
            row("Expression", 3, payload().vector((0, 1)).vector((0, 1)), (1,)),
            row("Carrier", 4, payload().string("rows").u64(1).u64(0)),
            row("Carrier", 5, payload().string("columns").u64(1).u64(0)),
            row("UnivariateFamily", 6, payload(), (1, 3, 4, 2)))
        workflow = CWorkflow.create({"producer_namespace": 91, "producer": 92,
            "revision": 1, "evidence_namespace": 93, "evidence_cutoff": 0,
            "purpose": "parameter domain"}, records, (5, 2))
        with DataWriter(self.root / "domains", name="domains") as writer:
            reference = writer.native("family", workflow)
            writer.finish()
        with DataStore(self.root / "domains") as store:
            restored = store.native(reference)
            self.assertEqual(restored.coefficient_digest, workflow.coefficient_digest)
            self.assertEqual(restored.value(5).evaluate("rank", 1)["fields"]["rank"], (1,))
            self.assertEqual(restored.value(2).evaluate("expression", 1)["fields"]["value"], 1)
            for index, operation in ((5, "rank"), (2, "expression")):
                with self.assertRaises((ValueError, RuntimeError)):
                    restored.value(index).evaluate(operation, 0)

    def test_transport_retains_free_coefficients_and_identified_field_response(self):
        from rexgraph_cpp import NativeSource
        from rexgraph.cpp_graded_stalk_transport import (cpp_graded_stalk_transport_model,
            graded_stalk_transport_observation_batch, graded_stalk_transport_query)
        from rexgraph.model_runtime import train_model, infer_model
        source = NativeSource.from_declaration({"profile": "declared-role-source-v1", "namespace": 101,
            "semantic_scale": 0, "max_grade": 1, "role_aspects": ("boundary",),
            "types": ({"name": "measurement", "min_arity": 0, "max_arity": 0},),
            "witnesses": (10, 20), "relations": ({"grade": 1, "key": 1, "type": 0,
                "occurrences": ({"key": 1, "lower": 10, "negative": True, "magnitude": 1,
                    "exists": True, "aspects": (0,)}, {"key": 2, "lower": 20, "negative": False,
                    "magnitude": 1, "exists": True, "aspects": (0,)})},),
            "fields": (), "metrics": (), "carriers": ()})
        state = cpp_graded_stalk_transport_model(source, ((('q',), ('q',)), (('q',),)),
            ((), ((), ())), ((0, 0, 0), (1, 0, 0)), (0,), grade=1)
        batch = graded_stalk_transport_observation_batch(state, ((1,),), ((0, (2,)),))
        state = train_model(state, source, batch)
        output = infer_model(state, source, graded_stalk_transport_query(state, ((2,),)))
        self.assertEqual(state.payload["fit"]["fit"]["kernel_shape"], (2, 1))
        self.assertEqual(output.execution_receipt["boundary_identified"], ((True,), (False,)))
        with self.assertRaises(TypeError):
            output.tensor()
        with DataWriter(self.root / "owners", name="owners") as writer:
            owner_ref = writer.native("source", source)
            writer.finish()
        with DataStore(self.root / "owners") as owners:
            link = record_link(owners, owner_ref)
            bound = native_owner(owners, owner_ref)
            with DataWriter(self.root / "transport", name="transport") as writer:
                refs = tuple(writer.native(name, value, dependencies=(link,))
                             for name, value in (("model", state), ("response", output)))
                writer.finish()
            with DataStore(self.root / "transport") as store:
                model, response = (store.native(ref, dependencies=(owners,), owners=(bound,)) for ref in refs)
                self.assertEqual(model.coefficient_digest, state.coefficient_digest)
                self.assertEqual(response.coefficient_digest, output.coefficient_digest)
                self.assertEqual(model.payload["fit"]["fit"]["kernel_shape"], (2, 1))
                self.assertIs(model.bind(bound.value).source, bound.value)
                self.assertIs(response.source.bind(bound.value).source, bound.value)
                replayed = infer_model(model, bound.value, graded_stalk_transport_query(model, ((2,),)))
                self.assertEqual(replayed.coefficient_digest, response.coefficient_digest)
                bindings = NativeBindings()
                actual = bindings.bind("native", owners, owner_ref)
                with bindings.scope() as scope:
                    parameters = scope.parameters({"source": model.source})
                    self.assertIs(parameters["source"].source, actual.value)
                    repeated = infer_model(model, parameters["source"].source,
                        graded_stalk_transport_query(model, ((2,),)))
                    self.assertEqual(repeated.coefficient_digest, response.coefficient_digest)


if __name__ == "__main__":
    unittest.main()
