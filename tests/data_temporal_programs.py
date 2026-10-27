"""Retain whole function programs through native storage and replay."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/rcdb"))
from gadget_data import DataStore, DataWriter, capture_program, replay_program
from gadget_data._records import graph


class NativeTemporalPrograms(unittest.TestCase):
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



if __name__ == "__main__":
    unittest.main()
