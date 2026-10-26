"""Actual source reuse and current native dataset admission across execution paths."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/rcdb"))
from gadget_data import DataStore, DataWriter, NativeBindings
from gadget_data.query import install


class Bindings(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.addCleanup(self.directory.cleanup)
        install()

    def source(self, name="data"):
        with DataWriter(self.root / name, name=name) as writer:
            writer.axis("rows", ("a", "b", "c"), offset=0)
            writer.axis("columns", ("x", "y"), offset=0)
            writer.field("counts", axes=("rows", "columns"), shape=(3, 2), quantity="count", units="count")
            writer.dense("counts", np.array(((0, 2), (3, 4), (5, 6)), dtype=np.int32),
                         np.ones((3, 2), dtype=bool), offset=(0, 0))
            writer.finish()
        store = DataStore(self.root / name)
        self.addCleanup(store.close)
        return store, next(iter(store.manifest["fields"]["counts"]["blocks"]))

    @staticmethod
    def program():
        from rcql import Program, ProgramInput, ProgramStep, parse
        return Program("selected_counts", (ProgramStep("field",
            parse("FROM $counts RETURN GADGET_DATA_FIELD($columns)")),), (ProgramInput("columns"),))

    def test_capture_cold_replay_and_aliases_share_one_original_owner(self):
        store, reference = self.source()
        bindings = NativeBindings()
        with patch.object(store, "read", wraps=store.read) as reads:
            owner = bindings.bind("counts", store, reference)
            self.assertIs(bindings.bind("other", store, reference), owner)
            self.assertEqual(reads.call_count, 1)
            reads.reset_mock()
            bindings.capture(self.program(), self.root / "run", parameters={"columns": (1, 0)})
            output = bindings.replay(self.root / "run", parameters={"columns": (1, 0)})
            self.assertIs(output.values[0].source.source, owner.value)
            np.testing.assert_array_equal(output.values[0].values, ((2, 0), (4, 3), (6, 5)))
            self.assertEqual(reads.call_count, 0)
        bindings.remove("other")
        with self.assertRaises(ValueError):
            bindings.replay(self.root / "run", parameters={"columns": (0, 1)})

    def test_field_program_retains_the_complete_family_and_exact_source_pin(self):
        from gadget_data import field_program
        from fractions import Fraction
        store, reference = self.source()
        bindings = NativeBindings()
        owner = bindings.bind("counts", store, reference)
        program = field_program(reference, source="counts", relation_family=True)
        bindings.capture(program, self.root / "family", parameters={"columns": (0, 1)})
        field, family = bindings.replay(self.root / "family", parameters={"columns": (0, 1)}).values
        self.assertIs(field.source.source, owner.value)
        self.assertIs(family.source.source, owner.value)
        self.assertEqual(family.directions.shape, (3, 1))
        direction = family.directions.apply((Fraction(1),))
        self.assertTrue(any(direction))
        self.assertTrue(all(isinstance(value, Fraction) for value in direction))
        for column in ((0, 3, 5), (2, 4, 6)):
            self.assertEqual(sum(value * coefficient for value, coefficient in zip(column, direction)), 0)
        wrong = field_program(dict(reference, version=2), source="counts")
        with self.assertRaisesRegex(ValueError, "exact program dependency"):
            bindings.capture(wrong, self.root / "wrong", parameters={"columns": (0,)})
        self.assertFalse((self.root / "wrong").exists())

    def test_journal_reuse_and_cold_verification_keep_actual_source(self):
        store, reference = self.source()
        bindings = NativeBindings()
        bindings.bind("counts", store, reference)
        with patch.object(store, "read", side_effect=AssertionError("source was decoded again")):
            with bindings.journal(self.program(), self.root / "journal") as journal:
                first = journal.append("same", parameters={"columns": (0, 1)}, occurrence={"row": 0})
                second = journal.append("same", parameters={"columns": (0, 1)}, occurrence={"row": 1})
                self.assertEqual(first, second)
                journal.finish()
            self.assertEqual(bindings.verify_journal(self.root / "journal"), {"rows": 2, "requests": 1})

    def test_dependent_native_owner_reuses_the_admitted_parent_store(self):
        from gadget_data import record_link
        from gadget_data._records import graph
        from rcql import Program, ProgramStep, parse
        store, reference = self.source()
        bindings = NativeBindings()
        bindings.bind("counts", store, reference)
        bindings.bind("counts_alias", store, reference)
        with DataWriter(self.root / "context", name="context") as writer:
            child = writer.native("context", graph("context", ("original",)),
                                  dependencies=(record_link(store, reference),))
            writer.finish()
        with DataStore(self.root / "context") as context:
            with patch.object(store, "read", side_effect=AssertionError("parent source decoded again")):
                bindings.bind("context", context, child, dependencies=(store,))
                program = Program("declared_context", (ProgramStep("value", parse("FROM $context RETURN 3")),))
                bindings.capture(program, self.root / "dependent", parameters={})
                self.assertEqual(bindings.replay(self.root / "dependent", parameters={}).values, (3,))
            bindings.remove("counts")
            store.close()
            with self.assertRaisesRegex(ValueError, "store is closed"):
                bindings.replay(self.root / "dependent", parameters={})

    def test_revocation_prevents_journal_reuse_and_publication(self):
        store, reference = self.source()
        bindings = NativeBindings()
        bindings.bind("counts", store, reference)
        with bindings.journal(self.program(), self.root / "journal") as journal:
            journal.append("first", parameters={"columns": (0,)})
            bindings.remove("counts")
            with self.assertRaisesRegex(ValueError, "bindings changed"):
                journal.append("second", parameters={"columns": (0,)})
            with self.assertRaisesRegex(ValueError, "bindings changed"):
                journal.finish()
        self.assertFalse((self.root / "journal").exists())

    def test_closed_store_refuses_warm_sources_and_capture(self):
        store, reference = self.source()
        bindings = NativeBindings()
        bindings.bind("counts", store, reference)
        store.close()
        with self.assertRaisesRegex(ValueError, "store is closed"):
            bindings.capture(self.program(), self.root / "run", parameters={"columns": (0,)})
        self.assertFalse((self.root / "run").exists())

    def test_source_revocation_during_execution_prevents_capture_publication(self):
        from rcql import Executor
        store, reference = self.source()
        bindings = NativeBindings()
        bindings.bind("counts", store, reference)
        execute = Executor.execute_program
        def close_after_execution(executor, program, **options):
            result = execute(executor, program, **options)
            store.close()
            return result
        with patch.object(Executor, "execute_program", close_after_execution):
            with self.assertRaisesRegex(ValueError, "store is closed"):
                bindings.capture(self.program(), self.root / "run", parameters={"columns": (0,)})
        self.assertFalse((self.root / "run").exists())

    def test_permissions_aliases_and_current_record_metadata_are_checked(self):
        from rcql import Executor, parse
        store, reference = self.source()
        bindings = NativeBindings()
        with self.assertRaises(ValueError):
            bindings.bind("db", store, reference)
        for permission in ("files", "file_read", "train", "gadget_execute", "*"):
            with self.assertRaises(ValueError):
                bindings.bind("counts", store, reference, permissions=(permission,))
        bindings.bind("counts", store, reference, permissions=("identity",))
        with self.assertRaises(ValueError):
            bindings.bind("counts", store, reference)
        with bindings.scope() as scope:
            with self.assertRaises(PermissionError):
                Executor(sources=scope.sources()).execute(parse("FROM $counts RETURN GADGET_DATA_FIELD([0])"))
        bindings.remove("counts")
        bindings.bind("counts", store, reference)
        with patch.object(store.store, "history", return_value=()):
            with self.assertRaisesRegex(ValueError, "address or retained envelope"):
                with bindings.scope():
                    pass

    def test_detached_parameters_bind_exact_owner_and_reject_other_versions(self):
        from dataclasses import replace
        from rcql import Executor, parse
        store, reference = self.source()
        bindings = NativeBindings()
        owner = bindings.bind("counts", store, reference)
        with bindings.scope() as scope:
            field, = Executor(sources=scope.sources()).execute(parse("FROM $counts RETURN GADGET_DATA_FIELD([0])")).values
            detached = replace(field, source=field.source.detached())
            rebound = scope.parameters({"field": detached})["field"]
            self.assertIs(rebound.source.source, owner.value)
            wrong = replace(detached, source=replace(detached.source, version=2))
            with self.assertRaisesRegex(ValueError, "one original admitted owner"):
                scope.parameters({"field": wrong})
            wrong = replace(field, source=replace(field.source, version=2))
            with self.assertRaisesRegex(ValueError, "one original admitted owner"):
                scope.parameters({"field": wrong})
        bindings.bind("restricted", store, reference, permissions=("identity",))
        with bindings.scope() as scope:
            with self.assertRaises(PermissionError):
                scope.parameters({"field": detached})

    def test_mutated_actual_owner_is_refused_before_capture(self):
        store, reference = self.source()
        bindings = NativeBindings()
        owner = bindings.bind("counts", store, reference)
        owner.value.attach_metadata(1, 0, "changed", True)
        with self.assertRaisesRegex(ValueError, "native source changed"):
            bindings.capture(self.program(), self.root / "run", parameters={"columns": (0,)})
        self.assertFalse((self.root / "run").exists())

    def test_equal_native_states_in_distinct_stores_require_an_explicit_owner(self):
        from gadget_data._records import graph
        from rexgraph.tensor_field import FieldSource
        bindings = NativeBindings()
        value = graph("carrier", ("a", "b"))
        for name in ("one", "two"):
            with DataWriter(self.root / name, name=name) as writer:
                reference = writer.native("carrier", value)
                writer.finish()
            store = DataStore(self.root / name)
            self.addCleanup(store.close)
            bindings.bind(name, store, reference)
        reference = FieldSource(value).detached()
        with bindings.scope() as scope:
            with self.assertRaisesRegex(ValueError, "one original admitted owner"):
                scope.parameters({"source": reference})


if __name__ == "__main__":
    unittest.main()
