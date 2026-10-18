"""Native lookup, retained owners and executable program contracts."""
import argparse
from contextlib import ExitStack
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "rcdb"))
from gadget_data import (AxisIndex, DataStore, DataWriter, Limits, build_axis_index,
                         field_program, record_link)
from gadget_data._records import declaration, graph
from gadget_data.query import install
from gadget_data.workflow import restore_record


class NativeWorkflow(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix=self._testMethodName + "-", dir=ROOT))
        install()
        install()

    def dataset(self, name="primary", *, missing=False):
        values = np.array([[1, 0], [0, 1], [1, 1], [2, 3]], dtype=np.int64)
        path = self.root / (name + ".npy")
        np.save(path, values)
        original = np.load(path, mmap_mode="r")
        with DataWriter(self.root / name, name=name) as writer:
            writer.source("original", path)
            for offset in (0, 2):
                writer.axis("features", ("same", "other"), offset=offset,
                    context=({"symbol": "TP53", "flag": False}, {"symbol": "MUC2", "flag": 0}))
            writer.axis("profiles", ("p0", "p1"), offset=0)
            writer.field("counts", axes=("features", "profiles"), shape=values.shape,
                         quantity="read_count", units="count")
            support = np.ones(values.shape, bool)
            if missing:
                support[0, 1] = False
            writer.dense("counts", original, support, offset=(0, 0),
                         origins=({"source": "original", "selector": {"array": "npy", "offset": (0, 0)}},))
            writer.finish()
        return self.root / name

    def execute(self, store, *, family=True):
        from rcql import BoundSource, Executor, SourcePolicy
        reference = store.manifest["fields"]["counts"]["blocks"][0]
        program = field_program(reference, relation_family=family)
        executor = Executor(sources={"db": BoundSource(store.store,
            SourcePolicy.allow("identity", "records", "read"))}, params={"columns": [0, 1]})
        return reference, program, executor, executor.execute_program(program).values

    def test_literal_index_keeps_occurrences_types_and_carriers(self):
        with DataStore(self.dataset()) as first, DataStore(self.dataset("second")) as second:
            with patch.object(first.store, "get_version", wraps=first.store.get_version) as reads:
                build_axis_index(((first, "features"), (second, "features")), self.root / "index",
                                 context_roles=("symbol", "flag"))
                self.assertFalse(any(call.args[0].startswith("fields/") for call in reads.call_args_list))
            with AxisIndex(self.root / "index") as index:
                self.assertEqual(index.verify_index(datasets=(first, second)),
                                 {"locators": 4, "original_positions": 8})
                matches = list(index.find("TP53", datasets=(first, second)))
                self.assertEqual([m["coordinate"] for m in matches], [0, 2, 0, 2])
                self.assertEqual({m["source"]["dataset"]["dataset"] for m in matches},
                                 {first.identity, second.identity})
                self.assertEqual(len(list(index.find("same", datasets=(first, second)))), 4)
                self.assertEqual([m["coordinate"] for m in index.find(False, datasets=(first, second))], [0, 2, 0, 2])
                self.assertEqual([m["coordinate"] for m in index.find(0, datasets=(first, second))], [1, 3, 1, 3])

    def test_index_fallback_and_pinned_versions(self):
        with DataStore(self.dataset()) as source:
            build_axis_index(((source, "features"),), self.root / "index", context_roles=("symbol",))
            with AxisIndex(self.root / "index", limits=Limits(catalog_bytes=2048)) as index:
                first = index.manifest["objects"][0]
                replacement = graph("replacement", ("changed",))
                index.store.put(first["id"], replacement, analytics=False)
                self.assertEqual([m["coordinate"] for m in index.find("TP53", datasets=(source,))], [0, 2])
                index.verify_index(datasets=(source,))
                with self.assertRaises(KeyError):
                    list(index.find("TP53", datasets=()))

    def test_complete_exact_family_matches_independent_equations(self):
        with DataStore(self.dataset()) as source:
            reference, program, executor, (field, family) = self.execute(source)
            original = np.load(self.root / "primary.npy")
            np.testing.assert_array_equal(field.values, original)
            basis = np.asarray([[-1, -2], [-1, -3], [1, 0], [0, 1]], dtype=object)
            actual = np.full(basis.shape, Fraction(0), dtype=object)
            for i, j, coefficient in family.directions.entries:
                actual[i, j] = coefficient
            np.testing.assert_array_equal(actual, basis)
            np.testing.assert_array_equal(original.T @ basis, np.zeros((2, 2), dtype=object))
            self.assertEqual(family.direction_count, 2)
            self.assertEqual(field.source.record_id, reference["id"])
            self.assertEqual(field.source.version, reference["version"])
        with DataStore(self.dataset("foreign")) as foreign:
            from rcql import BoundSource, Executor, SourcePolicy
            with self.assertRaises(ValueError):
                Executor(sources={"db": BoundSource(foreign.store, SourcePolicy.allow("identity", "records", "read"))},
                         params={"columns": [0, 1]}).execute_program(program)

    def test_missing_support_refuses_undeclared_learning(self):
        with DataStore(self.dataset(missing=True)) as source:
            with self.assertRaises(ValueError):
                self.execute(source)

    def test_native_result_and_program_cold_restoration(self):
        from rexgraph.io import save_rcbd, load_rcbd
        primary = self.dataset()
        with DataStore(primary) as source:
            ref, program, executor, (field, family) = self.execute(source)
            link = record_link(source, ref)
            expected = (field.coefficient_digest, family.coefficient_digest, program.coefficient_digest)
            with DataWriter(self.root / "results", name="results") as writer:
                refs = [writer.native(name, value, dependencies=(link,)) for name, value in
                        (("field", field), ("family", family), ("program", program))]
                writer.finish()
        with DataStore(primary) as source, DataStore(self.root / "results") as results:
            for ref, digest in zip(refs, expected, strict=True):
                restored = results.native(ref, dependencies=(source,))
                self.assertEqual(restored.coefficient_digest, digest)
                with self.assertRaises(KeyError):
                    results.native(ref)
                path = self.root / (ref["id"].split("/")[-1] + ".rcbd")
                save_rcbd(path, results.read(ref))
                cold = restore_record(load_rcbd(path), datasets=(source,))
                self.assertEqual(cold.coefficient_digest, digest)
            replay = self.execute(source)[2].execute_program(results.native(refs[2], dependencies=(source,)))
            self.assertEqual(replay.values[1].coefficient_digest, expected[1])
            self.assertEqual(results.verify(), 4)
        code = """from pathlib import Path
from gadget_data import DataStore
from rcql import BoundSource, Executor, SourcePolicy
from gadget_data.query import install
install()
with DataStore(Path(__import__('sys').argv[1])) as source, DataStore(Path(__import__('sys').argv[2])) as results:
    refs = tuple(results.manifest['objects'])
    program = results.native(refs[2], dependencies=(source,))
    executor = Executor(sources={'db': BoundSource(source.store, SourcePolicy.allow('identity','records','read'))}, params={'columns':[0,1]})
    assert executor.execute_program(program).values[1].coefficient_digest == __import__('sys').argv[3]
"""
        subprocess.run([sys.executable, "-B", "-c", code, str(primary), str(self.root / "results"), expected[1]], check=True)

    def test_owner_links_and_live_mutation_are_checked(self):
        from rexgraph.tensor_field import TensorField
        with DataStore(self.dataset()) as source:
            ref, _, _, (field,) = self.execute(source, family=False)
            with DataWriter(self.root / "bad", name="bad") as writer:
                with self.assertRaises(ValueError):
                    writer.native("field", field)
                foreign = record_link(source, ref)
                foreign["record"]["version"] += 1
                with self.assertRaises(ValueError):
                    writer.native("field", field, dependencies=(foreign,))
            field.source.source.attach_metadata(1, 0, "changed", True)
            with DataWriter(self.root / "mutated", name="mutated") as writer:
                with self.assertRaises(ValueError):
                    writer.native("field", field, dependencies=(record_link(source, ref),))
            with DataWriter(self.root / "unscoped", name="unscoped") as writer:
                with self.assertRaises(ValueError):
                    writer.native("field", TensorField(field.space, field.values, field.axes))

    def test_shared_context_and_relocated_complete_export(self):
        from rexgraph.object_identity import object_digest
        from gadget_data.export import export_bundle
        context = graph("original assertions", ("assertion0", "assertion1"),
                        context=({"donor": "D1", "flag": False}, {"donor": "D1", "flag": 0}))
        before = object_digest(context)
        with DataWriter(self.root / "context", name="context") as writer:
            ctx = writer.native("donors", context)
            writer.finish()
        self.assertEqual(object_digest(context), before)
        with DataStore(self.root / "context") as contexts:
            link = record_link(contexts, ctx)
            with DataWriter(self.root / "measurements", name="measurements") as writer:
                for offset in (0, 1):
                    writer.axis("sample", ("s" + str(offset),), offset=offset, context_references=(link,))
                writer.axis("value", ("v",), offset=0)
                writer.field("count", axes=("sample", "value"), shape=(2, 1), quantity="q", units="u")
                writer.dense("count", np.array([[0], [2]], dtype=np.int64), np.ones((2, 1), bool), offset=(0, 0))
                writer.finish()
            with DataStore(self.root / "measurements") as source:
                source.materialize("count", self.root / "selected", rows=(1, 0), columns=(0,), tile=(1, 1))
                source.select("count", rows=(0, 1), columns=(0,)).materialize(self.root / "selected-field")
                for path in ("selected", "selected-field"):
                    with DataStore(self.root / path) as selected:
                        for block in selected.manifest["axes"]["sample"]:
                            restored, = selected.contexts(block, dependencies=(contexts,))
                            self.assertEqual(restored.nE, 2)
                with self.assertRaises(ValueError):
                    export_bundle(source, self.root / "incomplete")
                locations = export_bundle(source, self.root / "complete", dependencies=(contexts,))
                identities = source.identity, contexts.identity
        with ExitStack() as stack:
            source, contexts = [stack.enter_context(DataStore(self.root / "complete" / locations["datasets"][key]))
                                for key in identities]
            for block in source.manifest["axes"]["sample"]:
                owner, = source.contexts(block, dependencies=(contexts,))
                self.assertEqual(owner.nE, 2)
                self.assertEqual(object_digest(owner), ctx["state_digest"])
                self.assertEqual(declaration(source.read(block))["context_roles"], ())

    def test_cpp_publication_and_contributing_owner_cold_restore(self):
        from gadget_data.bridge import export_gadget
        from rexgraph.source_binding import CppSourceReference
        from rexgraph.tensor_field import TensorField
        from rexgraph.type_accession import CoordinateSpace
        from rexgraph.io import load_rcbd, save_rcbd
        import rexgraph_cpp as cpp
        with DataStore(self.dataset()) as source:
            ref, _, _, (input_field,) = self.execute(source, family=False)
            path = export_gadget(source.select("counts", rows=range(4), columns=range(2)),
                self.root / "owner.rfcxs", importer=IMPORTER, namespace=635000, max_bytes=4 * 1024**2)
            native = cpp.NativeSource(path.read_bytes())
            reference = CppSourceReference(native, "model", 3)
            field = TensorField(CoordinateSpace("native rows", native.identity["coordinates"][1]),
                [[Fraction(1, 3)], [-2], [0], [2**70]],
                (CoordinateSpace("profile", ("p0",)),), reference, 1, "chain", dependencies=(input_field.source,))
            with DataWriter(self.root / "owners", name="native owners") as writer:
                owner_ref = writer.native("source", native)
                writer.finish()
            with DataStore(self.root / "owners") as owners:
                with DataWriter(self.root / "results", name="native result") as writer:
                    result_ref = writer.native("field", field, dependencies=(record_link(owners, owner_ref), record_link(source, ref)))
                    writer.finish()
                with DataStore(self.root / "results") as results:
                    restored = results.native(result_ref, dependencies=(source, owners))
                    self.assertEqual(restored.coefficient_digest, field.coefficient_digest)
                    self.assertEqual(restored.source.source.archive, native.archive)
                    self.assertEqual(restored.source.record_id, "model")
                    self.assertEqual(restored.source.version, 3)
                    self.assertIsNotNone(restored.dependencies[0].source)
                    with self.assertRaises(KeyError):
                        results.native(result_ref, dependencies=(owners,))
                    save_rcbd(self.root / "cpp-result.rcbd", results.read(result_ref))
                    cold = restore_record(load_rcbd(self.root / "cpp-result.rcbd"), datasets=(source, owners))
                    self.assertEqual(cold.coefficient_digest, field.coefficient_digest)
                    cold.check_state()

    def test_forged_axis_address_and_lookup_budget_refuse(self):
        with DataStore(self.dataset()) as source:
            block = source.manifest["axes"]["features"][0]
            source._axis_metadata(block)
            forged = dict(block, offset=10)
            with self.assertRaises(ValueError):
                source._axis_metadata(forged)
            build_axis_index(((source, "features"),), self.root / "index", context_roles=("symbol",))
            with AxisIndex(self.root / "index", limits=Limits(catalog_bytes=32)) as index:
                with self.assertRaises(MemoryError):
                    list(index.find("TP53", datasets=(source,)))

    def test_cli_lookup_uses_explicit_sources(self):
        import json
        first, second = self.dataset(), self.dataset("second")
        command = [sys.executable, "-B", "-m", "gadget_data"]
        subprocess.run([*command, "build-index", str(self.root / "index"),
            "--source", str(first), "--axis", "features", "--source", str(second), "--axis", "features",
            "--context-role", "symbol"], check=True, capture_output=True, text=True)
        result = subprocess.run([*command, "lookup", str(self.root / "index"), '"TP53"',
            "--source", str(first), "--source", str(second)], check=True, capture_output=True, text=True)
        self.assertEqual([json.loads(line)["coordinate"] for line in result.stdout.splitlines()], [0, 2, 0, 2])

    def test_program_reflection_transformation_and_composition(self):
        from rcql import BoundSource, Executor, NameRelation, SourcePolicy, bind, parse
        policy = SourcePolicy.allow("identity", "records", "read")
        with DataStore(self.dataset()) as source:
            ref, program, executor, values = self.execute(source)
            definition = program.to_record()
            reflective = Executor(sources={"definition": BoundSource(definition, policy)}, params={
                "sources": {"db": bind("db", source.store, policy)}, "parameters": {"columns": [0, 1]}})
            query = ('FROM $definition LET self=PROGRAM_READ() '
                     'LET change=TRANSFORM_PROGRAM(self,"step",["family","retained_family"]) '
                     'RETURN PROGRAM_TOPOLOGY(self,$sources,$parameters), '
                     'TRANSFORM_PROGRAM_COMPILE(change,self,$sources,$parameters)')
            topology, changed = reflective.execute(parse(query)).values
            self.assertGreater(topology.nE, 0)
            self.assertEqual(executor.execute_program(changed).values[1].coefficient_digest, values[1].coefficient_digest)
            composed = Executor(sources={"db": BoundSource(source.store, policy)}, params={
                "record": ref["id"], "version": ref["version"],
                "accession": NameRelation.operator("GADGET_DATA_FIELD"),
                "learning": NameRelation.operator("GADGET_RELATION_FAMILY"), "arguments": {"columns": [0, 1]}})
            result = composed.execute(parse('FROM RCDB_VERSION($db,$record,$version) '
                'LET change=TRANSFORM_COMPOSE($accession,$learning,["field"]) '
                'LET operation=TRANSFORM_COMPILE(change,$accession,$arguments) '
                'RETURN NAME_APPLY(operation,$arguments)')).values[0]
            self.assertEqual(result.coefficient_digest, values[1].coefficient_digest)

    def test_retained_field_pure_rcql_accession(self):
        from rcql import BoundSource, Executor, Program, ProgramInput, ProgramStep, OutputRef, SourcePolicy, parse
        with DataStore(self.dataset()) as source:
            ref, _, _, (field,) = self.execute(source, family=False)
            with DataWriter(self.root / "results", name="results") as writer:
                retained = writer.native("field", field, dependencies=(record_link(source, ref),))
                writer.finish()
            with DataStore(self.root / "results") as results:
                program = Program("retained", (ProgramStep("field", parse(
                    'FROM RCDB_VERSION($db,"objects/field",1) RETURN GADGET_RETAINED_FIELD($owners)')),),
                    (ProgramInput("owners"),), (("field", OutputRef("field")),))
                executor = Executor(sources={"db": BoundSource(results.store,
                    SourcePolicy.allow("identity", "records", "read"))}, params={"owners": [field.source]})
                result = executor.execute_program(program).values[0]
                self.assertEqual(result.coefficient_digest, field.coefficient_digest)
                executor.params["owners"] = []
                with self.assertRaises(ValueError):
                    executor.execute_program(program)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--importer", type=Path, required=True)
    args, extra = parser.parse_known_args()
    ROOT = args.root.resolve(); ROOT.mkdir(parents=True, exist_ok=True)
    IMPORTER = args.importer.resolve(strict=True)
    unittest.main(argv=[sys.argv[0], *extra])
