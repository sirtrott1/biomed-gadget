"""Original input, selective native I/O and cold materialization contracts."""
import argparse
from contextlib import closing
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "rcdb"))
from gadget_data import DataStore, DataWriter, Limits, SourceResolver
from gadget_data.readers import import_coordinate, import_dense, import_records, import_measurements


class DataIO(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix=self._testMethodName + "-", dir=ROOT))

    def test_selected_axis_metadata_obeys_selection_admission(self):
        labels = ("a" * 8192, "b" * 8192)
        context = ({"detail": "x" * 8192}, {"detail": "y" * 8192})
        with DataWriter(self.root / "data", name="metadata") as writer:
            writer.axis("rows", labels, offset=0, context=context)
            writer.axis("columns", ("value",), offset=0)
            writer.field("value", axes=("rows", "columns"), shape=(2, 1), quantity="q", units="u")
            for row in range(2):
                writer.dense("value", np.asarray([[row + 3]], dtype=np.int64), np.ones((1, 1), bool), offset=(row, 0))
            writer.finish()
        with DataStore(self.root / "data", limits=Limits(selection_bytes=512)) as store:
            with self.assertRaises(MemoryError):
                store.axis("rows", (0, 1))
            with self.assertRaises(MemoryError):
                store.select("value", rows=(0, 1), columns=(0,))
            with patch.object(store.store, "get_version", wraps=store.store.get_version) as read:
                with self.assertRaises(MemoryError):
                    store.provenance("value", rows=(0, 1), columns=(0,))
                self.assertFalse(any(call.args[0].startswith("fields/") for call in read.call_args_list))
        with DataStore(self.root / "data", limits=Limits(selection_bytes=64 * 1024)) as store:
            selected = store.select("value", rows=(0, 1), columns=(0,))
            self.assertEqual(selected.row_labels, labels)
            self.assertEqual(selected.row_context, context)
            np.testing.assert_array_equal(selected.values, [[3], [4]])

    def dense(self, *, missing=False):
        original = np.arange(72, dtype=np.int64).reshape(12, 6)
        if missing:
            original = original.astype(np.float32)
            original[2, 1] = np.nan
            original[0, 0] = -0.0
        np.save(self.root / "original.npy", original)
        import_dense(self.root / "original.npy", self.root / "data", reader="npy", name="original",
            quantity="measured_count", units="count", missing="nan" if missing else None,
            row_labels=("same" for _ in range(12)), tile=(3, 2))
        return original

    def test_catalog_and_selection_open_only_needed_blocks(self):
        original = self.dense()
        with DataStore(self.root / "data") as store:
            with patch.object(store.store, "get_version", wraps=store.store.get_version) as read:
                catalog = store.catalog()
                self.assertEqual(read.call_count, 0)
                provenance = store.provenance("measurements", rows=(8, 6), columns=(3, 2))
                self.assertTrue(all(call.args[0].startswith("catalog/") for call in read.call_args_list))
                self.assertEqual(provenance["field_records"][0]["origins"][0]["source"], "original")
                selected = store.select("measurements", rows=(8, 6), columns=(3, 2))
                np.testing.assert_array_equal(selected.values, original[np.ix_((8, 6), (3, 2))])
                self.assertEqual(selected.row_labels, ("same", "same"))
                fields = [call.args[0] for call in read.call_args_list if call.args[0].startswith("fields/")]
                self.assertEqual(fields, ["fields/measurements/7"])
                self.assertEqual(selected.tensor().values.tolist(), [[51, 50], [39, 38]])
                self.assertEqual(catalog["fields"]["measurements"]["shape"], (12, 6))
            whole = np.empty_like(original)
            for selected in store.iter_select("measurements"):
                whole[np.ix_(selected.rows, selected.columns)] = selected.values
            np.testing.assert_array_equal(whole, original)

    def test_zero_absence_ieee_bits_and_exact_math(self):
        original = self.dense(missing=True)
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=(0, 2), columns=(0, 1))
            self.assertTrue(np.signbit(selected.values[0, 0]))
            self.assertTrue(selected.observed[0, 0])
            self.assertFalse(selected.observed[1, 1])
            with self.assertRaises(ValueError):
                selected.tensor()
            from rexgraph.value import Absent
            self.assertIs(selected.exact().values()[1, 1], Absent)
            store.materialize("measurements", self.root / "view", rows=(0, 2), columns=(0, 1), tile=(1, 1))
        with DataStore(self.root / "view") as view:
            selected = view.select("measurements", rows=(0, 1), columns=(0, 1))
            self.assertTrue(np.signbit(selected.values[0, 0]))
            self.assertFalse(selected.observed[1, 1])
            axis = view.read(view.manifest["axes"]["rows"][1])
            from gadget_data._records import declaration
            self.assertEqual(declaration(axis)["mapping"]["coordinates"], (2,))

    def test_bigint_exact_and_independent_carriers(self):
        from rexgraph.value import Absent
        big = 2**190 + 33
        values = [[big, Fraction(-7, 13)], [0, Absent]]
        for name in ("left", "right"):
            with DataWriter(self.root / name, name=name) as writer:
                writer.axis("rows", ["a", "a"], offset=0, context=[{"flag": False}, {"flag": 0}])
                writer.axis("columns", ["x", "y"], offset=0)
                writer.field("value", axes=("rows", "columns"), shape=(2, 2), layout="exact", quantity="q", units="u")
                writer.exact("value", values, offset=(0, 0))
                writer.finish()
        with DataStore(self.root / "left") as left, DataStore(self.root / "right") as right:
            selected = left.select("value", rows=(1, 0), columns=(1, 0))
            self.assertEqual(selected.values[1, 1], big)
            self.assertIs(selected.values[0, 0], Absent)
            self.assertEqual(selected.row_context, ({"flag": 0}, {"flag": False}))
            with self.assertRaises(ValueError):
                left.read(right.manifest["fields"]["value"]["blocks"][0])
            a = left.select("value", rows=(0,), columns=(0,)).tensor()
            b = right.select("value", rows=(0,), columns=(0,)).tensor()
            self.assertFalse(a.source.matches(b.source))

    def test_sparse_occurrences_keep_duplicate_zero_and_original_order(self):
        path = self.root / "original.mtx"
        path.write_text("%%MatrixMarket matrix coordinate integer general\n% original order\n4 5 5\n1 2 0\n1 2 3\n4 5 8\n1 2 -1\n2 4 0\n")
        import_coordinate(path, self.root / "data", name="sparse", quantity="q", units="u",
                          implicit="observed_zero", block_entries=2)
        with DataStore(self.root / "data") as store:
            selected = store.occurrences("measurements", rows=(0,), columns=(1,))
            np.testing.assert_array_equal(selected.indices, [[0, 1]] * 3)
            np.testing.assert_array_equal(selected.values, [0, 3, -1])
            np.testing.assert_array_equal(selected.original_entries, [0, 1, 3])
            self.assertEqual(selected.implicit, "observed_zero")
            streamed = tuple(store.iter_occurrences("measurements"))
            self.assertEqual(sum(chunk.values.size for chunk in streamed), 5)
            store.materialize("measurements", self.root / "sparse-view", rows=(3, 0), columns=(4, 1), tile=(1, 1))
            with self.assertRaises(ValueError):
                store.select("measurements", rows=(0,), columns=(1,))
        with DataStore(self.root / "sparse-view") as view:
            selected = view.occurrences("measurements", rows=(0, 1), columns=(0, 1))
            np.testing.assert_array_equal(selected.values, [0, 3, 8, -1])
            np.testing.assert_array_equal(selected.original_entries, [0, 1, 2, 3])
            np.testing.assert_array_equal(selected.indices, [[1, 1], [1, 1], [0, 0], [1, 1]])

    def test_native_declaration_retains_rows_and_ordered_arity(self):
        from rexgraph import RelationSpec
        from rexgraph.io.declaration import DatasetDeclaration
        from rexgraph.io.records import RecordSchema, RecordField
        path = self.root / "rows.jsonl"
        path.write_text('\n'.join(json.dumps({"participants": ["a", "b", "c", "d"], "name": label})
                                   for label in ("first", "second")) + '\n')
        specification = DatasetDeclaration("jsonl", RecordSchema((RecordField("participants", "sequence"),
            RecordField("name", "string"))), RelationSpec("row", ("participants",), attributes=("name",)),
            options={"batch_size": 1})
        import_records(path, self.root / "data", declaration=specification, name="rows")
        with DataStore(self.root / "data") as store:
            self.assertEqual(len(store.manifest["records"]), 2)
            self.assertEqual(store.read(store.manifest["declarations"][0]).digest, specification.digest)
            for index, ref in enumerate(store.manifest["records"]):
                owner = store.read(ref)
                self.assertEqual(tuple(owner.relations.support_idx), (0, 1, 2, 3))
                self.assertEqual(owner.nE, 1)
                self.assertEqual(owner.relations.attributes[1][0]["name"], "first" if index == 0 else "second")

    def test_declared_measurement_table_exact_values_context_and_sources(self):
        from rexgraph import RelationSpec, ValueRules
        from rexgraph.io.declaration import DatasetDeclaration
        from rexgraph.io.records import RecordSchema, RecordField
        path = self.root / "table.jsonl"
        path.write_text('\n'.join(json.dumps(row) for row in (
            {"row": "same", "subject": "person", "flag": False, "a": "0", "b": "0.1"},
            {"row": "same", "subject": "person", "flag": True, "a": "123456789012345678901234567890", "b": ""})) + '\n')
        specification = DatasetDeclaration("jsonl", RecordSchema((RecordField("row"), RecordField("subject"),
            RecordField("flag", "bool"), RecordField("a", "number", ValueRules("decimal_exact", ("",))),
            RecordField("b", "number", ValueRules("decimal_exact", ("",))))),
            RelationSpec("pair", ("row", "subject")), options={"batch_size": 1})
        import_measurements(path, self.root / "data", declaration=specification, name="table",
            measurements=("a", "b"), label="row", context=("flag",), quantity="measurement", units="u",
            index_labels=True, indexed_context=("flag",))
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=(0, 1), columns=(0, 1))
            self.assertEqual(selected.values[0, 1], Fraction(1, 10))
            self.assertEqual(selected.values[1, 0], 123456789012345678901234567890)
            self.assertFalse(selected.observed[1, 1])
            self.assertEqual(selected.row_labels, ("same", "same"))
            self.assertEqual(selected.row_context, ({"flag": False}, {"flag": True}))
            self.assertEqual(len(store.manifest["records"]), 2)
            self.assertEqual([match["coordinate"] for match in store.find_axis("records", label="same", context={"flag": False})], [0])
            self.assertEqual(list(store.find_axis("records", context={"flag": 0})), [])
            selected.row_context[0]["flag"] = True
            self.assertEqual([match["coordinate"] for match in store.find_axis("records", context={"flag": False})], [0])
            with self.assertRaises(ValueError):
                selected.tensor()

    def test_sparse_decimal_and_bigint_are_exact(self):
        path = self.root / "original.mtx"
        path.write_text("%%MatrixMarket matrix coordinate real general\n2 3 2\n1 2 0.10000000000000000000000000000001\n2 3 -0.25\n")
        import_coordinate(path, self.root / "data", name="sparse", quantity="q", units="u", implicit="missing", block_entries=1)
        with DataStore(self.root / "data") as store:
            selected = store.occurrences("measurements", rows=(0, 1), columns=(1, 2))
            self.assertEqual(tuple(selected.values), (Fraction("0.10000000000000000000000000000001"), Fraction(-1, 4)))
            self.assertEqual(selected.implicit, "missing")

    def test_long_form_keeps_complete_groups_across_reader_batches(self):
        from rexgraph import RelationSpec
        from rexgraph.io.declaration import DatasetDeclaration
        from rexgraph.io.records import RecordSchema, RecordField
        path = self.root / "groups.csv"
        path.write_text("group,participant\ng,a\ng,b\ng,c\nh,d\nh,e\n")
        specification = DatasetDeclaration("csv", RecordSchema((RecordField("group"), RecordField("participant"))),
            RelationSpec("long", ("participant",), identity="group"), options={"batch_size": 1})
        import_records(path, self.root / "data", name="groups", declaration=specification, grouping="contiguous")
        with DataStore(self.root / "data") as store:
            owners = [store.read(reference) for reference in store.manifest["records"]]
            self.assertEqual([owner.nE for owner in owners], [1, 1])
            self.assertEqual([len(owner.relations.support_idx) for owner in owners], [3, 2])
            self.assertEqual([record["fields"]["participant"] for record in owners[0].get_metadata(1, 0, "input_records")], ["a", "b", "c"])
        path.write_text("group,participant\ng,a\ng,b\nh,c\nh,d\ng,e\ng,f\n")
        with self.assertRaises(ValueError):
            import_records(path, self.root / "failed", name="groups", declaration=specification, grouping="contiguous")
        self.assertFalse((self.root / "failed").exists())

    def test_typed_axis_labels_and_unsigned_endian_values(self):
        original = np.array([[0, 2**64 - 1]], dtype=">u8")
        np.save(self.root / "original.npy", original)
        import_dense(self.root / "original.npy", self.root / "data", reader="npy", name="unsigned",
                     row_labels=(Fraction(2, 3),), column_labels=(b"a", b"b"), quantity="q", units="u")
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=(0,), columns=(0, 1))
            self.assertEqual(selected.values.dtype, original.dtype)
            self.assertEqual(selected.values.tobytes(), original.tobytes())
            self.assertEqual(int(selected.tensor().values[0, 1]), 2**64 - 1)
            from gadget_data.bridge import export_gadget
            output = export_gadget(selected, self.root / "typed.rfcxs", importer=IMPORTER,
                                   namespace=904000, max_bytes=16 * 1024**2)
        record = json.loads(subprocess.run([str(CONSUMER), "--read", str(output)],
                             check=True, capture_output=True, text=True).stdout)
        from rexgraph.value_codec import unpack_value
        identities = unpack_value(bytes.fromhex(record["coordinates"]["identities_native"]["payload_hex"]))
        self.assertEqual(identities["row_labels"], (Fraction(2, 3),))
        self.assertEqual(identities["column_labels"], (b"a", b"b"))
        self.assertEqual(record["values"], ["0", str(2**64 - 1)])

    def test_refusals_before_decode_and_pinned_versions(self):
        self.dense()
        with DataStore(self.root / "data") as store:
            block = store.manifest["fields"]["measurements"]["blocks"][0]
            original = store.read(block)
            store.store.put("fields/measurements/0", original, meta={"changed": True}, analytics=False)
            self.assertEqual(store.store.get_record(block["id"]).version, 2)
            np.testing.assert_array_equal(store.select("measurements", rows=(0,), columns=(0,)).values, [[0]])
            store.limits = Limits(decode_bytes=1)
            with patch.object(store.store, "get_version", side_effect=AssertionError("decoded before admission")):
                with self.assertRaises(MemoryError):
                    store.read(block)
            store.limits = Limits(selection_bytes=1)
            with patch.object(store.store, "get_version", side_effect=AssertionError("decoded before admission")):
                with self.assertRaises(MemoryError):
                    store.select("measurements", rows=(0,), columns=(0,))
            substituted = dict(block, version=2)
            with self.assertRaises(ValueError):
                store.read(substituted)
            from gadget_data.cli import _coordinates
            store.limits = Limits()
            with self.assertRaises(MemoryError):
                store.select("measurements", rows=_coordinates("0:1000000000"), columns=(0,))

    def test_atomic_failure_original_mutation_and_fresh_destinations(self):
        path = self.root / "raw"
        path.write_bytes(b"original")
        with DataWriter(self.root / "failed", name="failure") as writer:
            writer.source("raw", path)
            path.write_bytes(b"changed")
            with self.assertRaises(ValueError):
                writer.finish()
        self.assertFalse((self.root / "failed").exists())
        path.write_bytes(b"original")
        with DataWriter(self.root / "restored", name="restored") as writer:
            writer.source("raw", path)
            path.write_bytes(b"changed!")
            path.write_bytes(b"original")
            with self.assertRaises(ValueError):
                writer.finish()
        self.assertFalse((self.root / "restored").exists())
        with DataWriter(self.root / "gap", name="gap") as writer:
            writer.axis("rows", [0, 1], offset=0)
            writer.axis("columns", [0], offset=0)
            writer.field("value", axes=("rows", "columns"), shape=(2, 1), quantity="q", units="u")
            writer.dense("value", np.array([[0]]), np.ones((1, 1), bool), offset=(0, 0))
            with self.assertRaises(ValueError):
                writer.finish()
        self.assertFalse((self.root / "gap").exists())
        (self.root / "exists").mkdir()
        (self.root / "exists" / "keep").write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            DataWriter(self.root / "exists", name="collision")
        self.assertEqual((self.root / "exists" / "keep").read_bytes(), b"keep")

    def test_relocation_original_resolution_and_portable_closure(self):
        from gadget_data.export import export_bundle
        self.dense()
        with DataStore(self.root / "data") as source:
            source.materialize("measurements", self.root / "view", rows=(11, 0), columns=(5, 1), tile=(1, 1))
            with DataStore(self.root / "view") as view:
                with self.assertRaises(ValueError):
                    export_bundle(view, self.root / "missing-parent")
                index = export_bundle(view, self.root / "bundle", dependencies=(source,),
                    resolvers={source.identity: SourceResolver(self.root)}, include_originals=True)
                self.assertEqual(index["root"], view.identity)
        relocated = self.root / "relocated"
        shutil.move(self.root / "bundle", relocated)
        self.assertTrue((relocated / "locations.json").is_file())
        with DataStore(relocated / index["datasets"][index["root"]]) as view:
            np.testing.assert_array_equal(view.select("measurements", rows=(0, 1), columns=(0, 1)).values, [[71, 67], [5, 1]])
        resolver = SourceResolver(paths={"original": self.root / "original.npy"})
        with DataStore(self.root / "data") as source:
            self.assertEqual(source.verify_sources(resolver)["original"], self.root / "original.npy")
            (self.root / "original.npy").write_bytes(b"wrong")
            with self.assertRaises(ValueError):
                source.verify_sources(resolver)

    def test_rcql_native_field_and_cold_cli(self):
        from rcql import BoundSource, SourcePolicy, run
        from gadget_data.query import install
        self.dense()
        install(); install()
        with DataStore(self.root / "data") as store:
            result = run("FROM RCDB_VERSION($db,$record,1) RETURN GADGET_DATA_FIELD($columns)",
                sources={"db": BoundSource(store.store, SourcePolicy.allow("read", "identity", "records"))},
                params={"record": "fields/measurements/0", "columns": (1, 0)}).values[0]
            self.assertEqual(result.values.tolist(), [[1, 0], [7, 6], [13, 12]])
            self.assertEqual(result.source.record_id, "fields/measurements/0")
            self.assertEqual(result.source.version, 1)
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "tools" / "rcdb"))
        command = [sys.executable, "-B", "-m", "gadget_data"]
        process = subprocess.run([*command, "explore", str(self.root / "data")], check=True, capture_output=True, text=True, env=env)
        self.assertEqual(json.loads(process.stdout)["name"], "original")
        subprocess.run([*command, "select", str(self.root / "data"), str(self.root / "selected"),
            "--field", "measurements", "--rows", "11,0", "--columns", "5,1"], check=True, capture_output=True, env=env)
        with DataStore(self.root / "selected") as selected:
            np.testing.assert_array_equal(selected.select("measurements", rows=(0, 1), columns=(0, 1)).values, [[71, 67], [5, 1]])

    def test_hdf5_uses_selected_hyperslabs(self):
        import h5py
        values = np.arange(240, dtype=np.float32).reshape(20, 12)
        with h5py.File(self.root / "original.h5", "w") as file:
            file.create_dataset("values", data=values, chunks=(5, 4))
            file["foreign"] = h5py.ExternalLink("absent.h5", "/values")
        import_dense(self.root / "original.h5", self.root / "data", reader="hdf5", dataset="values",
                     name="selected", quantity="q", units="u", tile=(5, 4))
        with DataStore(self.root / "data") as store:
            actual = np.empty_like(values)
            for selected in store.iter_select("measurements"):
                actual[np.ix_(selected.rows, selected.columns)] = selected.values
            np.testing.assert_array_equal(actual.view("uint32"), values.view("uint32"))
        with self.assertRaises(ValueError):
            import_dense(self.root / "original.h5", self.root / "foreign", reader="hdf5", dataset="foreign",
                         name="foreign", quantity="q", units="u")
        self.assertFalse((self.root / "foreign").exists())

    def test_gadget_native_export_actual_reader_and_field_selector(self):
        from gadget_data.bridge import export_gadget
        self.dense(missing=True)
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=(0, 2), columns=(0, 1))
            destination = export_gadget(selected, self.root / "field.rfcxs", importer=IMPORTER,
                namespace=900100, max_bytes=16 * 1024**2)
        read = json.loads(subprocess.run([str(CONSUMER), "--read", str(destination)],
            check=True, capture_output=True, text=True).stdout)
        self.assertEqual(read["values"], ["0", "1", "12", None])
        declaration = read["coordinates"]
        self.assertEqual(declaration["rows"], [0, 2])
        self.assertEqual(declaration["columns"], [0, 1])
        import rexgraph_cpp as cpp
        owner = cpp.NativeSource(destination.read_bytes(), max_bytes=16 * 1024**2)
        request = {"schema": "gadget.rcdb.field.query.v1", "source_namespace": "900101",
            "cutoff": owner.identity["sequence"], "tensor": 72, "input_carrier": None,
            "output_carrier": 71, "interpretation": "Original selected row and measurement coordinates with declared observed support"}
        (self.root / "request.json").write_text(json.dumps(request))
        subprocess.run([str(SELECTOR), "--input", str(destination), "--request", str(self.root / "request.json"),
            "--out", str(self.root / "selected.rfcxs")], check=True, capture_output=True)
        with self.assertRaises(FileExistsError):
            export_gadget(selected, destination, importer=IMPORTER, namespace=900100, max_bytes=16 * 1024**2)
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=(8, 6), columns=(3, 2))
            destination = export_gadget(selected, self.root / "matrix.rfcxs", importer=IMPORTER,
                namespace=900200, max_bytes=16 * 1024**2, source_rows=1)
            actual = json.loads(subprocess.run([str(CONSUMER), "--apply", str(destination)],
                check=True, capture_output=True, text=True).stdout)
            self.assertEqual(actual["values"], ["252", "192"])
            self.assertEqual(len(actual["coordinates"]["native_source"]["components"]), 2)
            original_run = subprocess.run

            def invalid_component(command, **kwargs):
                import hashlib
                path = Path(command[command.index("--components") + 1])
                manifest = json.loads(path.read_text())
                manifest["parts"][1]["offset"] = 0
                raw = json.dumps(manifest) + "\n"
                path.write_text(raw)
                control_path = Path(command[command.index("--coordinates") + 1])
                control = json.loads(control_path.read_text())
                control["native_source"]["components_sha256"] = hashlib.sha256(raw.encode()).hexdigest()
                control_path.write_text(json.dumps(control))
                return original_run(command, **kwargs)

            with patch("gadget_data.bridge.subprocess.run", side_effect=invalid_component):
                with self.assertRaises(subprocess.CalledProcessError):
                    export_gadget(selected, self.root / "overlap.rfcxs", importer=IMPORTER,
                        namespace=900300, max_bytes=16 * 1024**2, source_rows=1)
            self.assertFalse((self.root / "overlap.rfcxs").exists())

    def test_native_field_beyond_source_constructor_dimension_profile(self):
        from gadget_data.bridge import export_gadget
        original = np.zeros((1024, 1024), dtype=np.float32)
        original[0, 0], original[-1, -1] = 7, np.nan
        np.save(self.root / "original.npy", original)
        import_dense(self.root / "original.npy", self.root / "data", reader="npy", name="wide",
                     quantity="q", units="u", missing="nan")
        with DataStore(self.root / "data") as store:
            selected = store.select("measurements", rows=range(1024), columns=range(1024))
            output = export_gadget(selected, self.root / "wide.rfcxs", importer=IMPORTER,
                namespace=901000, max_bytes=32 * 1024**2, working_bytes=2 * 1024**3)
        result = json.loads(subprocess.run([str(CONSUMER), "--summary", str(output)],
            check=True, capture_output=True, text=True).stdout)
        self.assertEqual(result, {"coordinates": 1048576, "observed": 1048575, "nonzeros": 1, "first": "7", "last": None})
        import rexgraph_cpp as cpp
        owner = cpp.NativeSource(output.read_bytes(), max_bytes=32 * 1024**2,
                                 max_records=4 * 1048576, max_dimension=4 * 1048576, max_logical_size=4 * 1048576)
        self.assertEqual(owner.identity["namespace"], 901001)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--consumer", type=Path, required=True)
    parser.add_argument("--selector", type=Path, required=True)
    args, extra = parser.parse_known_args()
    ROOT = args.root.resolve(); ROOT.mkdir(parents=True, exist_ok=True)
    IMPORTER, CONSUMER, SELECTOR = (path.resolve() for path in (args.importer, args.consumer, args.selector))
    unittest.main(argv=[sys.argv[0], *extra])
