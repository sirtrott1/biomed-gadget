"""Complete request reuse preserves every original physical occurrence."""
from pathlib import Path
from fractions import Fraction as Q
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/rcdb"))
from gadget_data import (DataStore, DataWriter, ProgramJournal, field_program,
                         journal_occurrences, verify_journal)
from gadget_data.query import install


class Journal(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.addCleanup(self.directory.cleanup)
        install()
        with DataWriter(self.root / "source", name="source") as writer:
            writer.axis("rows", ("a", "b", "c", "d"), offset=0)
            writer.axis("measurements", ("x", "y"), offset=0)
            writer.field("field", axes=("rows", "measurements"), shape=(4, 2),
                         quantity="count", units="count")
            self.reference = writer.dense("field", np.array([[1, 0], [0, 1], [1, 1], [2, 3]]),
                                          np.ones((4, 2), bool), offset=(0, 0))
            writer.finish()

    def sources(self, source):
        from rcql import BoundSource, SourcePolicy
        return {"db": BoundSource(source.store, SourcePolicy.allow("read", "identity", "records"))}

    def test_one_thousand_rows_compute_one_complete_family(self):
        from rexgraph.section_calculus import SectionFamily
        with DataStore(self.root / "source") as source:
            program = field_program(self.reference, relation_family=True)
            sources = self.sources(source)
            with patch.object(SectionFamily, "solve", wraps=SectionFamily.solve) as solve:
                with ProgramJournal(self.root / "journal", program, sources=sources,
                                    contributors=((source, self.reference),)) as journal:
                    for i in range(1000):
                        journal.append("repeated identity", parameters={"columns": [0, 1]},
                                       occurrence={"sample": i, "observed": i % 2 == 0})
                    self.assertEqual((journal.rows, journal.requests), (1000, 1))
                    self.assertLessEqual(len(journal._pending), 128)
                    journal.finish()
                self.assertEqual(solve.call_count, 1)
            with DataStore(self.root / "journal") as stored:
                rows = list(journal_occurrences(stored, dependencies=(source,)))
                self.assertEqual([row["ordinal"] for row in rows], list(range(1000)))
                self.assertEqual([row["occurrence"]["sample"] for row in rows], list(range(1000)))
                self.assertEqual(len({row["request"]["state_digest"] for row in rows}), 1)
            with patch.object(SectionFamily, "solve", wraps=SectionFamily.solve) as solve:
                self.assertEqual(verify_journal(self.root / "journal", sources=sources, datasets=(source,)),
                                 {"rows": 1000, "requests": 1})
                self.assertEqual(solve.call_count, 1)

    def test_equal_particular_fields_do_not_merge_different_complete_families(self):
        with DataStore(self.root / "source") as source:
            with ProgramJournal(self.root / "journal", field_program(self.reference, relation_family=True),
                                sources=self.sources(source), contributors=((source, self.reference),)) as journal:
                a = journal.append("a", parameters={"columns": [0]})
                b = journal.append("b", parameters={"columns": [1]})
                self.assertNotEqual(a, b)
                self.assertEqual(journal.requests, 2)
                journal.finish()
            with DataStore(self.root / "journal") as stored:
                families = []
                for ref in (a, b):
                    request = stored.native(ref, dependencies=(source,)).get_metadata(1, 0, "gadget_journal_request")
                    families.append(stored.native(request["outputs"][1], dependencies=(source,)))
                for family in families:
                    self.assertEqual(family.particular.values.tolist(), [Q(0)] * 4)
                    self.assertEqual(family.direction_count, 3)
                    self.assertEqual({row for row, _, q in family.directions.entries if q}, set(range(4)))
                self.assertNotEqual(families[0].directions.entries, families[1].directions.entries)
            self.assertEqual(verify_journal(self.root / "journal", sources=self.sources(source), datasets=(source,)),
                             {"rows": 2, "requests": 2})

    def test_pending_physical_metadata_is_frozen_and_replay_needs_owners(self):
        with DataStore(self.root / "source") as source:
            metadata = {"observed": True, "original": [1, 2]}
            with ProgramJournal(self.root / "journal", field_program(self.reference),
                                sources=self.sources(source), contributors=((source, self.reference),)) as journal:
                reference = journal.append("original", parameters={"columns": [0]}, occurrence=metadata)
                metadata["original"][0] = 99
                reference["state_digest"] = "0" * 64
                journal.finish()
                with self.assertRaises(ValueError):
                    journal.append("later", parameters={"columns": [0]})
            with DataStore(self.root / "journal") as stored:
                row, = journal_occurrences(stored, dependencies=(source,))
                self.assertEqual(row["occurrence"]["original"], [1, 2])
                self.assertNotEqual(row["request"]["state_digest"], "0" * 64)
            with self.assertRaises(KeyError):
                verify_journal(self.root / "journal", sources=self.sources(source), datasets=())


if __name__ == "__main__":
    unittest.main()
