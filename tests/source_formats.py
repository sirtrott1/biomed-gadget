"""Compare source visitors with independent CSV, decimal and byte coordinates."""
import csv
from decimal import Decimal
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class SourceFormats(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def source(self, name, contents):
        path = self.root / name
        path.write_bytes(contents.encode())
        return path

    def run_reference(self, script, *args):
        result = subprocess.run([sys.executable, str(Path(__file__).with_name(script)), *map(str, args)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_bindingdb_occurrences_and_affinities(self):
        header = ['BindingDB Reactant_set_id', 'BindingDB MonomerID', 'Target Name',
                  'Target Source Organism According to Curator or DataSource',
                  'Ki (nM)', 'Kd (nM)', 'IC50 (nM)', 'EC50 (nM)']
        rows = [['9007199254740993', '7', 'reported construct', 'reported organism',
                 '0', '', '>= 0.000000000000000000000000001', '~1e100'],
                ['9007199254740993', '7', 'reported construct', 'reported organism',
                 '0', '', '>= 0.000000000000000000000000001', '~1e100'],
                ['3', '8', 'other construct', 'other organism', '', '-2', '1e-400', '']]
        source = self.source('binding.tsv', '\t'.join(header) + '\r\n' +
                             ''.join('\t'.join(row) + '\r\n' for row in rows))
        self.run_reference('bindingdb_bulk_reference.py', '--reader', BINDINGDB,
                           '--source', source, '--source-sha256',
                           hashlib.sha256(source.read_bytes()).hexdigest(),
                           '--out', self.root / 'binding-report.json')

    def test_pubchem_original_typed_values_and_links(self):
        assay = self.source('assay.json', '''{"PC_AssaySubmit":{"assay":{"descr":{
          "aid":{"id":100,"version":1},"name":"authored source","results":[
          {"tid":1,"name":"signal","type":1},{"tid":2,"name":"count","type":2},
          {"tid":3,"name":"flag","type":3},{"tid":4,"name":"annotation","type":4}]}},
          "data":[{"sid":100,"outcome":1,"data":[{"tid":1,"value":{"fval":-0.000e+00}},
          {"tid":2,"value":{"ival":9007199254740993}},{"tid":3,"value":{"bval":false}},
          {"tid":4,"value":{"sval":""}}]},{"sid":100,"data":[
          {"tid":1,"value":{"fval":1e309}},{"tid":2,"value":{"ival":-0}}]},
          {"sid":101,"data":[]}],"revoke":[101,101]}}''')
        concise = self.root / 'concise.csv'
        with concise.open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['AID', 'SID', 'CID', 'Target GeneID', 'Target Accession', 'Note'])
            writer.writerows([['100', '100', '100', '1', '', 'a\r\nb'],
                              ['100', '100', '200', '', 'NP_1.1', ''],
                              ['100', '100', '200', '', '', ''],
                              ['100', '101', '', '', '', '']])
        self.run_reference('pubchem_assay_reference.py', '--reader', PUBCHEM,
                           '--assay', assay, '--concise', concise, '--aid', '100',
                           '--assay-sha256', 'sha256:' + hashlib.sha256(assay.read_bytes()).hexdigest(),
                           '--concise-sha256', 'sha256:' + hashlib.sha256(concise.read_bytes()).hexdigest(),
                           '--out', self.root / 'pubchem-report')

    def test_gct_coordinates_missing_zero_and_repeated_donors(self):
        samples = ['GTEX-AA-0001', 'GTEX-AA-0002', 'GTEX-BB-0001']
        rows = [['ENSG0001.1', 'same label', '0', 'NA', '2'],
                ['ENSG0002.1', 'same label', '1', '3', '0'],
                ['ENSG0003.1', 'last', 'NA', '0', '9007199254740993']]
        source = self.source('matrix.gct', '#1.2\r\n3\t3\r\nName\tDescription\t' +
                             '\t'.join(samples) + '\r\n' +
                             ''.join('\t'.join(row) + '\r\n' for row in rows))
        attributes = self.source('samples.tsv', 'SAMPID\tSMTSD\tannotation\n' +
                                 ''.join(sample + '\tdeclared tissue\tNA\n' for sample in samples))
        subjects = self.source('subjects.tsv', 'SUBJID\tSEX\tAGE\nGTEX-AA\t1\t20-29\nGTEX-BB\t2\tNA\n')
        result = subprocess.run([GCT, '--dump', str(source), 'gct12', 'gene', 'integer_counts',
                                 str(attributes), str(subjects)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        actual = list(csv.reader(io.StringIO(result.stdout), delimiter='\t', quoting=csv.QUOTE_NONE))
        self.assertEqual(actual[0], ['H', '3', '3', *samples])
        occurrences = [row for row in actual if row[0] == 'R']
        self.assertEqual(len(occurrences), len(rows))
        lines = source.read_bytes().splitlines(keepends=True)
        offset = sum(map(len, lines[:3]))
        for index, (reported, original) in enumerate(zip(occurrences, rows, strict=True)):
            self.assertEqual(list(map(int, reported[1:5])), [index, index + 4, offset, len(lines[index + 3])])
            self.assertEqual(reported[5:7], original[:2])
            self.assertEqual(len(reported[7:]), len(samples))
            for token, expected in zip(reported[7:], original[2:], strict=True):
                support, value = token.split(':', 1)
                self.assertEqual(value, expected)
                self.assertEqual(support, '0' if expected == 'NA' else '1')
                if support == '1':
                    self.assertEqual(Decimal(value), Decimal(expected))
            offset += len(lines[index + 3])
        self.assertEqual(offset, source.stat().st_size)
        self.assertEqual(actual[-1], ['Z', '3', '3', '7', '2', '3'])
        joined = [row for row in actual if row[0] == 'S']
        self.assertEqual([row[1] for row in joined], samples)
        self.assertEqual([row[2] for row in joined], ['GTEX-AA', 'GTEX-AA', 'GTEX-BB'])
        self.assertEqual([row[4] for row in joined], ['1', '1', '0'])


if __name__ == '__main__':
    BINDINGDB, PUBCHEM, GCT = map(str, map(Path, sys.argv[1:4]))
    unittest.main(argv=sys.argv[:1])
