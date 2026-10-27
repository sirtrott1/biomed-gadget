"""Reject local paths and unpublished Git objects without displaying their contents."""
import importlib.util
import gzip
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('publication', Path(__file__).resolve().parents[1] / 'tools/check_publication.py')
publication = importlib.util.module_from_spec(spec)
spec.loader.exec_module(publication)


class PublicationContracts(unittest.TestCase):
    def test_attribution_and_commit_trailers_are_rejected(self):
        values = (b'co' + b'dex', b'chat' + b'gpt', b'Co' + b'-author: Example', b'Co' + b'-authored-by: Example',
                  b'Signed' + b'-off-by: Example')
        for value in values:
            with self.subTest(value=value):
                self.assertGreater(publication.content_findings(value, 'README.md'), 0)
                self.assertGreater(publication.content_findings(b'value = "' + value + b'"', 'src/value.py'), 0)

    def test_private_document_references_and_development_labels_are_rejected(self):
        reference = b'internal' + b'/plans/report.md'
        label = b'Checkpoint ' + b'C' + b'03'
        for value in (reference, label):
            with self.subTest(value=value):
                self.assertGreater(publication.content_findings(value, 'README.md'), 0)
                self.assertGreater(publication.content_findings(b'# ' + value, 'src/method.py'), 0)
                self.assertGreater(publication.content_findings(b'// ' + value, 'src/method.cpp'), 0)
        self.assertGreater(publication.content_findings(b'path = "' + reference + b'"', 'src/method.py'), 0)
        self.assertGreater(publication.path_findings('examples/private/notes.md'), 0)

    def test_comment_and_document_style_is_checked(self):
        words = (b'lever' + b'age', b'seam' + b'lessly', b'world' + b'-class', b'power' + b'ful')
        for value in (*words, b'-' * 8, b'=' * 8, '\u2014'.encode()):
            with self.subTest(value=value):
                self.assertGreater(publication.content_findings(value, 'README.md'), 0)
                self.assertGreater(publication.content_findings(b'# ' + value, 'src/method.py'), 0)
                self.assertGreater(publication.content_findings(b'"""' + value + b'"""', 'src/method.py'), 0)
                self.assertGreater(publication.content_findings(b'/* ' + value + b' */', 'src/method.cpp'), 0)

    def test_scientific_identifiers_and_markdown_tables_remain_valid(self):
        values = ((b'Transcripts A1 and A2 retain their source coordinates.', 'README.md'),
                  ('BED3\u20136 describes the supported file columns.'.encode(), 'README.md'),
                  (b'| Field | Type |\n| --- | --- |\n| count | integer |\n', 'README.md'),
                  (b'# C0 and C1 name the declared chain spaces.\n', 'src/method.py'),
                  (b'// The B1 boundary retains its original orientation.\n', 'src/method.cpp'),
                  (b'const char *url = "https://example.invalid";', 'src/method.cpp'))
        for value, name in values:
            with self.subTest(value=value):
                self.assertEqual(publication.content_findings(value, name), 0)

    def test_removed_historical_attribution_and_message_labels_are_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def git(*args):
                return subprocess.check_output(['git', '-C', str(root), *args], stderr=subprocess.PIPE)
            git('init', '--quiet', '--initial-branch=main')
            git('config', 'user.name', 'Publication fixture')
            git('config', 'user.email', 'fixture@example.invalid')
            path = root / 'README.md'
            path.write_bytes(b'co' + b'dex\n')
            git('add', 'README.md')
            git('commit', '--quiet', '-m', 'Describe the library')
            path.write_text('# Library\n')
            git('commit', '--quiet', '-am', 'Describe library inputs')
            self.assertEqual(publication.check(root), 0)
            self.assertGreater(publication.check(root, history=True), 0)
            git('checkout', '--quiet', '--orphan', 'clean')
            git('add', 'README.md')
            git('commit', '--quiet', '-m', 'Retain library inputs')
            git('branch', '-D', 'main')
            self.assertEqual(publication.check(root, history=True), 0)
            git('commit', '--quiet', '--allow-empty', '-m', 'Update checkpoint ' + 'E' + '1')
            self.assertGreater(publication.check(root, history=True), 0)

    def test_historical_assistant_identity_is_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def git(*args):
                subprocess.run(['git', '-C', str(root), *args], check=True, capture_output=True)
            git('init', '--quiet')
            git('config', 'user.name', 'co' + 'dex')
            git('config', 'user.email', 'fixture@example.invalid')
            (root / 'README.md').write_text('# Library\n')
            git('add', 'README.md')
            git('commit', '--quiet', '-m', 'Describe the library')
            self.assertGreater(publication.check(root, history=True), 0)

    def test_source_and_historical_content(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def git(*args):
                subprocess.run(['git', '-C', str(root), *args], check=True, capture_output=True)
            git('init', '--quiet')
            git('config', 'user.name', 'Publication fixture')
            git('config', 'user.email', 'fixture@example.invalid')
            path = root / 'README.md'
            path.write_text('# Library\n')
            git('add', 'README.md')
            git('commit', '--quiet', '-m', 'Add library guide')
            self.assertEqual(publication.check(root, history=True), 0)
            path.write_text('/'.join(('', 'home', 'example', 'source.tsv')))
            self.assertGreater(publication.check(root), 0)
            git('commit', '--quiet', '-am', 'Update library guide')
            path.write_text('# Library\n')
            git('commit', '--quiet', '-am', 'Describe library inputs')
            self.assertEqual(publication.check(root), 0)
            self.assertGreater(publication.check(root, history=True), 0)

    def test_publication_paths_and_project_credit(self):
        self.assertGreater(publication.path_findings('internal/notes.txt'), 0)
        self.assertGreater(publication.path_findings('.github/README.md'), 0)
        self.assertGreater(publication.path_findings('docs/build-review.json'), 0)
        self.assertGreater(publication.path_findings('tests/fixtures/result.rfcxs.gz'), 0)
        self.assertEqual(publication.path_findings('examples/native-data/README.md'), 0)
        self.assertEqual(publication.content_findings(b'Copyright 2026 biomed-gadget\n', 'NOTICE'), 0)
        self.assertGreater(publication.content_findings(b'Copyright 2026 Example Person\n', 'NOTICE'), 0)
        self.assertGreater(publication.content_findings(b'Measurements from the development host.', 'README.md'), 0)

    def test_compressed_content(self):
        local = '/'.join(('', 'home', 'example', 'source.tsv')).encode()
        self.assertGreater(publication.content_findings(gzip.compress(local)), 0)
        self.assertEqual(publication.content_findings(gzip.compress(b'declared test input')), 0)
        self.assertGreater(publication.content_findings(gzip.compress(b'a' * ((16 << 20) + 1))), 0)

    def test_generated_files_and_native_data(self):
        generated = ('AGENTS.md', 'experiments/notes.md', '.llama-source/CMakeLists.txt',
                     'tools/rcdb/dist/package.whl', 'tools/rcdb/build/lib/module.py',
                     'examples/counts.native/catalog', 'tools/store.sqlite3',
                     'tools/selected.rcbd', 'tools/measurements.rcdb', 'tools/counts.npy')
        for name in generated:
            with self.subTest(name=name):
                self.assertGreater(publication.path_findings(name), 0)
        for name in ('examples/measurements.csv', 'examples/scientific/coverage.rfcxs',
                     'src/native_record.cpp'):
            with self.subTest(name=name):
                self.assertEqual(publication.path_findings(name), 0)

    def test_generated_outputs_are_ignored(self):
        root = Path(__file__).resolve().parents[1]
        names = ('counts.native/catalog', 'measurements.rcdb', 'selected.rcbd',
                 'result.rfcxs.gz', 'cache.rcblob', 'store.sqlite3', 'store.sqlite3-wal',
                 'counts.npy', '.llama-source/README.md', 'tools/rcdb/dist/package.whl')
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run(['git', 'init', '--quiet', directory], check=True, capture_output=True)
            (Path(directory) / '.gitignore').write_bytes((root / '.gitignore').read_bytes())
            result = subprocess.check_output(['git', '-C', directory, 'check-ignore',
                                               '--no-index', '--stdin', '-z'],
                                              input=b'\0'.join(name.encode() for name in names) + b'\0')
            self.assertEqual(set(result.decode().strip('\0').split('\0')), set(names))

    def test_excluded_data_is_rejected_without_reading_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(['git', 'init', '--quiet', directory], check=True, capture_output=True)
            (root / 'measurements.npy').write_bytes(b'excluded data')
            with patch.object(Path, 'read_bytes', side_effect=AssertionError('Excluded data was read')):
                self.assertGreater(publication.check(root), 0)


if __name__ == '__main__':
    unittest.main()
