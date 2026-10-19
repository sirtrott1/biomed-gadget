"""Reject local paths and unpublished Git objects without displaying their contents."""
import importlib.util
import gzip
from pathlib import Path
import subprocess
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('publication', Path(__file__).resolve().parents[1] / 'tools/check_publication.py')
publication = importlib.util.module_from_spec(spec)
spec.loader.exec_module(publication)


class PublicationContracts(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
