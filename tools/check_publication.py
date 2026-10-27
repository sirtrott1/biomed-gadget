#!/usr/bin/env python3
"""Check source and Git history for material that does not belong in a publication."""
import argparse
import ast
import io
from pathlib import Path
import re
import subprocess
import sys
import tokenize
import zlib


LOCAL_PATH = re.compile(rb'/(?:home|Users)/[^/\s"\x27`<>]+/|~/(?:projects|dev)/')
WORK_NOTES = re.compile(rb'\b(?:worktree|handoff|session [ABC])\b', re.I)
DEVELOPMENT_REPORT = re.compile(rb'qualification\s+workstation|development\s+host|'
                                rb'(?:CPU|protocol|planning) qualification (?:on|used)', re.I)
ATTRIBUTION = re.compile(rb'\b(?:co[d]ex|chat[g]pt|co[- ]?author(?:ed)?(?:[- ]?by)?|signed[- ]?off[- ]?by)\b', re.I)
PRIVATE_REFERENCE = re.compile(rb'\b(?:internal|private|experiments)[/\\][^\s\x22\x27)`<>]+', re.I)
DEVELOPMENT_LABEL = re.compile(rb'\b(?:batch|checkpoint|milestone|phase|handoff|session|implementation|review|report)'
                               rb'\s*(?:[:#/-]\s*)?[A-FQ][0-9]{1,3}\b', re.I)
PROMOTIONAL = re.compile(rb'\b(?:leverage|delve|seamless(?:ly)?|empower(?:ing)?|supercharg\w*|'
                         rb'game.chang\w*|world.class|cutting.edge|state.of.the.art|groundbreaking|'
                         rb'revolutionary|best.in.class|transformative|effortlessly|holistic|synerg\w*|'
                         rb'harness|foster|comprehensive|elegant|powerful)\b', re.I)
DIVIDER = re.compile(rb'^[ \t]*(?:(?://|#|/\*|\*|<!--)[ \t]*)?[-=_*]{3,}(?:[ \t].*)?$', re.M)
PRIVATE_ROOTS = {'internal', 'private', 'research', 'verification', 'migration',
                 'data', 'models', 'results', 'archives', 'releases', 'experiments',
                 'build', 'dist', '.venv', '.rfc-source', '.rfc-download', '.llama-source'}
GENERATED_SUFFIXES = {'.whl', '.egg', '.o', '.obj', '.a', '.so', '.dylib', '.dll', '.exe',
                      '.rcdb', '.rcbd', '.rcblob', '.sqlite', '.sqlite3', '.npy', '.npz',
                      '.h5', '.hdf5', '.h5ad', '.db', '.zip', '.bundle'}


def prose(data, name):
    """Select documentation and comments while preserving technical identifiers."""
    suffix = Path(name).suffix.lower()
    if not name or suffix in ('.md', '.rst', '.txt', '.cff'):
        return data
    if suffix == '.py':
        text = data.decode('utf-8', errors='replace')
        parts = []
        try:
            parts.extend(token.string for token in tokenize.generate_tokens(io.StringIO(text).readline)
                         if token.type == tokenize.COMMENT)
            tree = ast.parse(text)
            parts.extend(value for node in ast.walk(tree)
                         if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                         and (value := ast.get_docstring(node)))
        except (SyntaxError, tokenize.TokenError, IndentationError):
            parts.extend(re.findall(r'(?m)^\s*#.*$', text))
        return '\n'.join(parts).encode()
    if suffix in ('.cpp', '.hpp', '.c', '.h'):
        tokens = re.findall(rb'"(?:\\.|[^"\\])*"|\x27(?:\\.|[^\x27\\])*\x27|//[^\n]*|/\*.*?\*/',
                            data, re.S)
        return b'\n'.join(token for token in tokens if token.startswith((b'//', b'/*')))
    return b'\n'.join(re.findall(rb'(?m)^\s*#.*$', data))


def path_findings(name):
    path = Path(name)
    if (path.parts and path.parts[0] in PRIVATE_ROOTS or
            any(part.lower() in {'internal', 'private', 'experiments', 'verification', 'migration'}
                for part in path.parts)):
        return 1
    if (path.name == 'AGENTS.md' or any(part in ('__pycache__', '.cache', 'build', 'dist', 'CMakeFiles')
                                      or part.endswith(('.native', '.rcdb', '.rcbd', '.egg-info'))
                                      for part in path.parts)
            or any(suffix in GENERATED_SUFFIXES for suffix in path.suffixes)):
        return 1
    if path.parent == Path('docs') and path.suffix == '.json':
        return 1
    if path.parent == Path('.github') and path.name.lower().startswith('readme.'):
        return 1
    if name.startswith('tests/fixtures/') and any(suffix in ('.rfcxs', '.rcbd', '.rcdb') for suffix in path.suffixes):
        return 1
    return 0


def content_findings(data, name=''):
    if data.startswith(b'\x1f\x8b'):
        decoder = zlib.decompressobj(31)
        try:
            expanded = decoder.decompress(data, (16 << 20) + 1)
        except zlib.error:
            return 1
        if not decoder.eof or decoder.unused_data or len(expanded) > 16 << 20:
            return 1
        return content_findings(expanded, name.removesuffix('.gz'))
    count = int(LOCAL_PATH.search(data) is not None)
    count += int(ATTRIBUTION.search(data) is not None)
    text = prose(data, name)
    # Test source can contain synthetic exclusion inputs. Its comments and
    # documentation follow the same publication rules as other source files.
    references = text if name.startswith('tests/') else data
    if name != '.gitignore':
        count += int(PRIVATE_REFERENCE.search(references) is not None)
    count += int(DEVELOPMENT_LABEL.search(text) is not None)
    count += int(PROMOTIONAL.search(text) is not None)
    count += int(DIVIDER.search(text) is not None)
    count += int('\u2014' in text.decode('utf-8', errors='replace'))
    if name.endswith(('.md', '.rst')):
        count += int(DEVELOPMENT_REPORT.search(data) is not None)
    if name in ('NOTICE', 'tools/rcdb/NOTICE'):
        notices = [line for line in data.splitlines() if line.startswith(b'Copyright ')]
        count += int(not notices or any(not re.fullmatch(rb'Copyright \d{4}(?:-\d{4})? biomed-gadget', line)
                                       for line in notices))
    if name == 'CITATION.cff':
        count += int(b'  - name: biomed-gadget\n' not in data or b'given-names:' in data or b'family-names:' in data)
    return count


def check(root, git_dir=None, history=False):
    git = (['git', '--git-dir=' + str(git_dir), '--work-tree=' + str(root)]
           if git_dir else ['git', '-C', str(root)])
    count = 0
    if not history:
        names = subprocess.check_output(git + ['ls-files', '-z', '--cached', '--others', '--exclude-standard']).decode().split('\0')
        for name in names:
            path = root / name
            if not name or not path.is_file():
                continue
            excluded = path_findings(name)
            count += excluded
            if not excluded:
                count += content_findings(path.read_bytes(), name)
        return count
    names = subprocess.check_output(git + ['log', '--all', '--root', '-m', '--format=', '--name-only', '-z']).decode().split('\0')
    count += sum(path_findings(name.strip()) for name in set(names) if name.strip())
    if count:
        return count
    objects = subprocess.check_output(git + ['rev-list', '--objects', '--all']).splitlines()
    reader = subprocess.Popen(git + ['cat-file', '--batch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for entry in objects:
        oid, _, name = entry.partition(b' ')
        reader.stdin.write(oid + b'\n')
        reader.stdin.flush()
        header = reader.stdout.readline().split()
        data = reader.stdout.read(int(header[2]))
        assert reader.stdout.read(1) == b'\n'
        if header[1] == b'blob':
            count += content_findings(data, name.decode())
        elif header[1] == b'commit':
            metadata, message = data.split(b'\n\n', 1)
            count += int(ATTRIBUTION.search(metadata) is not None)
            count += content_findings(message) + int(WORK_NOTES.search(message) is not None)
    reader.stdin.close()
    reader.stdout.close()
    if reader.wait():
        raise RuntimeError('Cannot read selected Git history')
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--git-dir', type=Path)
    parser.add_argument('--history', action='store_true')
    args = parser.parse_args()
    count = check(Path(__file__).resolve().parents[1], args.git_dir, args.history)
    if count:
        print(f'Publication check failed: {count} findings. Review the source privately.', file=sys.stderr)
        return 1
    print('PASS: publication source' + (' and history' if args.history else ''))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
