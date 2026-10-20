#!/usr/bin/env python3
"""Check source and Git history for material that does not belong in a publication."""
import argparse
from pathlib import Path
import re
import subprocess
import sys
import zlib


LOCAL_PATH = re.compile(rb'/(?:home|Users)/[^/\s"\x27`<>]+/|~/(?:projects|dev)/')
WORK_NOTES = re.compile(rb'\b(?:worktree|handoff|session [ABC])\b', re.I)
DEVELOPMENT_REPORT = re.compile(rb'qualification\s+workstation|development\s+host|'
                                rb'(?:CPU|protocol|planning) qualification (?:on|used)', re.I)
PRIVATE_ROOTS = {'internal', 'private', 'research', 'verification', 'migration',
                 'data', 'models', 'results', 'archives', 'releases', 'experiments',
                 'build', 'dist', '.venv', '.rfc-source', '.rfc-download', '.llama-source'}
GENERATED_SUFFIXES = {'.whl', '.egg', '.o', '.obj', '.a', '.so', '.dylib', '.dll', '.exe',
                      '.rcdb', '.rcbd', '.rcblob', '.sqlite', '.sqlite3', '.npy', '.npz',
                      '.h5', '.hdf5', '.h5ad', '.db', '.zip', '.bundle'}


def path_findings(name):
    path = Path(name)
    if path.parts and path.parts[0] in PRIVATE_ROOTS:
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
            message = data.split(b'\n\n', 1)[1]
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
