#!/usr/bin/env python3
"""Export public GADGET source with available Git history and verify its bytes."""
import argparse
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('publication', ROOT / 'tools/check_publication.py')
publication = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(publication)
PUBLIC_DIRECTORIES = ('.github', 'apps', 'cmake', 'examples', 'include', 'schemas', 'src', 'tests', 'tools')
PUBLIC_FILES = ('.clang-format', '.gitignore', 'CMakeLists.txt', 'CMakePresets.json', 'README.md', 'DEPENDENCIES.md', 'LICENSE', 'NOTICE', 'CITATION.cff')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1024 * 1024): h.update(chunk)
    return h.hexdigest()

def verify_snapshot(selection, files):
    head = subprocess.run([*selection, 'rev-parse', '--verify', 'HEAD'],
                          capture_output=True, text=True)
    if head.returncode:
        return
    entries = subprocess.check_output([*selection, 'ls-tree', '-rz',
                                       head.stdout.strip()]).split(b'\0')
    tracked = {}
    for entry in entries:
        if entry:
            attributes, name = entry.split(b'\t', 1)
            mode, kind, oid = attributes.decode().split()
            if kind != 'blob':
                raise ValueError('Source history contains an unsupported external object')
            tracked[name.decode()] = (mode, oid)
    selected = {path.relative_to(ROOT).as_posix() for path in files}
    if selected != tracked.keys():
        raise ValueError('Source selection differs from the selected commit; commit public files '
                         'before exporting history, or use --without-history')
    algorithm = subprocess.check_output([*selection, 'rev-parse', '--show-object-format'],
                                        text=True).strip()
    for path in files:
        mode, expected = tracked[path.relative_to(ROOT).as_posix()]
        actual = hashlib.new(algorithm, b'blob ' + str(path.stat().st_size).encode() + b'\0')
        with path.open('rb') as stream:
            while chunk := stream.read(1024 * 1024):
                actual.update(chunk)
        executable = bool(path.stat().st_mode & 0o100)
        if (actual.hexdigest() != expected or mode not in ('100644', '100755')
                or executable != (mode == '100755')):
            raise ValueError('Source differs from the selected commit; commit public changes '
                             'before exporting history, or use --without-history')

def object_inventory(git_dir):
    command = ['git', '--git-dir', str(git_dir), 'cat-file', '--batch-all-objects',
               '--batch-check=%(objectname) %(objecttype) %(objectsize)']
    h = hashlib.sha256()
    with subprocess.Popen(command, stdout=subprocess.PIPE) as process:
        while chunk := process.stdout.read(1024 * 1024):
            h.update(chunk)
        process.stdout.close()
        if process.wait():
            raise RuntimeError('Cannot inventory the complete Git object store')
    return h.digest()

def history(directory, git_dir=None):
    command = ['git', '--git-dir', str(git_dir)] if git_dir else ['git', '-C', str(ROOT)]
    head = subprocess.run([*command, 'rev-parse', '--verify', 'HEAD'],
                          capture_output=True, text=True)
    if head.returncode:
        if git_dir:
            raise ValueError('Selected Git history has no readable HEAD')
        print('Git history is unavailable; exporting a source snapshot')
        return None
    original = subprocess.check_output([*command, 'rev-parse', '--absolute-git-dir'], text=True).strip()
    target = directory / 'history.git'
    subprocess.run(['git', 'init', '--quiet', '--bare', str(target)], check=True)
    refs = ['for-each-ref', '--format=%(refname) %(objectname)', 'refs/heads/', 'refs/tags/']
    selected = subprocess.check_output([*command, *refs], text=True).splitlines()
    refspecs = [f'{name}:{name}' for name, _ in (line.split() for line in selected)]
    subprocess.run(['git', '--git-dir', str(target), 'fetch', '--quiet', '--no-tags',
                    '--no-write-fetch-head', original, *refspecs, head.stdout.strip()], check=True)
    symbolic = subprocess.run([*command, 'symbolic-ref', '--quiet', 'HEAD'],
                              capture_output=True, text=True)
    if symbolic.returncode == 0:
        subprocess.run(['git', '--git-dir', str(target), 'symbolic-ref', 'HEAD',
                        symbolic.stdout.strip()], check=True)
    else:
        subprocess.run(['git', '--git-dir', str(target), 'update-ref', '--no-deref',
                        'HEAD', head.stdout.strip()], check=True)
    # Ref journals, temporary packs and unreachable objects stay outside exports.
    inventory = object_inventory(target)
    subprocess.run(['git', '--git-dir', str(target), '-c', 'pack.threads=2',
                    'repack', '-a', '-d'], check=True,
                   stdout=subprocess.DEVNULL)
    for temporary in (target / 'objects/pack').glob('tmp_pack_*'):
        if temporary.is_file():
            temporary.unlink()
    assert object_inventory(target) == inventory, 'Selected Git objects changed'
    # Transport configuration has no machine local remotes or executable hooks.
    (target / 'config').write_text('[core]\n\trepositoryformatversion = 0\n'
                                  '\tfilemode = true\n\tbare = false\n'
                                  '\tlogallrefupdates = true\n')
    shutil.rmtree(target / 'hooks', ignore_errors=True)
    assert not (target / 'objects/info/alternates').exists(), 'External Git objects'
    assert subprocess.check_output(['git', '--git-dir', str(target), 'rev-parse',
                                    'HEAD'], text=True).strip() == head.stdout.strip(), 'Git HEAD changed'
    subprocess.run(['git', '--git-dir', str(target), '--work-tree', str(ROOT),
                    'read-tree', 'HEAD'], check=True)
    subprocess.run(['git', '--git-dir', str(target), 'fsck', '--full'], check=True,
                   stdout=subprocess.DEVNULL)
    if publication.check(ROOT, target, history=True):
        raise ValueError('Selected history failed the publication check')
    assert subprocess.check_output([*command, *refs]) == \
           subprocess.check_output(['git', '--git-dir', str(target), *refs]), 'Git refs changed'
    count = subprocess.check_output(['git', '--git-dir', str(target), 'rev-list',
                                     '--all', '--count'], text=True).strip()
    print(f'Git history: {count} commits; HEAD {head.stdout.strip()}')
    return target

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--build', type=Path)
    parser.add_argument('--rexgraph', '--sdk', dest='rexgraph', type=Path,
                        help='External RexGraphCpp repository or installed prefix')
    parser.add_argument('--consumed-only', action='store_true',
                        help='Verify the canonical headers actually used by this build')
    parser.add_argument('--without-history', action='store_true',
                        help='Export a source snapshot without Git history')
    parser.add_argument('--git-dir', type=Path,
                        help='Use an explicit Git metadata directory for history and public file selection')
    args = parser.parse_args()
    if args.out.exists(): parser.error('The output archive must be new')
    if args.git_dir and args.without_history: parser.error('--git-dir and --without-history are separate export modes')
    if args.git_dir: args.git_dir = args.git_dir.resolve(strict=True)
    if bool(args.build) != bool(args.rexgraph): parser.error('Supply --build and --rexgraph together')
    if args.consumed_only and not args.build: parser.error('--consumed-only requires --build')
    if args.build:
        command = [sys.executable, str(ROOT / 'tools/check_dependency.py'),
                   '--build', str(args.build), '--rexgraph', str(args.rexgraph)]
        if args.consumed_only: command.append('--consumed-only')
        subprocess.run(command, check=True)
    files = [ROOT / name for name in PUBLIC_FILES if (ROOT / name).is_file()]
    for directory in PUBLIC_DIRECTORIES:
        if (ROOT / directory).exists():
            files.extend(path for path in (ROOT / directory).rglob('*') if path.is_file()
                         and not any(part in ('__pycache__', '.git', '.cache') for part in path.parts)
                         and path.name != 'AGENTS.md' and path.suffix not in ('.pyc', '.pyo'))
    selection = (['git', '--git-dir', str(args.git_dir), '--work-tree', str(ROOT)]
                 if args.git_dir else ['git', '-C', str(ROOT)])
    candidates = subprocess.run([*selection, 'ls-files', '-z',
                                 '--cached', '--others', '--exclude-standard'],
                                capture_output=True, text=True)
    if candidates.returncode == 0:
        eligible = set(candidates.stdout.split('\0'))
        tracked = subprocess.run([*selection, 'ls-tree', '-rz', '--name-only', 'HEAD'],
                                 capture_output=True, text=True)
        if tracked.returncode == 0:
            eligible.update(tracked.stdout.split('\0'))
        files = [path for path in files if path.relative_to(ROOT).as_posix() in eligible]
    elif args.git_dir:
            raise ValueError('Selected Git repository cannot enumerate the public source')
    assert all(not path.is_symlink() for path in files), 'Public source symlinks require an explicit export declaration'
    files = sorted(set(files))
    if any(publication.path_findings(path.relative_to(ROOT).as_posix()) or
           publication.content_findings(path.read_bytes(), path.relative_to(ROOT).as_posix())
           for path in files):
        raise ValueError('Selected source failed the publication check')
    if not args.without_history:
        verify_snapshot(selection, files)
    inventory = {'gadget/' + path.relative_to(ROOT).as_posix(): digest(path) for path in files}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    staging = ROOT / 'build'
    staging.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='source-export-', dir=staging) as directory:
        git = None if args.without_history else history(Path(directory), args.git_dir)
        histories = sorted(path for path in git.rglob('*') if path.is_file()) if git else []
        assert all(not path.is_symlink() for path in histories), 'Unsupported Git symlink'
        inventory.update({'gadget/.git/' + path.relative_to(git).as_posix(): digest(path)
                          for path in histories})
        with zipfile.ZipFile(args.out, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for path in files: archive.write(path, 'gadget/' + path.relative_to(ROOT).as_posix())
            for path in histories: archive.write(path, 'gadget/.git/' + path.relative_to(git).as_posix())
            if git:
                archive.write(git / 'refs', 'gadget/.git/refs/')
            archive.writestr('SHA256SUMS', ''.join(f'{sha}  {name}\n' for name, sha in inventory.items()))
    with zipfile.ZipFile(args.out) as archive:
        for name, expected in inventory.items():
            h = hashlib.sha256()
            with archive.open(name) as stream:
                while chunk := stream.read(1024 * 1024): h.update(chunk)
            assert h.hexdigest() == expected, name
    print(f'Verified {len(files)} public files: {args.out}')
    print(f'SHA256 {digest(args.out)}')

if __name__ == '__main__': main()
