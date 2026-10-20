"""Verify public source selection excludes ignored files inside source directories."""
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

root = Path(sys.argv[1]).resolve()
shutil.rmtree(root, ignore_errors=True)
project = root / 'project'
(project / 'tools').mkdir(parents=True)
source = Path(__file__).resolve().parents[1]
shutil.copy2(source / 'tools/export_source.py', project / 'tools/export_source.py')
shutil.copy2(source / 'tools/check_publication.py', project / 'tools/check_publication.py')
(project / '.gitignore').write_text('tools/generated/\n*.rfcxs\n')
(project / 'README.md').write_text('# Export fixture\n')
(project / 'tools/method.py').write_text('value = 1\n')
(project / 'tools/method.py').chmod(0o755)
(project / 'tools/generated').mkdir()
(project / 'tools/generated/private.txt').write_text('Excluded generated output\n')
(project / 'tools/result.rfcxs').write_bytes(b'Excluded binary output')
(project / 'internal').mkdir()
(project / 'internal/notes.txt').write_text('Excluded workspace notes\n')
subprocess.run(['git', 'init', '--quiet', str(project)], check=True)
subprocess.run(['git', '-C', str(project), 'add', '.gitignore', 'README.md'], check=True)
archive = root / 'source.zip'
subprocess.run([sys.executable, str(project / 'tools/export_source.py'), '--out',
                str(archive), '--without-history'], check=True, capture_output=True, text=True)
with zipfile.ZipFile(archive) as output:
    assert set(output.namelist()) == {'gadget/.gitignore', 'gadget/README.md',
        'gadget/tools/method.py', 'gadget/tools/export_source.py',
        'gadget/tools/check_publication.py', 'SHA256SUMS'}
    assert output.read('gadget/tools/method.py') == (project / 'tools/method.py').read_bytes()
    assert (output.getinfo('gadget/tools/method.py').external_attr >> 16) & 0o777 == 0o755
print('Public source export contracts passed')

metadata = root / 'selected.git'
subprocess.run(['git', 'init', '--quiet', '--bare', str(metadata)], check=True)
selected = ['git', '--git-dir', str(metadata), '--work-tree', str(project)]
subprocess.run([*selected, 'add', '.gitignore', 'README.md', 'tools/method.py',
                'tools/export_source.py', 'tools/check_publication.py'], check=True)
subprocess.run([*selected, '-c', 'user.name=Export contract', '-c', 'user.email=export@example.invalid',
                'commit', '--quiet', '-m', 'Retain public export fixture'], check=True)
expected_head = subprocess.check_output([*selected, 'rev-parse', 'HEAD']).strip()
orphan_bytes = b'Unreferenced export fixture bytes\n'
orphan = subprocess.check_output([*selected, 'hash-object', '-w', '--stdin'], input=orphan_bytes).strip()
private_tree = subprocess.check_output([*selected, 'mktree'],
    input=b'100644 blob ' + orphan + b'\tlocal.txt\n').strip()
private_commit = subprocess.check_output([*selected, '-c', 'user.name=Export contract',
    '-c', 'user.email=export@example.invalid', 'commit-tree', private_tree.decode()],
    input=b'Local fixture state\n').strip()
subprocess.run([*selected, 'update-ref', 'refs/stash', private_commit.decode()], check=True)
temporary_pack = metadata / 'objects/pack/tmp_pack_fixture'
temporary_pack.write_bytes(b'Unfinished local pack fixture')
subprocess.run([*selected, 'branch', 'retained-fixture'], check=True)
private_note = 'fetch /private/workspaces/source-fixture'
journal = metadata / 'logs/refs/heads/retained-fixture'
journal.parent.mkdir(parents=True, exist_ok=True)
journal_bytes = (f'{"0" * len(expected_head)} {expected_head.decode()} '
                 f'Export contract <export@example.invalid> 1 +0000\t{private_note}\n').encode()
journal.write_bytes(journal_bytes)
full = root / 'history.zip'
subprocess.run([sys.executable, str(project / 'tools/export_source.py'), '--out', str(full),
                '--git-dir', str(metadata)], check=True, capture_output=True, text=True)
extracted = root / 'extracted'
with zipfile.ZipFile(full) as output:
    assert 'gadget/.git/HEAD' in output.namelist()
    assert not any('tmp_pack_' in name for name in output.namelist())
    assert not any(name.startswith('gadget/.git/logs/') for name in output.namelist())
    assert not any('generated' in name or 'internal/' in name or 'result.rfcxs' in name for name in output.namelist())
    output.extractall(extracted)
    for member in output.infolist():
        if not member.is_dir():
            (extracted / member.filename).chmod((member.external_attr >> 16) & 0o777)
received = extracted / 'gadget'
assert subprocess.check_output(['git', '-C', str(received), 'rev-parse', 'HEAD']).strip() == expected_head
assert subprocess.check_output(['git', '-C', str(received), 'rev-parse', 'retained-fixture']).strip() == expected_head
for object_id in (orphan, private_commit):
    assert subprocess.run(['git', '-C', str(received), 'cat-file', '-e', object_id.decode()],
                          capture_output=True).returncode != 0
assert subprocess.run(['git', '-C', str(received), 'rev-parse', '--verify', 'refs/stash'],
                      capture_output=True).returncode != 0
subprocess.run(['git', '-C', str(received), 'fsck', '--full', '--strict'], check=True, capture_output=True)
assert subprocess.check_output(['git', '-C', str(received), 'status', '--porcelain']) == b''
assert (received / 'tools/method.py').read_bytes() == (project / 'tools/method.py').read_bytes()
assert (received / 'tools/method.py').stat().st_mode & 0o777 == 0o755
assert temporary_pack.read_bytes() == b'Unfinished local pack fixture'
assert journal.read_bytes() == journal_bytes
assert subprocess.check_output([*selected, 'cat-file', 'blob', orphan.decode()]) == orphan_bytes
assert subprocess.check_output([*selected, 'rev-parse', 'refs/stash']).strip() == private_commit
before = full.read_bytes()
refused = subprocess.run([sys.executable, str(project / 'tools/export_source.py'), '--out', str(root / 'conflict.zip'),
    '--git-dir', str(metadata), '--without-history'], capture_output=True)
assert refused.returncode and not (root / 'conflict.zip').exists() and full.read_bytes() == before
print('Published branches, extracted index and exclusion of local refs and objects passed')

# A bare repository can have a valid commit and no index. Its tracked examples
# must survive source selection even when their extensions are ignored.
(project / 'examples').mkdir()
example = project / 'examples/reads.fastq'
example.write_bytes(b'@read\nA\n+\n!\n')
with (project / '.gitignore').open('a') as stream:
    stream.write('*.fastq\n')
subprocess.run([*selected, 'add', '.gitignore'], check=True)
subprocess.run([*selected, 'add', '-f', 'examples/reads.fastq'], check=True)
subprocess.run([*selected, '-c', 'user.name=Export contract',
                '-c', 'user.email=export@example.invalid', 'commit', '--quiet',
                '-m', 'Retain example input'], check=True)
updated_head = subprocess.check_output([*selected, 'rev-parse', 'HEAD']).strip()
without_index = root / 'without-index.git'
subprocess.run(['git', 'init', '--quiet', '--bare', '--initial-branch=main', str(without_index)], check=True)
subprocess.run(['git', '--git-dir', str(without_index), 'fetch', '--quiet', '--no-tags',
                '--no-write-fetch-head', str(metadata), updated_head.decode() + ':refs/heads/main'], check=True)
assert not (without_index / 'index').exists()
complete = root / 'complete.zip'
subprocess.run([sys.executable, str(project / 'tools/export_source.py'), '--out', str(complete),
                '--git-dir', str(without_index)], check=True, capture_output=True, text=True)
restored = root / 'complete'
with zipfile.ZipFile(complete) as output:
    assert output.read('gadget/examples/reads.fastq') == example.read_bytes()
    assert not any('internal/' in name or 'generated/' in name for name in output.namelist())
    output.extractall(restored)
    for member in output.infolist():
        if not member.is_dir():
            (restored / member.filename).chmod((member.external_attr >> 16) & 0o777)
subprocess.run(['git', '-C', str(restored / 'gadget'), 'fsck', '--full', '--strict'],
               check=True, capture_output=True)
assert subprocess.check_output(['git', '-C', str(restored / 'gadget'), 'status', '--porcelain']) == b''
assert not (without_index / 'index').exists()
print('Tracked example inputs survive exports from Git metadata without an index')

def refuse_inconsistent_source(label):
    destination = root / (label + '.zip')
    result = subprocess.run([sys.executable, str(project / 'tools/export_source.py'),
        '--out', str(destination), '--git-dir', str(without_index)], capture_output=True)
    assert result.returncode and not destination.exists(), label

readme = project / 'README.md'
original = readme.read_bytes()
readme.write_bytes(b'Uncommitted source text\n')
refuse_inconsistent_source('modified-source')
snapshot = root / 'snapshot.zip'
subprocess.run([sys.executable, str(project / 'tools/export_source.py'), '--out', str(snapshot),
                '--without-history'], check=True, capture_output=True)
with zipfile.ZipFile(snapshot) as output:
    assert output.read('gadget/README.md') == readme.read_bytes()
    assert not any(name.startswith('gadget/.git/') for name in output.namelist())
readme.write_bytes(original)
untracked = project / 'tools/untracked.py'
untracked.write_bytes(b'value = 2\n')
refuse_inconsistent_source('untracked-source')
untracked.unlink()
readme.unlink()
refuse_inconsistent_source('missing-source')
readme.write_bytes(original)
executable = project / 'tools/method.py'
executable.chmod(0o644)
refuse_inconsistent_source('changed-mode')
executable.chmod(0o755)
assert not (without_index / 'index').exists()
print('History exports refuse modified, untracked, missing and mode changed source; snapshots stay explicit')
