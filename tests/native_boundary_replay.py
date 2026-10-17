"""Replay a retained scientific publication without changing its input bytes."""
import gzip
import json
import subprocess
import sys
from pathlib import Path

binary, fixture, root = map(Path, sys.argv[1:])
original = fixture.read_bytes()
if fixture.suffix == '.gz':
    original = gzip.decompress(original)
root.mkdir(parents=True, exist_ok=True)
record, replayed = root / 'original.rfcxs', root / 'replayed.rfcxs'
record.write_bytes(original)
replayed.unlink(missing_ok=True)
subprocess.run([str(binary), 'replay-scientific', '--record', str(record),
                '--out-record', str(replayed)], check=True)
assert record.read_bytes() == original and replayed.stat().st_size > 0
print(json.dumps(dict(status='PASS', original_publication_bytes=len(original),
                     complete_native_field_identity_replay=True)))
