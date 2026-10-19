"""Exact byte stream publication, real output failure, concurrency and stop checks."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

publisher, contracts, root = (Path(value).resolve() for value in sys.argv[1:4])
root.mkdir(parents=True)
subprocess.run([str(contracts), str(root / "native")], check=True)
checks = []


def clean():
    assert not list(root.glob(".rfc-stage-*"))


def call(path, data, count=None, success=True, **kwargs):
    command = [str(publisher), "--out", str(path), "--bytes", str(len(data) if count is None else count)]
    result = subprocess.run(command, input=data, capture_output=True, **kwargs)
    assert (result.returncode == 0) == success, (command, result.stderr)
    if success:
        assert path.read_bytes() == data
    checks.append({"input_bytes": len(data), "declared_count": str(len(data) if count is None else count),
                   "exit_code": result.returncode, "success": success})
    return result


for index, size in enumerate([0, 1, 255, 1048575, 1048576, 1048577, 2097281]):
    raw = (bytes(range(256)) * (size // 256 + 1))[:size]
    path = root / f"complete space β {index}"
    call(path, raw)
    call(path, b"replacement", success=False)
    assert path.read_bytes() == raw
    clean()

for index, count in enumerate(["", "-1", "+1", "false", "1x", str(2**128)]):
    path = root / f"invalid-count-{index}"
    call(path, b"input", count, success=False)
    assert not path.exists()
    clean()

retry = root / "retry"
for raw, count in [(b"short", 9), (b"trailing", 4), (b"extra", 0)]:
    call(retry, raw, count, success=False)
    assert not retry.exists()
    clean()
call(retry, b"complete after refused inputs")
directory = root / "directory"
directory.mkdir()
call(directory, b"input", success=False)
assert directory.is_dir()
alias = root / "alias"
alias.symlink_to(retry)
call(alias, b"replacement", success=False)
assert alias.is_symlink() and retry.read_bytes() == b"complete after refused inputs"
dangling = root / "dangling"
dangling.symlink_to(root / "absent")
call(dangling, b"replacement", success=False)
assert dangling.is_symlink() and not (root / "absent").exists()
clean()

if os.name == "posix":
    import resource

    def limited():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (1024, 1024))

    path = root / "limited"
    raw = bytes(range(256)) * 8192
    result = call(path, raw, success=False, preexec_fn=limited)
    assert any(message in result.stderr for message in (b"file write failed", b"staged output flush failed", b"staged output close failed")), result.stderr
    assert not path.exists()
    clean()
    call(path, raw)

winner = root / "conditional"
payloads = [bytes([index]) * (1048576 + index) for index in range(1, 5)]
with ThreadPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(lambda raw: subprocess.run(
        [str(publisher), "--out", str(winner), "--bytes", str(len(raw))],
        input=raw, capture_output=True), payloads))
successful = [index for index, result in enumerate(results) if result.returncode == 0]
assert len(successful) == 1 and winner.read_bytes() == payloads[successful[0]]
clean()

# A stopped process cannot run destructors; its unpublished owned stage is retained.
stopped = root / "stopped"
with (root / "stopped.stderr").open("w") as stderr:
    child = subprocess.Popen([str(publisher), "--out", str(stopped), "--bytes", "2097152"],
                             stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=stderr)
    try:
        child.stdin.write(b"x" * 1048576 + b"partial")
        child.stdin.flush()
        deadline = time.monotonic() + 10
        # A stream may retain part of the first chunk in its output buffer.
        # Observe a real staged write without requiring the complete chunk on disk.
        def staged_bytes_present():
            for path in root.glob(".rfc-stage-*"):
                try:
                    if (path / "image").stat().st_size > 0:
                        return True
                except FileNotFoundError:
                    continue
            return False

        while not staged_bytes_present():
            assert child.poll() is None
            assert time.monotonic() < deadline, "owned writer did not enter staging"
            time.sleep(.01)
        child.kill()
        child.wait(timeout=10)
    finally:
        child.stdin.close()
        if child.poll() is None:
            child.kill()
            child.wait(timeout=10)
assert not stopped.exists()
orphans = [str(path.name) for path in root.glob(".rfc-stage-*")]
assert orphans
call(stopped, b"complete retry after process stop")
assert [str(path.name) for path in root.glob(".rfc-stage-*")] == orphans
report = {"status": "PASS", "checks": checks, "conditional_writers": 4, "complete_winners": 1,
          "real_filesystem_failure_and_retry": os.name == "posix",
          "stopped_writer_destination_absent": True, "retry_after_stop_complete": True,
          "owned_staging_residue_after_abrupt_stop": orphans,
          "scope": "Existing SDK staged exclusive local file publication; trusted parents/hard-link support. No fsync durability or automatic orphan recovery claim."}
(root / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")
print("PASS: exact bytes, real filesystem failure/retry, one conditional winner and stopped-writer namespace preservation")
