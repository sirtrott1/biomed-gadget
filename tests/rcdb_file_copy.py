"""Complete native file copies retain captured owners and exclusive outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tracemalloc
from unittest.mock import patch


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gadget", type=Path)
    parser.add_argument("publisher", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("root", type=Path)
    parser.add_argument("--whole", action="store_true")
    parser.add_argument("--refused", type=Path)
    parser.add_argument("--reopen", type=Path)
    args = parser.parse_args()
    binary, publisher, source, root = (path.resolve() for path in
                                      (args.gadget, args.publisher, args.source, args.root))
    spec = importlib.util.spec_from_file_location(
        "gadget_rcdb", Path(__file__).resolve().parents[1]/"tools/rcdb/gadget_rcdb.py")
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    before, size = digest(source), source.stat().st_size
    if args.reopen:
        output = args.reopen.resolve()
        adapter.validate_publication_file(output, gadget=binary)
        assert output.stat().st_size == size and digest(output) == before
        print("Fresh-process complete native file replay PASS")
        return
    root.mkdir(parents=True)
    options = dict(gadget=binary, publisher=publisher, max_bytes=size)
    checks = []

    def refuses(name, action, error=ValueError, contains=None):
        try:
            action()
        except error as refused:
            if contains:
                assert contains in str(refused), str(refused)
            checks.append({"case": name, "refused": True, "message": str(refused)})
        else:
            raise AssertionError("Expected whole-file copy refusal: "+name)

    copied = root/"complete space β.rfcxs"
    tracemalloc.start()
    assert adapter.copy_publication_file(source, copied, **options) == copied
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert digest(copied) == before and copied.stat().st_size == size
    if args.whole:
        assert size > 64 * 1024 * 1024 and peak < 8 * 1024 * 1024, peak
    for value in (0, -1, True, "100"):
        refuses("invalid-budget-"+str(value), lambda: adapter.copy_publication_file(
            source, root/"invalid-budget.rfcxs", **{**options, "max_bytes": value}))
    refuses("measured-byte-budget", lambda: adapter.copy_publication_file(
        source, root/"over-budget.rfcxs", **{**options, "max_bytes": size-1}),
        contains="byte budget")
    assert not (root/"invalid-budget.rfcxs").exists() and not (root/"over-budget.rfcxs").exists()
    refuses("existing-output", lambda: adapter.copy_publication_file(source, copied, **options), FileExistsError)
    assert digest(copied) == before
    command = [sys.executable, "-I", "-B", str(Path(__file__).resolve()), str(binary), str(publisher),
               str(source), str(root), "--reopen", str(copied)]
    process = subprocess.run(command, capture_output=True, text=True, timeout=300)
    assert process.returncode == 0, (process.stdout, process.stderr)
    print(process.stdout.strip(), flush=True)
    if not args.whole:
        raw = source.read_bytes()
        directory = root/"directory"
        directory.mkdir()
        refuses("directory-source", lambda: adapter.copy_publication_file(directory, root/"directory-copy", **options))
        for name, target in (("existing-link", copied), ("dangling-link", root/"absent")):
            alias = root/name
            alias.symlink_to(target)
            refuses(name, lambda: adapter.copy_publication_file(source, alias, **options), FileExistsError)
            assert alias.is_symlink() and alias.readlink() == target
        hard = root/"source-hardlink"
        hard.hardlink_to(copied)
        refuses("source-hardlink", lambda: adapter.copy_publication_file(source, hard, **options), FileExistsError)
        hard.unlink()
        if hasattr(os, "mkfifo"):
            fifo = root/"input.fifo"
            os.mkfifo(fifo)
            refuses("fifo-source", lambda: adapter.copy_publication_file(fifo, root/"fifo-copy", **options))
        invalid = root/"invalid.rfcxs"
        invalid.write_bytes(raw[:-1])
        refuses("truncated-native", lambda: adapter.copy_publication_file(invalid, root/"invalid-copy", **options))
        assert not (root/"invalid-copy").exists()
        ordinary = subprocess.run
        changing = root/"changing.rfcxs"
        for mutation in ("grow", "truncate", "same-size"):
            changing.write_bytes(raw)
            changed = []

            def mutate_after_capture(command, **keywords):
                result = ordinary(command, **keywords)
                if (Path(command[0]).resolve() == publisher and Path(command[2]).name == "source.rfcxs"
                        and not changed):
                    changing.write_bytes(raw+b"extra" if mutation == "grow" else raw[:-1]
                                         if mutation == "truncate" else raw[:-1]+bytes([raw[-1] ^ 1]))
                    changed.append(mutation)
                return result

            with patch.object(subprocess, "run", mutate_after_capture):
                refuses("source-"+mutation, lambda: adapter.copy_publication_file(
                    changing, root/(mutation+"-copy"), **options), contains="changed while capturing")
            assert changed == [mutation] and not (root/(mutation+"-copy")).exists()
        control = root/"different.json"
        replacement = root/"different.rfcxs"
        control.write_text('{"purpose":"different admitted native source"}')
        ordinary([str(binary), "native-record", "--input", str(control), "--out", str(replacement)], check=True)
        replaced = []

        def substitute_capture(command, **keywords):
            if len(command) > 1 and command[1] == "native-record":
                captured = Path(command[command.index("--input")+1])
                captured.write_bytes(replacement.read_bytes())
                replaced.append(str(captured))
            return ordinary(command, **keywords)

        with patch.object(subprocess, "run", substitute_capture):
            refuses("admission-substitution", lambda: adapter.copy_publication_file(
                source, root/"substituted-copy", **options), contains="changed during native admission")
        assert len(replaced) == 1 and not (root/"substituted-copy").exists()
        changing.write_bytes(raw)
        replaced.clear()

        def replace_locator(command, **keywords):
            if len(command) > 1 and command[1] == "native-record" and not replaced:
                changing.write_bytes(replacement.read_bytes())
                replaced.append(str(changing))
            return ordinary(command, **keywords)

        pinned = root/"captured-owner.rfcxs"
        with patch.object(subprocess, "run", replace_locator):
            adapter.copy_publication_file(changing, pinned, **options)
        assert replaced == [str(changing)] and digest(pinned) == before and digest(changing) != before
        checks.append({"case": "fixed-owner-after-locator-change", "complete_original_bytes_equal": True})
        winner = root/"conditional-winner.rfcxs"
        with ThreadPoolExecutor(max_workers=2) as workers:
            def publish():
                try:
                    adapter.copy_publication_file(source, winner, **options)
                except FileExistsError:
                    return False
                return True
            outcomes = list(workers.map(lambda _: publish(), range(2)))
        assert sum(outcomes) == 1 and digest(winner) == before
        checks.append({"case": "conditional-writers", "writers": 2, "complete_winners": 1})
        failed = root/"failed-publication.rfcxs"

        def limit_output():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))

        def write_failure(command, **keywords):
            if Path(command[0]).resolve() == publisher and Path(command[2]) == failed:
                keywords["preexec_fn"] = limit_output
            return ordinary(command, **keywords)

        with patch.object(subprocess, "run", write_failure):
            refuses("real-output-byte-exhaustion", lambda: adapter.copy_publication_file(
                source, failed, **options), subprocess.CalledProcessError)
        assert not failed.exists() and not list(root.glob(".rfc-stage-*"))
        adapter.copy_publication_file(source, failed, **options)
        assert digest(failed) == before
        checks.append({"case": "retry-after-output-exhaustion", "complete_original_bytes_equal": True})
    if args.refused:
        foreign = args.refused.resolve()
        original = digest(foreign)
        refuses("independent-native-profile", lambda: adapter.copy_publication_file(
            foreign, root/"independent-copy.rfcxs", **{**options, "max_bytes": foreign.stat().st_size}))
        assert digest(foreign) == original and not (root/"independent-copy.rfcxs").exists()
    assert digest(source) == before and digest(copied) == before
    receipt = {"status": "PASS", "source_bytes": size, "source_sha256": before,
               "python_peak_traced_allocation_bytes": peak, "checks": checks,
               "complete_original_bytes_equal": True, "fresh_process_native_replay": True,
               "scope": "Existing exact exclusive file writer and canonical native admission; complete byte ownership, not a new codec, scientific calculation or fsync durability claim"}
    (root/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print("Complete native file copy, captured owner, bounded Python allocation and output preservation PASS", flush=True)


if __name__ == "__main__":
    main()
