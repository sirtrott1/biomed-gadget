#!/usr/bin/env python3
"""Measure one Linux workload independently of compilation; preserve its exit code."""
import argparse
import json
import resource
import subprocess
import sys
import time
from contextlib import ExitStack
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--stdout", type=Path)
    parser.add_argument("--stderr", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if sys.platform != "linux" or not command:
        parser.error("requires Linux and a command after --")
    with ExitStack() as files:
        stdout = files.enter_context(args.stdout.open("wb")) if args.stdout else None
        stderr = files.enter_context(args.stderr.open("wb")) if args.stderr else None
        started = time.perf_counter()
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        image_peak, samples = 0, 0
        while child.poll() is None:
            try:
                for line in Path(f"/proc/{child.pid}/status").read_text().splitlines():
                    if line.startswith("VmHWM:"):
                        image_peak = max(image_peak, int(line.split()[1]))
                        samples += 1
                        break
            except (OSError, ValueError):
                pass
            time.sleep(0.005)
        returncode = child.wait()
        elapsed = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    report = dict(command=command, exit_code=returncode, elapsed_seconds=elapsed,
                  peak_rss_kib=usage.ru_maxrss, user_seconds=usage.ru_utime,
                  system_seconds=usage.ru_stime, scope="single workload child; no compilation",
                  process_image_sampled_peak_rss_kib=image_peak if samples else None,
                  process_image_samples=samples,
                  memory_scope="VmHWM of the direct child after exec, sampled every 5 ms; observed lower bound. RUSAGE_CHILDREN peak_rss_kib is retained separately and may include a pre-exec inherited floor.")
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))
    return returncode if returncode >= 0 else 128 - returncode


if __name__ == "__main__":
    sys.exit(main())
