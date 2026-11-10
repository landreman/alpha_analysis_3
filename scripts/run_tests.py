"""Run a pytest tier with the repository's total wall-time budget."""

import os
import signal
import subprocess
import sys
import time


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] not in ("fast", "full"):
        raise SystemExit("usage: run_tests.py {fast|full} [pytest arguments]")
    tier = sys.argv[1]
    budget = 120 if tier == "fast" else 300
    command = [sys.executable, "-m", "pytest", "--durations=20"]
    if tier == "fast":
        command += ["-m", "not slow"]
    command += sys.argv[2:]
    start = time.monotonic()
    process = subprocess.Popen(command, start_new_session=True)
    try:
        code = process.wait(timeout=budget)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        print(f"{tier} tier exceeded its {budget}s wall-time budget", file=sys.stderr)
        code = 124
    print(f"{tier} tier wall time: {time.monotonic() - start:.2f}s / {budget}s")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
