#!/usr/bin/env python3
"""Benchmark runner for Qelect SSLE protocol.

Runs docker containers with different thread counts, pins each to dedicated
CPU cores, saves results to timestamped files, then analyzes and outputs CSV.

## Design constraints & assumptions

### numcores ≤ group_size (enforced in MPE_docker_version.cpp)
If NUM_THREADS > group_size, that group_size is SKIPPED (continue). Division
by zero prevention. For 16 threads, group_size 2/4/8 are skipped.

### Core pinning — contiguous allocation only
Each docker container gets ``--cpuset-cpus=<start>-<end>``, a single
contiguous range. No two containers share a core. Jobs wait (block) when
no large-enough contiguous free block exists in the pool. Default ``-c 64``
covers one logical CPU per physical core on a 64-core HT Xeon, avoiding
hyper-thread interference.

### Session isolation
Each ``--skip-run`` false run creates a ``YYYYMMDD_HHMMSS/`` subdirectory
under ``--data-dir``.  Old results can never leak into new analysis.

### File naming
``result_<N>threads_run<M>.txt`` — simplified because the session directory
already carries the timestamp.

### NUM_THREADS via environment variable (no rebuild)
``MPE_docker_version.cpp`` reads ``getenv("NUM_THREADS")`` at runtime.
Defaults to 1 if unset. Changing thread count never requires a Docker
rebuild, only changing the C++ source does.

### ETA estimation
First run shows ``estimating...`` (no data). After at least one run
completes, ETA = avg_run_duration × remaining / active_parallelism.
Refreshed every 2 seconds via background ticker thread.
Only actually-running containers (not waiting-for-cores) count as active.

### Ctrl+C behaviour
1. ``_aborted`` event set → ``CorePool.allocate()`` stops granting cores
2. All tracked docker containers are ``docker stop``'d
3. ``os._exit(130)`` force-terminates the process

### CSV output units & columns
- ``benchmark_detail.csv``: scheme, party_count, total_time_ms (one row per
  measurement, no run column)
- ``benchmark_summary.csv``: scheme, party_count, runs, mean/std/min/max
- ``benchmark_mean.csv``: scheme, party_count, total_time_mean_ms (compact)
- ``threads`` is embedded in scheme name: ``Qelect`` (1 thread) or
  ``Qelect (N threads)`` (N > 1)
- Time values are **milliseconds** (raw us ÷ 1000.0)
"""

import argparse
import csv
import os
import re
import signal
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Global state for signal handling
# ---------------------------------------------------------------------------

_running_containers: set[str] = set()
_running_lock = threading.Lock()
_aborted = threading.Event()


# ---------------------------------------------------------------------------
# Core pool — contiguous allocation with blocking
# ---------------------------------------------------------------------------


class CorePool:
    """Manages a set of CPU cores, allocating contiguous intervals only."""

    def __init__(self, total: int):
        self.total = total
        self._free: list[tuple[int, int]] = [(0, total - 1)]
        self._lock = threading.Lock()
        self._cond = threading.Condition(self._lock)

    def allocate(self, n: int) -> tuple[int, int]:
        """Block until *n* contiguous cores are available; return (start, end).

        Raises SystemExit if Ctrl+C aborted.
        """
        with self._cond:
            while True:
                if _aborted.is_set():
                    sys.exit(1)
                for i, (lo, hi) in enumerate(self._free):
                    if hi - lo + 1 >= n:
                        a, b = lo, lo + n - 1
                        if b == hi:
                            del self._free[i]
                        else:
                            self._free[i] = (b + 1, hi)
                        return (a, b)
                self._cond.wait(timeout=1)

    def free(self, start: int, end: int) -> None:
        """Return an interval to the pool."""
        with self._cond:
            self._free.append((start, end))
            self._free.sort()
            merged: list[list[int]] = []
            for lo, hi in self._free:
                if merged and lo <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], hi)
                else:
                    merged.append([lo, hi])
            self._free = [(lo, hi) for lo, hi in merged]
            self._cond.notify_all()


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

_BLOCK_RE = re.compile(
    r"Total number of parties:\s+(\d+)\s+"
    r"Preprocessed time\s+:\s+\d+\s+us\.\s+"
    r"Total time\s+:\s+(\d+)\s+us\."
)

_FILENAME_RE = re.compile(r"result_(\d+)threads_run(\d+)\.txt$")


def parse_output(filepath: str) -> list[tuple[int, int]]:
    """Parse a result file, return [(party_count, total_time_us), ...]."""
    with open(filepath) as f:
        return [(int(m[0]), int(m[1])) for m in _BLOCK_RE.findall(f.read())]


def parse_filename(filename: str) -> tuple[int, int] | None:
    """Extract (threads, run_num) from 'result_Nthreads_runM.txt'."""
    m = _FILENAME_RE.search(filename)
    return (int(m.group(1)), int(m.group(2))) if m else None


# ---------------------------------------------------------------------------
# Progress display  (multi-line ANSI, auto-refresh ticker)
# ---------------------------------------------------------------------------


class Progress:
    """Thread-safe multi-line progress display with ETA.

    Tracks *waiting* (queued for cores) and *running* (docker started)
    separately.  ETA uses only running count for active-parallelism.
    """

    def __init__(self, total: int):
        self.total = total
        self.completed = 0
        self.start_time = time.time()
        self._lock = threading.Lock()
        self._run_times: list[float] = []
        self._waiting: dict[int, str] = {}  # job_id -> desc
        self._running: dict[int, str] = {}  # job_id -> desc
        self._done = False
        self._prev_lines = 0
        self._ticker = threading.Thread(target=self._tick, daemon=True)
        self._ticker.start()

    # -- public API ----------------------------------------------------------

    def set_waiting(self, job_id: int, desc: str):
        with self._lock:
            self._waiting[job_id] = desc
            self._refresh()

    def set_running(self, job_id: int, desc: str):
        with self._lock:
            self._waiting.pop(job_id, None)
            self._running[job_id] = desc
            self._refresh()

    def clear(self, job_id: int):
        with self._lock:
            self._waiting.pop(job_id, None)
            self._running.pop(job_id, None)

    def add_completed(self, elapsed: float):
        with self._lock:
            self.completed += 1
            self._run_times.append(elapsed)

    def finish(self):
        with self._lock:
            self._done = True
        if self._prev_lines > 0:
            print(f"\033[{self._prev_lines}A\033[J", file=sys.stderr, flush=True)
        elapsed = time.time() - self.start_time
        print(
            f"All {self.total} jobs finished in {_fmt_duration(elapsed)}.",
            file=sys.stderr,
        )

    # -- internal ------------------------------------------------------------

    def _tick(self):
        while True:
            time.sleep(2)
            with self._lock:
                if self._done:
                    return
                if self._running or self._waiting:
                    self._refresh()

    def _refresh(self):
        """Redraw the multi-line progress display."""
        elapsed = time.time() - self.start_time
        done, total = self.completed, self.total
        n_waiting, n_running = len(self._waiting), len(self._running)

        if self._run_times:
            avg = sum(self._run_times) / len(self._run_times)
            remaining = (total - done) * avg / n_running if n_running else 0
            eta = _fmt_duration(remaining)
        else:
            eta = "estimating..."

        lines = [
            f"[{done}/{total}]  running: {n_running}  waiting: {n_waiting}  "
            f"elapsed={_fmt_duration(elapsed)}  ETA={eta}"
        ]
        for jid in sorted(self._running):
            lines.append(f"  {self._running[jid]}")
        for jid in sorted(self._waiting):
            lines.append(f"  {self._waiting[jid]}")

        out = "\033[K" + "\n\033[K".join(lines)
        if self._prev_lines > 0:
            out = f"\033[{self._prev_lines}A" + out
        self._prev_lines = len(lines)
        print(out, file=sys.stderr, flush=True)


def _fmt_duration(seconds: float) -> str:
    h, m = divmod(int(seconds), 3600)
    m, s = divmod(m, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


# ---------------------------------------------------------------------------
# Signal handling
# ---------------------------------------------------------------------------


def _on_interrupt(_signum, _frame):
    """Ctrl+C: abort all waiters, stop containers, force exit."""
    _aborted.set()
    with _running_lock:
        names = list(_running_containers)
    if names:
        print(f"\nStopping {len(names)} docker container(s)...", file=sys.stderr)
        for name in names:
            subprocess.run(
                ["docker", "stop", name],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    os._exit(130)


signal.signal(signal.SIGINT, _on_interrupt)


# ---------------------------------------------------------------------------
# Single benchmark job
# ---------------------------------------------------------------------------


def run_one(
    threads: int,
    run_num: int,
    runs_total: int,
    data_dir: str,
    job_id: int,
    pool: CorePool,
    progress: Progress,
) -> dict:
    """Execute one docker run and return metadata dict."""
    progress.set_waiting(
        job_id, f"t={threads}  run={run_num}/{runs_total}  cores=[waiting]"
    )

    core_start, core_end = pool.allocate(threads)

    desc = f"t={threads}  run={run_num}/{runs_total}  cores=[{core_start}-{core_end}]"
    progress.set_running(job_id, desc)

    filename = f"result_{threads}threads_run{run_num}.txt"
    filepath = os.path.join(data_dir, filename)

    session_ts = os.path.basename(data_dir)
    container = f"qelect-{session_ts}-t{threads}-r{run_num}"
    with _running_lock:
        _running_containers.add(container)

    t0 = time.time()
    error = None

    try:
        with open(filepath, "w") as out:
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--name",
                    container,
                    "--cpuset-cpus",
                    f"{core_start}-{core_end}",
                    "-e",
                    f"NUM_THREADS={threads}",
                    "qelect_project",
                ],
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )
    except FileNotFoundError:
        error = "docker: command not found"
    except subprocess.CalledProcessError as e:
        error = f"exit={e.returncode}"
    except Exception as e:
        error = repr(e)
    finally:
        with _running_lock:
            _running_containers.discard(container)
        pool.free(core_start, core_end)

    elapsed = time.time() - t0
    progress.clear(job_id)
    progress.add_completed(elapsed)
    if error:
        print(f"  FAILED: {desc} - {error}", file=sys.stderr)

    return {
        "threads": threads,
        "run": run_num,
        "filepath": filepath,
        "cores": f"{core_start}-{core_end}",
        "elapsed": elapsed,
        "error": error,
    }


# ---------------------------------------------------------------------------
# Run phase — orchestrate all jobs
# ---------------------------------------------------------------------------


def run_benchmarks(args) -> tuple[list[dict], str]:
    """Launch all benchmark jobs; return (results, session_dir)."""
    thread_counts = [int(t.strip()) for t in args.threads.split(",")]

    # Build job list (larger threads first for better core packing)
    jobs: list[tuple[int, int]] = []
    for t in thread_counts:
        for r in range(1, args.runs + 1):
            jobs.append((t, r))
    jobs.sort(key=lambda x: (-x[0], x[1]))

    total = len(jobs)
    session_dir = os.path.join(args.data_dir, datetime.now().strftime("%Y%m%d_%H%M%S"))
    print(
        f"Benchmark: {len(thread_counts)} thread configs x {args.runs} runs "
        f"= {total} jobs on {args.total_cores} cores",
        file=sys.stderr,
    )
    print(f"Session dir: {session_dir}", file=sys.stderr, flush=True)

    os.makedirs(session_dir, exist_ok=True)

    pool = CorePool(args.total_cores)
    progress = Progress(total)
    results: list[dict] = []
    results_lock = threading.Lock()

    def _worker(idx: int, t: int, r: int):
        res = run_one(t, r, args.runs, session_dir, idx, pool, progress)
        with results_lock:
            results.append(res)

    wall_start = time.time()
    max_workers = min(args.total_cores, total)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for idx, (t, r) in enumerate(jobs):
            futures.append(executor.submit(_worker, idx, t, r))
        for f in futures:
            f.result()

    wall_elapsed = time.time() - wall_start
    progress.finish()
    print(f"Total wall-clock time: {_fmt_duration(wall_elapsed)}.", file=sys.stderr)

    results.sort(key=lambda x: (x["threads"], x["run"]))
    return results, session_dir


# ---------------------------------------------------------------------------
# Analyze phase — parse files and write CSVs
# ---------------------------------------------------------------------------


def _make_scheme_name(threads: int) -> str:
    return "Qelect" if threads == 1 else f"Qelect ({threads} threads)"


def _write_detail_csv(prefix: str, rows: list[dict]) -> str:
    path = f"{prefix}_detail.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, ["scheme", "party_count", "total_time_ms"])
        w.writeheader()
        w.writerows(rows)
    return path


def _write_summary_csv(
    prefix: str, data: dict[tuple[str, int], dict[int, list[float]]]
) -> str:
    path = f"{prefix}_summary.csv"
    fields = [
        "scheme",
        "party_count",
        "runs",
        "total_time_mean_ms",
        "total_time_std_ms",
        "total_time_min_ms",
        "total_time_max_ms",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        for (scheme, _), pc_map in sorted(data.items()):
            for pc in sorted(pc_map):
                vals = pc_map[pc]
                w.writerow(
                    {
                        "scheme": scheme,
                        "party_count": pc,
                        "runs": len(vals),
                        "total_time_mean_ms": f"{statistics.mean(vals):.3f}",
                        "total_time_std_ms": (
                            f"{statistics.stdev(vals):.3f}"
                            if len(vals) > 1
                            else "0.000"
                        ),
                        "total_time_min_ms": f"{min(vals):.3f}",
                        "total_time_max_ms": f"{max(vals):.3f}",
                    }
                )
    return path


def _write_mean_csv(
    prefix: str, data: dict[tuple[str, int], dict[int, list[float]]]
) -> str:
    path = f"{prefix}_mean.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, ["scheme", "party_count", "total_time_mean_ms"])
        w.writeheader()
        for (scheme, _), pc_map in sorted(data.items()):
            for pc in sorted(pc_map):
                vals = pc_map[pc]
                w.writerow(
                    {
                        "scheme": scheme,
                        "party_count": pc,
                        "total_time_mean_ms": f"{statistics.mean(vals):.3f}",
                    }
                )
    return path


def analyze_results(args) -> None:
    """Parse saved result files, write all CSV variants."""
    data_dir = Path(args.data_dir)
    if not data_dir.is_dir():
        sys.exit(f"Error: data dir '{args.data_dir}' not found.")

    files = sorted(data_dir.glob("result_*threads_run*.txt"))
    if not files:
        sys.exit(f"Error: no result_*threads_run*.txt files in {args.data_dir}")

    print(f"Found {len(files)} result file(s) in {args.data_dir}", file=sys.stderr)

    # Parse all files → [(threads, run, [(party_count, total_time_us), ...])]
    parsed: list[tuple[int, int, list[tuple[int, int]]]] = []
    for fp in files:
        info = parse_filename(fp.name)
        if info is None:
            print(f"  Warning: cannot parse '{fp.name}', skipping", file=sys.stderr)
            continue
        threads, run_num = info
        try:
            timing = parse_output(str(fp))
        except Exception as e:
            print(f"  Warning: failed to parse '{fp.name}': {e}", file=sys.stderr)
            continue
        if not timing:
            print(f"  Warning: no data in '{fp.name}', skipping", file=sys.stderr)
            continue
        parsed.append((threads, run_num, timing))
        print(
            f"  {fp.name}: threads={threads} run={run_num}, "
            f"{len(timing)} party counts",
            file=sys.stderr,
        )

    if not parsed:
        sys.exit("Error: no valid data to analyze.")

    parsed.sort(key=lambda x: (x[0], x[1]))

    # Build detail rows and summary aggregation
    detail_rows: list[dict] = []
    summary: dict[tuple[str, int], dict[int, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for threads, _run, timing in parsed:
        scheme = _make_scheme_name(threads)
        for party_count, total_time_us in timing:
            total_time_ms = total_time_us / 1000.0
            detail_rows.append(
                {
                    "scheme": scheme,
                    "party_count": party_count,
                    "total_time_ms": f"{total_time_ms:.3f}",
                }
            )
            summary[(scheme, threads)][party_count].append(total_time_ms)

    summary_plain = {k: dict(v) for k, v in summary.items()}

    # Write CSV files into data_dir
    csv_prefix = os.path.join(args.data_dir, args.output)

    p = _write_detail_csv(csv_prefix, detail_rows)
    print(f"Wrote {p} ({len(detail_rows)} rows)", file=sys.stderr)

    p = _write_summary_csv(csv_prefix, summary_plain)
    print(f"Wrote {p}", file=sys.stderr)

    p = _write_mean_csv(csv_prefix, summary_plain)
    print(f"Wrote {p}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(
        description="Qelect benchmark runner & analyzer",
        epilog="examples:\n"
        "  python3 benchmark.py -t 1 -r 1 -p 1     # quick test\n"
        "  python3 benchmark.py                     # full run (5x5)\n"
        "  python3 benchmark.py -t 1,2,8 -r 10      # custom config\n"
        "  python3 benchmark.py --skip-run          # analyze only",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    g = p.add_argument_group("run options")
    g.add_argument(
        "-t",
        "--threads",
        default="1,2,4,8,16",
        help="comma-separated thread counts (default: 1,2,4,8,16)",
    )
    g.add_argument(
        "-r",
        "--runs",
        type=int,
        default=5,
        help="number of runs per thread config (default: 5)",
    )
    g.add_argument(
        "-c",
        "--total-cores",
        type=int,
        default=64,
        help="total CPU cores available (default: 64)",
    )
    g.add_argument(
        "-p",
        "--parallel",
        default="auto",
        help="max concurrent containers, 'auto' or integer " "(default: auto)",
    )
    g.add_argument(
        "-d",
        "--data-dir",
        default="./data",
        help="directory for result files (default: ./data)",
    )
    g.add_argument(
        "--skip-run",
        action="store_true",
        help="skip docker runs, only analyze existing files",
    )

    g = p.add_argument_group("output options")
    g.add_argument(
        "-o",
        "--output",
        default="benchmark",
        help="CSV output prefix, produces <prefix>_detail.csv, "
        "<prefix>_summary.csv, <prefix>_mean.csv "
        "(default: benchmark)",
    )

    args = p.parse_args()

    # -- parameter validation ------------------------------------------------

    if args.runs < 1:
        sys.exit("Error: --runs must be >= 1")
    if args.total_cores < 1:
        sys.exit("Error: --total-cores must be >= 1")

    try:
        thread_counts = [int(t.strip()) for t in args.threads.split(",")]
    except ValueError:
        sys.exit("Error: --threads must be comma-separated positive integers")

    if not thread_counts:
        sys.exit("Error: --threads must not be empty")
    if any(t < 1 for t in thread_counts):
        sys.exit("Error: all thread counts must be >= 1")
    if any(t > args.total_cores for t in thread_counts):
        sys.exit(
            "Error: each thread count must be <= --total-cores " f"({args.total_cores})"
        )

    if args.parallel != "auto":
        try:
            args.parallel = int(args.parallel)
            if args.parallel < 1:
                raise ValueError
        except ValueError:
            sys.exit("Error: --parallel must be 'auto' or a positive integer")

    # -- run / analyze -------------------------------------------------------

    if not args.skip_run:
        results, session_dir = run_benchmarks(args)
        args.data_dir = session_dir
        failed = [r for r in results if r["error"]]
        if failed:
            for r in failed:
                print(
                    f"FAILED: threads={r['threads']} run={r['run']} "
                    f"file={r['filepath']}: {r['error']}",
                    file=sys.stderr,
                )
            sys.exit(f"Error: {len(failed)} benchmark job(s) failed")

    analyze_results(args)


if __name__ == "__main__":
    main()
