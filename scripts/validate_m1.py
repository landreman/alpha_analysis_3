"""Sequential, memory-guarded M1 matrix; full numerical evidence outside CI.

Run from the repository: python scripts/validate_m1.py --output outputs/m1.
The JSON includes failed cases, provenance, timings and active-budget snapshots.
"""

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import threading
import time

import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.conventions import LAMBDA_LEVELS, default_source, uniform_source
from alpha_analysis_3.diagnostics import plot_geometry, raw_transition_counts
from alpha_analysis_3.ledger import Geometry, pitch_ledger
from alpha_analysis_3.wells import wells_at_pitch

ROOT = Path(__file__).resolve().parents[1]


class Watchdog:
    """Abort above 2 GB RSS; record active 60/180/300 s snapshots."""

    def __init__(self):
        self.start = time.monotonic()
        self.peak = 0.0
        self.snapshots = []
        self.stage = "starting"
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        pending = [60, 180, 300]
        while not self.stop.wait(0.2):
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            rss *= 1e-9 if sys.platform == "darwin" else 1024e-9
            self.peak = max(self.peak, rss)
            elapsed = time.monotonic() - self.start
            if pending and elapsed >= pending[0]:
                row = dict(seconds=pending.pop(0), stage=self.stage, rss_gb=rss)
                self.snapshots.append(row)
                print("SNAPSHOT", json.dumps(row), flush=True)
            if rss > 2.0:
                print("ABORT: resident memory exceeded 2 GB", flush=True)
                os._exit(137)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join()


def test_ledger_matrix(rows):
    """DESIGN M1 offline acceptance; retain every failure in the JSON."""
    failures = []
    for row in rows:
        tolerance = 0.01 if row["file"] in ("DMercFail", "d23p4") else 0.06
        for name, result in row["sources"].items():
            ratio = result.get("ratio", float("nan"))
            if not np.isfinite(ratio) or abs(ratio - 1) > tolerance:
                failures.append((row["file"], row["lambda_n"], name, ratio, tolerance))
    return failures


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / "data/manifest.json").read_text())
    rows = []
    report = dict(
        code_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT / "src/alpha_analysis_3").glob("*.py"))
        },
        python=platform.python_version(),
        numpy=np.__version__,
        workers=1,
        controls=vars(args),
        source="h(rho)=(1-rho^2)^4; rho=sqrt(s); h=1 alongside",
        extrema_scope="s in [0,1], 65 radial x 64² angular seeds, exact-series optimization; numerical estimates",
        rows=rows,
        runs=[],
    )
    for entry in manifest["files"]:
        if args.file and entry["id"] != args.file:
            continue
        path = ROOT / "data" / entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        with Watchdog() as guard:

            def progress(stage):
                guard.stage = stage

            start = time.monotonic()
            # Four extra graded inner surfaces, as required by the design.
            s = np.unique(
                np.r_[np.geomspace(0.001, 0.02, 5), np.linspace(0.02, 1, args.surfaces)]
            )
            geometry = Geometry(
                BoozerField(path),
                s,
                nt=args.table,
                nphi=args.table,
                density=args.nodes / (2 * np.pi),
                periods=args.periods,
                samples_per_period=args.samples,
                progress=progress,
            )
            setup = time.monotonic() - start
            print(entry["id"], "geometry", setup, flush=True)
            for level in LAMBDA_LEVELS:
                guard.stage = f'{entry["id"]} lambda={level}'
                t = time.monotonic()
                b = geometry.Bmin + level * (geometry.Bmax - geometry.Bmin)
                wells = [
                    wells_at_pitch(
                        l.surface,
                        l.points[:, 0],
                        l.points[:, 1],
                        l.H,
                        l.chain,
                        b,
                        order=args.order,
                    )
                    for l in geometry.layers
                ]
                results = {
                    name: pitch_ledger(
                        geometry,
                        b,
                        h,
                        volume_nphi=args.volume_phi,
                        order=args.order,
                        wells=wells,
                    )
                    for name, h in [
                        ("uniform", uniform_source),
                        ("peaked", default_source),
                    ]
                }
                counts = [
                    raw_transition_counts(l, w, b)
                    for l, w in zip(geometry.layers, wells)
                ]
                row = dict(
                    file=entry["id"],
                    data_sha256=entry["sha256"],
                    lambda_n=level,
                    b=b,
                    Bmin=geometry.Bmin,
                    Bmax=geometry.Bmax,
                    sources=results,
                    raw_T_Gamma=counts,
                    count_mismatches=sum(a != c for a, c in counts),
                    seconds=time.monotonic() - t,
                    failed_probes=[
                        dict(
                            s=l.surface.s,
                            node=int(j),
                            status=int(w.status[j]),
                            theta=float(l.points[j, 0]),
                            phi=float(l.points[j, 1]),
                            H=float(l.H[j]),
                            first_min=float(l.chain.min_H[j, 0]),
                            z_in=float(w.z_in[j]),
                            z_out=float(w.z_out[j]),
                        )
                        for l, w in zip(geometry.layers, wells)
                        for j in np.flatnonzero(w.status >= 2)
                    ],
                )
                rows.append(row)
                print(json.dumps(row), flush=True)
                (output / "matrix.json").write_text(json.dumps(report, indent=2))
            plot_geometry(geometry, output / f'{entry["id"]}-geometry.png', b)
            run = dict(
                file=entry["id"],
                setup_seconds=setup,
                total_seconds=time.monotonic() - start,
                peak_rss_gb=guard.peak,
                snapshots=guard.snapshots,
                volumes={
                    name: geometry.volume(h)
                    for name, h in [
                        ("uniform", uniform_source),
                        ("peaked", default_source),
                    ]
                },
                surfaces=geometry.s.tolist(),
                extrema=geometry.extrema,
                nodes=sum(len(l.H) for l in geometry.layers),
                unpartnered_nodes=sum(int(np.sum(p < 0)) for p in geometry.partners),
            )
            report["runs"].append(run)
            print("RUN", json.dumps(run), flush=True)
            del geometry, wells
            gc.collect()
        report["failures"] = test_ledger_matrix(rows)
        (output / "matrix.json").write_text(json.dumps(report, indent=2))
    return bool(report["failures"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/m1")
    parser.add_argument(
        "--file", choices=["PCA", "TURBO", "DMercFail", "d23p4", "n3are"]
    )
    parser.add_argument("--surfaces", type=int, default=24)
    parser.add_argument("--table", type=int, default=512)
    parser.add_argument("--nodes", type=int, default=128)
    parser.add_argument("--periods", type=int, default=48)
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--volume-phi", type=int, default=512)
    parser.add_argument("--order", type=int, default=64)
    raise SystemExit(run(parser.parse_args()))
