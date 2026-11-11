"""Sequential M2 transition census with resource and provenance records.

This is offline evidence, not a unit-test matrix. Failed surfaces and every
unresolved event remain in the output; a count miss makes the run nonzero.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.tables import Surface, global_extrema
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.transitions import surface_transitions
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    report = dict(
        revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (ROOT / "src/alpha_analysis_3").glob("*.py")
        },
        controls=vars(args),
        python=platform.python_version(),
        workers=1,
        thread_controls={
            k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
        },
        runs=[],
        rows=[],
        scope="Transition geometry only; no cell ledger or f claim",
    )

    def save():
        temporary = target.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(target)

    entries = json.loads((ROOT / "data/manifest.json").read_text())["files"]
    for entry in entries:
        if args.file and entry["id"] != args.file:
            continue
        with Watchdog() as guard:
            start = time.monotonic()
            path = ROOT / "data" / entry["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
            field = BoozerField(path)
            extrema = global_extrema(field)
            radii = np.unique(
                np.r_[np.geomspace(0.001, 0.02, 5), np.linspace(0.02, 1, args.surfaces)]
            )
            if args.radius is not None:
                radii = np.array([args.radius])
            for s in radii:
                sf = Surface(field, s, args.table, args.table)
                arcs = place_nodes(sf, ridge_arcs(sf), args.nodes / (2 * np.pi))
                for level in args.levels:
                    guard.stage = f"{entry['id']} s={s:.8g} lambda={level}"
                    t = time.monotonic()
                    b = extrema["min"]["B"] + level * (
                        extrema["max"]["B"] - extrema["min"]["B"]
                    )
                    result = surface_transitions(sf, arcs, b, args.factor, args.order)
                    row = dict(
                        file=entry["id"],
                        data_sha256=entry["sha256"],
                        s=float(s),
                        lambda_n=level,
                        b=b,
                        extrema=extrema,
                        transitions=result.transitions,
                        gamma_max=result.gamma_max,
                        resolved_events=len(result.events),
                        unknown_edges=len(result.unknown_edges),
                        factor=result.factor,
                        samples=result.samples,
                        failures=result.failures,
                        unknown_geometry=[
                            dict(
                                left=a.point.tolist(),
                                right=c.point.tolist(),
                                statuses=[a.status, c.status],
                                K=[a.K, c.K],
                            )
                            for a, c in result.unknown_edges
                        ],
                        events=[
                            dict(
                                point=e.point.tolist(),
                                partner=e.partner.tolist(),
                                parent=e.parent,
                                child_1=e.child_1,
                                child_3=e.child_3,
                                width=e.width,
                            )
                            for e in result.events
                        ],
                        seconds=time.monotonic() - t,
                    )
                    report["rows"].append(row)
                    print(
                        guard.stage,
                        "T/G",
                        result.transitions,
                        result.gamma_max,
                        "events",
                        len(result.events),
                        "unknown",
                        len(result.unknown_edges),
                        "factor",
                        result.factor,
                        "seconds",
                        round(row["seconds"], 3),
                        flush=True,
                    )
                    save()
                del sf
            report["runs"].append(
                dict(
                    file=entry["id"],
                    seconds=time.monotonic() - start,
                    rss_gb=guard.peak,
                    snapshots=guard.snapshots,
                )
            )
            save()
    failures = [
        r
        for r in report["rows"]
        if r["transitions"] != r["gamma_max"]
        or r["resolved_events"] != r["transitions"]
        or r["unknown_edges"]
    ]
    print("Unresolved surfaces:", len(failures), flush=True)
    return bool(failures)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="outputs/m2/transitions.json")
    parser.add_argument("--file")
    parser.add_argument("--radius", type=float)
    parser.add_argument("--surfaces", type=int, default=24)
    parser.add_argument("--nodes", type=int, default=128)
    parser.add_argument("--table", type=int, default=512)
    parser.add_argument("--factor", type=int, default=16)
    parser.add_argument("--order", type=int, default=32)
    parser.add_argument(
        "--levels", type=float, nargs="+", default=[0.05, 0.1, 0.5, 0.8, 0.9, 0.95]
    )
    raise SystemExit(run(parser.parse_args()))
