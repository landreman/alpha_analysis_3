"""M2 strip coverage, population closure and real-input cut-cell probes.

Run equilibria sequentially under a 2 GB watchdog. Reports preserve incomplete
weights and failures; coverage is not filled from the independent reference.
"""

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.cells import (
    SurfaceCache,
    coverage_ledger,
    cut_cells,
    geometric_frontier,
    insert_surfaces,
    stitch_ridges,
    polish_folds,
)
from alpha_analysis_3.conventions import default_source, uniform_source
from alpha_analysis_3.ledger import (
    Geometry,
    measure_pitch,
    pitch_ledger,
    volume_surface,
)
from alpha_analysis_3.diagnostics import plot_transition_cells
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
        workers=1,
        source="h(sqrt(s))=(1-s)^4; h=1 alongside",
        rows=[],
        runs=[],
        scope="M2 coverage and cut cells; no reachability or f",
    )

    def save():
        temp = target.with_suffix(".tmp")
        temp.write_text(json.dumps(report, indent=2) + "\n")
        temp.replace(target)

    entries = json.loads((ROOT / "data/manifest.json").read_text())["files"]
    for entry in entries:
        if args.file and entry["id"] != args.file:
            continue
        path = ROOT / "data" / entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        with Watchdog() as guard:
            start = time.monotonic()
            field = BoozerField(path)
            radii = np.unique(
                np.r_[np.geomspace(0.001, 0.02, 5), np.linspace(0.02, 1, args.surfaces)]
            )
            density = args.nodes / (2 * np.pi)

            def progress(stage):
                guard.stage = f"{entry['id']} {stage}"

            geometry = Geometry(
                field,
                radii,
                nt=args.table,
                nphi=args.table,
                density=density,
                progress=progress,
                chain_budget_bytes=args.chain_budget_mib * 1024**2,
            )
            fold_records = polish_folds(geometry)
            insertions = []
            active_fields = [
                geometry.Bmin + v * (geometry.Bmax - geometry.Bmin) for v in args.levels
            ]
            while len(geometry.s) < args.max_surfaces:
                indicators = geometric_frontier(geometry, density, active_fields)
                # Open-fold tails need angular refinement, not repeated
                # insertion into the same radial intervals.
                candidates = []
                for row in indicators:
                    partner = row["partner_frontier_measure"]
                    if partner > 0:
                        candidates.append(
                            (partner, (row["s_bottom"] + row["s_top"]) / 2)
                        )
                if not candidates or sum(
                    c[0] for c in candidates
                ) < args.geometric_tolerance * sum(
                    r["total_measure"] for r in indicators
                ):
                    break
                candidates.sort(reverse=True)
                selected = [
                    s
                    for _, s in candidates[
                        : min(8, args.max_surfaces - len(geometry.s))
                    ]
                ]
                insertions.extend(selected)
                insert_surfaces(
                    geometry, field, selected, density, args.max_surfaces, progress
                )
                fold_records.extend(polish_folds(geometry))
            setup = time.monotonic() - start
            for level in args.levels:
                progress(f"lambda={level}")
                t = time.monotonic()
                b = geometry.Bmin + level * (geometry.Bmax - geometry.Bmin)
                measure = measure_pitch(
                    geometry,
                    b,
                    transverse_tolerance=args.control,
                    acceptance_tolerance=0.01,
                    retain_edges=True,
                    clip_folds=True,
                )
                references = []
                for layer in geometry.layers:
                    references.append(
                        volume_surface(
                            layer.surface,
                            b,
                            nphi=args.volume_phi,
                            theta_factor=args.volume_theta_factor,
                        )
                    )
                sources = {}
                for name, h in [
                    ("uniform", uniform_source),
                    ("peaked", default_source),
                ]:
                    baseline = pitch_ledger(
                        geometry, b, h, measure=measure, volume_values=references
                    )
                    coverage = coverage_ledger(geometry, measure, h, density)
                    coverage["volume_weight"] = baseline["volume_weight"]
                    coverage["ratio"] = (
                        coverage["cell_weight"] / baseline["volume_weight"]
                        if baseline["volume_weight"] > 0
                        else np.nan
                    )
                    coverage["resolved_plus_lower_ratio"] = (
                        (
                            coverage["resolved_weight"]
                            + coverage["capped_weight_lower_bound"]
                        )
                        / baseline["volume_weight"]
                        if baseline["volume_weight"] > 0
                        else np.nan
                    )
                    coverage["refinement_error"] = baseline["refinement_relative_error"]
                    sources[name] = coverage
                cut_probe = None
                if args.cut_probe:
                    # Exercise all cell classes on the radial strip carrying
                    # the greatest known population. Full graph weights are M3.
                    weights = measure.weights[:, 0] + measure.weights[:, 1]
                    i = int(
                        np.argmax(np.diff(geometry.s) * (weights[:-1] + weights[1:]))
                    )
                    a, c = geometry.layers[i : i + 2]
                    mesh = stitch_ridges(
                        a.arcs,
                        c.arcs,
                        geometry.s[i],
                        geometry.s[i + 1],
                        a.surface.iota,
                        c.surface.iota,
                        a.surface.nfp,
                        density,
                    )
                    probe = cut_cells(
                        SurfaceCache(field, args.table),
                        mesh.face_points,
                        b,
                        max_depth=args.depth,
                        max_cells=20000,
                        frontier_faces=mesh.unknown_faces,
                    )
                    kinds = {
                        kind: sum(c.kind == kind for c in probe.cells)
                        for kind in ("regular", "bank", "unknown", "inactive")
                    }
                    cut_probe = dict(
                        s_bottom=float(geometry.s[i]),
                        s_top=float(geometry.s[i + 1]),
                        kinds=kinds,
                        evaluations=probe.evaluations,
                        transition_segments=len(probe.transitions),
                    )
                    if args.plots and level == 0.8:
                        directory = Path(args.plots)
                        directory.mkdir(parents=True, exist_ok=True)
                        plot_transition_cells(
                            probe,
                            field,
                            directory / f"{entry['id']}-cells.png",
                            f"{entry['id']}, lambda_n={level}",
                        )
                row = dict(
                    file=entry["id"],
                    data_sha256=entry["sha256"],
                    lambda_n=level,
                    b=b,
                    extrema=geometry.extrema,
                    sources=sources,
                    capped_nodes=measure.capped_nodes,
                    failed_nodes=measure.failed_nodes,
                    cut_probe=cut_probe,
                    seconds=time.monotonic() - t,
                )
                report["rows"].append(row)
                save()
                print(
                    entry["id"],
                    level,
                    "radii",
                    len(geometry.s),
                    "sources",
                    {
                        k: (v["ratio"], v["frontier_fraction"], v["refinement_error"])
                        for k, v in sources.items()
                    },
                    "cut",
                    cut_probe,
                    "seconds",
                    row["seconds"],
                    flush=True,
                )
            report["runs"].append(
                dict(
                    file=entry["id"],
                    seconds=time.monotonic() - start,
                    setup_seconds=setup,
                    insertions=insertions,
                    fold_records=fold_records,
                    radii=geometry.s.tolist(),
                    peak_rss_gb=guard.peak,
                    retained_chain_bytes=geometry.chain_budget.retained_bytes,
                    snapshots=guard.snapshots,
                )
            )
            save()
            del geometry, measure
            gc.collect()
    failures = []
    for row in report["rows"]:
        for name, v in row["sources"].items():
            ratio = v["ratio"]
            closed = (
                (np.isfinite(ratio) and abs(ratio - 1) <= 0.01)
                if not row["capped_nodes"]
                else v["resolved_plus_lower_ratio"] <= 1.01
            )
            if (
                not closed
                or row["failed_nodes"]
                or not np.isfinite(v["frontier_fraction"])
                or v["frontier_fraction"] > 0.005
                or v["refinement_error"] > 0.01
            ):
                failures.append((row["file"], row["lambda_n"], name))
    report["failures"] = failures
    save()
    print("Failures:", failures, flush=True)
    return bool(failures)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--file")
    p.add_argument("--output", default="outputs/m2/coverage.json")
    p.add_argument("--nodes", type=int, default=128)
    p.add_argument("--surfaces", type=int, default=24)
    p.add_argument("--max-surfaces", type=int, default=96)
    p.add_argument("--chain-budget-mib", type=int, default=768)
    p.add_argument("--table", type=int, default=512)
    p.add_argument("--volume-phi", type=int, default=1024)
    p.add_argument("--volume-theta-factor", type=int, default=4)
    p.add_argument("--control", type=float, default=0.0025)
    p.add_argument("--geometric-tolerance", type=float, default=0.0005)
    p.add_argument("--cut-probe", action="store_true")
    p.add_argument("--plots")
    p.add_argument("--depth", type=int, default=3)
    p.add_argument(
        "--levels", type=float, nargs="+", default=[0.05, 0.1, 0.5, 0.8, 0.9, 0.95]
    )
    raise SystemExit(run(p.parse_args()))
