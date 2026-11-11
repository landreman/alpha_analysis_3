"""Thirty real-input cut-cell strip probes, selected from coverage reports.

This exercises M2 cell construction at every required pitch. Full-slice
population closure is in validate_m2_coverage; loss weights and f remain M3.
Only one equilibrium/strip is processed at a time under the RSS watchdog.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace

import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.cells import SurfaceCache, cut_cells, polish_folds, stitch_ridges
from alpha_analysis_3.conventions import LAMBDA_LEVELS
from alpha_analysis_3.diagnostics import plot_transition_cells
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.tables import Surface
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    selected = {}
    sources = {}
    for name in args.reports:
        report = json.loads(Path(name).read_text())
        sources[name] = hashlib.sha256(Path(name).read_bytes()).hexdigest()
        for row in report["rows"]:
            selected[row["file"], row["lambda_n"]] = (row, report["controls"], name)
    manifest = json.loads((ROOT / "data/manifest.json").read_text())["files"]
    expected = {(entry["id"], level) for entry in manifest for level in LAMBDA_LEVELS}
    assert expected <= set(
        selected
    ), f"missing required cases: {expected-set(selected)}"
    report = dict(
        revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (ROOT / "src/alpha_analysis_3").glob("*.py")
        },
        controls=vars(args),
        input_report_hashes=sources,
        workers=1,
        rows=[],
        runs=[],
        scope="one highest-resolved-population radial strip per required case; no f or full cell-loss weights",
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        temporary = output.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(output)

    save()
    for entry in manifest:
        path = ROOT / "data" / entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        with Watchdog() as guard:
            start = time.monotonic()
            field = BoozerField(path)
            for level in LAMBDA_LEVELS:
                t = time.monotonic()
                source, controls, source_name = selected[entry["id"], level]
                b = source["b"]
                strip = max(
                    source["sources"]["uniform"]["strips"],
                    key=lambda r: r["resolved_weight"],
                )
                radii = [strip["s_bottom"], strip["s_top"]]
                guard.stage = f"{entry['id']} lambda={level}, strip={radii}"
                layers = []
                density = controls["nodes"] / (2 * np.pi)
                for s in radii:
                    sf = Surface(field, s, controls["table"], controls["table"])
                    arcs = place_nodes(sf, ridge_arcs(sf), density)
                    layers.append(SimpleNamespace(surface=sf, arcs=arcs))
                folds = polish_folds(SimpleNamespace(layers=layers))
                a, c = layers
                mesh = stitch_ridges(
                    a.arcs,
                    c.arcs,
                    *radii,
                    a.surface.iota,
                    c.surface.iota,
                    field.nfp,
                    density,
                )
                sampler = SurfaceCache(field, controls["table"])
                sampler.cache[radii[0]] = a.surface
                sampler.cache[radii[1]] = c.surface
                result = cut_cells(
                    sampler,
                    mesh.face_points,
                    b,
                    max_depth=args.depth,
                    order=args.order,
                    max_cells=args.max_cells,
                    frontier_faces=mesh.unknown_faces,
                )
                counts = Counter(cell.kind for cell in result.cells)
                reasons = Counter(
                    cell.reason for cell in result.cells if cell.kind == "unknown"
                )
                areas = Counter()
                for cell in result.cells:
                    p = cell.points
                    alpha = (
                        p[:, 1] - np.asarray(field.iota(p[:, 0])) * p[:, 2] / field.nfp
                    )
                    area = (
                        abs(
                            np.sum(
                                p[:, 0] * np.roll(alpha, -1)
                                - alpha * np.roll(p[:, 0], -1)
                            )
                        )
                        / 2
                    )
                    areas[cell.kind] += float(area)
                residuals = [
                    abs(e.parent - e.child_1 - e.child_3) / e.parent
                    for segment in result.transitions
                    for e in segment
                ]
                row = dict(
                    file=entry["id"],
                    data_sha256=entry["sha256"],
                    lambda_n=level,
                    b=b,
                    extrema=source["extrema"],
                    source_report=source_name,
                    nodes=controls["nodes"],
                    radii=radii,
                    counts=dict(counts),
                    areas=dict(areas),
                    unknown_reasons=dict(reasons),
                    fold_records=folds,
                    transitions=len(result.transitions),
                    evaluations=result.evaluations,
                    maximum_relative_additivity_residual=max(residuals, default=None),
                    maximum_relative_barrier_height_residual=max(
                        (
                            abs(e.height_residual) / b
                            for segment in result.transitions
                            for e in segment
                        ),
                        default=None,
                    ),
                    seconds=time.monotonic() - t,
                )
                report["rows"].append(row)
                save()
                print(
                    entry["id"],
                    level,
                    dict(counts),
                    "segments",
                    len(result.transitions),
                    "seconds",
                    row["seconds"],
                    flush=True,
                )
                if args.plots and level == 0.8:
                    directory = Path(args.plots)
                    directory.mkdir(parents=True, exist_ok=True)
                    plot_transition_cells(
                        result,
                        field,
                        directory / f"{entry['id']}-cells.png",
                        f"{entry['id']}, lambda_n={level}",
                    )
                del sampler, result, layers, mesh
            report["runs"].append(
                dict(
                    file=entry["id"],
                    seconds=time.monotonic() - start,
                    peak_rss_gb=guard.peak,
                    snapshots=guard.snapshots,
                )
            )
            save()
    failed = [
        (r["file"], r["lambda_n"])
        for r in report["rows"]
        if r["counts"].get("regular", 0) + r["counts"].get("bank", 0) == 0
    ]
    report["failures"] = failed
    save()
    return bool(failed)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reports", nargs="+", required=True)
    p.add_argument("--output", default="outputs/m2/cell-probes.json")
    p.add_argument("--plots")
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--order", type=int, default=32)
    p.add_argument("--max-cells", type=int, default=20000)
    raise SystemExit(run(p.parse_args()))
