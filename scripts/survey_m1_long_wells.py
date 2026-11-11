"""ADR 0001 long-well survey on production geometry, one equilibrium at a time.

Default window matches the researcher's diagnostic survey: lambda_n=.40:.01:.99.
All integrated values are scoped to this window, not the complete pitch integral
or f. Residual weight is a partial-K lower bound, never a filled ledger gap.
"""

import argparse
import gc
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
from scipy.integrate import trapezoid

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.conventions import default_source, uniform_source, source_at_s
from alpha_analysis_3.ledger import (
    Geometry,
    measure_pitch,
    volume_surface,
    slice_refinement_error,
)
from alpha_analysis_3.wells import ACTIVE
from alpha_analysis_3.conventions import LAMBDA_LEVELS
from alpha_analysis_3.diagnostics import raw_transition_counts, plot_geometry
from validate_m1 import Watchdog, test_ledger_matrix
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / "data/manifest.json").read_text())
    report = dict(
        code_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in list(sorted((ROOT / "src/alpha_analysis_3").glob("*.py")))
            + [Path(__file__).resolve()]
        },
        python=platform.python_version(),
        numpy=np.__version__,
        controls=vars(args),
        workers=1,
        error_scope="Numerical window integral of resolved weights and partial-K lower bounds. "
        "No upper bound on missing tails or unmeasured pitches; no f or connectivity.",
        runs=[],
        matrix_rows=[],
    )
    sources = dict(uniform=uniform_source, peaked=default_source)
    selected = manifest["files"]
    if args.start_file:
        selected = selected[
            next(i for i, e in enumerate(selected) if e["id"] == args.start_file) :
        ]
    for entry in selected:
        if args.file and entry["id"] != args.file:
            continue
        path = ROOT / "data" / entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        rows = []
        with Watchdog() as guard:
            start = time.monotonic()
            s = np.unique(
                np.r_[np.geomspace(0.001, 0.02, 5), np.linspace(0.02, 1, args.surfaces)]
            )
            g = Geometry(
                BoozerField(path),
                s,
                nt=args.table,
                nphi=args.table,
                density=args.nodes / (2 * np.pi),
                periods=64,
                progress=lambda stage: setattr(
                    guard, "stage", f'{entry["id"]} {stage}'
                ),
            )
            setup = time.monotonic() - start
            levels = np.linspace(args.lambda_min, args.lambda_max, args.levels)
            if args.matrix:
                levels = np.unique(np.round(np.r_[levels, LAMBDA_LEVELS], 12))
            # Layer-major evaluation keeps only this surface's temporary
            # derivative tables, reusing them across pitches. Kernels, pitch
            # order within each layer, source weights and controls are identical.
            nb, ns = len(levels), len(g.layers)
            weights = np.zeros((nb, ns, 5))
            volumes = np.full((nb, ns), np.nan)
            frontier = np.zeros((nb, ns))
            stats = [[] for _ in levels]
            counts = [[] for _ in levels]
            caps = np.zeros(nb, dtype=int)
            failures = caps.copy()
            seconds = np.zeros(nb)
            diagnostic = np.array(
                [
                    args.matrix
                    and np.any(np.isclose(x, LAMBDA_LEVELS, rtol=0, atol=1e-12))
                    for x in levels
                ]
            )
            for i, layer in enumerate(g.layers):
                t_layer = time.monotonic()
                layer.surface.restore_derivatives()
                try:
                    for j, level in enumerate(levels):
                        guard.stage = (
                            f'{entry["id"]} s={g.s[i]:.6f}, lambda={level:.3f}'
                        )
                        t = time.monotonic()
                        b = g.Bmin + level * (g.Bmax - g.Bmin)
                        measure = measure_pitch(
                            SimpleNamespace(layers=[layer]),
                            b,
                            max_periods=256,
                            transverse_tolerance=(
                                0.01 if entry["id"] in ("DMercFail", "d23p4") else 0.06
                            ),
                        )
                        weights[j, i] = measure.weights[0]
                        caps[j] += measure.capped_nodes
                        failures[j] += measure.failed_nodes
                        stats[j].append(measure.diagnostics[0])
                        if diagnostic[j]:
                            volumes[j, i] = volume_surface(
                                layer.surface, b, args.volume_phi
                            )
                            counts[j].append(
                                raw_transition_counts(layer, measure.wells[0], b)
                            )
                            missing = (
                                g.partners[i] < 0
                                if i < len(g.partners)
                                else np.zeros(len(layer.H), bool)
                            )
                            good = missing & (measure.wells[0].status == ACTIVE)
                            frontier[j, i] = float(
                                layer.measure[good] @ measure.wells[0].K[good]
                            )
                        seconds[j] += time.monotonic() - t
                finally:
                    layer.surface.release_derivatives()
                print(
                    "SURFACE",
                    entry["id"],
                    i,
                    float(g.s[i]),
                    time.monotonic() - t_layer,
                    flush=True,
                )
                np.savez_compressed(
                    output / f'{entry["id"]}-checkpoint.npz',
                    s=g.s,
                    levels=levels,
                    weights=weights,
                    volume=volumes,
                    completed_surfaces=i + 1,
                )
            for j, level in enumerate(levels):
                b = g.Bmin + level * (g.Bmax - g.Bmin)
                diagnostics = stats[j]
                source_weights = {
                    name: trapezoid(
                        source_at_s(h, g.s)[:, None] * weights[j], g.s, axis=0
                    ).tolist()
                    for name, h in sources.items()
                }
                tolerance = 0.01 if entry["id"] in ("DMercFail", "d23p4") else 0.06
                errors = {
                    name: slice_refinement_error(g.s, weights[j], diagnostics, h)
                    for name, h in sources.items()
                }
                converged = all(x[1] < tolerance for x in errors.values())
                row = dict(
                    lambda_n=float(level),
                    b=b,
                    seconds=float(seconds[j]),
                    capped_nodes=int(caps[j]),
                    failed_nodes=int(failures[j]),
                    converged=converged,
                    all_local_refinements_converged=all(
                        d["converged"] for d in diagnostics
                    ),
                    refinement_errors={
                        name: dict(absolute=x[0], relative=x[1])
                        for name, x in errors.items()
                    },
                    maximum_periods=max(d["maximum_periods"] for d in diagnostics),
                    maximum_factor=max(d["factor"] for d in diagnostics),
                    extension_seconds=sum(
                        d["base_extension"]["extension_seconds"]
                        + d["extension_seconds"]
                        for d in diagnostics
                    ),
                    refinement_seconds=sum(
                        d["refinement_seconds"] for d in diagnostics
                    ),
                    nonconverged_surfaces=[
                        dict(
                            s=float(g.s[i]),
                            factor=d["factor"],
                            change=d["relative_change"],
                        )
                        for i, d in enumerate(diagnostics)
                        if not d["converged"]
                    ],
                    source_weights=source_weights,
                )
                rows.append(row)
                if diagnostic[j]:
                    results = {}
                    for name, h in sources.items():
                        value = source_weights[name]
                        volume = float(trapezoid(source_at_s(h, g.s) * volumes[j], g.s))
                        resolved, lower = value[:2]
                        total = resolved if caps[j] + failures[j] == 0 else float("nan")
                        results[name] = dict(
                            ridge_weight=total,
                            resolved_weight=resolved,
                            capped_weight_lower_bound=lower,
                            resolved_plus_lower=resolved + lower,
                            gap_to_volume=volume - resolved - lower,
                            volume_weight=volume,
                            ratio=total / volume if volume > 0 else float("nan"),
                            capped_nodes=int(caps[j]),
                            failed_nodes=int(failures[j]),
                            refinement_converged=errors[name][1] < tolerance,
                            refinement_error=errors[name][0],
                            refinement_relative_error=errors[name][1],
                            all_local_refinements_converged=all(
                                d["converged"] for d in diagnostics
                            ),
                            frontier_weight=float(
                                trapezoid(source_at_s(h, g.s) * frontier[j], g.s)
                            ),
                        )
                    report["matrix_rows"].append(
                        dict(
                            file=entry["id"],
                            data_sha256=entry["sha256"],
                            lambda_n=float(level),
                            b=b,
                            Bmin=g.Bmin,
                            Bmax=g.Bmax,
                            seconds=float(seconds[j]),
                            sources=results,
                            raw_T_Gamma=counts[j],
                            count_mismatches=sum(a != c for a, c in counts[j]),
                            adaptive=diagnostics,
                        )
                    )
                print(entry["id"], json.dumps(row), flush=True)
            partial = dict(
                file=entry["id"],
                data_sha256=entry["sha256"],
                setup_seconds=setup,
                controls=vars(args),
                s=g.s.tolist(),
                extrema=g.extrema,
                rows=rows,
                source_hashes=report["source_hashes"],
                code_revision=report["code_revision"],
                nodes=sum(len(l.H) for l in g.layers),
                unpartnered_nodes=sum(int(np.sum(p < 0)) for p in g.partners),
            )
            window = [
                row
                for row in rows
                if args.lambda_min - 1e-12 <= row["lambda_n"] <= args.lambda_max + 1e-12
            ]
            b = np.array([row["b"] for row in window])
            summary = {}
            for name, h in sources.items():
                weights = np.array([row["source_weights"][name] for row in window])
                integrated = trapezoid(weights / b[:, None] ** 2, b, axis=0)
                volume = g.volume(h)
                summary[name] = dict(
                    V_h=volume,
                    resolved_window_integral=float(integrated[0]),
                    capped_window_lower_bound=float(integrated[1]),
                    capped_lower_fraction_Vh=float(integrated[1] / volume),
                    capped_lower_fraction_total_birth=float(
                        integrated[1] / (2 * volume)
                    ),
                    budget_fraction_Vh=0.005,
                    budget_status=(
                        "lower_bound_exceeds_budget"
                        if integrated[1] > 0.005 * volume
                        else "not_certified_by_a_lower_bound_or_partial_pitch_window"
                    ),
                    resolved_tail_fraction={
                        str(limit): float(integrated[j] / integrated[0])
                        for j, limit in enumerate((48, 64, 128), 2)
                    },
                )
            partial.update(
                total_seconds=time.monotonic() - start,
                peak_rss_gb=guard.peak,
                snapshots=guard.snapshots,
                summary=summary,
                maximum_periods=max(row["maximum_periods"] for row in rows),
                capped_slices=sum(row["capped_nodes"] > 0 for row in rows),
                failed_slices=sum(row["failed_nodes"] > 0 for row in rows),
                nonconverged_slices=sum(not row["converged"] for row in rows),
                local_ceiling_slices=sum(
                    bool(row["nonconverged_surfaces"]) for row in rows
                ),
                extension_seconds=sum(row["extension_seconds"] for row in rows),
                refinement_seconds=sum(row["refinement_seconds"] for row in rows),
            )
            report["runs"].append(partial)
            report["matrix_failures"] = test_ledger_matrix(report["matrix_rows"])
            plot_geometry(
                g,
                output / f'{entry["id"]}-geometry.png',
                g.Bmin + 0.95 * (g.Bmax - g.Bmin),
            )
            (output / f'{entry["id"]}.json').write_text(json.dumps(partial, indent=2))
            (output / "survey.json").write_text(json.dumps(report, indent=2))
            print(
                "RUN",
                entry["id"],
                json.dumps({k: v for k, v in partial.items() if k != "rows"}),
                flush=True,
            )
            del g, measure
            gc.collect()
    return int(
        bool(report.get("matrix_failures"))
        or any(
            r["failed_slices"]
            or r["nonconverged_slices"]
            or any(
                x["budget_status"] == "lower_bound_exceeds_budget"
                for x in r["summary"].values()
            )
            for r in report["runs"]
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/m1-long-survey")
    parser.add_argument(
        "--start-file", choices=["PCA", "TURBO", "DMercFail", "d23p4", "n3are"]
    )
    parser.add_argument(
        "--file", choices=["PCA", "TURBO", "DMercFail", "d23p4", "n3are"]
    )
    parser.add_argument("--surfaces", type=int, default=48)
    parser.add_argument("--nodes", type=int, default=256)
    parser.add_argument("--table", type=int, default=512)
    parser.add_argument("--levels", type=int, default=60)
    parser.add_argument("--lambda-min", type=float, default=0.4)
    parser.add_argument("--lambda-max", type=float, default=0.99)
    parser.add_argument("--matrix", action="store_true")
    parser.add_argument("--volume-phi", type=int, default=1024)
    raise SystemExit(run(parser.parse_args()))
