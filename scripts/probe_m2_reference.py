"""Separate angular-reference error from ridge weight on a fixed PCA grid.

The positive strip ledger determines per-surface ridge weights from the
inactive first strip. Those weights stay fixed while the independent angular
integral changes theta and phi resolution separately. No ledger gap is filled.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import minimize

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.conventions import default_source, uniform_source, source_at_s
from alpha_analysis_3.ledger import volume_surface
from alpha_analysis_3.tables import Surface
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    input_path = Path(args.report)
    source = json.loads(input_path.read_text())
    row = next(
        r
        for r in source["rows"]
        if r["file"] == args.file and r["lambda_n"] == args.level
    )
    strips = row["sources"]["uniform"]["strips"]
    radii = np.array([strips[0]["s_bottom"]] + [s["s_top"] for s in strips])
    assert (
        strips[0]["resolved_weight"] == 0.0
    ), "Need an inactive first strip to anchor reconstruction"
    ridge = [0.0]
    for strip, ds in zip(strips, np.diff(radii)):
        ridge.append(2 * strip["resolved_weight"] / ds - ridge[-1])
    ridge = np.asarray(ridge)
    assert np.min(ridge) > -1e-8
    entry = next(
        e
        for e in json.loads((ROOT / "data/manifest.json").read_text())["files"]
        if e["id"] == args.file
    )
    path = ROOT / "data" / entry["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    resolutions = [(1024, 1), (2048, 1), (2048, 2), (2048, 4), (2048, 8), (4096, 4)]
    result = dict(
        revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (ROOT / "src/alpha_analysis_3").glob("*.py")
        },
        input_report_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        controls=vars(args),
        file=entry["id"],
        data_sha256=entry["sha256"],
        b=row["b"],
        extrema=row["extrema"],
        rows=[],
        sources={},
        workers=1,
        scope="Fixed ridge weights and fixed radial grid; independent angular-reference refinement only",
    )
    start = time.monotonic()
    with Watchdog() as guard:
        field = BoozerField(path)
        values = []
        for s, R in list(zip(radii, ridge)) + [(s, None) for s in args.extra_radius]:
            sf = Surface(
                field, s, source["controls"]["table"], source["controls"]["table"]
            )
            guard.stage = f"{args.file} lambda={args.level} reference s={s}"
            maximum = None
            if abs(sf.Bmax - row["b"]) < 2 * sf.extremum_grid_error:
                i, j = np.unravel_index(np.argmax(sf.B), sf.B.shape)
                point = np.array([i * 2 * np.pi / sf.nt, j * 2 * np.pi / sf.nphi])
                opt = minimize(
                    lambda x: -float(sf.exact(x[0], x[1])),
                    point,
                    method="BFGS",
                    options={"gtol": 1e-10},
                )
                maximum = -float(opt.fun)
            angular = [
                volume_surface(sf, row["b"], nphi=nphi, theta_factor=factor)
                for nphi, factor in resolutions
            ]
            if R is not None:
                values.append(angular)
            result["rows"].append(
                dict(
                    s=float(s),
                    ridge_resolved_per_ds=float(R) if R is not None else None,
                    volume_per_ds=angular,
                    grid_Bmax=sf.Bmax,
                    grid_extremum_error=sf.extremum_grid_error,
                    polished_maximum_estimate=maximum,
                )
            )
            print(float(s), R, angular, flush=True)
        values = np.asarray(values)
        for name, h in [("uniform", uniform_source), ("peaked", default_source)]:
            weight = source_at_s(h, radii)
            known = float(trapezoid(weight * ridge, radii))
            np.testing.assert_allclose(
                known, row["sources"][name]["resolved_weight"], rtol=1e-10
            )
            references = trapezoid(weight[:, None] * values, radii, axis=0)
            result["sources"][name] = dict(
                ridge_resolved_weight=known,
                reference_weights=references.tolist(),
                ratios=(known / references).tolist(),
                relative_theta_change=float(
                    abs(references[4] - references[3]) / references[4]
                ),
                relative_phi_change=float(
                    abs(references[5] - references[3]) / references[5]
                ),
            )
        result.update(
            seconds=time.monotonic() - start,
            peak_rss_gb=guard.peak,
            snapshots=guard.snapshots,
            resolutions=[dict(nphi=n, theta_factor=f) for n, f in resolutions],
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result["false_zero_surfaces"] = [
        r["s"]
        for r in result["rows"]
        if r["polished_maximum_estimate"] is not None
        and r["grid_Bmax"] <= row["b"] < r["polished_maximum_estimate"]
    ]
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["sources"], flush=True)
    return bool(result["false_zero_surfaces"])


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--report", required=True)
    p.add_argument("--file", default="PCA")
    p.add_argument("--level", type=float, default=0.95)
    p.add_argument("--extra-radius", type=float, nargs="*", default=[])
    p.add_argument("--output", default="outputs/m2/reference-refinement.json")
    raise SystemExit(run(p.parse_args()))
