"""Independent counterexample to zero-weight fold endpoints and T-only gluing.

This is diagnostic evidence, not a proposed replacement for the state-space
contract. The oracle uses an explicit analytic field and independent roots and
quadrature. Production chains are compared without repairing their statuses.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.cells import polish_folds
from alpha_analysis_3.chains import trace_chains
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.tables import Surface
from alpha_analysis_3.wells import wells_at_pitch, ACTIVE
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


class FoldField:
    nfp = 1
    xm = np.array([0, 0, 1, 1, 0])
    xn = np.array([0, -1, -1, 1, -2])

    def bmnc(self, s):
        return np.array([6, 0, 0, 0, 0])

    def bmns(self, s):
        return np.array([0, -2, 0.25, -0.25, 1])

    def G(self, s):
        return 3 * (1 + 0.2 * s)

    def I(self, s):
        return 0.0

    def iota(self, s):
        return 0.0


def oracle(theta, s=0.5):
    """Finite physical well, irrespective of which maximum represents it."""
    amplitude = -2 + 0.5 * np.cos(theta)

    def B(phi):
        return 6 + amplitude * np.sin(phi) + np.sin(2 * phi)

    left = brentq(lambda p: B(p) - 5, 0.001, 2 * np.pi / 3, xtol=1e-14)
    right = brentq(lambda p: B(p) - 5, 2 * np.pi / 3, np.pi, xtol=1e-14)
    integrals, errors = [], []
    for power in (0.5, -0.5):

        def integrand(t):
            half = (right - left) / 2
            phi = (right + left) / 2 + half * np.sin(t)
            return (
                3
                * (1 + 0.2 * s)
                / B(phi)
                * (1 - B(phi) / 5) ** power
                * half
                * np.cos(t)
            )

        value, error = quad(integrand, -np.pi / 2, np.pi / 2, epsabs=1e-9, epsrel=1e-9)
        integrals.append(value)
        errors.append(error)
    return dict(z_in=left, z_out=right, J=integrals[0], K=integrals[1], errors=errors)


def maxima(theta):
    amplitude = -2 + 0.5 * np.cos(theta)
    discriminant = np.sqrt(amplitude**2 + 32)
    points = [-np.arccos((-amplitude - discriminant) / 8)]
    open_cosine = (-amplitude + discriminant) / 8
    if open_cosine < 1:
        points.append(np.arccos(open_cosine))
    return np.array(points)


def run(args):
    result = dict(
        revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_hashes={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__)]
            + list((ROOT / "src/alpha_analysis_3").glob("*.py"))
        },
        formula="B=6+(-2+0.5*cos(theta))*sin(phi)+sin(2*phi); C=3*(1+0.2*s); iota=0; b=5",
        fold=dict(theta=float(np.pi / 2), phi=0, H=6, D3B=-6, b=5),
        limit=oracle(np.pi / 2),
        analytic_transition_count=0,
        analytic_gamma_max_count=0,
        count_scope="All maxima remain above b in this neighbourhood; no H=b crossing or interior barrier event",
        rows=[],
        constant_J=[],
        real_fold_limits=[],
        scope="Diagnostic counterexample only; no corrected population, coverage, connectivity or f claim",
    )
    with Watchdog() as guard:
        guard.stage = "analytic fold identity and chain-resolution counterexample"
        sf = Surface(FoldField(), 0.5, 256, 256)
        for offset in (-0.01, -0.001, -0.0001, 0.0001, 0.001, 0.01):
            theta = np.pi / 2 + offset
            phi = maxima(theta)
            th = np.full(len(phi), theta)
            row = dict(
                offset=offset, oracle=oracle(theta), maxima_phi=phi.tolist(), scans=[]
            )
            for samples in (128, 512, 4096):
                chain = trace_chains(sf, th, phi, periods=2, samples_per_period=samples)
                well = wells_at_pitch(sf, th, phi, sf.b(th, phi), chain, 5)
                row["scans"].append(
                    dict(
                        samples_per_period=samples,
                        status=well.status.tolist(),
                        first_minimum=chain.min_H[:, 0].tolist(),
                        J=[float(v) if np.isfinite(v) else None for v in well.J],
                        K=[float(v) if np.isfinite(v) else None for v in well.K],
                    )
                )
            assert row["scans"][-1]["status"] == ([0, 1] if offset < 0 else [1])
            result["rows"].append(row)
        target = result["limit"]["J"]
        for s in (0.4, 0.5, 0.6):
            theta = brentq(
                lambda t: oracle(t, s)["J"] - target,
                np.pi / 2 - 0.5,
                np.pi / 2 + 0.5,
            )
            representative = (
                "fold boundary"
                if abs(theta - np.pi / 2) < 1e-10
                else "open" if theta < np.pi / 2 else "closed"
            )
            result["constant_J"].append(
                dict(
                    s=s,
                    theta=theta,
                    representative=representative,
                    **oracle(theta, s),
                )
            )
        if args.real_report:
            source = json.loads(Path(args.real_report).read_text())
            entry = next(
                e
                for e in json.loads((ROOT / "data/manifest.json").read_text())["files"]
                if e["id"] == args.real_file
            )
            path = ROOT / "data" / entry["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
            guard.stage = f"{args.real_file} one-surface finite fold limits"
            sf = Surface(BoozerField(path), args.real_radius, 512, 512)
            arcs = place_nodes(sf, ridge_arcs(sf), 256 / (2 * np.pi))
            records = polish_folds(
                SimpleNamespace(layers=[SimpleNamespace(surface=sf, arcs=arcs)])
            )
            for j, arc in enumerate(arcs):
                if arc.closed:
                    continue
                for end in (0, 1):
                    if not arc.folds_resolved[end]:
                        continue
                    p = arc.endpoints[end]
                    H, D3 = float(sf.exact(*p)), float(sf.exact(*p, 3))
                    if sf.sign * D3 >= 0:
                        continue
                    th, ph = np.array([p[0]]), np.array([p[1]])
                    chain = trace_chains(sf, th, ph, periods=2, samples_per_period=512)
                    for row in source["rows"]:
                        if row["file"] != args.real_file:
                            continue
                        well = wells_at_pitch(
                            sf, th, ph, np.array([H]), chain, row["b"]
                        )
                        if well.status[0] == ACTIVE:
                            result["real_fold_limits"].append(
                                dict(
                                    file=args.real_file,
                                    data_sha256=entry["sha256"],
                                    s=args.real_radius,
                                    arc=j,
                                    end=end,
                                    point=p.tolist(),
                                    H=H,
                                    D3B=D3,
                                    lambda_n=row["lambda_n"],
                                    b=row["b"],
                                    J=float(well.J[0]),
                                    K=float(well.K[0]),
                                    scope="One-sided limit evaluated at the polished degenerate point; not a strict ridge state or a corrected ledger",
                                )
                            )
            result["real_probe"] = dict(
                fold_records=records,
                report_sha256=hashlib.sha256(
                    Path(args.real_report).read_bytes()
                ).hexdigest(),
            )
        result.update(peak_rss_gb=guard.peak, snapshots=guard.snapshots)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        phi = np.linspace(-2.5, 3.2, 1600)
        for offset in (-0.1, 0, 0.1):
            theta = np.pi / 2 + offset
            axes[0].plot(
                phi,
                6 + (-2 + 0.5 * np.cos(theta)) * np.sin(phi) + np.sin(2 * phi),
                label=f"theta - pi/2 = {offset:g}",
            )
        axes[0].axhline(5, color="black", ls="--", label="b = 5")
        axes[0].plot(0, 6, "ko")
        axes[0].set(
            xlabel="phi (unwrapped)",
            ylabel="B",
            title="The well survives a fold above b",
        )
        axes[0].legend(fontsize=8)
        radial = np.linspace(0.4, 0.6, 31)
        theta = [
            brentq(
                lambda t: oracle(t, s)["J"] - target, np.pi / 2 - 0.5, np.pi / 2 + 0.5
            )
            for s in radial
        ]
        axes[1].plot(radial, np.array(theta) - np.pi / 2)
        axes[1].axhline(0, color="black", ls="--")
        axes[1].text(0.405, 0.015, "closed maximum represents the well")
        axes[1].text(0.405, -0.035, "open maximum represents the same well")
        axes[1].set(
            xlabel="s",
            ylabel="theta - pi/2",
            title=f"One constant-J curve, J = {target:.8f}",
        )
        plot = Path(args.plot)
        plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot, dpi=150)
        plt.close(fig)
    print(
        json.dumps(
            dict(
                limit=result["limit"],
                constant_J=result["constant_J"],
                real_finite_limits=len(result["real_fold_limits"]),
                peak_rss_gb=result["peak_rss_gb"],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="outputs/m2/fold-identity.json")
    parser.add_argument("--plot")
    parser.add_argument("--real-report")
    parser.add_argument("--real-file", default="n3are")
    parser.add_argument("--real-radius", type=float, default=0.7922826086956523)
    run(parser.parse_args())
