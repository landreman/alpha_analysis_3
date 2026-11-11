"""Plot dyadic quadrature nodes around an originally capped PCA well.

This is a single-surface diagnostic using the production tables/chains/wells;
it neither computes f nor fills unknown weights from a reference integral.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.chains import trace_chains
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.tables import Surface
from alpha_analysis_3.wells import ACTIVE, CAPPED, FAILED, adaptive_wells_at_pitch
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def main(output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][0]
    s, b = 0.8295652173913044, 11.81493985662058
    target = np.array([15.982749366888166, 13.057817282260391])
    with Watchdog():
        sf = Surface(BoozerField(ROOT / "data" / entry["file"]), s)
        arcs = place_nodes(sf, ridge_arcs(sf), 128 / (2 * np.pi))
        arc, index = min(
            (
                (a, int(np.argmin(np.linalg.norm(a.points - target, axis=1))))
                for a in arcs
            ),
            key=lambda pair: np.linalg.norm(pair[0].points[pair[1]] - target),
        )
        if index == 0 or index == len(arc.H) - 1:
            raise ValueError("diagnostic target requires two interior adjacent edges")
        endpoints = arc.points[index - 1 : index + 2]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
        axes[0].plot(
            endpoints[:, 0] - target[0], endpoints[:, 1] - target[1], "k-", lw=1
        )
        for factor in (1, 2, 4, 8, 16):
            guesses = np.concatenate(
                [
                    a[None, :]
                    + (c - a)[None, :] * np.linspace(0, 1, factor + 1)[:, None]
                    for a, c in zip(endpoints[:-1], endpoints[1:])
                ]
            )
            for _ in range(12):
                d1, d2 = sf.db(guesses[:, 0], guesses[:, 1]), sf.d2b(
                    guesses[:, 0], guesses[:, 1]
                )
                step = np.divide(
                    d1, d2, out=np.zeros_like(d1), where=np.abs(d2) > 1e-14
                )
                guesses -= np.clip(step, -0.02, 0.02)[:, None] * np.array(
                    [sf.iota, sf.nfp]
                )
            assert np.all(sf.d2b(guesses[:, 0], guesses[:, 1]) < 0)
            H = sf.b(guesses[:, 0], guesses[:, 1])
            ch = trace_chains(sf, guesses[:, 0], guesses[:, 1], periods=64)
            well, extended, _ = adaptive_wells_at_pitch(
                sf, guesses[:, 0], guesses[:, 1], H, ch, b
            )
            assert not np.any(well.status == FAILED)
            if factor == 1:
                center = int(np.argmin(np.linalg.norm(guesses - target, axis=1)))
                exit_index = well.exit_index[center]
                periods = np.abs(extended.max_z[center, : exit_index + 1]) / sf.L
                heights = extended.max_H[center, : exit_index + 1]
                axes[2].plot(periods, heights, ".-", ms=2, lw=0.7)
                axes[2].axhline(b, color="k", ls="--", label="b")
                axes[2].axvline(64, color="r", ls=":", label="base cap")
                axes[2].plot(periods[-1], heights[-1], "go", label="first exit")
            alpha = guesses[:, 0] - sf.iota * guesses[:, 1] / sf.nfp
            alpha -= target[0] - sf.iota * target[1] / sf.nfp
            active = well.status == ACTIVE
            axes[1].plot(
                alpha[active], well.K[active], ".-", ms=4, lw=0.7, label=f"{factor}x"
            )
            cap = well.status == CAPPED
            if np.any(cap):
                axes[1].scatter(
                    alpha[cap],
                    well.partial_K[cap],
                    marker="^",
                    label=f"{factor}x capped lower bound",
                )
            if factor == 16:
                axes[0].scatter(
                    guesses[:, 0] - target[0],
                    guesses[:, 1] - target[1],
                    c=alpha,
                    cmap="viridis",
                    s=15,
                )
        axes[0].plot(0, 0, "ro", label="original 115-period node")
        axes[0].set(
            xlabel="theta offset [rad]",
            ylabel="phi offset [rad]",
            title="Refined ridge edges, PCA",
        )
        axes[0].legend(fontsize=8)
        axes[1].set(
            xlabel="alpha offset [rad]",
            ylabel="K [length units]",
            title=f"Transverse well time, s={s:.3f}, lambda=0.95",
        )
        axes[1].legend(fontsize=8)
        axes[2].set(
            xlabel="field periods from ridge node",
            ylabel="maximum B [input units]",
            title="Original long well: first return above b",
        )
        axes[2].legend(fontsize=8)
        fig.savefig(output, dpi=140)
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/m1-adaptive-strip.png")
    main(parser.parse_args().output)
