"""Diagnostic counts and plots for ridge geometry and chains (DESIGN §9 M1).

Raw T counts are probes of the interior-barrier test, not a certification of
located transitions. Bisection, subdivision and residual unknown cells are M2.
"""

import numpy as np

from .wells import ACTIVE


def raw_transition_counts(layer, well, b):
    """Report (raw T, Gamma_max) per surface before M2 subdivision."""
    transitions = gamma = offset = 0
    sign = layer.surface.sign
    for arc in layer.arcs:
        n = len(arc.H)
        for k in range(n if arc.closed else n - 1):
            a, c = offset + k, offset + (k + 1) % n
            gamma += int((layer.H[a] > b) != (layer.H[c] > b))
            if well.status[a] != ACTIVE or well.status[c] != ACTIVE:
                continue
            za = layer.chain.max_z[a, well.exit_index[a]] * sign
            zc = layer.chain.max_z[c, well.exit_index[c]] * sign
            short, long = (a, c) if za < zc else (c, a)
            zs = min(za, zc)
            maxima = layer.chain.max_z[long] * sign
            before = maxima[: well.exit_index[long]]
            if len(before):
                previous = np.r_[0, maxima[:-1]]
                following = np.r_[maxima[1:], np.nan]
                spacing = np.fmin(maxima - previous, following - maxima)
                near = np.abs(before - zs) < 0.25 * spacing[: len(before)]
                transitions += int(np.any(near))
        offset += n
    return transitions, gamma


def plot_geometry(geometry, path, b=None):
    """Ridges on the angular torus plus a signed-zeta extremum chain."""
    import matplotlib.pyplot as plt

    layer = geometry.layers[len(geometry.layers) // 2]
    sf = layer.surface
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axes[0].imshow(
        sf.B.T,
        origin="lower",
        extent=(0, 2 * np.pi, 0, 2 * np.pi),
        aspect="auto",
        cmap="viridis",
    )
    for arc in layer.arcs:
        pts = arc.points.copy() % (2 * np.pi)
        jumps = np.flatnonzero(np.linalg.norm(np.diff(pts, axis=0), axis=1) > np.pi) + 1
        for piece in np.split(pts, jumps):
            axes[0].plot(piece[:, 0], piece[:, 1], "r.-", ms=1, lw=0.8)
    axes[0].set(
        xlabel="theta [rad]",
        ylabel="phi = N_fp zeta [rad]",
        title=f"Ridge maxima, s={sf.s:.3g}",
    )
    if len(layer.points):
        theta, phi = layer.points[0]
        z = np.linspace(0, sf.sign * 6 * sf.L, 2400)
        axes[1].plot(z, sf.b(theta + sf.iota * z, phi + sf.nfp * z), lw=1)
        for zs, hs, marker in (
            (layer.chain.max_z[0], layer.chain.max_H[0], "ro"),
            (layer.chain.min_z[0], layer.chain.min_H[0], "bo"),
        ):
            keep = np.isfinite(zs) & (np.abs(zs) <= 6 * sf.L)
            axes[1].plot(zs[keep], hs[keep], marker, ms=3)
        if b is not None:
            axes[1].axhline(b, color="k", ls="--", label="bounce field b")
            axes[1].legend()
    axes[1].set(
        xlabel="signed zeta offset [rad]",
        ylabel="B [input units]",
        title=f"Forward chain, sign(C)={sf.sign:+g}",
    )
    fig.savefig(path, dpi=140)
    plt.close(fig)
