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


def plot_transition_cells(result, field, path, title="Transition-cut ridge cells"):
    """Show cell banks, the child-3 ridge point and independently checked J.

    Diagnostics only: polygons use unwrapped theta, while ridge points on the
    angular torus are reduced modulo 2pi. No reachability is inferred here.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    colours = {
        "regular": "#94c5e0",
        "bank": "#48a86b",
        "unknown": "#d86b5a",
        "inactive": "#dddddd",
    }
    polygons = [c.points[:, [0, 1]] for c in result.cells]
    axes[0].add_collection(
        PolyCollection(
            polygons,
            facecolors=[colours[c.kind] for c in result.cells],
            edgecolors="#333333",
            linewidths=0.15,
        )
    )
    axes[0].autoscale_view()
    axes[0].set(xlabel="s = rho²", ylabel="unwrapped theta [rad]", title=title)
    axes[0].legend(
        handles=[Patch(color=v, label=k) for k, v in colours.items()], fontsize=8
    )
    events = [e for segment in result.transitions for e in segment]
    for e in events:
        x = e.point[-2:] % (2 * np.pi)
        y = e.partner[-2:] % (2 * np.pi)
        axes[1].plot(x[0], x[1], "ko", ms=3)
        axes[1].plot(y[0], y[1], "rx", ms=4)
        if np.linalg.norm(x - y) < np.pi:
            axes[1].plot([x[0], y[0]], [x[1], y[1]], "k--", lw=0.5)
    axes[1].set(
        xlabel="theta mod 2pi [rad]",
        ylabel="phi mod 2pi [rad]",
        title="Transition (black), child 3 (red)",
        xlim=(0, 2 * np.pi),
        ylim=(0, 2 * np.pi),
    )
    for k, e in enumerate(events):
        axes[2].plot(k, e.parent, "ko", ms=3)
        axes[2].plot(k, e.child_1, "b.", ms=3)
        axes[2].plot(k, e.child_3, "r.", ms=3)
        axes[2].plot(k, e.child_1 + e.child_3, "g+", ms=5)
    axes[2].set(
        xlabel="Located transition endpoint",
        ylabel="J [input length units]",
        title="Parent (black), children (blue/red), sum (+)",
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)
