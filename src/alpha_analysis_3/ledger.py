"""M1 ridge-strip ledger and independent volume quadrature (DESIGN §§2.3,6).

Weights use s=rho², radians and the input B/current units. The ledger integrates
K |d alpha| ds over all states, without a connectivity classification. M2 will
replace strips with fully covered transition-cut cells. Radial matching is
stored but never used to discard birth weight here.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import trapezoid

from .chains import trace_chains
from .conventions import source_at_s
from .ridge import place_nodes, radial_partners, ridge_arcs
from .tables import Surface, TWO_PI, global_extrema
from .wells import ACTIVE, CAPPED, FAILED, wells_at_pitch


@dataclass
class Layer:
    surface: Surface
    arcs: list
    chain: object
    points: np.ndarray
    H: np.ndarray
    measure: np.ndarray


def node_measure(arcs):
    """Dual edge measure at each node, including each closed curve's seam."""
    measures = []
    for arc in arcs:
        if arc.measure is not None:
            measures.append(arc.measure)
        elif arc.closed:
            measures.append((arc.dalpha + np.roll(arc.dalpha, 1)) / 2)
        else:
            measures.append((np.r_[0, arc.dalpha] + np.r_[arc.dalpha, 0]) / 2)
    return np.concatenate(measures) if measures else np.empty(0)


class Geometry:
    """Pitch-independent layers and one-direction integer radial partners."""

    def __init__(
        self,
        field,
        s_grid,
        nt=512,
        nphi=512,
        density=128 / TWO_PI,
        periods=48,
        samples_per_period=128,
        progress=None,
    ):
        self.s = np.asarray(s_grid, float)
        if (
            len(self.s) < 2
            or np.any(np.diff(self.s) <= 0)
            or self.s[0] <= 0
            or self.s[-1] > 1
        ):
            raise ValueError("s grid must increase within (0,1]")
        self.layers = []
        for s in self.s:
            sf = Surface(field, s, nt, nphi)
            arcs = place_nodes(sf, ridge_arcs(sf), density)
            points = np.vstack([a.points for a in arcs]) if arcs else np.empty((0, 2))
            heights = np.concatenate([a.H for a in arcs]) if arcs else np.empty(0)
            ch = trace_chains(
                sf, points[:, 0], points[:, 1], periods, samples_per_period
            )
            sf.release_derivatives()
            self.layers.append(Layer(sf, arcs, ch, points, heights, node_measure(arcs)))
            if progress:
                progress(f"surface s={s:.6g}, nodes={len(points)}")
        self.partners = [
            radial_partners(a.arcs, b.arcs)
            for a, b in zip(self.layers, self.layers[1:])
        ]
        self.extrema = global_extrema(field)
        self.Bmin = self.extrema["min"]["B"]
        self.Bmax = self.extrema["max"]["B"]

    def volume(self, h):
        """V_h in DESIGN §2.3; denominator of f is 2 V_h."""
        vals = [
            abs(l.surface.C) * TWO_PI * l.surface.L * np.mean(1 / l.surface.B**2)
            for l in self.layers
        ]
        return float(trapezoid(source_at_s(h, self.s) * vals, self.s))


def volume_surface(surface, b, nphi=512, order=64):
    """Independent angular integral of |C|/(B sqrt(1-B/b)) [B<b].

    Integrate theta at fixed phi, locating ALL crossings on the table grid;
    a sine substitution removes the square-root endpoint singularities.
    This uses no ridges, chains, well K, partners or ridge measures. A separate
    phi grid and order can be refined independently of the ridge ledger.
    """
    if surface.Bmax <= b:
        return 0.0
    phi = (np.arange(nphi) + 0.5) * TWO_PI / nphi
    theta = np.arange(surface.nt + 1) * TWO_PI / surface.nt
    values = surface.b(theta[:, None], phi[None, :])
    crossing = (values[:-1] < b) != (values[1:] < b)
    rows, columns = np.nonzero(crossing)
    lo, hi = theta[rows].copy(), theta[rows + 1].copy()
    low = values[rows, columns] < b
    for _ in range(42):
        mid = (lo + hi) / 2
        move = (surface.b(mid, phi[columns]) < b) == low
        lo, hi = np.where(move, mid, lo), np.where(move, hi, mid)
    roots = (lo + hi) / 2
    # Merge periodic below-b intervals across theta=0. Fully below-b theta
    # lines have no singular endpoints and use periodic trapezoidal quadrature.
    left, right, ph = [], [], []
    total = 0.0
    for j in range(nphi):
        r = np.sort(roots[columns == j])
        if not len(r):
            if values[0, j] < b:
                B = values[:-1, j]
                total += float(np.mean(1 / (B * np.sqrt(1 - B / b))) * TWO_PI)
            continue
        boundaries = np.r_[r, r[0] + TWO_PI]
        for a, c in zip(boundaries[:-1], boundaries[1:]):
            if surface.b((a + c) / 2, phi[j]) < b:
                left.append(a)
                right.append(c)
                ph.append(phi[j])
    t, w = np.polynomial.legendre.leggauss(order)
    t, w = t * np.pi / 2, w * np.pi / 2
    for start in range(0, len(left), 128):
        l = np.asarray(left[start : start + 128])[:, None]
        r = np.asarray(right[start : start + 128])[:, None]
        p = np.asarray(ph[start : start + 128])[:, None]
        half = (r - l) / 2
        th = (l + r) / 2 + half * np.sin(t)
        B = surface.b(th, p)
        g = 1 - B / b
        if np.any(g <= 0):
            raise FloatingPointError("volume quadrature crossed a bounce boundary")
        total += float(np.sum(w * half * np.cos(t) / (B * np.sqrt(g))))
    return total * TWO_PI / nphi * abs(surface.C) / surface.nfp


def pitch_ledger(geometry, b, h, volume_nphi=512, order=64, wells=None):
    """Return ridge/volume weights and explicit unresolved diagnostics.

    A capped or failed well makes the ridge total incomplete (NaN), rather
    than silently assigning zero weight. resolved_weight remains diagnostic.
    Partner-frontier weight is reported separately and stays in the total.
    """
    ridge, volume, frontier = [], [], []
    capped = failed = 0
    for i, l in enumerate(geometry.layers):
        well = (
            wells[i]
            if wells is not None
            else wells_at_pitch(
                l.surface, l.points[:, 0], l.points[:, 1], l.H, l.chain, b, order
            )
        )
        good = well.status == ACTIVE
        capped += int(np.sum(well.status == CAPPED))
        failed += int(np.sum(well.status == FAILED))
        ridge.append(float(np.sum(l.measure[good] * well.K[good])))
        volume.append(volume_surface(l.surface, b, volume_nphi, order))
        missing = (
            geometry.partners[i] < 0
            if i < len(geometry.partners)
            else np.zeros(len(l.H), bool)
        )
        frontier.append(
            float(np.sum(l.measure[good & missing] * well.K[good & missing]))
        )
    source = source_at_s(h, geometry.s)
    resolved = float(trapezoid(source * ridge, geometry.s))
    independent = float(trapezoid(source * volume, geometry.s))
    total = resolved if capped + failed == 0 else np.nan
    return dict(
        ridge_weight=total,
        resolved_weight=resolved,
        volume_weight=independent,
        ratio=total / independent if independent > 0 else np.nan,
        frontier_weight=float(trapezoid(source * frontier, geometry.s)),
        capped_nodes=capped,
        failed_nodes=failed,
    )
