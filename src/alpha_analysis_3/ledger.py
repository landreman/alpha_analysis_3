"""M1 ridge-strip ledger and independent volume quadrature (DESIGN §§2.3,6).

Weights use s=rho², radians and the input B/current units. The ledger integrates
K |d alpha| ds over all states, without a connectivity classification. M2 will
replace strips with fully covered transition-cut cells. Radial matching is
stored but never used to discard birth weight here.
"""

from dataclasses import dataclass
import time

import numpy as np
from scipy.integrate import trapezoid

from .chains import ChainBudget, trace_chains
from .conventions import source_at_s
from .ridge import place_nodes, radial_partners, ridge_arcs
from .tables import Surface, TWO_PI, global_extrema
from .wells import (
    ACTIVE,
    CAPPED,
    FAILED,
    INACTIVE,
    wells_at_pitch,
    adaptive_wells_at_pitch,
)


@dataclass
class Layer:
    surface: Surface
    arcs: list
    chain: object
    points: np.ndarray
    H: np.ndarray
    measure: np.ndarray
    chain_budget: object = None


@dataclass
class PitchMeasure:
    """Per-surface weights, before source reweighting; DESIGN §6 / ADR 0001.

    Columns are resolved, capped lower bound, resolved tails >48/>64/>128
    periods. A capped lower bound is not an estimate of the full missing weight.
    """

    weights: np.ndarray
    wells: list
    diagnostics: list
    capped_nodes: int
    failed_nodes: int
    converged: bool
    tolerance: float = 0.01
    edge_weights: list | None = None


def _point_weights(surface, well, chain):
    result = np.zeros((len(well.status), 5))
    active = well.status == ACTIVE
    result[active, 0] = well.K[active]
    capped = well.status == CAPPED
    result[capped, 1] = well.partial_K[capped]
    indices = np.flatnonzero(active)
    periods = np.zeros(len(well.status))
    periods[indices] = (
        np.abs(chain.max_z[indices, well.exit_index[indices]]) / surface.L
    )
    for column, limit in enumerate((48, 64, 128), 2):
        ids = np.flatnonzero(active & (periods > limit))
        result[ids, column] = well.K[ids]
    return result, periods


def _relative_strip_change(previous, current, known):
    """Sum magnitudes before summing strips; opposite errors cannot cancel."""
    difference = np.sum((current - previous)[:, :2], axis=1)
    return float(np.sum(np.abs(difference))) / max(abs(known), 1e-300)


def _strip_integrals(widths, values, factor, rule):
    """Integrate sampled K d alpha with positive dyadic quadrature weights."""
    if rule == "trapezoid":
        return (
            widths[:, None]
            * (values[:, 0] / 2 + values[:, -1] / 2 + values[:, 1:-1].sum(axis=1))
            / factor
        )
    return (
        widths[:, None]
        * (
            values[:, 0]
            + values[:, -1]
            + 4 * values[:, 1:-1:2].sum(axis=1)
            + 2 * values[:, 2:-1:2].sum(axis=1)
        )
        / (3 * factor)
    )


def slice_refinement_error(s, weights, diagnostics, h):
    """Source-weighted full-slice estimate from absolute edge changes.

    Magnitudes are summed before applying source and radial weights. Opposite
    strips or surfaces cannot cancel. Local convergence, including a 16x
    ceiling, remains separately reported. This estimates changes in resolved
    plus partial-lower weight, not the unknown full weight of capped states.
    """
    known = np.sum(weights[:, :2], axis=1)
    changes = np.array([d.get("relative_change", 0.0) for d in diagnostics]) * known
    source = source_at_s(h, s)
    total = float(trapezoid(source * known, s))
    error = float(trapezoid(source * changes, s))
    return error, error / max(abs(total), 1e-300)


def _refine_layer(
    layer,
    well,
    b,
    tolerance,
    max_factor,
    max_periods,
    order,
    strip_rule,
    edge_output=None,
    clip_folds=False,
):
    """Dyadic positive quadrature on long-well and activity-boundary edges.

    Every base node remains in the quadrature. Convergence is the change in
    selected strip weight divided by the complete surface's known weight;
    applying the slice tolerance on every surface also bounds its source-
    weighted radial sum. No volume-reference values enter this calculation.
    """
    sf = layer.surface
    budget = layer.chain_budget
    point, periods = _point_weights(sf, well, layer.chain)
    base = layer.measure @ point
    long = (periods > 64) | (well.status == CAPPED)
    edges = []
    offset = 0
    for arc in layer.arcs:
        n = len(arc.H)
        for k in range(n if arc.closed else n - 1):
            a, c = offset + k, offset + (k + 1) % n
            output_id = -1
            if edge_output is not None:
                output_id = len(edge_output)
                edge_output.append((a, c, arc.dalpha[k] * (point[a] + point[c]) / 2))
            boundary = (well.status[a] != well.status[c]) or (
                (layer.H[a] > b) != (layer.H[c] > b)
            )
            if not (long[a] or long[c] or boundary):
                continue
            right = arc.points[(k + 1) % n].copy()
            if arc.closed and k == n - 1:
                right += TWO_PI * np.array(arc.wind)
            width = arc.dalpha[k]
            edges.append(
                (
                    arc.points[k],
                    right,
                    width,
                    np.stack((point[a], point[c])),
                    np.array([well.status[a], well.status[c]]),
                    np.array([a, c]),
                    width * (point[a] + point[c]) / 2,
                    output_id,
                )
            )
        if not arc.closed and arc.endpoints is not None:
            for node, fold in ((0, 0), (n - 1, 1)):
                idx = offset + node
                endpoint = arc.endpoints[fold]
                alpha = endpoint[0] - sf.iota * endpoint[1] / sf.nfp
                width = abs(alpha - arc.alpha[node])
                output_id = -1
                if edge_output is not None:
                    output_id = len(edge_output)
                    edge_output.append((-1, idx, width * point[idx]))
                resolved_fold = (
                    arc.folds_resolved is not None and arc.folds_resolved[fold]
                )
                if not (long[idx] or (clip_folds and resolved_fold)):
                    continue
                # At a fold m0=H0, so no strictly active starting well exists.
                edges.append(
                    (
                        endpoint,
                        arc.points[node],
                        width,
                        np.stack((np.zeros(5), point[idx])),
                        np.array([INACTIVE, well.status[idx]]),
                        np.array([-1, idx]),
                        width * point[idx],
                        output_id,
                    )
                )
        offset += n
    stats = dict(
        long_nodes=int(long.sum()),
        marked_edges=len(edges),
        factor=1,
        converged=True,
        relative_change=0.0,
        levels=[],
        extension_seconds=0.0,
        refinement_seconds=0.0,
        maximum_periods=int(np.max(layer.chain.periods, initial=0)),
        strip_rule=strip_rule,
    )
    caps = set(np.flatnonzero(well.status == CAPPED).tolist())
    failed = int(np.sum(well.status == FAILED))
    if not edges:
        return base, stats, len(caps), failed
    original_edges = np.array([e[6] for e in edges])
    original = original_edges.sum(axis=0)
    unmarked = base - original
    widths = np.array([e[2] for e in edges])
    left = np.array([e[0] for e in edges])
    right = np.array([e[1] for e in edges])
    values = np.array([e[3] for e in edges])
    statuses = np.array([e[4] for e in edges])
    previous = original_edges
    start = time.monotonic()
    had_derivatives = sf.DB is not None
    sf.restore_derivatives()
    try:
        factor = 2
        stats["converged"] = False
        while factor <= max_factor:
            fractions = np.arange(1, factor, 2) / factor
            guesses = (
                left[:, None, :] + (right - left)[:, None, :] * fractions[None, :, None]
            )
            nodes = guesses.reshape(-1, 2).copy()
            for _ in range(12):
                d1, d2 = sf.db(nodes[:, 0], nodes[:, 1]), sf.d2b(
                    nodes[:, 0], nodes[:, 1]
                )
                step = np.divide(
                    d1, d2, out=np.zeros_like(d1), where=np.abs(d2) > 1e-14
                )
                nodes -= np.clip(step, -0.02, 0.02)[:, None] * np.array(
                    [sf.iota, sf.nfp]
                )
            H = sf.b(nodes[:, 0], nodes[:, 1])
            valid = (sf.d2b(nodes[:, 0], nodes[:, 1]) < 0) & (
                np.abs(sf.db(nodes[:, 0], nodes[:, 1])) < 1e-9 * np.maximum(H, 1)
            )
            # Auxiliary quadrature points are pitch-specific. H<=b is already
            # inactive. For the others stop once the first exit is known,
            # doubling a short initial scan through 64 and then to the same
            # 256-period ceiling. Base geometry chains still start at 64.
            ids = np.flatnonzero(valid & (H > b))
            new_values = np.zeros((len(nodes), 5))
            new_status = np.full(len(nodes), INACTIVE, dtype=np.int8)
            new_status[~valid] = FAILED
            extension = dict(maximum_periods=0, extension_seconds=0.0)
            if len(ids):
                initial = min(2, int(np.min(layer.chain.periods)), max_periods)
                ch = trace_chains(
                    sf,
                    nodes[ids, 0],
                    nodes[ids, 1],
                    initial,
                    layer.chain.samples_per_period,
                    allocation_check=budget.check_new if budget is not None else None,
                )
                new, ch, extension = adaptive_wells_at_pitch(
                    sf,
                    nodes[ids, 0],
                    nodes[ids, 1],
                    H[ids],
                    ch,
                    b,
                    order,
                    max_periods,
                    allocation_check=(
                        budget.check_extension if budget is not None else None
                    ),
                )
                new_values[ids] = _point_weights(sf, new, ch)[0]
                new_status[ids] = new.status
                del ch, new
            refined = np.empty((len(edges), factor + 1, 5))
            refined[:, ::2] = values
            refined[:, 1::2] = new_values.reshape(len(edges), len(fractions), 5)
            state = np.empty((len(edges), factor + 1), dtype=np.int8)
            state[:, ::2] = statuses
            state[:, 1::2] = new_status.reshape(len(edges), len(fractions))
            values, statuses = refined, state
            per_edge = _strip_integrals(widths, values, factor, strip_rule)
            strip = per_edge.sum(axis=0)
            known = float(np.sum((unmarked + strip)[:2]))
            change = _relative_strip_change(previous, per_edge, known)
            stats["levels"].append(
                dict(
                    factor=factor, strip_weights=strip.tolist(), relative_change=change
                )
            )
            stats.update(
                factor=factor,
                relative_change=change,
                maximum_periods=max(
                    stats["maximum_periods"], extension["maximum_periods"]
                ),
            )
            stats["extension_seconds"] += extension["extension_seconds"]
            # Compare two subdivided levels, retaining all original samples.
            if factor >= 4 and change < tolerance and not np.any(statuses == FAILED):
                stats["converged"] = True
                break
            previous = per_edge
            factor *= 2
    finally:
        if not had_derivatives:
            sf.release_derivatives()
    stats["refinement_seconds"] = time.monotonic() - start
    internal_caps = int(np.sum(statuses[:, 1:-1] == CAPPED))
    failed += int(np.sum(statuses[:, 1:-1] == FAILED))
    if edge_output is not None:
        for edge, value in zip(edges, per_edge):
            a, c = edge[5]
            edge_output[edge[7]] = (int(a), int(c), value.copy())
    return unmarked + strip, stats, len(caps) + internal_caps, failed


def measure_pitch(
    geometry,
    b,
    order=64,
    max_periods=256,
    transverse_tolerance=0.01,
    max_factor=16,
    refine=True,
    acceptance_tolerance=None,
    strip_rule="trapezoid",
    retain_edges=False,
    clip_folds=False,
):
    """Evaluate reusable ridge weights with ADR 0001 adaptive policy.

    The independent volume integral is deliberately absent. Source profiles
    and volume-reference resolution can be changed without recomputing wells.
    transverse_tolerance controls local stopping; acceptance_tolerance, when
    supplied, keeps the source-weighted slice gate fixed under tighter control.
    """
    if max_factor not in (4, 8, 16):
        raise ValueError("transverse refinement factor must be 4, 8 or 16")
    if retain_edges and not refine:
        raise ValueError("retained cell-edge weights require refinement")
    if strip_rule not in ("trapezoid", "simpson"):
        raise ValueError("strip rule must be trapezoid or simpson")
    if not 0 < transverse_tolerance < 1:
        raise ValueError("transverse tolerance must be between zero and one")
    accepted = (
        transverse_tolerance if acceptance_tolerance is None else acceptance_tolerance
    )
    if not transverse_tolerance <= accepted < 1:
        raise ValueError("quadrature control must be at least as strict as acceptance")
    weights = []
    wells = []
    diagnostics = []
    edge_weights = [] if retain_edges else None
    capped = failed = 0
    for layer in geometry.layers:
        budget = layer.chain_budget
        well, layer.chain, extension = adaptive_wells_at_pitch(
            layer.surface,
            layer.points[:, 0],
            layer.points[:, 1],
            layer.H,
            layer.chain,
            b,
            order,
            max_periods,
            allocation_check=budget.check_extension if budget is not None else None,
        )
        if budget is not None:
            budget.track(layer, layer.chain)
        if refine:
            edge_output = [] if retain_edges else None
            value, stats, cap, bad = _refine_layer(
                layer,
                well,
                b,
                transverse_tolerance,
                max_factor,
                max_periods,
                order,
                strip_rule,
                edge_output,
                clip_folds,
            )
            if retain_edges:
                edge_weights.append(edge_output)
        else:
            value = layer.measure @ _point_weights(layer.surface, well, layer.chain)[0]
            stats = dict(
                factor=1,
                converged=True,
                long_nodes=int(np.sum(well.status == CAPPED)),
                marked_edges=0,
                levels=[],
                refinement_seconds=0.0,
                extension_seconds=0.0,
                maximum_periods=extension["maximum_periods"],
            )
            cap = int(np.sum(well.status == CAPPED))
            bad = int(np.sum(well.status == FAILED))
        stats["base_extension"] = extension
        weights.append(value)
        wells.append(well)
        diagnostics.append(stats)
        capped += cap
        failed += bad
    return PitchMeasure(
        np.asarray(weights),
        wells,
        diagnostics,
        capped,
        failed,
        all(s["converged"] for s in diagnostics),
        accepted,
        edge_weights,
    )


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
        periods=64,
        samples_per_period=128,
        progress=None,
        chain_budget_bytes=768 * 1024**2,
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
        self.chain_budget = ChainBudget(chain_budget_bytes)
        for s in self.s:
            sf = Surface(field, s, nt, nphi)
            arcs = place_nodes(sf, ridge_arcs(sf), density)
            points = np.vstack([a.points for a in arcs]) if arcs else np.empty((0, 2))
            heights = np.concatenate([a.H for a in arcs]) if arcs else np.empty(0)
            ch = trace_chains(
                sf,
                points[:, 0],
                points[:, 1],
                periods,
                samples_per_period,
                allocation_check=self.chain_budget.check_new,
            )
            sf.release_derivatives()
            layer = Layer(
                sf, arcs, ch, points, heights, node_measure(arcs), self.chain_budget
            )
            self.layers.append(layer)
            self.chain_budget.track(layer, ch)
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


def volume_surface(surface, b, nphi=512, order=64, theta_factor=1):
    """Independent angular integral of |C|/(B sqrt(1-B/b)) [B<b].

    Integrate theta at fixed phi, locating ALL crossings on the table grid;
    a sine substitution removes the square-root endpoint singularities.
    This uses no ridges, chains, well K, partners or ridge measures. A separate
    phi grid and order can be refined independently of the ridge ledger.
    """
    if theta_factor not in (1, 2, 4, 8, 16):
        raise ValueError("volume theta refinement must be dyadic and <=16")
    if surface.Bmin - surface.extremum_grid_error >= b:
        return 0.0
    if surface.Bmax <= b:
        return 0.0
    phi = (np.arange(nphi) + 0.5) * TWO_PI / nphi
    count = surface.nt * theta_factor
    theta = np.arange(count + 1) * TWO_PI / count
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


def pitch_ledger(
    geometry,
    b,
    h,
    volume_nphi=512,
    order=64,
    wells=None,
    measure=None,
    volume_values=None,
    max_periods=256,
    transverse_tolerance=0.01,
    max_factor=16,
    refine=True,
    acceptance_tolerance=None,
    strip_rule="trapezoid",
):
    """Return ridge/volume weights and explicit unresolved diagnostics.

    A capped or failed well makes the ridge total incomplete (NaN), rather
    than silently assigning zero weight. resolved_weight remains diagnostic.
    Partner-frontier weight is reported separately and stays in the total.
    Without precomputed wells or a measure, apply the accepted adaptive policy.
    """
    if measure is None and wells is None:
        measure = measure_pitch(
            geometry,
            b,
            order=order,
            max_periods=max_periods,
            transverse_tolerance=transverse_tolerance,
            max_factor=max_factor,
            refine=refine,
            acceptance_tolerance=acceptance_tolerance,
            strip_rule=strip_rule,
        )
    ridge, lower, volume, frontier = [], [], [], []
    capped = failed = 0
    for i, l in enumerate(geometry.layers):
        well = (
            measure.wells[i]
            if measure is not None
            else (
                wells[i]
                if wells is not None
                else wells_at_pitch(
                    l.surface, l.points[:, 0], l.points[:, 1], l.H, l.chain, b, order
                )
            )
        )
        good = well.status == ACTIVE
        capped += int(np.sum(well.status == CAPPED))
        failed += int(np.sum(well.status == FAILED))
        ridge.append(
            float(measure.weights[i, 0])
            if measure is not None
            else float(np.sum(l.measure[good] * well.K[good]))
        )
        cap = well.status == CAPPED
        lower.append(
            float(measure.weights[i, 1])
            if measure is not None
            else float(np.sum(l.measure[cap] * well.partial_K[cap]))
        )
        volume.append(
            volume_values[i]
            if volume_values is not None
            else volume_surface(l.surface, b, volume_nphi, order)
        )
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
    capped_lower = float(trapezoid(source * lower, geometry.s))
    if measure is not None:
        capped, failed = measure.capped_nodes, measure.failed_nodes
    total = resolved if capped + failed == 0 else np.nan
    return dict(
        ridge_weight=total,
        resolved_weight=resolved,
        capped_weight_lower_bound=capped_lower,
        resolved_plus_lower=resolved + capped_lower,
        gap_to_volume=independent - resolved - capped_lower,
        refinement_error=(
            slice_refinement_error(geometry.s, measure.weights, measure.diagnostics, h)[
                0
            ]
            if measure is not None
            else 0.0
        ),
        refinement_relative_error=(
            slice_refinement_error(geometry.s, measure.weights, measure.diagnostics, h)[
                1
            ]
            if measure is not None
            else 0.0
        ),
        refinement_converged=(
            slice_refinement_error(geometry.s, measure.weights, measure.diagnostics, h)[
                1
            ]
            < measure.tolerance
            if measure is not None
            else True
        ),
        all_local_refinements_converged=(
            measure.converged if measure is not None else True
        ),
        tail_weights=(
            {
                str(limit): float(trapezoid(source * measure.weights[:, j], geometry.s))
                for j, limit in enumerate((48, 64, 128), 2)
            }
            if measure is not None
            else {}
        ),
        volume_weight=independent,
        ratio=total / independent if independent > 0 else np.nan,
        frontier_weight=float(trapezoid(source * frontier, geometry.s)),
        capped_nodes=capped,
        failed_nodes=failed,
    )
