"""Interior-barrier transitions on the ridge surface (DESIGN §§3.4,4.2).

Coordinates are unwrapped (theta, phi=N_fp*zeta) radians. Chain offsets are
signed zeta radians; J and K have the input length units. Capped and failed
states remain unknown. This module does not infer connectivity from actions.
"""

from dataclasses import dataclass

import numpy as np

from .chains import trace_chains
from .wells import ACTIVE, INACTIVE, FAILED, adaptive_wells_at_pitch


@dataclass
class State:
    point: np.ndarray
    status: int
    J: float
    K: float
    z_in: float
    z_out: float
    exit_index: int
    max_z: np.ndarray
    max_H: np.ndarray
    H: float = np.nan


@dataclass
class Transition:
    resolved: bool
    point: np.ndarray
    partner: np.ndarray
    parent: float
    child_1: float
    child_3: float
    width: float
    reason: str
    height_residual: float = np.nan


def sample_states(surface, points, b, order=64, max_periods=256):
    """Project ridge guesses along field lines and evaluate bounded wells.

    Projection preserves alpha; a failed projection is an unknown state rather
    than a state on the wrong extremum. Temporary chains are
    bounded by the same 256-period ceiling as the production geometry.
    """
    points = np.asarray(points, float).reshape(-1, 2).copy()
    if len(points) > 4096:
        raise MemoryError("transition sample batch exceeds 4096 nodes")
    # Exact-series Newton polishing avoids table derivative reconstruction at
    # each bisection. Evaluations are bounded batches in Surface.exact.
    for _ in range(20):
        d1 = surface.exact(points[:, 0], points[:, 1], 1)
        d2 = surface.exact(points[:, 0], points[:, 1], 2)
        step = np.divide(d1, d2, out=np.zeros_like(d1), where=np.abs(d2) > 1e-14)
        points -= np.clip(step, -0.02, 0.02)[:, None] * np.array(
            [surface.iota, surface.nfp]
        )
        if np.max(np.abs(d1), initial=0) < 1e-12:
            break
    good = (surface.exact(points[:, 0], points[:, 1], 2) < 0) & (
        np.abs(surface.exact(points[:, 0], points[:, 1], 1)) <= 1e-8
    )
    ids = np.flatnonzero(good)
    H = surface.b(points[:, 0], points[:, 1])
    result = [
        State(
            p,
            FAILED,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1,
            np.array([np.nan]),
            np.array([np.nan]),
            float(h),
        )
        for p, h in zip(points, H)
    ]
    if not len(ids):
        return result
    selected = points[ids]
    had = surface.DB is not None
    surface.restore_derivatives()
    try:
        chain = trace_chains(surface, selected[:, 0], selected[:, 1], periods=2)
        well, chain, _ = adaptive_wells_at_pitch(
            surface,
            selected[:, 0],
            selected[:, 1],
            H[ids],
            chain,
            b,
            order=order,
            max_periods=max_periods,
        )
    finally:
        if not had:
            surface.release_derivatives()
    for i, idx in enumerate(ids):
        result[idx] = State(
            points[idx],
            int(well.status[i]),
            float(well.J[i]),
            float(well.K[i]),
            float(well.z_in[i]),
            float(well.z_out[i]),
            int(well.exit_index[i]),
            chain.max_z[i, : max(1, int(well.exit_index[i]) + 1)].copy(),
            chain.max_H[i, : max(1, int(well.exit_index[i]) + 1)].copy(),
            float(H[idx]),
        )
    return result


def interior_barrier(a, b, fraction=0.25):
    """Return (long-side index, interior index), or None for a regular edge.

    The shorter well's exit must occur inside the longer well's chain,
    within a fraction of its local maximum spacing. A bounce-point jump or
    an exit-index change alone is insufficient (DESIGN §4.2.3).
    """
    if a.status != ACTIVE or b.status != ACTIVE:
        return None
    long_side = int(abs(b.z_out) > abs(a.z_out))
    short, long = (a, b) if long_side else (b, a)
    n = long.exit_index
    if n <= 0:
        return None
    z = long.max_z[: n + 1]
    previous = np.r_[0.0, z[:-1]]
    following = z[1:]
    gap = np.minimum(np.abs(z[:-1] - previous[:-1]), np.abs(following - z[:-1]))
    candidate = np.flatnonzero(
        np.abs(z[:-1] - short.max_z[short.exit_index]) < fraction * gap
    )
    if not len(candidate):
        return None
    k = int(candidate[np.argmin(np.abs(z[candidate] - short.max_z[short.exit_index]))])
    return long_side, k


def locate_transition(surface, left, right, b, tolerance=1e-8, order=64):
    """Bisect an interior-barrier edge and independently evaluate child 3.

    Banks are evaluated on their own sides, never by splitting or assigning
    J_parent-J_child1 to child 3. Additivity is a measured residual, so a
    non-additive result is unresolved. Width is Euclidean angular radians.
    """
    radial = len(left) == 3

    def sample(points):
        if not radial:
            return sample_states(surface, points, b, order)
        result = []
        for p in points:
            state = sample_states(surface(float(p[0])), [p[1:]], b, order)[0]
            state.point = np.r_[p[0], state.point]
            result.append(state)
        return result

    a, c = sample([left, right])
    failure = lambda why, width: Transition(
        False,
        (a.point + c.point) / 2,
        np.full_like(a.point, np.nan),
        np.nan,
        np.nan,
        np.nan,
        width,
        why,
    )
    found = interior_barrier(a, c)
    if found is None:
        return failure(
            "no interior-barrier bracket", float(np.linalg.norm(c.point - a.point))
        )
    for _ in range(48):
        width = float(np.linalg.norm(c.point - a.point))
        if width <= tolerance:
            break
        advanced = False
        # At an exactly sampled separatrix K diverges although J is finite.
        # Keep two finite banks by trying off-centre points; never coerce that
        # failed K to a regular state or zero weight.
        for fraction in (0.5, 0.375, 0.625):
            mid = sample([a.point + fraction * (c.point - a.point)])[0]
            if mid.status != ACTIVE:
                continue
            first, second = interior_barrier(a, mid), interior_barrier(mid, c)
            if first is not None and second is None:
                c = mid
                advanced = True
                break
            if second is not None and first is None:
                a = mid
                advanced = True
                break
        if not advanced:
            return failure("multiple or unresolved transitions in edge", width)
    found = interior_barrier(a, c)
    if found is None:
        return failure("lost interior-barrier bracket", width)
    side, k = found
    long, short = (c, a) if side else (a, c)
    # A symmetric root can coincide with a bank to machine precision. Move
    # the evaluation a fraction of the final bracket into the short bank;
    # event location/width still come from the original enclosing bracket.
    # Additivity must pass with these independently evaluated finite banks.
    finite_short = sample([short.point + 0.5 * (short.point - long.point)])[0]
    if (
        finite_short.status == ACTIVE
        and interior_barrier(long, finite_short) is not None
    ):
        short = finite_short
    # The short bank has H_k>b; start its other child at that very maximum.
    # Match the physical barrier using z, rather than assuming a common index
    # across folds of the chain.
    j = int(np.nanargmin(np.abs(short.max_z - long.max_z[k])))
    sf = surface(float(short.point[0])) if radial else surface
    flow = np.array([0, sf.iota, sf.nfp] if radial else [sf.iota, sf.nfp])
    y = short.point + short.max_z[j] * flow
    child = sample([y])[0]
    child_action = child.J if child.status == ACTIVE else _action_only(sf, child, b)
    if not np.isfinite(child_action):
        return failure("child 3 is unresolved", width)
    residual = abs(long.J - short.J - child_action)
    if residual > 2e-6 * max(abs(long.J), 1e-12):
        return failure("bank actions fail additivity", width)
    point = (a.point + c.point) / 2
    critical_sf = surface(float(point[0])) if radial else surface
    angles = point[-2:].copy()
    flow = np.array([critical_sf.iota, critical_sf.nfp])
    for _ in range(12):
        d1 = float(critical_sf.exact(*angles, 1))
        d2 = float(critical_sf.exact(*angles, 2))
        if abs(d1) < 1e-12:
            break
        if d2 >= 0:
            return failure("critical starting ridge is unresolved", width)
        angles -= np.clip(d1 / d2, -0.02, 0.02) * flow
    if np.linalg.norm(angles - point[-2:]) > max(width, tolerance):
        return failure("critical ridge projection exceeds location bracket", width)
    if (
        abs(float(critical_sf.exact(*angles, 1))) > 1e-8
        or float(critical_sf.exact(*angles, 2)) >= 0
    ):
        return failure("critical starting maximum did not polish", width)
    point = point.copy()
    point[-2:] = angles
    z = float(long.max_z[k])
    for _ in range(12):
        y = angles + z * flow
        d1 = float(critical_sf.exact(*y, 1))
        d2 = float(critical_sf.exact(*y, 2))
        if d2 >= 0:
            return failure("critical interior maximum is unresolved", width)
        if abs(d1) < 1e-12:
            break
        z -= np.clip(d1 / d2, -0.02, 0.02)
    y = angles + z * flow
    previous = 0.0 if k == 0 else long.max_z[k - 1]
    gap = min(abs(long.max_z[k] - previous), abs(long.max_z[k + 1] - long.max_z[k]))
    if abs(z - long.max_z[k]) > 0.25 * gap:
        return failure("critical partner changed maximum identity", width)
    if abs(float(critical_sf.exact(*y, 1))) > 1e-8:
        return failure("critical interior maximum did not polish", width)
    partner = np.r_[point[0], y] if radial else y
    return Transition(
        True,
        point,
        partner,
        long.J,
        short.J,
        child_action,
        width,
        "",
        float(critical_sf.exact(*y)) - b,
    )


def _action_only(surface, state, b, order=64):
    """Resolve finite J when K quadrature loses an endpoint to roundoff.

    This never repairs a state's status or K. Only a known finite well with
    resolved bounce points is eligible, and any negative radicand beyond
    floating-point roundoff makes its action unknown. At a marginal maximum
    J has a finite limit even though K diverges.
    """
    if state.exit_index < 0 or not np.isfinite(state.z_in + state.z_out):
        return np.nan
    count = max(1, int(np.ceil(abs(state.z_out - state.z_in) / surface.L * 16)))
    if count > 4096:
        raise MemoryError("action-only segment ceiling exceeded")
    breaks = np.linspace(state.z_in, state.z_out, count + 1)
    x, w = np.polynomial.legendre.leggauss(order)
    t = x * np.pi / 2
    total = 0.0
    for start in range(0, count, 128):
        left = breaks[:-1][start : start + 128, None]
        right = breaks[1:][start : start + 128, None]
        half = (right - left) / 2
        z = (right + left) / 2 + half * np.sin(t)
        B = surface.b(
            state.point[-2] + surface.iota * z, state.point[-1] + surface.nfp * z
        )
        g = 1 - B / b
        if np.min(g) < -64 * np.finfo(float).eps or not np.isfinite(g).all():
            return np.nan
        total += float(
            np.sum(
                w
                * np.pi
                / 2
                * abs(surface.C)
                / B
                * abs(half)
                * np.cos(t)
                * np.sqrt(np.maximum(g, 0))
            )
        )
    return total


def well_assignment(a, c, expected_variation):
    """Classify a ridge edge without silently accepting incompatible wells.

    expected_variation is a length from neighbouring regular J differences,
    not a tolerance on a transition. An unresolved edge needs subdivision or
    unknown incidence; it must never become an ordinary constant-J link.
    """
    if a.status != ACTIVE or c.status != ACTIVE:
        return "boundary" if a.status == 0 or c.status == 0 else "unknown"
    if interior_barrier(a, c) is not None:
        return "transition"
    if abs(a.exit_index - c.exit_index) >= 2:
        return "subdivide"
    if abs(a.J - c.J) > max(5 * expected_variation, 1e-10):
        return "unknown"
    return "regular"


@dataclass
class SurfaceTransitions:
    """Refined ridge-edge census; unresolved edges retain their endpoints."""

    events: list
    unknown_edges: list
    transitions: int
    gamma_max: int
    factor: int
    samples: int
    failures: list
    unknown_weight_per_ds: float
    resolved_unknown_weight_per_ds: float
    total_resolved_weight_per_ds: float


def surface_transitions(surface, arcs, b, max_factor=16, order=32, tolerance=1e-8):
    """DESIGN §6.2: refine a count mismatch, report every remaining edge.

    Uniform dyadic subdivision on a mismatched surface also tests edges whose
    endpoints hide an even number of crossings. Chain jumps and assignment
    failures are retained as unknown if they survive the finite ceiling.
    """
    if max_factor not in (1, 2, 4, 8, 16):
        raise ValueError("transition refinement factor must be dyadic and <=16")
    if not np.isfinite(b) or b <= 0:
        raise ValueError("bounce field must be positive and finite")
    if (
        b <= surface.Bmin - surface.extremum_grid_error
        or b >= surface.Bmax + surface.extremum_grid_error
    ):
        # Fourier grid-extremum bounds establish an empty active set; this
        # does not infer passing from a capped or finite field-line search.
        return SurfaceTransitions([], [], 0, 0, 1, 0, [], 0.0, 0.0, 0.0)
    segments = []
    folds = {}
    unresolved_folds = []

    def key(p):
        return tuple(np.round(np.asarray(p) % (2 * np.pi), 12))

    for arc in arcs:
        for k in range(len(arc.points) if arc.closed else len(arc.points) - 1):
            a, c = arc.points[k], arc.points[(k + 1) % len(arc.points)].copy()
            if arc.closed and k == len(arc.points) - 1:
                c += 2 * np.pi * np.array(arc.wind)
            segments.append((a, c))
        if not arc.closed and arc.folds_resolved is not None:
            for end, node in ((0, 0), (1, -1)):
                if arc.folds_resolved[end]:
                    p = arc.endpoints[end]
                    folds[key(p)] = float(surface.b(*p))
                    segments.append((p, arc.points[node]))
        if not arc.closed and arc.endpoints is not None:
            for end, node in ((0, 0), (1, -1)):
                if arc.folds_resolved is None or not arc.folds_resolved[end]:
                    unresolved_folds.append((arc.endpoints[end], arc.points[node]))
    factor = 1
    while True:
        edges = []
        for a, c in segments:
            p = a + (c - a) * np.arange(factor + 1)[:, None] / factor
            edges.extend(zip(p[:-1], p[1:]))
        states = []
        points = np.asarray(edges).reshape(-1, 2)
        for j in range(0, len(points), 1024):
            batch = points[j : j + 1024]
            regular = [i for i, p in enumerate(batch) if key(p) not in folds]
            sampled = sample_states(surface, batch[regular], b, order)
            iterator = iter(sampled)
            for p in batch:
                if key(p) in folds:
                    # A verified generic fold is not a trapped starting well:
                    # its first maximum and minimum coincide. This creates an
                    # explicit boundary state, never repairs a failed well.
                    states.append(
                        State(
                            p,
                            INACTIVE,
                            np.nan,
                            np.nan,
                            np.nan,
                            np.nan,
                            -1,
                            np.array([np.nan]),
                            np.array([np.nan]),
                            folds[key(p)],
                        )
                    )
                else:
                    states.append(next(iterator))
        pairs = list(zip(states[::2], states[1::2]))
        nt = sum(interior_barrier(a, c) is not None for a, c in pairs)
        ng = sum((a.H > b) != (c.H > b) for a, c in pairs)
        differences = [
            abs(a.J - c.J)
            for a, c in pairs
            if a.status == c.status == ACTIVE and interior_barrier(a, c) is None
        ]
        global_scale = float(np.median(differences)) if differences else 0.0
        changes = np.array(
            [
                (
                    abs(a.J - c.J)
                    if a.status == c.status == ACTIVE and interior_barrier(a, c) is None
                    else np.nan
                )
                for a, c in pairs
            ]
        )
        scales = []
        for j in range(len(pairs)):
            local = changes[max(0, j - 5) : j + 6]
            local = local[np.isfinite(local)]
            scales.append(
                max(float(np.median(local)) if len(local) else 0.0, global_scale)
            )
        bad = [
            (a, c)
            for (a, c), scale in zip(pairs, scales)
            if well_assignment(a, c, scale) in ("subdivide", "unknown")
        ]
        events, unknown, failures = [], list(bad), []
        for a, c in pairs:
            if interior_barrier(a, c) is not None:
                event = locate_transition(
                    surface, a.point, c.point, b, tolerance, order
                )
                if event.resolved:
                    events.append(event)
                else:
                    unknown.append((a, c))
                    failures.append(
                        dict(
                            reason=event.reason,
                            left=a.point.tolist(),
                            right=c.point.tolist(),
                        )
                    )
        if factor == max_factor or (nt == ng and not bad and not failures):
            break
        factor *= 2
    # A residual count miss gives no location for the missing transition.
    # Account the entire surface as unresolved, rather than choose an edge.
    if nt != ng:
        unknown = pairs
    for fold, node in unresolved_folds:
        state = State(
            fold,
            FAILED,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1,
            np.array([np.nan]),
            np.array([np.nan]),
            float(surface.b(*fold)),
        )
        other = sample_states(surface, [node], b, order)[0]
        unknown.append((state, other))
        failures.append(
            dict(
                reason="unresolved fold endpoint",
                left=fold.tolist(),
                right=node.tolist(),
            )
        )

    def population(edges):
        total = 0.0
        complete = True
        for a, c in edges:
            width = abs(
                c.point[0]
                - a.point[0]
                - surface.iota * (c.point[1] - a.point[1]) / surface.nfp
            )
            for state in (a, c):
                if state.status == ACTIVE and np.isfinite(state.K):
                    total += width * state.K / 2
                elif state.status != INACTIVE:
                    complete = False
        return total, complete

    unknown_resolved, complete = population(unknown)
    total, _ = population(pairs)
    return SurfaceTransitions(
        events,
        unknown,
        nt,
        ng,
        factor,
        len(states),
        failures,
        unknown_resolved if complete else np.nan,
        unknown_resolved,
        total,
    )
