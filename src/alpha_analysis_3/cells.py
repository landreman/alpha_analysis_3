"""Bidirectional ridge-strip coverage (DESIGN §§4.1.5,4.2.4).

The core uses integer node IDs and NumPy arrays. Angles are unwrapped
(theta,phi=N_fp*zeta) radians; area is |ds wedge d alpha|. Unequal node
counts are stitched by triangles, so top nodes missed by forward partners
are retained. Polished folds close open-arc end strips; unresolved folds and
unmatched arcs remain explicit frontiers.
"""

from dataclasses import dataclass, replace
from collections import OrderedDict

import numpy as np

from .ridge import radial_partners
from .tables import TWO_PI
from .tables import Surface
from .transitions import (
    sample_states,
    interior_barrier,
    locate_transition,
    well_assignment,
)
from .wells import ACTIVE, INACTIVE


@dataclass
class Partners:
    forward: np.ndarray
    reverse: np.ndarray
    uncovered_top: np.ndarray


@dataclass
class StripMesh:
    vertices: np.ndarray
    faces: np.ndarray
    face_points: np.ndarray
    measure: np.ndarray
    frontier_bottom: np.ndarray
    frontier_top: np.ndarray
    fold_gaps: list
    unknown_faces: np.ndarray
    vertex_layers: np.ndarray
    vertex_nodes: np.ndarray


def match_layers(bottom, top, density=128 / TWO_PI):
    """Both nearest-node maps with threshold scaled by linear node density.

    A reverse partner is not evidence of cell coverage: uncovered_top records
    nodes not reached by the forward map, which stitching must repair.
    """
    if not np.isfinite(density) or density <= 0:
        raise ValueError("ridge density must be positive")
    tolerance = 0.4 * (128 / TWO_PI) / density
    forward = radial_partners(bottom, top, tolerance)
    reverse = radial_partners(top, bottom, tolerance)
    covered = np.zeros(len(reverse), bool)
    covered[forward[forward >= 0]] = True
    return Partners(forward, reverse, ~covered)


def _arc_ids(arcs):
    return np.repeat(np.arange(len(arcs)), [len(a.points) for a in arcs])


def _sequence(arc, start=0, reverse=False):
    n = len(arc.points)
    if arc.closed:
        order = np.r_[np.arange(start, n), np.arange(start), start]
        p = arc.points[order].copy()
        p[n - start :] += TWO_PI * np.array(arc.wind)
    else:
        order = np.arange(n)
        if reverse:
            order = order[::-1]
        p = arc.points[order].copy()
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    return order, p, distance / distance[-1]


def stitch_ridges(
    bottom, top, s_bottom, s_top, iota_bottom, iota_top, nfp=1, density=128 / TWO_PI
):
    """Triangulate mutually partnered arcs, covering both boundary polylines.

    Every retained ridge edge appears exactly once on its strip boundary.
    Reconnecting or ambiguous arcs remain frontier; callers must insert
    surfaces or account their population as unknown. No action-based matching
    is used. Local unwrapped face coordinates preserve the closed seam term.
    """
    if s_top <= s_bottom or nfp < 1:
        raise ValueError("surfaces must increase and N_fp must be positive")
    match = match_layers(bottom, top, density)
    b_ids, t_ids = _arc_ids(bottom), _arc_ids(top)
    nb, nt = len(b_ids), len(t_ids)
    bp = np.vstack([a.points for a in bottom]) if nb else np.empty((0, 2))
    tp = np.vstack([a.points for a in top]) if nt else np.empty((0, 2))
    vertices = np.vstack(
        (
            np.column_stack((np.full(nb, s_bottom), bp)),
            np.column_stack((np.full(nt, s_top), tp)),
        )
    )
    vertices = list(vertices)
    vertex_layers = [0] * nb + [1] * nt
    vertex_nodes = list(range(nb)) + list(range(nt))
    fb, ft = np.ones(nb, bool), np.ones(nt, bool)
    offsets_b = np.r_[0, np.cumsum([len(a.points) for a in bottom])]
    offsets_t = np.r_[0, np.cumsum([len(a.points) for a in top])]
    faces, coords, gaps, allowances = [], [], [], []
    for j, a in enumerate(bottom):
        votes = match.forward[offsets_b[j] : offsets_b[j + 1]]
        valid = votes[votes >= 0]
        if not len(valid):
            continue
        candidates, counts = np.unique(t_ids[valid], return_counts=True)
        k = int(candidates[np.argmax(counts)])
        if np.max(counts) <= len(votes) / 2:
            continue
        reverse = match.reverse[offsets_t[k] : offsets_t[k + 1]]
        good = reverse >= 0
        if np.sum(good & (b_ids[np.maximum(reverse, 0)] == j)) <= len(reverse) / 2:
            continue
        c = top[k]
        if (
            a.closed != c.closed
            or a.wind != c.wind
            or not ft[offsets_t[k] : offsets_t[k + 1]].all()
        ):
            continue
        start = 0
        reverse_top = False
        if a.closed:
            distance = np.linalg.norm(
                (c.points - a.points[0] + np.pi) % TWO_PI - np.pi, axis=1
            )
            start = int(np.argmin(distance))
        else:

            def distance(x, y):
                return np.linalg.norm((x - y + np.pi) % TWO_PI - np.pi)

            reverse_top = distance(a.points[0], c.points[-1]) + distance(
                a.points[-1], c.points[0]
            ) < distance(a.points[0], c.points[0]) + distance(
                a.points[-1], c.points[-1]
            )
        bo, bpts, bt = _sequence(a)
        to, tpts, tt = _sequence(c, start, reverse_top)
        spacing = max(
            np.max(np.linalg.norm(np.diff(bpts, axis=0), axis=1)),
            np.max(np.linalg.norm(np.diff(tpts, axis=0), axis=1)),
        )
        tpts += TWO_PI * np.round((bpts[0] - tpts[0]) / TWO_PI)
        ib = it = 0
        while ib < len(bo) - 1 or it < len(to) - 1:
            take_bottom = it == len(to) - 1 or (
                ib < len(bo) - 1 and bt[ib + 1] <= tt[it + 1]
            )
            if take_bottom:
                ids = [
                    offsets_b[j] + bo[ib],
                    offsets_b[j] + bo[ib + 1],
                    nb + offsets_t[k] + to[it],
                ]
                p = [
                    [s_bottom, *bpts[ib]],
                    [s_bottom, *bpts[ib + 1]],
                    [s_top, *tpts[it]],
                ]
                ib += 1
            else:
                ids = [
                    offsets_b[j] + bo[ib],
                    nb + offsets_t[k] + to[it + 1],
                    nb + offsets_t[k] + to[it],
                ]
                p = [[s_bottom, *bpts[ib]], [s_top, *tpts[it + 1]], [s_top, *tpts[it]]]
                it += 1
            faces.append(ids)
            coords.append(p)
            allowances.append(spacing)
        fb[offsets_b[j] : offsets_b[j + 1]] = (
            match.forward[offsets_b[j] : offsets_b[j + 1]] < 0
        )
        ft[offsets_t[k] : offsets_t[k + 1]] = (
            match.reverse[offsets_t[k] : offsets_t[k + 1]] < 0
        )
        if not a.closed:
            # Nodes lie at half-spacing from folds. These omitted end strips
            # are explicit unknown area, not silently covered triangles.
            resolved = (
                a.folds_resolved is not None
                and c.folds_resolved is not None
                and np.all(a.folds_resolved)
                and np.all(c.folds_resolved)
            )
            if resolved:
                for end, node in ((0, 0), (1, -1)):
                    top_end = 1 - end if reverse_top else end
                    bottom_fold = a.endpoints[end].copy()
                    top_fold = c.endpoints[top_end].copy()
                    bottom_fold += TWO_PI * np.round(
                        (bpts[node] - bottom_fold) / TWO_PI
                    )
                    top_fold += TWO_PI * np.round((tpts[node] - top_fold) / TWO_PI)
                    bi, ti = len(vertices), len(vertices) + 1
                    vertices.extend(([s_bottom, *bottom_fold], [s_top, *top_fold]))
                    vertex_layers.extend((0, 1))
                    vertex_nodes.extend((-1, -1))
                    bn, tn = offsets_b[j] + bo[node], nb + offsets_t[k] + to[node]
                    quad_ids = [bi, bn, tn, ti]
                    quad = np.array(
                        [
                            [s_bottom, *bottom_fold],
                            [s_bottom, *bpts[node]],
                            [s_top, *tpts[node]],
                            [s_top, *top_fold],
                        ]
                    )
                    for ids in ((0, 1, 2), (0, 2, 3)):
                        faces.append([quad_ids[idx] for idx in ids])
                        coords.append(quad[list(ids)])
                        allowances.append(spacing)
                continue
            for arc, off, s, iota in (
                (a, offsets_b[j], s_bottom, iota_bottom),
                (c, nb + offsets_t[k], s_top, iota_top),
            ):
                if arc.endpoints is None:
                    raise ValueError("open ridge must retain its fold endpoints")
                for end, node in ((0, 0), (1, len(arc.points) - 1)):
                    d = arc.points[node] - arc.endpoints[end]
                    area = (s_top - s_bottom) * abs(d[0] - iota * d[1] / nfp) / 2
                    gaps.append((off + node, area))
    faces = np.asarray(faces, dtype=int).reshape(-1, 3)
    coords = np.asarray(coords, dtype=float).reshape(-1, 3, 3)
    fraction = (coords[:, :, 0] - s_bottom) / (s_top - s_bottom)
    iota = iota_bottom + (iota_top - iota_bottom) * fraction
    alpha = coords[:, :, 1] - iota * coords[:, :, 2] / nfp
    ds = coords[:, 1:, 0] - coords[:, :1, 0]
    da = alpha[:, 1:] - alpha[:, :1]
    area = np.abs(ds[:, 0] * da[:, 1] - ds[:, 1] * da[:, 0]) / 2
    unknown = np.zeros(len(faces), bool)
    threshold = 0.4 * (128 / TWO_PI) / density
    for j, (ids, p) in enumerate(zip(faces, coords)):
        for k in range(3):
            nxt = (k + 1) % 3
            if p[k, 0] != p[nxt, 0]:
                distance = np.linalg.norm(
                    (p[k, 1:] - p[nxt, 1:] + np.pi) % TWO_PI - np.pi
                )
                if distance > threshold + allowances[j]:
                    unknown[j] = True
                    for idx in ids:
                        node = vertex_nodes[idx]
                        if node >= 0:
                            (fb if vertex_layers[idx] == 0 else ft)[node] = True
    return StripMesh(
        np.asarray(vertices).reshape(-1, 3),
        faces,
        coords,
        area,
        fb,
        ft,
        gaps,
        unknown,
        np.asarray(vertex_layers),
        np.asarray(vertex_nodes),
    )


class SurfaceCache:
    """Bound temporary radial tables for transition and cell subdivision."""

    def __init__(self, field, table=512, capacity=2):
        if capacity not in (1, 2, 3, 4):
            raise ValueError("temporary surface capacity must be between 1 and 4")
        self.field, self.table, self.capacity = field, table, capacity
        self.cache = OrderedDict()

    def __call__(self, s):
        if s not in self.cache:
            if len(self.cache) == self.capacity:
                self.cache.popitem(last=False)
            self.cache[s] = Surface(self.field, s, self.table, self.table)
        self.cache.move_to_end(s)
        return self.cache[s]


@dataclass
class Cell:
    """One regular, bank, inactive or unknown polygon in (s,theta,phi)."""

    points: np.ndarray
    J: np.ndarray
    K: np.ndarray
    kind: str
    depth: int
    reason: str = ""
    transition: int = -1
    parent_face: int = -1
    vertex_ids: np.ndarray | None = None
    cut_edges: np.ndarray | None = None


@dataclass
class CutCells:
    cells: list
    transitions: list
    evaluations: int
    vertex_points: np.ndarray | None = None


def cut_cells(
    surface_at,
    faces,
    b,
    max_depth=4,
    order=32,
    expected_variation=None,
    max_cells=100000,
    frontier_faces=None,
):
    """Split two-cut triangles into banks; subdivide one-/three-cut cells.

    faces are unwrapped (s,theta,phi) triples. All unresolved children are
    retained as unknown polygons at the finite depth/cell ceiling. Transition
    endpoints carry bank actions and NaN K: M3 must integrate interior K nodes,
    rather than quadrature at the logarithmic singularity itself.
    """
    if not 0 <= max_depth <= 8 or not 1 <= max_cells <= 100000:
        raise ValueError("cell subdivision controls exceed finite bounds")
    if len(faces) > max_cells:
        raise MemoryError("initial face count exceeds cell budget")
    cache = {}
    cells, events = [], []
    event_cache = {}
    evaluations = 0
    cache_bytes = 0

    def sample(p):
        nonlocal evaluations, cache_bytes
        key = tuple(p)
        if key not in cache:
            sf = surface_at(float(p[0]))
            state = sample_states(sf, [p[1:]], b, order)[0]
            state.point = np.r_[p[0], state.point]
            # Bound transient state storage independently of subdivision.
            size = state.max_z.nbytes + state.max_H.nbytes + state.point.nbytes + 128
            if len(cache) >= 16384 or cache_bytes + size > 64 * 1024**2:
                cache.clear()
                cache_bytes = 0
            cache[key] = state
            cache_bytes += size
            evaluations += 1
        return cache[key]

    def retain(p, states, kind, depth, parent_face, why="", transition=-1):
        cells.append(
            Cell(
                np.asarray(p),
                np.array([x.J for x in states]),
                np.array([x.K for x in states]),
                kind,
                depth,
                why,
                transition,
                parent_face,
                cut_edges=np.full(len(p), kind == "unknown", bool),
            )
        )

    faces = np.asarray(faces, float).reshape(-1, 3, 3)
    frontier_faces = (
        np.zeros(len(faces), bool)
        if frontier_faces is None
        else np.asarray(frontier_faces, bool)
    )
    if frontier_faces.shape != (len(faces),):
        raise ValueError("frontier mask must identify every input face")
    vertices = np.unique(faces.reshape(-1, 3), axis=0)
    # Evaluate boundary vertices in surface batches before edge bisections
    # start evicting temporary tables. Chain storage stays bounded per batch.
    for s in np.unique(vertices[:, 0]):
        selected = vertices[vertices[:, 0] == s]
        sf = surface_at(float(s))
        for start in range(0, len(selected), 512):
            points = selected[start : start + 512]
            states = sample_states(sf, points[:, 1:], b, order)
            for p, state in zip(points, states):
                state.point = np.r_[s, state.point]
                size = (
                    state.max_z.nbytes + state.max_H.nbytes + state.point.nbytes + 128
                )
                if len(cache) >= 16384 or cache_bytes + size > 64 * 1024**2:
                    cache.clear()
                    cache_bytes = 0
                cache[tuple(p)] = state
                cache_bytes += size
                evaluations += 1
    local = {}
    changes = []
    if expected_variation is None:
        for p in faces:
            states = [sample(x) for x in p]
            for j in range(3):
                a, c = states[j], states[(j + 1) % 3]
                if a.status == c.status == ACTIVE and interior_barrier(a, c) is None:
                    change = abs(a.J - c.J)
                    changes.append(change)
                    for idx in (j, (j + 1) % 3):
                        local.setdefault(tuple(p[idx]), []).append(change)
    scale = float(np.median(changes)) if changes else 0.0
    local = {key: float(np.median(value)) for key, value in local.items()}
    pending = [(np.asarray(p, float), 0, i) for i, p in enumerate(faces)]
    while pending:
        p, depth, parent_face = pending.pop()
        states = [sample(x) for x in p]
        if frontier_faces[parent_face]:
            retain(
                p, states, "unknown", depth, parent_face, "frontier or uncovered strip"
            )
            continue
        classifications = [
            well_assignment(
                states[j],
                states[(j + 1) % 3],
                (
                    expected_variation
                    if expected_variation is not None
                    else max(
                        local.get(tuple(p[j]), scale / 2**depth),
                        local.get(tuple(p[(j + 1) % 3]), scale / 2**depth),
                    )
                ),
            )
            for j in range(3)
        ]
        cuts = [j for j, k in enumerate(classifications) if k == "transition"]
        if all(x.status == INACTIVE for x in states):
            retain(p, states, "inactive", depth, parent_face)
            continue
        if not cuts and all(k == "regular" for k in classifications):
            retain(p, states, "regular", depth, parent_face)
            continue
        if (
            len(cuts) == 2
            and all(x.status == ACTIVE for x in states)
            and len(pending) + len(cells) + 2 <= max_cells
        ):
            located = {}
            for j in cuts:

                def key(point):
                    return tuple(np.round(np.r_[point[0], point[1:] % TWO_PI], 12))

                edge_key = tuple(sorted((key(p[j]), key(p[(j + 1) % 3]))))
                if edge_key not in event_cache:
                    event_cache[edge_key] = locate_transition(
                        surface_at, p[j], p[(j + 1) % 3], b, order=order
                    )
                event = event_cache[edge_key]
                midpoint = (p[j] + p[(j + 1) % 3]) / 2
                shift = np.r_[
                    0.0, TWO_PI * np.round((midpoint[1:] - event.point[1:]) / TWO_PI)
                ]
                located[j] = replace(
                    event, point=event.point + shift, partner=event.partner + shift
                )
            if all(e.resolved for e in located.values()):
                # Uncut edges connect vertices on one bank. With two cuts in
                # a triangle one vertex is isolated, the other two agree.
                lone = next(j for j in range(3) if j in cuts and (j - 1) % 3 in cuts)
                reference = interior_barrier(states[cuts[0]], states[(cuts[0] + 1) % 3])
                long_vertex = (cuts[0] + reference[0]) % 3
                lone_long = long_vertex == lone
                tag = len(events)
                events.append(tuple(located[j] for j in cuts))
                for lone_bank in (True, False):
                    bank_points, actions, times = [], [], []
                    for j in range(3):
                        if (j == lone) == lone_bank:
                            bank_points.append(p[j])
                            actions.append(states[j].J)
                            times.append(states[j].K)
                        if j in located:
                            e = located[j]
                            bank_points.append(e.point)
                            long_bank = lone_long if lone_bank else not lone_long
                            actions.append(e.parent if long_bank else e.child_1)
                            times.append(np.nan)
                    cells.append(
                        Cell(
                            np.array(bank_points),
                            np.array(actions),
                            np.array(times),
                            "bank",
                            depth,
                            "",
                            tag,
                            parent_face,
                            cut_edges=np.array(
                                [
                                    all(
                                        any(
                                            np.array_equal(v, e.point)
                                            for e in located.values()
                                        )
                                        for v in (
                                            bank_points[k],
                                            bank_points[(k + 1) % len(bank_points)],
                                        )
                                    )
                                    for k in range(len(bank_points))
                                ]
                            ),
                        )
                    )
                continue
        if depth == max_depth or len(pending) + len(cells) + 4 > max_cells:
            limit = "subdivision ceiling" if depth == max_depth else "cell budget"
            retain(
                p,
                states,
                "unknown",
                depth,
                parent_face,
                f"{limit}; {len(cuts)} cuts; edges {','.join(classifications)}",
            )
            continue
        a, c, d = p
        ac, cd, da = (a + c) / 2, (c + d) / 2, (d + a) / 2
        pending.extend(
            (np.array(q), depth + 1, parent_face)
            for q in ((a, ac, da), (ac, c, cd), (da, cd, d), (ac, cd, da))
        )
    registry = {}
    vertices = []
    for cell in cells:
        ids = []
        for p in cell.points:
            canonical = np.r_[p[0], np.mod(p[1:], TWO_PI)]
            canonical[1:][np.isclose(canonical[1:], TWO_PI, rtol=0, atol=1e-12)] = 0.0
            key = tuple(np.round(canonical, 12))
            if key not in registry:
                registry[key] = len(vertices)
                vertices.append(canonical)
            ids.append(registry[key])
        cell.vertex_ids = np.asarray(ids, int)
    return CutCells(cells, events, evaluations, np.asarray(vertices).reshape(-1, 3))


def coverage_ledger(geometry, measure, h, density=128 / TWO_PI):
    """Partition the population ledger onto all strip boundary triangles.

    M2 uses positive boundary quadrature, assigning Delta s/2 times each
    surface edge integral to its unique incident triangle. Thus every surface
    edge contributes once from each adjacent radial strip, including the
    seam. Unmeshed frontiers retain their own quadrature contribution. M3
    replaces this low-order population rule with clipped interior cell nodes.

    No independent volume or ledger gap is used to construct weights. Full
    weights stay NaN in the presence of caps/failures; partial capped weights
    are separately labelled lower bounds.
    """
    from .conventions import source_at_s

    if measure.edge_weights is None:
        raise ValueError("coverage ledger needs retained ridge-edge quadrature")
    rows = []
    resolved = partial = frontier = frontier_partial = 0.0
    full_frontier_known = True
    for i, (bottom, top) in enumerate(zip(geometry.layers, geometry.layers[1:])):
        mesh = stitch_ridges(
            bottom.arcs,
            top.arcs,
            geometry.s[i],
            geometry.s[i + 1],
            bottom.surface.iota,
            top.surface.iota,
            bottom.surface.nfp,
            density,
        )
        ownership = {}
        for face in mesh.faces:
            for a, c in zip(face, np.roll(face, -1)):
                if mesh.vertex_layers[a] == mesh.vertex_layers[c]:
                    key = (
                        int(mesh.vertex_layers[a]),
                        *sorted((int(mesh.vertex_nodes[a]), int(mesh.vertex_nodes[c]))),
                    )
                    ownership[key] = ownership.get(key, 0) + 1
        ds = geometry.s[i + 1] - geometry.s[i]
        strip_frontier = strip_resolved = 0.0
        for side, mask in enumerate((mesh.frontier_bottom, mesh.frontier_top)):
            index = i + side
            factor = ds / 2 * float(source_at_s(h, geometry.s[index]))
            for a, c, value in measure.edge_weights[index]:
                key = (side, *sorted((a, c)))
                count = ownership.get(key, 0)
                if count > 1:
                    raise ValueError("ridge boundary edge belongs to multiple cells")
                unknown = count == 0 or (a >= 0 and mask[a]) or mask[c]
                weight = factor * value[0]
                lower = factor * value[1]
                resolved += weight
                partial += lower
                strip_resolved += weight
                if unknown:
                    frontier += weight
                    frontier_partial += lower
                    strip_frontier += weight
                    if (
                        lower > 0
                        or (a >= 0 and measure.wells[index].status[a] > ACTIVE)
                        or measure.wells[index].status[c] > ACTIVE
                    ):
                        full_frontier_known = False
        rows.append(
            dict(
                s_bottom=float(geometry.s[i]),
                s_top=float(geometry.s[i + 1]),
                cells=len(mesh.faces),
                unpartnered_bottom=int(mesh.frontier_bottom.sum()),
                unpartnered_top=int(mesh.frontier_top.sum()),
                uncovered_after_stitch=int(np.sum(mesh.frontier_top)),
                resolved_weight=strip_resolved,
                frontier_resolved_weight=strip_frontier,
            )
        )
    total = (
        resolved if not measure.capped_nodes and not measure.failed_nodes else np.nan
    )
    return dict(
        cell_weight=total,
        resolved_weight=resolved,
        capped_weight_lower_bound=partial,
        frontier_weight=frontier if full_frontier_known else np.nan,
        frontier_resolved_weight=frontier,
        frontier_partial_lower_bound=frontier_partial,
        frontier_fraction=(
            frontier / resolved if full_frontier_known and resolved > 0 else np.nan
        ),
        strips=rows,
    )


def insert_surfaces(
    geometry, field, radii, density=128 / TWO_PI, max_surfaces=128, progress=None
):
    """Insert pitch-independent layers without replacing existing work.

    Chains share the geometry's global allocation budget. The finite surface
    ceiling is a cost guard: any residual frontier remains unknown.
    """
    from .chains import trace_chains
    from .ledger import Layer, node_measure
    from .ridge import ridge_arcs, place_nodes

    requested = sorted(set(float(s) for s in radii) - set(geometry.s.tolist()))
    if len(requested) + len(geometry.s) > max_surfaces:
        raise ValueError("adaptive surface ceiling exceeded")
    reference = geometry.layers[0]
    layers = list(geometry.layers)
    for s in requested:
        if not geometry.s[0] < s < geometry.s[-1]:
            raise ValueError("inserted surface must lie inside existing radial range")
        sf = Surface(field, s, reference.surface.nt, reference.surface.nphi)
        arcs = place_nodes(sf, ridge_arcs(sf), density)
        points = np.vstack([a.points for a in arcs]) if arcs else np.empty((0, 2))
        heights = np.concatenate([a.H for a in arcs]) if arcs else np.empty(0)
        chain = trace_chains(
            sf,
            points[:, 0],
            points[:, 1],
            int(np.min(reference.chain.periods)),
            reference.chain.samples_per_period,
            allocation_check=geometry.chain_budget.check_new,
        )
        sf.release_derivatives()
        layer = Layer(
            sf, arcs, chain, points, heights, node_measure(arcs), geometry.chain_budget
        )
        geometry.chain_budget.track(layer, chain)
        layers.append(layer)
        if progress:
            progress(f"inserted surface s={s:.9g}, nodes={len(points)}")
    layers.sort(key=lambda l: l.surface.s)
    geometry.layers = layers
    geometry.s = np.array([l.surface.s for l in layers])
    geometry.partners = [
        radial_partners(a.arcs, c.arcs) for a, c in zip(layers, layers[1:])
    ]


def geometric_frontier(geometry, density=128 / TWO_PI, active_fields=None):
    """Return unknown transverse area by strip to guide surface insertion.

    This is a geometric refinement indicator, never a population-weight or
    connectivity error bound. Open-fold end strips are retained explicitly.
    """
    result = []
    activity = []
    for layer in geometry.layers:
        if active_fields is None:
            activity.append(np.ones(len(layer.H), bool))
        else:
            b = np.asarray(active_fields, float)
            m = layer.chain.min_H[:, 0]
            activity.append(
                np.any(
                    (layer.H[:, None] > b)
                    & ((m[:, None] < b) | ~np.isfinite(m[:, None])),
                    axis=1,
                )
            )
    for i, (a, c) in enumerate(zip(geometry.layers, geometry.layers[1:])):
        mesh = stitch_ridges(
            a.arcs,
            c.arcs,
            geometry.s[i],
            geometry.s[i + 1],
            a.surface.iota,
            c.surface.iota,
            a.surface.nfp,
            density,
        )
        ds = geometry.s[i + 1] - geometry.s[i]
        unknown = (
            ds
            / 2
            * (
                a.measure[mesh.frontier_bottom & activity[i]].sum()
                + c.measure[mesh.frontier_top & activity[i + 1]].sum()
            )
        )
        partner_unknown = unknown
        nb = len(a.points)
        fold_unknown = sum(
            area
            for node, area in mesh.fold_gaps
            if (activity[i][node] if node < nb else activity[i + 1][node - nb])
        )
        unknown += fold_unknown
        result.append(
            dict(
                s_bottom=float(geometry.s[i]),
                s_top=float(geometry.s[i + 1]),
                frontier_measure=float(unknown),
                partner_frontier_measure=float(partner_unknown),
                fold_frontier_measure=float(fold_unknown),
                total_measure=float(
                    ds
                    / 2
                    * (a.measure[activity[i]].sum() + c.measure[activity[i + 1]].sum())
                ),
                census_change=len(a.arcs) != len(c.arcs),
            )
        )
    return result


def polish_folds(geometry):
    """Locate generic ridge folds on the exact series, DESIGN §4.1.3.

    Solve D_parallel B=D_parallel²B=0 in angular coordinates. A nonzero third
    derivative identifies a generic max/min annihilation. Degenerate or
    failed roots retain unresolved end strips, never fabricated coverage.
    """
    from scipy.optimize import root
    from .ledger import node_measure

    records = []
    for layer in geometry.layers:
        sf = layer.surface
        for j, arc in enumerate(layer.arcs):
            if arc.closed or arc.endpoints is None:
                continue
            if arc.folds_resolved is not None and np.all(arc.folds_resolved):
                continue
            arc.folds_resolved = np.zeros(2, bool)
            for end, guess in enumerate(arc.endpoints.copy()):
                solution = root(
                    lambda p: np.array([sf.exact(p[0], p[1], k) for k in (1, 2)]),
                    guess,
                    tol=1e-11,
                )
                point = solution.x
                residual = np.array([sf.exact(point[0], point[1], k) for k in (1, 2)])
                scale = max(sf.Bmax, 1.0)
                resolved = (
                    np.max(np.abs(residual)) < 1e-9 * scale
                    and np.linalg.norm(point - guess) < 0.05
                    and abs(float(sf.exact(point[0], point[1], 3))) > 1e-9 * scale
                )
                if resolved:
                    arc.endpoints[end] = point
                    arc.folds_resolved[end] = True
                records.append(
                    dict(
                        s=sf.s,
                        arc=j,
                        end=end,
                        resolved=bool(resolved),
                        residual=residual.tolist(),
                        point=point.tolist(),
                    )
                )
            alpha = arc.endpoints[:, 0] - sf.iota * arc.endpoints[:, 1] / sf.nfp
            faces = np.r_[alpha[0], (arc.alpha[:-1] + arc.alpha[1:]) / 2, alpha[1]]
            arc.measure = np.abs(np.diff(faces))
        layer.measure = node_measure(layer.arcs)
    return records
