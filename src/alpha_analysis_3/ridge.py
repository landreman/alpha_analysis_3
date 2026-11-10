"""Concave ridge arcs, winding and transverse measure (DESIGN §§3.1,3.5,4.1).

Angles are (theta,phi=N_fp*zeta), in radians, unwrapped along each arc.
Periodic contour stitching follows the prototype; uniform arc length, fold
endpoints, contractible loops and seam orientation are handled explicitly.
"""

from dataclasses import dataclass

import contourpy
import numpy as np

from .tables import TWO_PI


@dataclass
class Ridge:
    points: np.ndarray
    closed: bool
    wind: tuple[int, int]
    H: np.ndarray | None = None
    alpha: np.ndarray | None = None
    dalpha: np.ndarray | None = None
    measure: np.ndarray | None = None


def _same(a, b, tol=1e-7):
    return np.linalg.norm((a - b + np.pi) % TWO_PI - np.pi) < tol


def _stitch(lines):
    pieces = [a.copy() for a in lines if len(a) >= 2]
    curves = []
    while pieces:
        curve = pieces.pop()
        while not _same(curve[0], curve[-1]):
            found = False
            for j, piece in enumerate(pieces):
                if _same(curve[-1], piece[-1]):
                    piece = piece[::-1]
                if _same(curve[-1], piece[0]):
                    shift = TWO_PI * np.round((curve[-1] - piece[0]) / TWO_PI)
                    curve = np.vstack((curve, (piece + shift)[1:]))
                    pieces.pop(j)
                    found = True
                    break
            if not found:
                break
        curves.append(curve)
    return curves


def ridge_arcs(surface):
    """Extract D_parallel B=0, D_parallel²B<0; retain folds as endpoints."""
    grid = np.pad(surface.DB, ((0, 1), (0, 1)), mode="wrap")
    # Shift a contour exactly on the box edge by half a grid cell. Otherwise
    # contourpy can omit a maximum that lies identically on that boundary.
    offset = np.pi / surface.nt
    theta = np.arange(surface.nt + 1) * TWO_PI / surface.nt + offset
    phi = np.linspace(0, TWO_PI, surface.nphi + 1) + np.pi / surface.nphi
    grid = surface.db(theta[:, None], phi[None, :])
    generator = contourpy.contour_generator(x=phi, y=theta, z=grid)
    curves = _stitch([line[:, ::-1] for line in generator.lines(0)])
    arcs = []
    for curve in curves:
        closed = _same(curve[0], curve[-1])
        wind = (
            tuple(np.rint((curve[-1] - curve[0]) / TWO_PI).astype(int))
            if closed
            else (0, 0)
        )
        d2 = surface.d2b(curve[:, 0], curve[:, 1])
        if np.all(d2 < 0) and closed:
            if wind < (0, 0):
                curve = curve[::-1]
                wind = tuple(-v for v in wind)
            arcs.append(Ridge(curve, True, wind))
            continue
        if not np.any(d2 < 0):
            continue
        if closed:
            # Start at a minimum branch, keeping the unwrapped winding shift.
            j = int(np.argmax(d2))
            curve = np.vstack((curve[j:-1], curve[: j + 1] + TWO_PI * np.array(wind)))
            d2 = surface.d2b(curve[:, 0], curve[:, 1])
        current = []
        for j in range(len(curve) - 1):
            a, b = curve[j : j + 2]
            da, db = d2[j : j + 2]
            if da * db < 0:
                fold = a + (b - a) * (-da) / (db - da)
                if da < 0:
                    current.extend((a, fold))
                    if len(current) > 2:
                        arcs.append(Ridge(np.array(current), False, (0, 0)))
                    current = []
                else:
                    current = [fold]
            elif da < 0:
                current.append(a)
        if d2[-1] < 0:
            current.append(curve[-1])
        if len(current) > 2:
            arcs.append(Ridge(np.array(current), False, (0, 0)))
    return arcs


def place_nodes(surface, arcs, density=128 / TWO_PI):
    """Uniform arc-length nodes, heights and |d alpha|, with closing seam."""
    result = []
    for arc in arcs:
        points = arc.points
        distances = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
        length = distances[-1]
        if length < 1e-10:
            continue
        count = max(8, int(np.ceil(density * length)))
        u = (np.arange(count) + (0 if arc.closed else 0.5)) * length / count
        nodes = np.column_stack(
            [np.interp(u, distances, points[:, j]) for j in range(2)]
        )
        # Project along the field line onto the maximum. This keeps alpha fixed
        # and avoids contour-grid offsets in the starting chain extremum.
        for _ in range(3):
            derivative = surface.db(nodes[:, 0], nodes[:, 1])
            curvature = surface.d2b(nodes[:, 0], nodes[:, 1])
            step = np.divide(
                derivative,
                curvature,
                out=np.zeros_like(derivative),
                where=np.abs(curvature) > 1e-12,
            )
            step = np.clip(step, -0.02, 0.02)
            nodes -= step[:, None] * np.array([surface.iota, surface.nfp])
        alpha = nodes[:, 0] - surface.iota * nodes[:, 1] / surface.nfp
        delta = np.diff(alpha)
        if arc.closed:
            p, q = arc.wind
            seam = TWO_PI * p - surface.iota * surface.L * q
            delta = np.r_[delta, alpha[0] + seam - alpha[-1]]
        if arc.closed:
            measure = (np.abs(delta) + np.roll(np.abs(delta), 1)) / 2
        else:
            endpoints = (
                points[[0, -1], 0] - surface.iota * points[[0, -1], 1] / surface.nfp
            )
            faces = np.r_[endpoints[0], (alpha[:-1] + alpha[1:]) / 2, endpoints[1]]
            measure = np.abs(np.diff(faces))
        result.append(
            Ridge(
                nodes,
                arc.closed,
                arc.wind,
                surface.b(nodes[:, 0], nodes[:, 1]),
                alpha,
                np.abs(delta),
                measure,
            )
        )
    return result


def radial_partners(bottom, top, tolerance=0.4):
    """One-direction nearest torus nodes, integer IDs or -1 (M1 scope)."""
    a = np.vstack([r.points for r in bottom]) if bottom else np.empty((0, 2))
    b = np.vstack([r.points for r in top]) if top else np.empty((0, 2))
    result = np.full(len(a), -1, dtype=int)
    if len(b):
        for start in range(0, len(a), 128):
            d = (a[start : start + 128, None] - b[None] + np.pi) % TWO_PI - np.pi
            distance = np.linalg.norm(d, axis=2)
            ids = np.argmin(distance, axis=1)
            good = distance[np.arange(len(ids)), ids] < tolerance
            result[start : start + len(ids)] = np.where(good, ids, -1)
    return result
