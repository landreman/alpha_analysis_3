"""Pitch activity, exits, bounce points, J and K (DESIGN §§2.2,3.3,4.2).

J and K are lengths; signed zeta offsets follow sign(C). Failed or capped
states have NaN action/time, never zero. Endpoint sine quadrature follows the
prototype, with composition at every extremum to resolve long, rippled wells.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

INACTIVE, ACTIVE, CAPPED, FAILED = 0, 1, 2, 3


@dataclass
class Wells:
    status: np.ndarray
    exit_index: np.ndarray
    z_in: np.ndarray
    z_out: np.ndarray
    J: np.ndarray
    K: np.ndarray


@lru_cache(maxsize=8)
def _gauss(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return x * np.pi / 2, w * np.pi / 2


def wells_at_pitch(surface, theta, phi, heights, chain, b, order=64):
    """First above-b return; capped/failed wells retain NaN J and K.

    Work is vectorized over nodes and bounded batches of extremum segments.
    """
    if not np.isfinite(b) or b <= 0:
        raise ValueError("bounce field must be positive and finite")
    theta, phi, heights = np.broadcast_arrays(theta, phi, heights)
    count = len(theta)
    status = np.full(count, INACTIVE, dtype=np.int8)
    exits = np.full(count, -1, dtype=int)
    zi, zo, J, K = (np.full(count, np.nan) for _ in range(4))
    candidate = heights > b
    missing = candidate & ~np.isfinite(chain.min_H[:, 0])
    status[missing] = CAPPED
    active = candidate & (chain.min_H[:, 0] < b)
    above = chain.max_H > b
    has_exit = above.any(axis=1)
    exits[active & has_exit] = np.argmax(above, axis=1)[active & has_exit]
    status[active & ~has_exit] = CAPPED
    indices = np.flatnonzero(active & has_exit)
    if not len(indices):
        return Wells(status, exits, zi, zo, J, K)
    ne = exits[indices]
    exit_z = chain.max_z[indices, ne]
    preceding = np.where(
        np.isfinite(chain.min_z[indices])
        & (surface.sign * chain.min_z[indices] < surface.sign * exit_z[:, None]),
        surface.sign * chain.min_z[indices],
        -np.inf,
    )
    last_min = surface.sign * np.max(preceding, axis=1)
    good = np.isfinite(last_min)
    status[indices[~good]] = FAILED
    indices, exit_z, last_min = indices[good], exit_z[good], last_min[good]
    th, ph = theta[indices], phi[indices]

    def bisect(lo, hi, falling):
        for _ in range(48):
            mid = (lo + hi) / 2
            high = surface.b(th + surface.iota * mid, ph + surface.nfp * mid) > b
            move = high if falling else ~high
            lo, hi = np.where(move, mid, lo), np.where(move, hi, mid)
        return (lo + hi) / 2

    zi[indices] = bisect(np.zeros(len(indices)), chain.min_z[indices, 0].copy(), True)
    zo[indices] = bisect(last_min.copy(), exit_z.copy(), False)
    owners, lefts, rights = [], [], []
    for j in indices:
        extrema = np.r_[chain.max_z[j], chain.min_z[j]]
        inner = extrema[
            np.isfinite(extrema)
            & (surface.sign * extrema > surface.sign * zi[j])
            & (surface.sign * extrema < surface.sign * zo[j])
        ]
        inner = inner[np.argsort(surface.sign * inner)]
        breaks = np.r_[zi[j], inner, zo[j]]
        owners.extend([j] * (len(breaks) - 1))
        lefts.extend(breaks[:-1])
        rights.extend(breaks[1:])
    owners = np.asarray(owners, dtype=int)
    lefts, rights = np.asarray(lefts), np.asarray(rights)
    sums_j, sums_k = np.zeros(count), np.zeros(count)
    failed = np.zeros(count, dtype=bool)
    t, w = _gauss(order)
    for start in range(0, len(owners), 128):
        ids = owners[start : start + 128]
        l, r = lefts[start : start + 128, None], rights[start : start + 128, None]
        half = (r - l) / 2
        z = (r + l) / 2 + half * np.sin(t)
        B = surface.b(
            theta[ids, None] + surface.iota * z, phi[ids, None] + surface.nfp * z
        )
        g = 1 - B / b
        bad = np.any((g <= 0) | ~np.isfinite(g), axis=1)
        failed[ids[bad]] = True
        safe = np.where(g > 0, g, np.nan)
        dl = abs(surface.C) / B * np.abs(half) * np.cos(t)
        np.add.at(sums_j, ids, np.sum(w * dl * np.sqrt(safe), axis=1))
        np.add.at(sums_k, ids, np.sum(w * dl / np.sqrt(safe), axis=1))
    status[indices] = np.where(failed[indices], FAILED, ACTIVE)
    valid = indices[~failed[indices]]
    J[valid], K[valid] = sums_j[valid], sums_k[valid]
    return Wells(status, exits, zi, zo, J, K)
