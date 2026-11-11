"""Pitch activity, exits, bounce points, J and K (DESIGN §§2.2,3.3,4.2).

J and K are lengths; signed zeta offsets follow sign(C). Failed or capped
states have NaN action/time, never zero. Endpoint sine quadrature follows the
prototype, with composition at every extremum to resolve long, rippled wells.
"""

from dataclasses import dataclass
from functools import lru_cache
import time

from .chains import trace_chains, chain_bytes, replace_chain_rows

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
    partial_K: np.ndarray


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
    # A cap is an unknown full K, with a separately accumulated finite lower
    # bound. A bound of zero when no below-b point was observed is not zero K.
    partial = np.zeros(count)
    limits = surface.sign * surface.L * chain.periods
    if len(indices):
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
        active_ids = indices[good]
        active_exit, active_min = exit_z[good], last_min[good]
    else:
        active_ids = np.empty(0, dtype=int)
        active_exit = active_min = np.empty(0)

    def bisect(ids, lo, hi, falling):
        th, ph = theta[ids], phi[ids]
        for _ in range(48):
            mid = (lo + hi) / 2
            high = surface.b(th + surface.iota * mid, ph + surface.nfp * mid) > b
            move = high if falling else ~high
            lo, hi = np.where(move, mid, lo), np.where(move, hi, mid)
        return (lo + hi) / 2

    zi[active_ids] = bisect(
        active_ids, np.zeros(len(active_ids)), chain.min_z[active_ids, 0].copy(), True
    )
    zo[active_ids] = bisect(active_ids, active_min.copy(), active_exit.copy(), False)
    ends = zo.copy()
    cap_ids = np.flatnonzero(status == CAPPED)
    if len(cap_ids):
        end_B = surface.b(
            theta[cap_ids] + surface.iota * limits[cap_ids],
            phi[cap_ids] + surface.nfp * limits[cap_ids],
        )
        entered = (chain.min_H[cap_ids, 0] < b) | (end_B < b)
        cap_ids = cap_ids[entered]
        end_B = end_B[entered]
        incoming = np.where(
            chain.min_H[cap_ids, 0] < b, chain.min_z[cap_ids, 0], limits[cap_ids]
        )
        zi[cap_ids] = bisect(cap_ids, np.zeros(len(cap_ids)), incoming, True)
        ends[cap_ids] = limits[cap_ids]
        # An outgoing bounce can precede a not-yet-sampled exit maximum.
        # It bounds the observed integral but does not change the capped status.
        returned = end_B >= b
        ids = cap_ids[returned]
        if len(ids):
            last = surface.sign * np.nanmax(surface.sign * chain.min_z[ids], axis=1)
            ends[ids] = bisect(ids, last, limits[ids].copy(), False)
    indices = np.r_[active_ids, cap_ids]
    owners, lefts, rights = [], [], []
    for j in indices:
        extrema = np.r_[chain.max_z[j], chain.min_z[j]]
        inner = extrema[
            np.isfinite(extrema)
            & (surface.sign * extrema > surface.sign * zi[j])
            & (surface.sign * extrema < surface.sign * ends[j])
        ]
        inner = inner[np.argsort(surface.sign * inner)]
        breaks = np.r_[zi[j], inner, ends[j]]
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
    status[indices[failed[indices]]] = FAILED
    status[active_ids[~failed[active_ids]]] = ACTIVE
    valid = active_ids[~failed[active_ids]]
    J[valid], K[valid] = sums_j[valid], sums_k[valid]
    partial[indices] = sums_k[indices]
    partial[status == FAILED] = np.nan
    return Wells(status, exits, zi, zo, J, K, partial)


def adaptive_wells_at_pitch(
    surface,
    theta,
    phi,
    heights,
    chain,
    b,
    order=64,
    max_periods=256,
    allocation_check=None,
):
    """ADR 0001: double only capped chains up to a finite ceiling.

    Returns (wells, updated chains, extension diagnostics). Original arrays are
    preserved. Nodes still capped retain NaN full J/K and partial_K as a lower
    bound. The ceiling is a cost guard, never a passing classification.
    """
    if not 1 <= max_periods <= 256:
        raise ValueError("adaptive ceiling must be between 1 and 256 periods")
    theta, phi, heights = np.broadcast_arrays(theta, phi, heights)
    well = wells_at_pitch(surface, theta, phi, heights, chain, b, order)
    first = int(np.sum(well.status == CAPPED))
    extended = set()
    start = time.monotonic()
    had_derivatives = surface.DB is not None
    try:
        while True:
            ids = np.flatnonzero(
                (well.status == CAPPED) & (chain.periods < max_periods)
            )
            if not len(ids):
                break
            surface.restore_derivatives()
            targets = np.minimum(2 * chain.periods[ids], max_periods)
            for target in np.unique(targets):
                selected = ids[targets == target]
                traced = trace_chains(
                    surface,
                    theta[selected],
                    phi[selected],
                    int(target),
                    chain.samples_per_period,
                    allocation_check=(
                        (lambda size: allocation_check(chain, size))
                        if allocation_check is not None
                        else None
                    ),
                )
                chain = replace_chain_rows(
                    chain,
                    selected,
                    traced,
                    allocation_check=(
                        (
                            lambda size: allocation_check(
                                chain, size + chain_bytes(traced)
                            )
                        )
                        if allocation_check is not None
                        else None
                    ),
                )
                new = wells_at_pitch(
                    surface,
                    theta[selected],
                    phi[selected],
                    heights[selected],
                    traced,
                    b,
                    order,
                )
                for name in (
                    "status",
                    "exit_index",
                    "z_in",
                    "z_out",
                    "J",
                    "K",
                    "partial_K",
                ):
                    getattr(well, name)[selected] = getattr(new, name)
                extended.update(selected.tolist())
                del traced, new
    finally:
        if not had_derivatives:
            surface.release_derivatives()
    stats = dict(
        initial_capped_nodes=first,
        extended_nodes=len(extended),
        maximum_periods=int(np.max(chain.periods, initial=0)),
        extension_seconds=time.monotonic() - start,
    )
    return well, chain, stats
