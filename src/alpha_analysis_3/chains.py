"""Bounded along-line extremum chains (DESIGN §§3.2,4.1).

Offsets max_z/min_z are signed zeta radians from the ridge node, increasing
in physical distance along sign(C). Heights have the input B units. The prototype's bounded scans are adapted to bracket roots of the derivative,
which retains shallow first minima missed by field-value differences.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Chains:
    max_z: np.ndarray
    max_H: np.ndarray
    min_z: np.ndarray
    min_H: np.ndarray
    limit: float
    samples_per_period: int


def trace_chains(surface, theta, phi, periods=48, samples_per_period=128):
    """Trace up to K_max periods, with bounded temporary arrays (32 lines)."""
    if (
        not 1 <= periods <= 256
        or not 16 <= samples_per_period <= 4096
        or periods * samples_per_period > 65536
    ):
        raise ValueError("chain controls exceed the bounded sampling range")
    theta, phi = np.broadcast_arrays(theta, phi)
    if theta.ndim != 1:
        raise ValueError("chain nodes must be one-dimensional")
    dz = surface.L / samples_per_period
    distance = np.arange(periods * samples_per_period + 1) * dz
    maxima, minima = [], []
    for start in range(0, len(theta), 32):
        th, ph = theta[start : start + 32], phi[start : start + 32]
        z = surface.sign * distance
        derivative = surface.sign * surface.db(
            th[:, None] + surface.iota * z, ph[:, None] + surface.nfp * z
        )
        # A ridge is falling immediately forward, including a very shallow
        # max/min pair wholly inside the first sampling interval.
        derivative[:, 0] = -np.finfo(float).tiny
        for is_max, output in ((True, maxima), (False, minima)):
            mask = (
                (derivative[:, :-1] > 0) & (derivative[:, 1:] <= 0)
                if is_max
                else (derivative[:, :-1] < 0) & (derivative[:, 1:] >= 0)
            )
            rows, columns = np.nonzero(mask)
            lo, hi = distance[columns].copy(), distance[columns + 1].copy()
            for _ in range(34):
                mid = (lo + hi) / 2
                value = surface.sign * surface.db(
                    th[rows] + surface.iota * surface.sign * mid,
                    ph[rows] + surface.nfp * surface.sign * mid,
                )
                move = value > 0 if is_max else value < 0
                lo, hi = np.where(move, mid, lo), np.where(move, hi, mid)
            offsets = surface.sign * (lo + hi) / 2
            heights = surface.b(
                th[rows] + surface.iota * offsets, ph[rows] + surface.nfp * offsets
            )
            for row in range(len(th)):
                take = rows == row
                output.append((offsets[take], heights[take]))

    def pack(rows):
        width = max((len(r[0]) for r in rows), default=0)
        offsets = np.full((len(theta), max(1, width)), np.nan)
        heights = offsets.copy()
        for j, (x, h) in enumerate(rows):
            offsets[j, : len(x)], heights[j, : len(x)] = x, h
        return offsets, heights

    mz, mh = pack(maxima)
    nz, nh = pack(minima)
    return Chains(
        mz, mh, nz, nh, float(surface.sign * distance[-1]), samples_per_period
    )
