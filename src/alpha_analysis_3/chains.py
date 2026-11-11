"""Bounded along-line extremum chains (DESIGN §§3.2,4.1).

Offsets max_z/min_z are signed zeta radians from the ridge node, increasing
in physical distance along sign(C). Heights have the input B units. The prototype's bounded scans are adapted to bracket roots of the derivative,
which retains shallow first minima missed by field-value differences.
"""

from dataclasses import dataclass
import weakref

import numpy as np

CHAIN_ARRAYS = ("max_z", "max_H", "min_z", "min_H", "periods")


def chain_bytes(chain):
    """Owned chain-array bytes, including NaN padding and per-row limits."""
    return sum(getattr(chain, name).nbytes for name in CHAIN_ARRAYS)


class ChainBudget:
    """Bound chain storage across all layers, including replacement peaks.

    Weak references track the retained geometry without extending its lifetime.
    The 768 MiB default leaves room for field tables and sampling temporaries;
    the driver's independent 2 GB RSS watchdog remains mandatory.
    """

    def __init__(self, limit_bytes=768 * 1024**2):
        if limit_bytes <= 0:
            raise ValueError("chain budget must be positive")
        self.limit_bytes = int(limit_bytes)
        self._retained = {}

    def track(self, owner, chain):
        self._retained[id(owner)] = weakref.ref(chain)

    def _live(self):
        return {
            id(ch): ch for ref in self._retained.values() if (ch := ref()) is not None
        }

    @property
    def retained_bytes(self):
        return sum(chain_bytes(ch) for ch in self._live().values())

    def check_new(self, allocation_bytes):
        if self.retained_bytes + allocation_bytes > self.limit_bytes:
            raise MemoryError("global chain allocation budget exceeded")

    def check_extension(self, current, allocation_bytes):
        extra = 0 if id(current) in self._live() else chain_bytes(current)
        self.check_new(extra + allocation_bytes)


@dataclass
class Chains:
    max_z: np.ndarray
    max_H: np.ndarray
    min_z: np.ndarray
    min_H: np.ndarray
    limit: float
    samples_per_period: int
    periods: np.ndarray


def trace_chains(
    surface, theta, phi, periods=64, samples_per_period=128, allocation_check=None
):
    """Trace up to K_max periods, with bounded temporary arrays (32 lines)."""
    if (
        not 1 <= periods <= 256
        or not 16 <= samples_per_period <= 4096
        or periods * samples_per_period > 1_048_576
    ):
        raise ValueError("chain controls exceed the bounded sampling range")
    theta, phi = np.broadcast_arrays(theta, phi)
    if theta.ndim != 1:
        raise ValueError("chain nodes must be one-dimensional")
    dz = surface.L / samples_per_period
    distance = np.arange(periods * samples_per_period + 1) * dz
    maxima, minima = [], []
    batch = max(1, min(32, 2_097_152 // len(distance)))
    for start in range(0, len(theta), batch):
        th, ph = theta[start : start + batch], phi[start : start + batch]
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

    widths = [
        max(1, max((len(r[0]) for r in rows), default=0)) for rows in (maxima, minima)
    ]
    if allocation_check is not None:
        allocation_check(
            2 * len(theta) * sum(widths) * np.dtype(float).itemsize
            + len(theta) * np.dtype(int).itemsize
        )

    def pack(rows):
        width = max((len(r[0]) for r in rows), default=0)
        if len(theta) * max(1, width) > 4_000_000:
            raise MemoryError("extremum table exceeds the bounded chain allocation")
        offsets = np.full((len(theta), max(1, width)), np.nan)
        heights = offsets.copy()
        for j, (x, h) in enumerate(rows):
            offsets[j, : len(x)], heights[j, : len(x)] = x, h
        return offsets, heights

    mz, mh = pack(maxima)
    nz, nh = pack(minima)
    return Chains(
        mz,
        mh,
        nz,
        nh,
        float(surface.sign * distance[-1]),
        samples_per_period,
        np.full(len(theta), periods, dtype=int),
    )


def chain_rows(chain, ids):
    """Copy selected integer node rows, retaining each row's traced limit."""
    return Chains(
        *(
            getattr(chain, name)[ids].copy()
            for name in ("max_z", "max_H", "min_z", "min_H")
        ),
        chain.limit,
        chain.samples_per_period,
        chain.periods[ids].copy(),
    )


def replace_chain_rows(chain, ids, extended, allocation_check=None):
    """Replace only retraced rows; pad other rows with NaN, never fake extrema."""
    allocation_bytes = chain.periods.nbytes
    for name in CHAIN_ARRAYS[:-1]:
        old, new = getattr(chain, name), getattr(extended, name)
        elements = len(old) * max(old.shape[1], new.shape[1])
        if elements > 4_000_000:
            raise MemoryError("extended table exceeds the bounded chain allocation")
        allocation_bytes += elements * old.dtype.itemsize
    if allocation_check is not None:
        allocation_check(allocation_bytes)
    arrays = []
    for name in ("max_z", "max_H", "min_z", "min_H"):
        old, new = getattr(chain, name), getattr(extended, name)
        width = max(old.shape[1], new.shape[1])
        result = np.full((len(old), width), np.nan)
        result[:, : old.shape[1]] = old
        result[ids] = np.nan
        result[ids, : new.shape[1]] = new
        arrays.append(result)
    periods = chain.periods.copy()
    periods[ids] = extended.periods
    limit = chain.limit if abs(chain.limit) > abs(extended.limit) else extended.limit
    return Chains(*arrays, limit, chain.samples_per_period, periods)
