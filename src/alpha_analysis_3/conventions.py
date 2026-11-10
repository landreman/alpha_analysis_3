"""Project-wide numerical conventions.

Implements the accepted decisions of docs/DESIGN.md section 1 and the conventions of
section 2. Every quantity here is dimensionless or in the units of the boozmn file
(B in tesla).
"""

from __future__ import annotations

import numpy as np

# Accuracy target: an error magnitude, not a full interval width (DESIGN section 1).
RELATIVE_TOLERANCE = 0.02
ABSOLUTE_TOLERANCE = 0.01

# Runtime and memory targets per equilibrium (DESIGN section 7).
RUNTIME_TARGET_SECONDS = 300.0
MEMORY_TARGET_GB = 2.0

# Diagnostic pitch levels (DESIGN section 1). lambda_n labels levels only; the conserved
# variable along a contour is b itself.
LAMBDA_LEVELS = (0.05, 0.1, 0.5, 0.8, 0.9, 0.95)


def numerical_error_target(f):
    """Allowed numerical error max(0.02 f, 0.01) for a loss fraction f.

    Accepts scalars or arrays. The floor applies at f = 0: returning zero without a
    calculation never satisfies the target.
    """
    return np.maximum(RELATIVE_TOLERANCE * np.abs(f), ABSOLUTE_TOLERANCE)


def pitch_from_lambda(lambda_n, b_min, b_max):
    """Bounce field b = B_min + lambda_n (B_max - B_min) with radially global extrema."""
    lambda_n = np.asarray(lambda_n, dtype=float)
    return b_min + lambda_n * (b_max - b_min)


def lambda_from_pitch(b, b_min, b_max):
    """Inverse of :func:`pitch_from_lambda`."""
    return (np.asarray(b, dtype=float) - b_min) / (b_max - b_min)


def default_source(rho):
    """Centrally peaked alpha source h(rho) = (1 - rho^2)^4, rho = sqrt(s)."""
    rho = np.asarray(rho, dtype=float)
    return (1.0 - rho**2) ** 4


def uniform_source(rho):
    """The comparison source h = 1, kept as a second column in every report."""
    return np.ones_like(np.asarray(rho, dtype=float))


def source_at_s(h, s):
    """Evaluate a source profile h(rho) at flux label s = rho^2 (DESIGN section 2.4)."""
    return h(np.sqrt(np.asarray(s, dtype=float)))


def field_line_direction(G, I, iota):
    """Physical direction of motion along the field line, sign(G + iota I).

    Raises if C = G + iota I vanishes anywhere, because the along-line arc length
    dl = |C| / B |dzeta| then degenerates (DESIGN section 2.1).
    """
    C = np.asarray(G, dtype=float) + np.asarray(iota, dtype=float) * np.asarray(
        I, dtype=float
    )
    if np.any(C == 0.0):
        raise ValueError(
            "G + iota I vanishes; the field-line parametrization degenerates"
        )
    return np.sign(C)
