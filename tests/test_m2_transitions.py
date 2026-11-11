"""Independent split/merge geometry and legacy-event expectations (DESIGN M2)."""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.tables import Surface
from alpha_analysis_3.transitions import (
    sample_states,
    interior_barrier,
    locate_transition,
    well_assignment,
    surface_transitions,
)
from alpha_analysis_3.ridge import ridge_arcs, place_nodes


class SplitField:
    """B=2+q-q*cos(phi)+.8*cos(2phi), q=-.1+.2*cos(theta)."""

    nfp = 1
    xm = np.array([0, 1, 0, 1, 1, 0])
    xn = np.array([0, 0, 1, 1, -1, 2])

    def bmnc(self, s):
        return np.array([1.9, 0.2, 0.1, -0.1, -0.1, 0.8])

    def bmns(self, s):
        return np.zeros(6)

    def iota(self, s):
        return 0.0

    def G(self, s):
        return getattr(self, "current", 3.0)

    def I(self, s):
        return 0.0


@pytest.mark.parametrize("sign", [-1, 1])
def test_synthetic_split_merge(sign):
    field = SplitField()
    field.current = sign * 3.0
    sf = Surface(field, 0.5, 64, 64)
    left, right = sample_states(sf, [[1.3, 0], [1.8, 0]], 2.6)
    assert interior_barrier(left, right) is not None
    event = locate_transition(sf, [1.3, 0], [1.8, 0], 2.6, tolerance=2e-9)
    assert event.resolved, event.reason
    assert event.point[0] == pytest.approx(np.pi / 2, abs=2e-9)
    assert event.width <= 2e-9
    assert event.parent == pytest.approx(event.child_1 + event.child_3, rel=2e-7)

    # Direct integral of the prescribed field, independent of tables/chains.
    def B(z):
        return 1.9 + 0.1 * np.cos(z) + 0.8 * np.cos(2 * z)

    incoming = brentq(lambda z: B(z) - 2.6, 0, 1)
    child = quad(
        lambda z: 3 / B(z) * np.sqrt(max(0, 1 - B(z) / 2.6)),
        incoming,
        np.pi,
        epsabs=1e-11,
        epsrel=1e-11,
    )[0]
    assert event.child_1 == pytest.approx(child, rel=2e-7)
    assert event.child_3 == pytest.approx(child, rel=2e-7)
    census = surface_transitions(
        sf, place_nodes(sf, ridge_arcs(sf), 16 / (2 * np.pi)), 2.6, max_factor=2
    )
    assert census.transitions == census.gamma_max == 2
    assert len(census.events) == 2


def test_interior_barrier_rejects_smooth_bounce_motion():
    sf = Surface(SplitField(), 0.5, 64, 64)
    states = sample_states(sf, [[1.0, 0], [1.4, 0]], 2.6)
    assert abs(states[0].z_out - states[1].z_out) > 0.01
    assert interior_barrier(*states) is None


def test_well_assignment_catches_inconsistent_edges():
    from dataclasses import replace

    sf = Surface(SplitField(), 0.5, 64, 64)
    a, c = sample_states(sf, [[0.8, 0], [0.9, 0]], 2.6)
    assert well_assignment(a, c, 0.1) == "regular"
    c = replace(c, J=19.0)
    a = replace(a, J=1e-19)
    assert well_assignment(a, c, 0.1) == "unknown"


def test_dmercfail_generic_event():
    root = Path(__file__).resolve().parents[1]
    record = json.loads((root / "data/legacy_probes.json").read_text())["generic_event"]
    entry = next(
        e
        for e in json.loads((root / "data/manifest.json").read_text())["files"]
        if e["id"] == "DMercFail"
    )
    sf = Surface(
        BoozerField(root / "data" / entry["file"]), record["parameter_s_alpha"][0]
    )
    # The start maximum is the one upstream of the old positive-zeta well.
    alpha = record["parameter_s_alpha"][1]

    def start(a):
        z = brentq(
            lambda z: float(sf.exact(a + sf.iota * z, sf.nfp * z, 1)), -1.6, -1.4
        )
        return [a + sf.iota * z, sf.nfp * z]

    event = locate_transition(
        sf, start(alpha - 0.005), start(alpha + 0.005), record["b"]
    )
    assert event.resolved, event.reason
    found = event.point[0] - sf.iota * event.point[1] / sf.nfp
    partner_alpha = event.partner[0] - sf.iota * event.partner[1] / sf.nfp
    assert partner_alpha == pytest.approx(found, abs=1e-12)
    assert abs(found - alpha) < 0.005
    assert event.parent == pytest.approx(record["action_estimates"]["parent"], rel=2e-5)
    assert event.parent == pytest.approx(event.child_1 + event.child_3, rel=2e-6)


def test_transition_count_matches_gamma_max():
    root = Path(__file__).resolve().parents[1]
    entry = json.loads((root / "data/manifest.json").read_text())["files"][2]
    b = json.loads((root / "data/legacy_probes.json").read_text())["generic_event"]["b"]
    sf = Surface(BoozerField(root / "data" / entry["file"]), 0.8)
    arcs = place_nodes(sf, ridge_arcs(sf), 64 / (2 * np.pi))
    from types import SimpleNamespace
    from alpha_analysis_3.cells import polish_folds

    polish_folds(SimpleNamespace(layers=[SimpleNamespace(surface=sf, arcs=arcs)]))
    result = surface_transitions(sf, arcs, b, max_factor=4)
    assert result.transitions == result.gamma_max == 2
    assert len(result.events) == 2 and not result.unknown_edges
