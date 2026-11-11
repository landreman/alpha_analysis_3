"""Independent M1 expectations: Fourier sums, topology and elliptic integrals."""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ellipk, elliprf, elliprj

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.tables import Surface
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.chains import trace_chains
from alpha_analysis_3.wells import wells_at_pitch, ACTIVE, CAPPED, INACTIVE
from alpha_analysis_3.wells import adaptive_wells_at_pitch
from alpha_analysis_3.ledger import (
    Geometry,
    pitch_ledger,
    volume_surface,
    measure_pitch,
)
from alpha_analysis_3.conventions import uniform_source, default_source

ROOT = Path(__file__).resolve().parents[1]


class FourierModel:
    """Small synthetic Fourier field, with no production evaluator in its oracle."""

    nfp = 1

    def __init__(self, m=(0, 1), n=(0, 0), c=(2.0, 0.5), sine=None, C=3.0, iota=0.7):
        self.xm, self.xn = np.array(m), np.array(n)
        self.c = np.array(c)
        self.sine = np.zeros_like(self.c) if sine is None else np.array(sine)
        self.current, self.transform = C, iota

    def bmnc(self, s):
        return self.c

    def bmns(self, s):
        return self.sine

    def G(self, s):
        return self.current

    def I(self, s):
        return 0.0

    def iota(self, s):
        return self.transform


@pytest.mark.parametrize(
    "entry",
    json.loads((ROOT / "data/manifest.json").read_text())["files"],
    ids=lambda e: e["id"],
)
def test_tables_match_exact_series(entry):
    field = BoozerField(ROOT / "data" / entry["file"])
    sf = Surface(field, 0.63, 512, 512)
    rng = np.random.default_rng(53)
    ij = rng.integers(0, 512, (2, 80))
    theta, phi = 2 * np.pi * ij / 512
    phase = theta[:, None] * field.xm - phi[:, None] * field.xn / field.nfp
    c, d = field.bmnc(sf.s), field.bmns(sf.s)
    k = sf.iota * field.xm - field.xn
    exact = [
        np.cos(phase) @ c + np.sin(phase) @ d,
        -np.sin(phase) @ (k * c) + np.cos(phase) @ (k * d),
        -np.cos(phase) @ (k * k * c) - np.sin(phase) @ (k * k * d),
    ]
    for table, expected in zip((sf.B, sf.DB, sf.D2B), exact):
        np.testing.assert_allclose(
            table[ij[0], ij[1]], expected, atol=3e-12, rtol=3e-13
        )
    theta, phi = rng.uniform(-5, 10, (2, 80))
    for name, method in (("B", sf.b), ("D_B", sf.db), ("D2_B", sf.d2b)):
        expected = getattr(field, name)(sf.s, theta, phi / sf.nfp)
        error = np.max(np.abs(method(theta, phi) - expected))
        assert error < 1e-7 * np.max(np.abs(expected))
    # Sine coefficients and an unrepresented mode cannot silently disappear.
    sine = Surface(FourierModel(sine=(0, 0.2)), 0.5, 64, 64)
    assert sine.b(np.pi / 2, 0) == pytest.approx(2.2, abs=1e-12)
    with pytest.raises(ValueError, match="resolution"):
        Surface(FourierModel(m=(0, 40)), 0.5, 64, 64)


def test_axisymmetric_ridge_is_one_closed_curve():
    sf = Surface(FourierModel(), 0.5, 64, 64)
    arcs = place_nodes(sf, ridge_arcs(sf), density=32 / (2 * np.pi))
    assert len(arcs) == 1
    a = arcs[0]
    assert a.closed and a.wind == (0, 1)
    np.testing.assert_allclose(a.H, 2.5, atol=1e-12)
    assert a.dalpha.sum() == pytest.approx(2 * np.pi * sf.iota, abs=1e-12)
    # The seam is one complete edge, equal to every other edge in this model.
    np.testing.assert_allclose(a.dalpha, 2 * np.pi * sf.iota / 32, atol=1e-12)


@pytest.mark.parametrize("sign", [-1, 1])
def test_action_and_time_of_cosine_well(sign):
    sf = Surface(FourierModel(C=sign * 3.0), 0.5, 512, 512)
    theta, phi = np.array([0.0]), np.array([0.3])
    chain = trace_chains(sf, theta, phi, periods=3)
    assert np.sign(chain.max_z[0, 0]) == sign
    b = 2.1
    well = wells_at_pitch(sf, theta, phi, np.array([2.5]), chain, b)
    # B = Bmin + 2a sin^2(t/2). Pi(n,m) expressed in Carlson integrals.
    a, minimum = 0.5, 1.5
    m, n = (b - minimum) / (2 * a), -(b - minimum) / minimum
    pi = elliprf(0, 1 - m, 1) + n / 3 * elliprj(0, 1 - m, 1, 1 - n)
    expected_k = 4 * 3 / abs(sf.iota) * np.sqrt(b / (2 * a)) / minimum * pi
    expected_j = expected_k - 4 * 3 / abs(sf.iota) / np.sqrt(2 * a * b) * ellipk(m)
    assert well.status[0] == ACTIVE
    assert well.K[0] == pytest.approx(expected_k, rel=1e-8)
    assert well.J[0] == pytest.approx(expected_j, rel=1e-8)


def test_well_exit_and_cap_statuses():
    # The slow modulation makes consecutive local maxima alternate in height.
    sf = Surface(
        FourierModel(m=(0, 0, 0), n=(0, -1, -5), c=(2, 0.3, 0.5)), 0.5, 128, 128
    )
    th, ph = np.zeros(3), np.array([0.0, 0.0, 0.0])
    chain = trace_chains(sf, th, ph, periods=4, samples_per_period=256)
    wells = wells_at_pitch(sf, th, ph, np.full(3, 2.8), chain, 2.65)
    assert wells.status[0] == ACTIVE and wells.exit_index[0] > 0
    first_above = np.flatnonzero(chain.max_H[0] > 2.65)[0]
    assert wells.exit_index[0] == first_above
    # A hidden above-b maximum must end the well even if a taller maximum follows.
    short = wells_at_pitch(sf, th, ph, np.full(3, 2.8), chain, 2.3)
    assert short.exit_index[0] == 0
    endpoint_chain = trace_chains(sf, th, ph, periods=1, samples_per_period=256)
    endpoint = wells_at_pitch(sf, th, ph, np.full(3, 2.8), endpoint_chain, 2.79)
    # A return exactly at the cap is a resolved exit, rather than an absent one.
    assert np.all(endpoint.status == ACTIVE)
    inactive = wells_at_pitch(sf, th, ph, np.full(3, 2.8), chain, 1.0)
    assert np.all(inactive.status == INACTIVE)

    slow = Surface(FourierModel(iota=0.2), 0.5, 64, 64)
    long_chain = trace_chains(slow, th[:1], ph[:1], periods=8)
    long = wells_at_pitch(slow, th[:1], ph[:1], np.array([2.5]), long_chain, 2.1)
    assert long.status[0] == ACTIVE
    assert abs(long.z_out[0]) > 3 * slow.L
    capped_chain = trace_chains(slow, th[:1], ph[:1], periods=1)
    capped = wells_at_pitch(slow, th[:1], ph[:1], np.array([2.5]), capped_chain, 2.1)
    assert capped.status[0] == CAPPED
    assert np.isnan(capped.K).all() and np.isnan(capped.J).all()


def test_cosine_volume_and_ledger_measure():
    geometry = Geometry(
        FourierModel(),
        [0.1, 0.5, 0.9],
        nt=64,
        nphi=64,
        density=32 / (2 * np.pi),
        periods=3,
    )
    # Integral_0^2pi dtheta/(B0+a cos(theta))².
    expected_volume = 0.8 * 3 * (2 * np.pi) ** 2 * 2 / (4 - 0.25) ** 1.5
    assert geometry.volume(uniform_source) == pytest.approx(expected_volume, rel=1e-12)
    for h in (uniform_source, default_source):
        result = pitch_ledger(geometry, 2.1, h, volume_nphi=32)
        assert result["ratio"] == pytest.approx(1, rel=2e-7)
    # Unknown well weight must invalidate the total, rather than become zero.
    short = Geometry(
        FourierModel(iota=0.2),
        [0.1, 0.9],
        nt=64,
        nphi=64,
        density=16 / (2 * np.pi),
        periods=1,
    )
    result = pitch_ledger(
        short, 2.1, uniform_source, volume_nphi=32, max_periods=1, refine=False
    )
    assert result["capped_nodes"] > 0
    assert np.isnan(result["ridge_weight"]) and np.isnan(result["ratio"])


def test_fold_arcs_and_radial_frontiers():
    from alpha_analysis_3.ridge import radial_partners

    # B=2+.6 cos(theta) cos(phi)+.1 cos(2phi); iota=0.
    sf = Surface(
        FourierModel(m=(0, 1, 1, 0), n=(0, 1, -1, -2), c=(2, 0.3, 0.3, 0.1), iota=0),
        0.5,
        128,
        128,
    )
    arcs = ridge_arcs(sf)
    assert len(arcs) == 2 and all(not a.closed for a in arcs)
    end = np.vstack([a.points[[0, -1]] for a in arcs])
    np.testing.assert_allclose(np.abs(np.cos(end[:, 0])), 2 / 3, atol=0.002)
    nodes = place_nodes(sf, arcs, density=32 / (2 * np.pi))
    for a in nodes:
        assert np.all(sf.d2b(a.points[1:-1, 0], a.points[1:-1, 1]) < 0)
    ids = radial_partners(nodes, nodes)
    np.testing.assert_array_equal(ids, np.arange(len(ids)))
    np.testing.assert_array_equal(radial_partners(nodes, []), -1)


def test_ledger_closes():
    """Fast production path on one manifest field, two global pitch levels."""
    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][2]
    geometry = Geometry(
        BoozerField(ROOT / "data" / entry["file"]),
        np.linspace(0.02, 1, 8),
        nt=128,
        nphi=128,
        density=64 / (2 * np.pi),
        periods=24,
    )
    for level in (0.1, 0.8):
        b = geometry.Bmin + level * (geometry.Bmax - geometry.Bmin)
        for h in (uniform_source, default_source):
            result = pitch_ledger(geometry, b, h, volume_nphi=128)
            assert abs(result["ratio"] - 1) < 0.02, result


def test_shallow_first_minimum_is_not_skipped():
    """PCA regression: a genuine minimum precedes the first scan sample."""
    from scipy.optimize import brentq

    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][0]
    field = BoozerField(ROOT / "data" / entry["file"])
    s, theta, phi = 0.7017391304347826, 2.128024127790614, 1.9521188388690494
    sf = Surface(field, s)
    exact = brentq(
        lambda z: float(field.D_B(s, theta + sf.iota * z, phi / sf.nfp + z)),
        0.002,
        0.012,
    )
    ch = trace_chains(sf, np.array([theta]), np.array([phi]), periods=2)
    assert ch.min_z[0, 0] == pytest.approx(exact, abs=1e-7)
    assert ch.min_z[0, 0] < sf.L / 128
    w = wells_at_pitch(
        sf,
        np.array([theta]),
        np.array([phi]),
        np.array([sf.b(theta, phi)]),
        ch,
        8.147184858800973,
    )
    assert w.status[0] == INACTIVE


def test_real_well_beyond_design_cap_stays_unknown():
    """PCA lambda=.8 needs 74 periods; a 64-period scan must not pass it."""
    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][0]
    sf = Surface(BoozerField(ROOT / "data" / entry["file"]), 0.44608695652173913)
    theta, phi = np.array([9.935229442252941]), np.array([6.903643850431844])
    b = 11.167688974652413
    ch = trace_chains(sf, theta, phi, periods=64)
    w = wells_at_pitch(sf, theta, phi, sf.b(theta, phi), ch, b)
    assert w.status[0] == CAPPED and w.exit_index[0] == -1
    assert np.isnan(w.J[0]) and np.isnan(w.K[0])
    # This extended scan is diagnostic evidence, not a default cap change.
    extended = trace_chains(sf, theta, phi, periods=128)
    w = wells_at_pitch(sf, theta, phi, sf.b(theta, phi), extended, b)
    assert w.status[0] == ACTIVE
    assert extended.max_z[0, w.exit_index[0]] / sf.L == pytest.approx(
        73.80031579, abs=1e-7
    )


def test_adaptive_extension_only_retraces_capped_nodes():
    sf = Surface(FourierModel(iota=0.2), 0.5, 64, 64)
    theta, phi = np.array([0.0, np.pi]), np.zeros(2)
    H = sf.b(theta, phi)
    ch = trace_chains(sf, theta, phi, periods=1)
    well, extended, stats = adaptive_wells_at_pitch(
        sf, theta, phi, H, ch, 2.1, max_periods=8
    )
    np.testing.assert_array_equal(well.status, [ACTIVE, INACTIVE])
    np.testing.assert_array_equal(extended.periods, [8, 1])
    assert stats["extended_nodes"] == 1
    assert ch.periods.tolist() == [1, 1]
    # Releasing extraction tables must not break a later extension.
    sf.release_derivatives()
    released, _, _ = adaptive_wells_at_pitch(sf, theta, phi, H, ch, 2.1, max_periods=8)
    np.testing.assert_allclose(released.K[0], well.K[0], rtol=1e-12)


def test_ceiling_keeps_unknown_with_partial_weight_lower_bound():
    sf = Surface(FourierModel(iota=0.2), 0.5, 128, 128)
    theta, phi = np.array([0.0]), np.array([0.3])
    H = sf.b(theta, phi)
    ch = trace_chains(sf, theta, phi, periods=1)
    a, ca, _ = adaptive_wells_at_pitch(sf, theta, phi, H, ch, 2.1, max_periods=2)
    b, cb, _ = adaptive_wells_at_pitch(sf, theta, phi, H, ch, 2.1, max_periods=3)
    full, _, _ = adaptive_wells_at_pitch(sf, theta, phi, H, ch, 2.1, max_periods=8)
    assert a.status[0] == b.status[0] == CAPPED
    assert ca.periods[0] == 2 and cb.periods[0] == 3
    assert np.isnan(a.J[0]) and np.isnan(a.K[0])
    assert 0 < a.partial_K[0] < b.partial_K[0] < full.K[0]


@pytest.mark.parametrize("strip_rule", ["trapezoid", "simpson"])
def test_long_strip_refinement_preserves_measure_and_reports_level(strip_rule):
    geometry = Geometry(
        FourierModel(iota=0.01),
        [0.1, 0.9],
        nt=64,
        nphi=64,
        density=8 / (2 * np.pi),
        periods=64,
    )
    # The survey may retain the current surface's temporary derivative tables;
    # the same analytic measure must hold with that cache already populated.
    geometry.layers[0].surface.restore_derivatives()
    measure = measure_pitch(
        geometry, 2.1, transverse_tolerance=0.01, strip_rule=strip_rule
    )
    result = pitch_ledger(
        geometry, 2.1, uniform_source, volume_nphi=32, measure=measure
    )
    assert measure.converged and measure.capped_nodes == measure.failed_nodes == 0
    assert result["ratio"] == pytest.approx(1, rel=2e-7)
    for stats in measure.diagnostics:
        assert stats["long_nodes"] == 8
        assert stats["factor"] == 4 and stats["converged"]
    assert result["tail_weights"]["64"] == pytest.approx(
        result["ridge_weight"], rel=1e-12
    )


def test_capped_ledger_is_incomplete_with_partial_lower_bound():
    geometry = Geometry(
        FourierModel(iota=0.001),
        [0.1, 0.9],
        nt=64,
        nphi=64,
        density=8 / (2 * np.pi),
        periods=64,
    )
    measure = measure_pitch(geometry, 2.1, transverse_tolerance=0.01)
    result = pitch_ledger(
        geometry, 2.1, uniform_source, volume_nphi=32, measure=measure
    )
    assert result["capped_nodes"] > 0 and np.isnan(result["ratio"])
    assert result["resolved_weight"] == 0
    assert 0 < result["capped_weight_lower_bound"] < result["volume_weight"]
    assert result["gap_to_volume"] == pytest.approx(
        result["volume_weight"] - result["capped_weight_lower_bound"]
    )
    assert all(
        np.isnan(w.K).all() and np.all(w.status == CAPPED) for w in measure.wells
    )


def test_opposite_strip_errors_do_not_cancel():
    from alpha_analysis_3.ledger import _relative_strip_change

    # Independent polynomials u² and 1-u² on [0,1] have integrals 1/3 and
    # 2/3. One-panel trapezoids give (.5,.5), two panels (.375,.625).
    # Their total is unchanged even though neither strip has converged.
    previous = np.array([[0.5, 0], [0.5, 0]])
    refined = np.array([[0.375, 0], [0.625, 0]])
    assert refined[:, 0].sum() == previous[:, 0].sum()
    assert _relative_strip_change(previous, refined, 1.0) == pytest.approx(0.25)


def test_dyadic_strip_quadrature_matches_polynomial_integrals():
    from alpha_analysis_3.ledger import _strip_integrals

    # Independent integrals on alpha in [0,2]: integral(1+(alpha/2)^3)=5/2,
    # integral((alpha/2)^2)=2/3. Both sampled densities are nonnegative.
    u = np.linspace(0, 1, 5)
    values = np.stack((1 + u**3, u**2), axis=-1)[None]
    result = _strip_integrals(np.array([2.0]), values, 4, "simpson")
    np.testing.assert_allclose(result[0], [2.5, 2 / 3], rtol=1e-14)


def test_slice_error_uses_absolute_changes_and_source_weights():
    from alpha_analysis_3.ledger import slice_refinement_error

    # Simpson or cancellation is not used: radial trapezoids have positive
    # weights (1/4,1/2,1/4). Two surface errors therefore add, rather than cancel.
    s = np.array([0.0, 0.5, 1.0])
    weights = np.zeros((3, 5))
    weights[:, 0] = [1.0, 2.0, 3.0]
    diagnostics = [
        dict(relative_change=0.1),
        dict(relative_change=0.2),
        dict(relative_change=0.1),
    ]
    error, relative = slice_refinement_error(s, weights, diagnostics, uniform_source)
    assert error == pytest.approx(0.3)
    assert relative == pytest.approx(0.15)
    # Source h(rho)=1-rho² gives radial weights (1/4,1/4,0).
    error, relative = slice_refinement_error(
        s, weights, diagnostics, lambda rho: 1 - rho**2
    )
    assert error == pytest.approx(0.125)
    assert relative == pytest.approx(1 / 6)


def test_marginal_surface_needs_tighter_quadrature_control():
    from types import SimpleNamespace
    from alpha_analysis_3.ledger import Layer, node_measure

    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][3]
    sf = Surface(BoozerField(ROOT / "data" / entry["file"]), 0.8540425531914895)
    arcs = place_nodes(sf, ridge_arcs(sf), density=512 / (2 * np.pi))
    points = np.vstack([a.points for a in arcs])
    H = np.concatenate([a.H for a in arcs])
    chain = trace_chains(sf, points[:, 0], points[:, 1], periods=64)
    layer = Layer(sf, arcs, chain, points, H, node_measure(arcs))
    # Independent angular volume, rather than another ridge discretization.
    b = 3.316712409154169
    reference = volume_surface(sf, b, nphi=2048, order=96)
    result = measure_pitch(
        SimpleNamespace(layers=[layer]),
        b,
        transverse_tolerance=0.0025,
        acceptance_tolerance=0.01,
    )
    assert result.capped_nodes == result.failed_nodes == 0
    assert abs(result.weights[0, 0] / reference - 1) < 0.02
