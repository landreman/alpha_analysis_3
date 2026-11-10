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
from alpha_analysis_3.ledger import Geometry, pitch_ledger, volume_surface
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


def test_tables_match_exact_series():
    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][0]
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
    result = pitch_ledger(short, 2.1, uniform_source, volume_nphi=32)
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
