"""Coverage of unequal ridge discretizations, including periodic seams."""

import numpy as np
import pytest

from alpha_analysis_3.cells import match_layers, stitch_ridges
from alpha_analysis_3.ridge import Ridge


def circle(n, shift=0):
    phi = np.arange(n) * 2 * np.pi / n + shift
    return Ridge(np.column_stack((np.zeros(n), phi)), True, (0, 1))


def test_bidirectional_coverage_includes_unreached_top_nodes():
    bottom, top = [circle(8)], [circle(19, 0.02)]
    match = match_layers(bottom, top, density=128 / (2 * np.pi))
    assert len(np.unique(match.forward)) < 19
    mesh = stitch_ridges(bottom, top, 0.2, 0.6, iota_bottom=0.7, iota_top=0.7)
    assert not mesh.frontier_bottom.any() and not mesh.frontier_top.any()
    assert set(mesh.faces.ravel()) == set(range(27))
    # Each boundary edge occurs once; interior triangle edges twice.
    edges = {}
    for face in mesh.faces:
        for a, b in zip(face, np.roll(face, -1)):
            key = tuple(sorted((a, b)))
            edges[key] = edges.get(key, 0) + 1
    assert sum(v == 1 for v in edges.values()) == 8 + 19
    assert all(v in (1, 2) for v in edges.values())
    assert mesh.measure.sum() == pytest.approx(0.4 * 0.7 * 2 * np.pi, rel=1e-12)


def test_unmatched_arc_is_explicit_frontier():
    bottom = [circle(8)]
    top = [
        circle(8),
        Ridge(np.column_stack((np.full(8, 2.0), np.linspace(0, 1, 8))), False, (0, 0)),
    ]
    mesh = stitch_ridges(bottom, top, 0.2, 0.6, iota_bottom=0.7, iota_top=0.7)
    assert np.all(mesh.frontier_top[8:])
    assert not mesh.frontier_top[:8].any()
    assert np.all(np.isfinite(mesh.measure))


def test_partner_threshold_scales_with_density():
    a, b = [circle(16)], [circle(16)]
    b[0].points[:, 0] += 0.15
    assert np.all(match_layers(a, b, density=128 / (2 * np.pi)).forward >= 0)
    assert np.all(match_layers(a, b, density=512 / (2 * np.pi)).forward == -1)


def test_open_arcs_align_orientation_and_include_polished_folds():
    from alpha_analysis_3.cells import polish_folds
    from alpha_analysis_3.ledger import Geometry
    from test_m1_geometry import FourierModel

    field = FourierModel(
        m=(0, 0, 1, 1, 0),
        n=(0, -1, -1, 1, -2),
        c=(6, 0, 0, 0, 0),
        sine=(0, -2, 0.25, -0.25, 1),
        iota=0,
    )
    g = Geometry(field, [0.2, 0.6], nt=64, nphi=64, density=32 / (2 * np.pi), periods=3)
    records = polish_folds(g)
    assert records and all(r["resolved"] for r in records)
    for layer in g.layers:
        for arc in layer.arcs:
            if not arc.closed:
                np.testing.assert_allclose(np.cos(arc.endpoints[:, 0]), 0, atol=1e-10)
                np.testing.assert_allclose(np.sin(arc.endpoints[:, 1]), 0, atol=1e-10)
    a, c = [l.arcs for l in g.layers]
    for arc in c:
        if not arc.closed:
            arc.points = arc.points[::-1].copy()
            arc.endpoints = arc.endpoints[::-1].copy()
    mesh = stitch_ridges(a, c, 0.2, 0.6, 0, 0, density=32 / (2 * np.pi))
    assert not mesh.frontier_bottom.any() and not mesh.frontier_top.any()
    assert not mesh.fold_gaps
    assert mesh.measure.sum() == pytest.approx(0.4 * 3 * np.pi, rel=1e-10)


def test_transition_cut_cells_keep_separate_banks():
    from alpha_analysis_3.cells import SurfaceCache, cut_cells
    from test_m2_transitions import SplitField

    sampler = SurfaceCache(SplitField(), table=64)
    face = np.array([[[0.3, 1.3, 0], [0.6, 1.3, 0], [0.3, 1.8, 0]]])
    result = cut_cells(sampler, face, 2.6, max_depth=2)
    assert len(result.cells) == 2
    assert all(c.kind == "bank" for c in result.cells)
    assert len(result.transitions) == 1
    assert all(c.parent_face == 0 and sum(c.cut_edges) == 1 for c in result.cells)
    assert len(set(result.cells[0].vertex_ids) & set(result.cells[1].vertex_ids)) == 2
    for event in result.transitions[0]:
        assert event.parent == pytest.approx(event.child_1 + event.child_3, rel=2e-6)

    def area(p):
        a = p[1:, :2] - p[:1, :2]
        # Shoelace works on the quad as well as the triangle bank.
        return (
            abs(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)))
            / 2
        )

    assert sum(area(c.points) for c in result.cells) == pytest.approx(
        area(face[0]), abs=1e-8
    )


def test_shared_transition_edge_has_identical_incidence():
    from alpha_analysis_3.cells import SurfaceCache, cut_cells
    from test_m2_transitions import SplitField

    sampler = SurfaceCache(SplitField(), table=64)
    points = np.array([[0.3, 1.3, 0], [0.6, 1.3, 0], [0.6, 1.8, 0], [0.3, 1.8, 0]])
    result = cut_cells(sampler, points[[[0, 1, 2], [0, 2, 3]]], 2.6)
    assert len(result.cells) == 4 and len(result.transitions) == 2
    first, second = result.transitions
    common = [
        a.point for a in first for c in second if np.array_equal(a.point, c.point)
    ]
    assert len(common) == 1
    shared = [
        c.vertex_ids[np.all(c.points == common[0], axis=1)][0] for c in result.cells
    ]
    assert len(set(shared)) == 1


def test_one_cut_cell_subdivides_and_retains_residual_unknown():
    from alpha_analysis_3.cells import SurfaceCache, cut_cells
    from test_m2_transitions import SplitField

    sampler = SurfaceCache(SplitField(), table=64)
    face = np.array([[[0.3, 1.3, 0], [0.6, 1.8, 0], [0.3, 1.8, np.pi]]])
    result = cut_cells(sampler, face, 2.6, max_depth=1)
    assert len(result.cells) >= 4
    assert any(c.kind == "unknown" for c in result.cells)
    assert all(c.depth == 1 for c in result.cells)


def test_three_cut_cell_subdivides_at_crossing_barriers():
    from alpha_analysis_3.cells import SurfaceCache, cut_cells
    from alpha_analysis_3.transitions import sample_states, interior_barrier
    from test_m2_transitions import SplitField

    class ThreeWells(SplitField):
        xm = np.array([0, 0, 0, 1, 1, 1, 1])
        xn = np.array([0, 3, 1, 1, -1, 2, -2])

        def bmnc(self, s):
            q = 0.2 - s
            return np.array([2 + q, 0.8, -q, 0.075, -0.075, -0.0375, 0.0375])

        def bmns(self, s):
            return np.zeros(7)

    sampler = SurfaceCache(ThreeWells(), table=64)
    points = np.array([[0.2, 0, 0], [0.4, np.arcsin(-2 / 3), 0], [0.5, 0, 0]])
    states = [sample_states(sampler(p[0]), [p[1:]], 2.6)[0] for p in points]
    assert [s.exit_index for s in states] == [0, 1, 2]
    assert all(
        interior_barrier(states[j], states[(j + 1) % 3]) is not None for j in range(3)
    )
    result = cut_cells(sampler, [points], 2.6, max_depth=1)
    assert len(result.cells) >= 4
    assert all(c.depth == 1 for c in result.cells)
    assert any(c.kind == "unknown" for c in result.cells)


def test_coverage_closes_ledger():
    from alpha_analysis_3.cells import coverage_ledger
    from alpha_analysis_3.ledger import Geometry, measure_pitch
    from alpha_analysis_3.conventions import uniform_source
    from test_m1_geometry import FourierModel

    g = Geometry(
        FourierModel(),
        [0.1, 0.5, 0.9],
        nt=128,
        nphi=128,
        density=16 / (2 * np.pi),
        periods=3,
    )
    measure = measure_pitch(g, 2.1, retain_edges=True)
    result = coverage_ledger(g, measure, uniform_source, density=16 / (2 * np.pi))
    expected = 0.8 * measure.weights[0, 0]
    assert result["cell_weight"] == pytest.approx(expected, rel=1e-12)
    # Independent direct line quadrature of the prescribed cosine field.
    from scipy.integrate import quad

    bounce = np.arccos((2.1 - 2) / 0.5)

    def integrand(theta):
        B = 2 + 0.5 * np.cos(theta)
        return 3 / 0.7 / B / np.sqrt(1 - B / 2.1)

    K = quad(integrand, bounce, 2 * np.pi - bounce, epsabs=1e-9, epsrel=1e-10)[0]
    assert result["cell_weight"] == pytest.approx(0.8 * 0.7 * 2 * np.pi * K, rel=1e-8)
    assert result["frontier_fraction"] == 0


def test_surface_insertion_repairs_displaced_ridges():
    from alpha_analysis_3.cells import insert_surfaces, geometric_frontier
    from alpha_analysis_3.ledger import Geometry
    from test_m1_geometry import FourierModel

    class MovingRidge(FourierModel):
        def bmnc(self, s):
            return np.array([2, 0.5 * np.cos(2 * s)])

        def bmns(self, s):
            return np.array([0, 0.5 * np.sin(2 * s)])

    f = MovingRidge()
    g = Geometry(f, [0.1, 0.9], nt=64, nphi=64, density=128 / (2 * np.pi), periods=3)
    assert geometric_frontier(g)[0]["frontier_measure"] > 0
    insert_surfaces(g, f, np.arange(0.2, 0.9, 0.1), density=128 / (2 * np.pi))
    assert sum(r["frontier_measure"] for r in geometric_frontier(g)) == 0
