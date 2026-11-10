"""M0: the accepted tolerance, pitch levels, source profile and direction sign."""

import numpy as np
import pytest

from alpha_analysis_3 import conventions as c


def test_error_target_and_pitch_conventions():
    # max(0.02 f, 0.01): the floor is an error magnitude at small f.
    assert c.numerical_error_target(0.0) == 0.01
    assert c.numerical_error_target(0.5) == pytest.approx(0.01)
    assert c.numerical_error_target(1.0) == pytest.approx(0.02)
    np.testing.assert_allclose(
        c.numerical_error_target(np.array([0.0, 0.25, 0.75])), [0.01, 0.01, 0.015]
    )

    assert c.LAMBDA_LEVELS == (0.05, 0.1, 0.5, 0.8, 0.9, 0.95)
    b = c.pitch_from_lambda(np.asarray(c.LAMBDA_LEVELS), 2.0, 4.0)
    np.testing.assert_allclose(b, [2.1, 2.2, 3.0, 3.6, 3.8, 3.9])
    np.testing.assert_allclose(c.lambda_from_pitch(b, 2.0, 4.0), c.LAMBDA_LEVELS)

    assert c.RUNTIME_TARGET_SECONDS == 300.0
    assert c.MEMORY_TARGET_GB == 2.0


def test_source_is_evaluated_at_rho_not_s():
    assert c.default_source(0.0) == 1.0
    assert c.default_source(1.0) == 0.0
    s = 0.25  # rho = 0.5
    assert c.source_at_s(c.default_source, s) == pytest.approx((1 - 0.25) ** 4)
    # Evaluating at s instead of rho would give (1 - 0.0625)^4; the two must differ.
    assert c.source_at_s(c.default_source, s) != pytest.approx((1 - s**2) ** 4)
    np.testing.assert_array_equal(c.uniform_source(np.linspace(0, 1, 5)), 1.0)


def test_field_line_direction_sign():
    assert c.field_line_direction(1.0, 0.0, 0.5) == 1.0
    assert c.field_line_direction(-1.0, 0.0, 0.5) == -1.0
    # G + iota I with I large and negative flips the direction.
    assert c.field_line_direction(1.0, -4.0, 0.5) == -1.0
    with pytest.raises(ValueError):
        c.field_line_direction(1.0, -2.0, 0.5)
