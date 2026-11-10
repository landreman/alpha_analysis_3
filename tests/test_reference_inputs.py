"""M0: the copied inputs are intact and the field reader's conventions hold.

Nothing here is evidence about f. See docs/DESIGN.md section 9, M0.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.conventions import LAMBDA_LEVELS, pitch_from_lambda

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((ROOT / "data/manifest.json").read_text())


@pytest.mark.parametrize("entry", MANIFEST["files"], ids=lambda row: row["id"])
def test_five_real_inputs_and_fourier_conventions(entry):
    path = ROOT / "data" / entry["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    assert path.stat().st_size == entry["bytes"]
    field = BoozerField.from_boozmn(path)
    assert field.nfp == entry["nfp"]

    # B at a few points equals the direct mode sum of the file's coefficients.
    s = float(field.s_half[len(field.s_half) // 2])
    theta = np.array([0.17, 1.1, 2.7])
    zeta = np.array([0.13, 0.23, 0.37])
    phase = np.outer(theta, field.xm) - np.outer(zeta, field.xn)
    expected = np.cos(phase) @ field.bmnc(s) + np.sin(phase) @ field.bmns(s)
    np.testing.assert_allclose(field.B(s, theta, zeta), expected, rtol=1e-13)
    # Periodicity over one field period in zeta.
    np.testing.assert_allclose(
        field.B(s, theta, zeta + 2 * np.pi / field.nfp), expected, rtol=1e-13
    )
    # The along-line derivative D_par B = iota dB/dtheta + dB/dzeta by finite differences.
    step = 1e-6
    iota = field.iota(s)
    derivative = (
        field.B(s, theta + iota * step, zeta + step)
        - field.B(s, theta - iota * step, zeta - step)
    ) / (2 * step)
    np.testing.assert_allclose(
        field.D_B(s, theta, zeta), derivative, atol=2e-8, rtol=1e-7
    )
    # The six diagnostic pitches are ordered and inside the legacy extrema.
    pitches = pitch_from_lambda(
        np.asarray(LAMBDA_LEVELS),
        entry["legacy_B_min_estimate"],
        entry["legacy_B_max_estimate"],
    )
    assert len(pitches) == 6 and np.all(np.diff(pitches) > 0)
    assert np.isfinite(field.B([0, 1], 0.2, 0.3)).all()


def test_reference_documents_unchanged():
    provenance = json.loads((ROOT / "docs/references/provenance.json").read_text())
    for entry in provenance["verified_unchanged"]:
        digest = hashlib.sha256((ROOT / entry["path"]).read_bytes()).hexdigest()
        assert digest == entry["sha256"], entry["path"]
