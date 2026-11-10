"""Reproduce the M1 cap obstruction without changing the production cap.

The two PCA nodes are preserved from the failed 30-case ledger. This is a
diagnostic experiment at 256 periods, not a new accepted production setting.
"""

import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.chains import trace_chains
from alpha_analysis_3.tables import Surface
from alpha_analysis_3.wells import wells_at_pitch
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]
PROBES = [
    dict(
        s=0.44608695652173913,
        theta=9.935229442252941,
        phi=6.903643850431844,
        b=11.167688974652413,
        lambda_n=0.8,
    ),
    dict(
        s=0.8295652173913044,
        theta=15.982749366888166,
        phi=13.057817282260391,
        b=11.81493985662058,
        lambda_n=0.95,
    ),
]


def main():
    entry = json.loads((ROOT / "data/manifest.json").read_text())["files"][0]
    field = BoozerField(ROOT / "data" / entry["file"])
    report = dict(
        data_sha256=entry["sha256"],
        source_hashes={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT / "src/alpha_analysis_3").glob("*.py"))
        },
        probes=[],
    )
    with Watchdog() as guard:
        start = time.monotonic()
        for probe in PROBES:
            sf = Surface(field, probe["s"])
            th, ph = np.array([probe["theta"]]), np.array([probe["phi"]])
            H = sf.b(th, ph)
            record = dict(probe, scans=[])
            for periods, samples in [(48, 128), (64, 128), (256, 128), (128, 512)]:
                guard.stage = (
                    f'PCA {probe["lambda_n"]}: {periods} periods, {samples} samples'
                )
                t = time.monotonic()
                chain = trace_chains(sf, th, ph, periods, samples)
                well = wells_at_pitch(sf, th, ph, H, chain, probe["b"])
                indices = np.flatnonzero(chain.max_H[0] > probe["b"])
                row = dict(
                    periods=periods,
                    samples_per_period=samples,
                    status=int(well.status[0]),
                    seconds=time.monotonic() - t,
                )
                if len(indices):
                    exit_index = int(indices[0])
                    exit_z = float(chain.max_z[0, exit_index])
                    # Independent exact-series roots at every candidate maximum
                    # up to and including the exit; no table interpolation here.
                    heights = []
                    for z in chain.max_z[0, : exit_index + 1]:
                        root = brentq(
                            lambda x: float(
                                field.D_B(
                                    probe["s"], th[0] + sf.iota * x, ph[0] / sf.nfp + x
                                )
                            ),
                            z - sf.L / samples / 2,
                            z + sf.L / samples / 2,
                            xtol=1e-13,
                        )
                        heights.append(
                            float(
                                field.B(
                                    probe["s"],
                                    th[0] + sf.iota * root,
                                    ph[0] / sf.nfp + root,
                                )
                            )
                        )
                    assert max(heights[:-1]) < probe["b"] < heights[-1]
                    coarse = wells_at_pitch(sf, th, ph, H, chain, probe["b"], order=32)
                    row.update(
                        exit_index=exit_index,
                        exit_periods=exit_z / sf.L,
                        exact_exit_height=heights[-1],
                        exact_largest_prior_height=max(heights[:-1]),
                        J=float(well.J[0]),
                        K=float(well.K[0]),
                        relative_J_order32=float(abs(coarse.J[0] / well.J[0] - 1)),
                        relative_K_order32=float(abs(coarse.K[0] / well.K[0] - 1)),
                    )
                record["scans"].append(row)
            report["probes"].append(record)
        report.update(
            total_seconds=time.monotonic() - start,
            peak_rss_gb=guard.peak,
            snapshots=guard.snapshots,
        )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
