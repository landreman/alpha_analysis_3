"""Independent angular refinement for the difficult production tail probes.

Reads measured ridge weights; recomputes only the independent volume reference.
No missing ridge weight is inferred from that reference.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.conventions import default_source, uniform_source, source_at_s
from alpha_analysis_3.ledger import volume_surface
from alpha_analysis_3.tables import Surface
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = {
        x["id"]: x
        for x in json.loads((ROOT / "data/manifest.json").read_text())["files"]
    }
    records = []
    for name, levels in [
        ("DMercFail", [0.97, 0.99]),
        ("d23p4", [0.9, 0.95, 0.98, 0.99]),
    ]:
        measured = json.loads(
            (ROOT / "outputs" / f"m1-dense-{name}" / f"{name}.json").read_text()
        )
        checkpoint = np.load(
            ROOT / "outputs" / f"m1-dense-{name}" / f"{name}-checkpoint.npz"
        )
        field = BoozerField(ROOT / "data" / manifest[name]["file"])
        s = np.array(measured["s"])
        by_level = {x["lambda_n"]: x for x in measured["rows"]}
        references = np.zeros((3, len(levels), len(s)))
        with Watchdog() as guard:
            for i, radius in enumerate(s):
                sf = Surface(field, radius)
                for j, level in enumerate(levels):
                    b = by_level[level]["b"]
                    for k, nphi in enumerate((1024, 2048, 4096)):
                        guard.stage = (
                            f"{name} s={radius:.5f}, lambda={level}, phi={nphi}"
                        )
                        references[k, j, i] = volume_surface(sf, b, nphi, order=96)
            for j, level in enumerate(levels):
                row = by_level[level]
                ids = np.flatnonzero(
                    np.isclose(checkpoint["levels"], level, rtol=0, atol=1e-12)
                )
                raw = checkpoint["weights"][ids[0]]
                for k, nphi in enumerate((1024, 2048, 4096)):
                    source_results = {}
                    for source, h in [
                        ("uniform", uniform_source),
                        ("peaked", default_source),
                    ]:
                        weights = source_at_s(h, s)
                        volume = float(trapezoid(weights * references[k, j], s))
                        resolved, lower = row["source_weights"][source][:2]
                        source_results[source] = dict(
                            resolved_weight=resolved,
                            capped_lower_bound=lower,
                            volume_weight=volume,
                            ratio=(
                                resolved / volume if row["capped_nodes"] == 0 else None
                            ),
                            resolved_plus_lower_ratio=(resolved + lower) / volume,
                            lower_bound_check=bool(resolved + lower <= volume * 1.01),
                        )
                    record = dict(
                        file=name,
                        lambda_n=level,
                        b=row["b"],
                        nphi=nphi,
                        order=96,
                        data_sha256=manifest[name]["sha256"],
                        capped_nodes=row["capped_nodes"],
                        sources=source_results,
                        per_surface=[
                            dict(
                                s=float(radius),
                                ridge=float(raw[i, 0]),
                                lower=float(raw[i, 1]),
                                volume=float(references[k, j, i]),
                            )
                            for i, radius in enumerate(s)
                        ],
                    )
                    records.append(record)
                    print(
                        json.dumps(
                            {x: y for x, y in record.items() if x != "per_surface"}
                        ),
                        flush=True,
                    )
            checkpoint.close()
    result = dict(
        source_hashes={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT / "src/alpha_analysis_3").glob("*.py"))
        },
        records=records,
        error_scope="Independent angular-refinement estimates, not enclosures. Unknown full K stays unavailable.",
    )
    (ROOT / "outputs/m1-tail-ledger-refinement.json").write_text(
        json.dumps(result, indent=2)
    )


if __name__ == "__main__":
    main()
