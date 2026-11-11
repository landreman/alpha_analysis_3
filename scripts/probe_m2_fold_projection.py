"""Diagnose fold-strip projections without assigning failed points a weight."""

import argparse
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np

from alpha_analysis_3.boozer_field import BoozerField
from alpha_analysis_3.cells import polish_folds
from alpha_analysis_3.ridge import ridge_arcs, place_nodes
from alpha_analysis_3.tables import Surface
from validate_m1 import Watchdog

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    source = json.loads(Path(args.report).read_text())
    row = next(r for r in source["rows"] if r["file"] == args.file)
    strips = row["sources"]["uniform"]["strips"]
    radii = [strips[0]["s_bottom"]] + [r["s_top"] for r in strips]
    entry = next(
        e
        for e in json.loads((ROOT / "data/manifest.json").read_text())["files"]
        if e["id"] == args.file
    )
    report = dict(
        revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        controls=vars(args),
        data_sha256=entry["sha256"],
        rows=[],
    )
    with Watchdog() as guard:
        field = BoozerField(ROOT / "data" / entry["file"])
        for s in radii:
            guard.stage = f"{args.file} fold projection s={s}"
            sf = Surface(
                field, s, source["controls"]["table"], source["controls"]["table"]
            )
            arcs = place_nodes(
                sf, ridge_arcs(sf), source["controls"]["nodes"] / (2 * np.pi)
            )
            polish_folds(
                SimpleNamespace(layers=[SimpleNamespace(surface=sf, arcs=arcs)])
            )
            for j, arc in enumerate(arcs):
                if arc.closed:
                    continue
                for end, node in ((0, 0), (1, -1)):
                    fold, p = arc.endpoints[end], arc.points[node]
                    guesses = fold + (p - fold) * np.arange(1, 16)[:, None] / 16
                    outcomes = {}
                    for mode in ("table", "exact"):
                        points = guesses.copy()
                        for _ in range(20 if mode == "exact" else 12):
                            d1 = (
                                sf.db(points[:, 0], points[:, 1])
                                if mode == "table"
                                else sf.exact(points[:, 0], points[:, 1], 1)
                            )
                            d2 = (
                                sf.d2b(points[:, 0], points[:, 1])
                                if mode == "table"
                                else sf.exact(points[:, 0], points[:, 1], 2)
                            )
                            step = np.divide(
                                d1, d2, out=np.zeros_like(d1), where=abs(d2) > 1e-14
                            )
                            points -= np.clip(step, -0.02, 0.02)[:, None] * np.array(
                                [sf.iota, sf.nfp]
                            )
                        d1 = sf.exact(points[:, 0], points[:, 1], 1)
                        d2 = sf.exact(points[:, 0], points[:, 1], 2)
                        valid = (abs(d1) < 1e-8) & (d2 < 0)
                        outcomes[mode] = dict(
                            valid=int(valid.sum()),
                            max_shift=float(
                                np.max(np.linalg.norm(points - guesses, axis=1))
                            ),
                            DB=d1.tolist(),
                            D2B=d2.tolist(),
                            points=points.tolist(),
                        )
                    if outcomes["table"]["valid"] < 15 or float(sf.exact(*p, 2)) >= 0:
                        record = dict(
                            s=s,
                            arc=j,
                            end=end,
                            fold=fold.tolist(),
                            node=p.tolist(),
                            node_D2=float(sf.exact(*p, 2)),
                            outcomes=outcomes,
                        )
                        report["rows"].append(record)
                        print(
                            s,
                            j,
                            end,
                            "node D2",
                            record["node_D2"],
                            {
                                k: (v["valid"], v["max_shift"])
                                for k, v in outcomes.items()
                            },
                            flush=True,
                        )
        report.update(peak_rss_gb=guard.peak, snapshots=guard.snapshots)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--report", required=True)
    p.add_argument("--file", default="n3are")
    p.add_argument("--output", default="outputs/m2/fold-projection.json")
    run(p.parse_args())
