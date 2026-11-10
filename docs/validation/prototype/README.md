# Prototype evidence (2026-11-09/10)

Runs of `prototype/proto_f2.py` on the five reference files, single process, one run at a
time, Apple M1 laptop with 16 GB. Base = 24 surfaces (s_min = 0.02), 128 nodes per 2π of
ridge length, 24 pitch slices, 4 action bins per median cell range; doubled = 48 / 256 /
48 / 4. Bins are capped dyadically at 128 per cell. Logs are in [logs/](logs/); the
command lines are recorded in [prototype/README.md](../../../prototype/README.md).

| file | run | ledger (h = 1) | f (h = 1) lower..upper | trapped fraction (h = 1) | f (peaked h) lower..upper | unknown weight (peaked) | wall time [s] |
|---|---|---|---|---|---|---|---|
| DMerc | base | 0.997–1.002 | 0.240..0.241 | 0.627 | 0.0381..0.0385 | 0.0004 | 10 |
| DMerc | doubled | 0.997–1.002 | 0.242..0.242 | 0.631 | 0.0402..0.0406 | 0.0003 | 64 |
| d23p4 | base | 0.965–0.999 | 0.256..0.261 | 0.399 | 0.1431..0.1468 | 0.0026 | 8 |
| d23p4 | doubled | 0.934–1.008 | 0.256..0.258 | 0.398 | 0.1415..0.1440 | 0.0011 | 53 |
| d23p4 | s_min = 0.001 (+4 inner surfaces) | 0.965–0.999 | 0.252..0.257 | 0.397 | 0.1353..0.1390 | 0.0025 | 9 |
| PCA | base | 0.927–1.005 | 0.153..0.218 | 0.388 | 0.0662..0.0968 | 0.0086 | 11 |
| PCA | doubled | 0.970–1.010 | 0.143..0.205 | 0.388 | 0.0551..0.0877 | 0.0030 | 63 |
| TURBO | base | 0.944–1.055 | 0.094..0.098 | 0.368 | 0.0036..0.0090 | 0.0054 | 8 |
| TURBO | doubled | 0.903–1.001 | 0.083..0.085 | 0.371 | 0.0024..0.0056 | 0.0032 | 46 |
| n3are | base | 0.982–1.047 | 0.087..0.252 | 0.286 | 0.0065..0.2033 | 0.0239 | 10 |
| n3are | doubled | 0.981–1.032 | 0.106..0.233 | 0.284 | 0.0173..0.1661 | 0.0093 | 52 |

"Ledger" is the ratio of summed cell weights to the independent volume integral of trapped
weight, range over the six λₙ levels. Previous-team trapped fractions (h = 1): DMercFail
0.630, d23p4 0.396, PCA 0.385, TURBO 0.367, n3are 0.278.

What the table shows and does not show:

- DMercFail and d23p4 change by ≤ 0.003 between base and doubled for both sources, with
  brackets ≤ 0.005 wide: converged at the target level, by the design's own estimate.
- PCA and n3are have small unknown weight but wide brackets because one-directional
  nearest-node matching leaves top nodes of reconnecting ridge arcs uncovered, and each
  uncovered node is a wildcard in the upper bound. The brackets narrow with resolution, so
  the width is a coverage artefact; closing it is milestone M2.
- TURBO moves 12 % between resolutions for the same reason.
- The capped graph and the uncapped graph agree to 0.0006 on the eight cases run both ways.
- No bisection, no polishing, no adaptivity, no subdivision of 1-/3-cut cells: the
  prototype is a lower bound on what the production code should achieve, not acceptance.

Θ maps (left panel lower bound, right panel upper bound, colour = reachable fraction of the
cell's J range, −0.2 = unknown; horizontal axis α mod 2π, vertical axis s):

- `d23p4_…_lam0.1.png`, `lam0.5.png`, `lam0.9.png`: W7-X at three pitches.
- `DMerc_…_lam0.9.png`: barely trapped DMercFail population with confined pockets.
- `PCA_…_lam0.5.png`: confined core for s < 0.7, upper bound inflated by uncovered nodes.
- `TURBO_…_lam0.8.png`: losses only for s > 0.8.
- `n3are_…_lam0.5.png`: the coverage problem, lower bound confined, upper bound lost.
