# M2: transitions and ridge-cell coverage

M2 is in progress. The acceptance row in STATUS remains unchecked. These
records retain the initial failures; they do not claim a completed milestone,
connectivity, Q(b), f, or the final accuracy/runtime target.

The baseline is `27f251d37a3ddd64784f5cead2d68a5563ba88a0`, with green
[GitHub Tests](https://github.com/landreman/alpha_analysis_3/actions/runs/43065702135).
The five exact manifest files, global extrema and six conserved bounce fields
are used throughout. Both sources are reported separately.

## Implementation and independent checks

- `transitions.py` requires the short exit to occur inside the long well's
  chain. It bisects the edge, evaluates the two finite banks and independently
  traces the other child from the interior maximum. The additivity residual
  must be at most 2e-6 relative; event-coordinate brackets are at most 1e-8.
  A finite-action fallback does not repair a failed time integral or its status.
- The synthetic two-well field has analytically known transition coordinates
  and independently integrated actions. The legacy DMercFail event is found at
  its recorded actual b, with its recorded parent action. Smooth bounce-point
  motion is not a transition. A three-well crossing forces cell subdivision.
- `cells.py` matches both directions, scales the partner distance with density,
  and triangulates strips with unequal node counts. Integer vertex IDs,
  parent-face IDs, cut-edge flags and shared transition locations are retained.
  One-/three-cut cells subdivide; residual cells stay unknown.
- Open arcs are oriented before stitching. Generic fold endpoints are polished
  on the exact Fourier series using D_parallel B = D_parallel² B = 0 and a
  nonzero third derivative. Resolved fold end strips are included as cells;
  failed or degenerate folds remain frontier. An analytic fold model recovers
  the full area 3 pi Delta s, including a reversed top arc.
- The M2 population check assigns each surface-edge integral once to its
  incident radial boundary triangle, with positive Delta s/2 weights. Unmeshed
  frontiers remain in the ledger. It uses no independent-volume gap to fill
  missing population. M3 still needs clipped cell weights with interior K
  nodes; this boundary rule is not a claim about within-cell loss weights.
- The independent cosine-field line integral checks cell-ledger normalization.
  Capped full weights remain unknown; partial-K contributions are lower bounds.

## Preserved first attempts

The [initial census](base-census-factor4.csv) evaluates all 840 combinations
of 28 surfaces and the 30 manifest/pitch cases. It uses 128 nodes per 2 pi,
512² tables, order 32 and up to fourfold edge refinement. Twenty-eight
surfaces retain count, assignment or localization uncertainty. The controls,
source hashes, resource records and raw-report hash are in
[initial-census-provenance.json](initial-census-provenance.json).

The [initial coverage sweep](coverage-initial.csv) uses up to 64 adaptive
surfaces and the same node density. All 30 cases and both sources are retained,
including every miss of the 1% ledger, 0.5% frontier or 1% transverse-change
criteria. The [provenance](coverage-initial-provenance.json) records the exact
controls and active snapshots. This attempt preceded open-arc orientation and
polished fold-strip repair. In particular, PCA lambda_n=0.05 has partnered
nodes but 1.42% peaked-source weight in the omitted fold end strips. Adding
radial surfaces cannot repair those angular gaps.

The initial PCA 0.95 transverse change is 6.62% for the peaked source. This
remains a failed estimate. Its population excess must be separated from
independent angular-reference error by refinement; it is not attributed to
missing coverage alone.

## Reproduction

Run from an installed clean environment, one equilibrium at a time:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/validate_m2.py \
  --output outputs/m2/transitions.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/validate_m2_coverage.py \
  --nodes 256 --volume-phi 2048 --cut-probe --plots outputs/m2/plots \
  --output outputs/m2/coverage.json
```

The drivers enforce a 2 GB resident-memory watchdog, retain active
60/180/300-second snapshots and exit nonzero on unresolved acceptance checks.
Full raw records and larger artifacts stay outside Git.
