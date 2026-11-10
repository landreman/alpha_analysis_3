# Design-phase prototype

`proto_f.py` and `proto_f2.py` are the exploratory implementation that produced the
evidence in [docs/validation/prototype](../docs/validation/prototype/README.md) on
2026-11-09/10. They are kept as a **reference**, not as production code:

- They are not installed, not formatted, not covered by the test suite, and never imported
  from `src/` or `tests/`.
- They use the production field reader (`alpha_analysis_3.boozer_field`) and the data
  files in `data/`; the import paths were adapted for this repository and nothing else was
  changed.
- Their one-directional partner matching, point-sampled geometry and lack of bisection,
  polishing and subdivision are known limitations that DESIGN §4 marks as production
  refinements. The dyadic bin cap, the memory watchdog and the interior-barrier transition
  test are the parts that were revised after measurements.

Run one case at a time (each run needs 0.4–0.9 GB; eleven in parallel once filled the
laptop's swap):

```bash
source .venv/bin/activate
cd prototype
MEM_CAP_GB=3 python proto_f2.py DMerc 24 128 24 0.02 0 4    # base resolution, ~10 s
MEM_CAP_GB=3 python proto_f2.py d23p4 48 256 48 0.02 0 4    # doubled, ~1 min
```

Arguments: file key (`DMerc`, `d23p4`, `PCA`, `TURBO`, `n3are`), number of surfaces, ridge
nodes per 2π, pitch slices, s_min, extra graded inner surfaces, bins per cell. Θ maps are
written to `prototype/out/`, which is ignored by Git. Environment variables
`MAX_BINS_PER_CELL` (default 128) and `MEM_CAP_GB` (default 3) control the graph cap and
the memory watchdog.

When porting a function into `src/`, cite this origin in the docstring and add the
independent test that DESIGN requires; the prototype's own agreement with itself is not
evidence.
