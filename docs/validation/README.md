# Numerical evidence

Commit concise reports with a reproducible driver, code/data/configuration identifiers,
source profile, actual b values, resolution controls, all outcomes including failures,
numerical scope, timings, worker counts and any active budget snapshots.

Store large arrays and generated plots under the ignored `outputs/` directory or in CI
artifacts; keep committed figures few and small. Evidence must be reproducible from this
repository alone.

- [bootstrap.md](bootstrap.md): M0 checks.
- [prototype/](prototype/README.md): the pre-repository prototype's five-file results,
  logs and Θ maps. Evidence for the design, not acceptance evidence for any milestone.

Full 30-case matrices and five-file convergence tables run outside the unit-test tiers
(fast under 120 s, full under 300 s).
