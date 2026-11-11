# Milestone status

[DESIGN.md section 9](DESIGN.md#9-milestones) defines the sequence. Only this
repository's M0–M5 milestones apply; the milestones of `landreman/alpha_analysis` and
`landreman/alpha_analysis_2` are neither prerequisites nor completion evidence.

| ID | Deliverable | Dependencies | Done | Evidence / PR |
| --- | --- | --- | --- | --- |
| M0 | Repository, scientific contract, data and workflow | None | [x] | [Bootstrap checks](validation/bootstrap.md); [Tests run 42825860916](https://github.com/landreman/alpha_analysis_3/actions/runs/42825860916) |
| M1 | Pitch-independent geometry and the per-pitch ledger | M0 | [x] | [Post-ADR evidence](validation/m1/adaptive.md); [PR #1](https://github.com/landreman/alpha_analysis_3/pull/1) |
| M2 | Transitions, cells and complete coverage | M1 | [ ] | [Draft PR #2](https://github.com/landreman/alpha_analysis_3/pull/2); [evidence and blocker](validation/m2/README.md); proposed [ADR 0002](adr/0002-fold-well-identity.md) |
| M3 | Reachability and the bracket: first complete coarse f | M2 | [ ] | |
| M4 | Outer integral, axis and convergence: accuracy acceptance | M3 | [ ] | |
| M5 | Speed and robustness within the runtime target | M4 | [ ] | |

**Next milestone: M2.** M1's 30 base-grid ledger cases meet the 6%/1% ratio
criteria for both sources. The same cases at declared production density meet
both ledger and transverse-refinement criteria. The coarse PCA 0.95 peaked
convergence miss remains recorded; production density resolves it. Accepted
[ADR 0001](adr/0001-long-well-cap-and-m1-ledger.md) is implemented: base chains
use 64 periods, only capped rows double through 128/256, and residual capped
wells retain unknown full J/K with a separate partial-K lower bound. All four
original capped nodes now exit. The [initial failed evidence](validation/m1/README.md)
is archived, not overwritten.

M2 is blocked on the fold-well identity decision in proposed ADR 0002. An
independent analytic counterexample and three n3are fold limits show that a
surviving well can have finite J/K when its representative maximum disappears.
Treating every fold endpoint as inactive is incorrect, and a narrow missed
max/min pair can duplicate a well without failing T=Gamma_max. All 30 candidate
cases and both sources are retained, but their apparently passing ratios are
provisional. The n3are sweep remains incomplete with failed geometry. Resolve
fold continuation and extremum certification, then revalidate the M1/M2
ledgers and complete the M2 census/cell gates. No M3 work has begun.

The denser n3are lambda_n=0.5 production case has two residual capped nodes at
256 periods. Both source ledgers remain incomplete and meet ADR 0001's
resolved-plus-lower-bound inequality; the other 58 production source checks
have finite ratios. No unknown full weight is replaced with that lower bound.
The partial lower bounds are 1.86% (uniform) and 5.53% (peaked) of that slice's
independent volume weight; they are not the full missing weight.

M2 must build bidirectional coverage, transitions and transition-cut cells.
The current one-direction partners, raw T/Gamma mismatches and frontier weights
are reported diagnostics, not complete coverage or connectivity. Auxiliary
population-quadrature nodes are not transition cells. Keep unknown weight explicit.
Scale the radial-partner threshold with node density when replacing the M1
fixed 0.4 rad diagnostic. The full `test_ledger_matrix` is an offline driver,
kept outside unit CI; the fast real-input ledger regression remains in CI.

M2/M4 must also explain ledger ratios above one, surface by surface. Missing
coverage alone cannot explain excess weight: PCA 0.95 peaked is 1.03344 at
production density, with a 1.66% transverse-change estimate. Separate angular
volume-reference error, radial interpolation near B_max(s)=b, ridge
reconnections/double counting, and cell coverage. The independent volume also
changes with the radial grid (n3are 0.5: about 301.5 to 307.1); angular ledger
closure is not evidence of radial convergence.

For M3/M4: the extra lambda_n=0.99 uniform-source tail probes in DMercFail and
d23p4 still miss 1% transverse convergence at 512 nodes and 16x subdivision.
Their residual capped weights are lower bounds, not certified totals. The new
0.5% of V_h budget applies from M3 onward and cannot be certified from these
lower bounds or the partial pitch window. Wildcard connectivity effects remain
unmeasured. Full f accuracy and the 300 s target remain unproven; some production
diagnostic sweeps exceed 300 s. No M2 implementation is included in PR #1.

Source: h(ρ) = (1 − ρ²)⁴ with h = 1 as a second column. Numerical error target
max(0.02 f, 0.01). Runtime target 300 s per equilibrium, memory under 2 GB.

The Claude review workflow needs the repository secret `CLAUDE_CODE_OAUTH_TOKEN` and
Claude GitHub App access; see [DEVELOPMENT.md](DEVELOPMENT.md).
