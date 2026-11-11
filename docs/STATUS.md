# Milestone status

[DESIGN.md section 9](DESIGN.md#9-milestones) defines the sequence. Only this
repository's M0–M5 milestones apply; the milestones of `landreman/alpha_analysis` and
`landreman/alpha_analysis_2` are neither prerequisites nor completion evidence.

| ID | Deliverable | Dependencies | Done | Evidence / PR |
| --- | --- | --- | --- | --- |
| M0 | Repository, scientific contract, data and workflow | None | [x] | [Bootstrap checks](validation/bootstrap.md); [Tests run 42825860916](https://github.com/landreman/alpha_analysis_3/actions/runs/42825860916) |
| M1 | Pitch-independent geometry and the per-pitch ledger | M0 | [x] | [Post-ADR evidence](validation/m1/adaptive.md); [PR #1](https://github.com/landreman/alpha_analysis_3/pull/1) |
| M2 | Transitions, cells and complete coverage | M1 | [ ] | |
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

The denser n3are lambda_n=0.5 production case has two residual capped nodes at
256 periods. Both source ledgers remain incomplete and meet ADR 0001's
resolved-plus-lower-bound inequality; the other 58 production source checks
have finite ratios. No unknown full weight is replaced with that lower bound.

M2 must build bidirectional coverage, transitions and transition-cut cells.
The current one-direction partners, raw T/Gamma mismatches and frontier weights
are reported diagnostics, not complete coverage or connectivity. Auxiliary
population-quadrature nodes are not transition cells. Keep unknown weight explicit.

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
