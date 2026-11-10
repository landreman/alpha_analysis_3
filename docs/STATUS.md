# Milestone status

[DESIGN.md section 9](DESIGN.md#9-milestones) defines the sequence. Only this
repository's M0–M5 milestones apply; the milestones of `landreman/alpha_analysis` and
`landreman/alpha_analysis_2` are neither prerequisites nor completion evidence.

| ID | Deliverable | Dependencies | Done | Evidence / PR |
| --- | --- | --- | --- | --- |
| M0 | Repository, scientific contract, data and workflow | None | [x] | [Bootstrap checks](validation/bootstrap.md); [Tests run 42825860916](https://github.com/landreman/alpha_analysis_3/actions/runs/42825860916) |
| M1 | Pitch-independent geometry and the per-pitch ledger | M0 | [ ] | |
| M2 | Transitions, cells and complete coverage | M1 | [ ] | |
| M3 | Reachability and the bracket: first complete coarse f | M2 | [ ] | |
| M4 | Outer integral, axis and convergence: accuracy acceptance | M3 | [ ] | |
| M5 | Speed and robustness within the runtime target | M4 | [ ] | |

**Next milestone: M1.** Start from the design, not from the prototype's structure: the
prototype (`prototype/`) is evidence and a reference for conventions. Its surface tables,
ridge extraction, chains and well quadrature are sound and may be ported function by
function where they meet the code standards; its one-directional partner matching, its
transition links and its graph are known to be incomplete (DESIGN §4 marks the production
refinements). Keep the ledger and the transition count as reported checks from the first
milestone on.

M1 has a draft implementation and [five-file evidence](validation/m1/README.md), but
remains incomplete. Confirmed PCA wells first return after about 74 and 115 field
periods, beyond the planned 48–64-period cap; their capped weight is not replaced
with zero. The 30-case ledger also has a d23p4 ratio just outside its 1% gate.
[Proposed ADR 0001](adr/0001-long-well-cap-and-m1-ledger.md) records the required
decision and follow-up checks. Do not start M2 until M1's gates pass.

Source: h(ρ) = (1 − ρ²)⁴ with h = 1 as a second column. Numerical error target
max(0.02 f, 0.01). Runtime target 300 s per equilibrium, memory under 2 GB.

The Claude review workflow needs the repository secret `CLAUDE_CODE_OAUTH_TOKEN` and
Claude GitHub App access; see [DEVELOPMENT.md](DEVELOPMENT.md).
