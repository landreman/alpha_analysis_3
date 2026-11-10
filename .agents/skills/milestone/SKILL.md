---
name: milestone
description: Implement one alpha_analysis_3 milestone from docs/DESIGN.md section 9 through real-input verification, STATUS, a ready PR and its Claude review. Use for an explicit milestone or next-milestone request; planning-only work does not start implementation.
---

# Implement one milestone

First make sure the latest version of the "main" branch has been pulled from github.
Read root `AGENTS.md`, `docs/STATUS.md`, the selected milestone in `docs/DESIGN.md`
section 9, and `.claude/commands/review-milestone.md`. Paths are repository-relative.
For "next", select the first incomplete **M0–M5** milestone whose dependencies are
complete. Do not import the earlier repositories' work queues. An explicit request for a
named milestone takes precedence when its dependencies are satisfied.

## Baseline and scope

Verify dependency evidence and green GitHub Tests on the baseline. Inspect `git status`
before branching; preserve other work and use an isolated checkout or worktree when
necessary. Default branch: `codex/<milestone>-<slug>`. Do not switch or pull over
uncommitted work. Implement one milestone unless the researcher asks for more.

The ridge-surface formulation, the convergence-plus-bracket error estimate, the source
h(ρ) = (1 − ρ²)⁴, the tolerance max(0.02 f, 0.01), the 300 s and 2 GB targets, and the base
dependencies are already authorized. This does not permit missing transitions, action
drift, unaccounted population, or a capped chain classified as passing or confined.

The prototype in `prototype/` is a reference. Port a function only when it meets the code
standards in AGENTS, say so in the PR, and never import `prototype/` from `src/` or `tests/`.

## Implement and verify

Derive a small set of discriminating tests from the milestone's named acceptance tests,
with analytic or independently computed expectations. Confirm the missing behavior fails
for the intended reason, implement it, and keep the fast and full test budgets. Run
`make check`; use `make smoke` for packaging or dependency changes. Apply the one or two
load-bearing physics mutations that DESIGN section 6 lists for the milestone, observe the
intended failures, and revert. Documentation-only changes need no fabricated mutation.

Exercise the milestone's machinery on the manifest's five files and the six global pitch
levels, with the refinement and independent checks the milestone requires. Keep full
matrices outside unit CI. Preserve all failed cases and record provenance, source, actual b,
resolution controls, unknown weight, timings, memory and the 60/180/300 s snapshots. Run
whole-equilibrium computations one at a time; the laptop has 16 GB. An all-unknown output,
a wide bracket, or a small local residual does not meet an acceptance gate. M3 is
deliberately the first complete coarse f: raise a feasibility problem there if coverage
or scaling remains poor.

Update the completed STATUS row in the same PR and add concise handoff notes. Write new
scientific decisions in DESIGN or an ADR, not in status history.

## PR, CI and adversarial review

A milestone implementation request authorizes pushing its branch, creating or updating its
draft PR, marking it ready after checks, and requesting and responding to the configured
Claude review. Do not ask again for those routine steps. It does not authorize merging or
starting another milestone.

In the PR text, lead with the concrete behavior and purpose, then map each
acceptance criterion to its named test/evidence. Include errors, unresolved scope,
timings, reverted mutations, and deviations. Explain problems in plain language;
if the algorithm looks unlikely to meet the f/runtime target, say so explicitly.

In the PR, include an explicit section with a human-readable summary of what this PR is about and for, and any problems that arose with achieving the goals for the milestone. For this part of the PR text, the text should not sound like AI slop. Do not assume the human reader knows the jargon - explain any concepts in simple terms that the reader needs to know. If the planned design seems unlikely to succeed, in terms of computing the metric f for the set of boozmn files in a reasonable wallclock time, such that we should re-think the overall design plan, point this out and explain why.

Wait for GitHub `Tests` on the final code revision and fix failures. When the milestone
gates pass, run `gh pr ready <number>` to trigger Claude Code Review; do not wait for
another researcher message. Read the resulting adversarial findings and address applicable
ones. Disagreements need evidence.

For a further review after substantive fixes, convert to draft, push, wait for Tests on
that revision, and mark ready again. Do not toggle readiness merely to repeat an unchanged
review. If Claude is running, keep checking until it completes or fails; do not return
immediately after requesting it. Missing credentials, workflow failures or an absent review
are reported separately and do not masquerade as a favorable review or invalidate unrelated
Tests. See `docs/DEVELOPMENT.md`.

Stop with the PR ready and **unmerged**. Summarize the result, validation, material limits,
PR URL and next milestone. Never merge without explicit instruction.

## New scientific blockers

Follow AGENTS' STOP process for a new physics or acceptance ambiguity, an unmet gate that
requires relaxation, a new base dependency or core boundary, or hidden uncertainty: record a
proposed ADR and leave the PR draft for a concrete researcher decision. Do not silently
choose a new metric, tolerance or input subset. Routine fixes and the decisions already
accepted in DESIGN do not require another approval.
