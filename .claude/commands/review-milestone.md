---
description: Adversarial scientific review of an alpha_analysis_3 milestone PR
argument-hint: "[PR number, or current branch]"
allowed-tools: Bash(git diff:*), Bash(git log:*), Bash(git status:*), Bash(gh pr:*), Bash(make*), Bash(pytest*), Bash(python*), Bash(black*), Read, Grep, Glob
---

Review $ARGUMENTS against AGENTS.md, docs/DESIGN.md, docs/STATUS.md and accepted ADRs.
Read the milestone criteria yourself; do not rely on the PR description. This repository
uses the ridge-surface formulation with numerical convergence estimates and an
unknown-weight bracket, not field certification. Source h(ρ) = (1 − ρ²)⁴ with h = 1
alongside; target error max(0.02 f, 0.01); runtime under 300 s and memory under 2 GB per
equilibrium. A correct unit-tested component is not evidence that f has converged or that
the required real cases work.

## Scope and completion

Check the selected M0–M5 milestone and its dependencies. Planning-only work does not
implement numerical milestones. An implementation PR may mark its own STATUS row complete
once its acceptance is demonstrated; it need not already be merged. No earlier repository's
gate applies.

## Scientific challenge

- Re-derive the measure and sign identities touched by the change: source at ρ not s,
  sign(G + ιI), dℓ = |C|/B |dζ|, the seam term 2πp − ιLq for closed ridge curves, the b⁻²
  pitch weight, |ds ∧ dα| = |∂_u θ_M − ι ∂_u ζ_M| ds du, the denominator.
- Check the ridge extraction: concave arcs only, winding numbers, periodic stitching, open
  arcs ending at folds, node density, and what happens where arcs reconnect between
  surfaces (unpartnered and uncovered nodes must become unknown weight).
- Check chains and wells: first chain maximum above b as the exit, hidden interior maxima,
  bounce points, the K_max cap as unknown and never as passing or confined, J and K
  quadrature with both endpoint singularities, long wells.
- Check transitions: the interior-barrier test, bisection, the per-surface T = Γ_max count
  as a reported check, additivity A_parent = A_child1 + A_child3, multi-jump edges, 1- and
  3-cut cells, and the well-assignment consistency check. Below-b folds are not transitions.
- Check reachability: links only on the shared edge's J range, whole-segment transition
  links including the child-3 cells at y, seeds at s = 1, the axis rule, dyadic bins that
  nest exactly, and that a Reeb edge cannot be traversed. Equal actions or overlapping
  ranges do not justify a link.
- Check that a small unknown cell can affect remote population, that the upper bound
  includes every wildcard incidence and the lower bound only justified paths, and that
  weights are not double counted.

## Tests that can fail

Identify the independent expectation behind each material claim. Spot-check the claimed
mutations and propose a missing test only for a credible wrong-physics path or an unmet
criterion. Avoid exhaustive coverage demands, tests of prose or implementation details, and
vacuous self-comparisons. A function ported from `prototype/` is not independent of the
prototype's evidence and must be disclosed. Never accept xfail or skip, looser tolerances,
or easier input subsets as the sole solution to a failure.

## Accuracy, feasibility and performance

Distinguish the N vs 2N estimate, the unknown-weight bracket, and any statistical interval.
Two agreeing grids or a tiny residual do not establish f accuracy; inspect the ledger ratio,
the count check and independent refinement directions. The 0.01 floor is an error magnitude,
not a full interval width or permission to return zero without a calculation.

Verify all five exact files and six global pitch levels in the milestone's evidence.
Distinguish early kernel probes from complete f acceptance. Check actual b, source, hashes,
all failures, cold setup, JIT and cache accounting, workers, memory, and 60/180/300 s
snapshots. A broad unknown result is not a runtime success. At M3 explicitly assess whether
the end-to-end strategy looks capable of the final goal; do not defer every difficulty.

Check fast under 120 s, full under 300 s, individual fast about 20 s and slow about 90 s.
Require a fast production regression for each scientific claim, but keep expensive full
matrices outside unit CI.

## Output

Give a short plain-language account of what the PR accomplishes and any important
obstacle. Return actionable findings with severity (blocking / should-fix / note),
location, problem and remedy; then a scoped **merge / fix first / needs ADR** verdict. If
there are no findings, say so. Do not invent criticism, describe failed CI as green, or
impose stronger mathematical certification than the agreed target. For infrastructure-only
PRs inspect usable workflows, contract consistency and setup evidence rather than demanding
unimplemented numerical capabilities.

Post the top-level review comment even with no findings. Add inline comments only for
confirmed code-specific issues. A merge recommendation does not itself authorize merging.

## Human-readable summary

In your review, include a human-readable summary of what this PR is about and for, and any
problems that arose with achieving the goals for the milestone. For this part of the review
text, the text should not sound like AI slop. Do not assume the human reader knows the
jargon - explain any concepts in simple terms that the reader needs to know. If the planned
design seems unlikely to succeed, in terms of computing the metric f for the set of boozmn
files in a reasonable wallclock time, such that we should re-think the overall design plan,
point this out and explain why.
