# Working in alpha_analysis_3

Read `docs/DESIGN.md`, then `docs/STATUS.md`. DESIGN governs the physics, the numerical
contract and the M0–M5 milestones. The original note in `docs/references/` defines the
metric; report any physics disagreement rather than silently choosing. `docs/adr/` records
new decisions. The earlier repositories `landreman/alpha_analysis` and
`landreman/alpha_analysis_2` are independent: do not edit them, inherit their work queues,
or require their gates. `docs/LESSONS.md` says what their records taught.

## Scope and autonomy

Use `.agents/skills/milestone/SKILL.md` for an explicit milestone or "next milestone"
request. Select the first incomplete milestone with satisfied dependencies. A planning or
bootstrap request does not start numerical milestones. Work through implementation,
verification, the required real-input evidence, STATUS, a draft PR, CI, readiness and the
configured Claude review without asking for routine confirmation. Never merge unless the
researcher explicitly asks. Preserve other work; use an isolated checkout or worktree when
needed.

Accepted and not reopened: the ridge-surface formulation; numerical convergence estimates
with an unknown-weight bracket instead of rigorous enclosures; source h(ρ) = (1 − ρ²)⁴ with
h = 1 alongside; error target max(0.02 f, 0.01); under 300 s and 2 GB per equilibrium;
NumPy, SciPy, Numba, Matplotlib and contourpy as base dependencies; multiple cores. A fast
answer with a wide bracket is not success.

## Environment and commands

On the researcher's Mac create a clean `.venv` from the existing interpreter:

```bash
/Users/mattland/opt/miniconda3/envs/20220806-03/bin/python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

Do not create another conda environment or use `--system-site-packages`. Do not use
`20250627-01-libE`; its Python 3.13 readline extension segfaults during pytest capture.
Other hosts and CI use clean Python 3.10–3.12 environments.

- `make test`: fast tier, under 120 s; individual fast tests about 20 s or less.
- `make test-full`: full tier, under 300 s; individual slow tests about 90 s or less.
- `make check`: Black plus the full tier; the readiness gate.
- `make format`: apply Black. `make smoke`: clean install and import check.

The runner enforces the tier totals. Reduce test cost before marking a test slow: lower
resolution, fewer surfaces, one file, shared fixtures with cached loaded fields. Every
scientific claim keeps a fast production-path test with an independent expectation.
Expensive 30-case matrices and five-file convergence tables run outside unit CI and leave
concise reproducible evidence under `docs/validation/`.

**Memory.** The laptop has 16 GB. Run whole-equilibrium computations one at a time, never
as a parallel batch of processes; give every long driver a resident-memory watchdog that
aborts instead of swapping; bound every data structure whose size is set by a
data-dependent ratio. A background job that fills swap costs the researcher hours.

## Code and scientific invariants

New code belongs in `src/alpha_analysis_3/` in the modules DESIGN §5 lists. Keep the core in
NumPy arrays and integer IDs; SciPy sparse connected components is the one graph primitive.
Numba kernels use `cache=True` and have NumPy or analytic references. Docstrings give
units, coordinates and the DESIGN section implemented. Format with Black.

Conserve b = B_bounce, never a radius-dependent normalized pitch. Use s = ρ² and evaluate
the source at ρ. Keep sign(G + ιI), dℓ = |C|/B |dζ|, the seam term 2πp − ιLq of closed
ridge curves, the b⁻² pitch weight, the interior-barrier transition test, the per-surface
T = Γ_max count, additivity A_parent = A_child1 + A_child3, edge-range bin links, whole-
segment transition links, dyadic bounded bins, and the lower/upper bracket.

Never turn a failed or capped well into passing, confined, zero action or zero weight; drop
unpartnered or uncovered nodes without accounting; treat a count mismatch as no
transition; link cells by proximity in J; allow ordinary action drift; or omit missing
transitions. Small unknown weight does not bound its global influence on connectivity.
Distinguish numerical estimates, the unknown-weight bracket and statistical intervals.

`prototype/` is a reference implementation from the design phase. Never import it from
`src/` or `tests/`; port functions only where they meet the standards above and say so.

## Required real inputs

`data/manifest.json` names five boozmn files (PCA, TURBO, DMercFail, d23p4, n3are) and
their hashes. All new numerical machinery must be exercised on all five, and fixed-pitch
machinery on λₙ = 0.05, 0.1, 0.5, 0.8, 0.9, 0.95 with

    b = B_min_global + λₙ (B_max_global − B_min_global),

extrema over all radii. Record the extrema used and the actual b. The 30 cases are
diagnostic probes that cover typical trapping classes, not the complete pitch integral.
Include the refinement and independent checks the milestone requires. Never replace hard
cases with easy ones or claim no trapped population because a finite search found none.

Record source, code and data identifiers, resolution controls, outcomes, error scope, failed
attempts, timings, worker counts, memory and active 60/180/300 s snapshots. Keep large
generated artifacts outside Git; keep concise reports and drivers.

## Definition of done and PRs

Complete the named DESIGN acceptance tests, the required real evidence, the one or two
load-bearing mutation checks, diagnostic plots for new geometry, and the local gates.
Update only the completed milestone in STATUS with necessary handoff notes. Push a
`codex/<milestone>-<slug>` branch and open a draft PR using the template: behavior,
acceptance evidence, mutations, runtime and error scope, limitations, stated plainly.

Every PR includes a human-readable summary of what the PR is about and for, and any
problems that arose, written for a reader who does not know the jargon; see the milestone
skill for the exact requirement.

Wait for GitHub `Tests` on the final revision. Then mark the PR ready to trigger Claude
Code Review. Address applicable findings; explain disagreements with evidence. If fixes
warrant another review, draft the PR, push, wait for green Tests, mark ready again. Wait
while a review is running; report a missing or failed review separately. Claude
availability is not a numerical dependency gate. Leave the finished PR ready and unmerged.
A documentation-only PR needs appropriate infrastructure checks, not invented physics
mutations.

## STOP conditions

For a new ambiguity that changes physics or acceptance, a gate that can only be met by
relaxing it, a new base dependency or core-boundary change, or a failure whose only
apparent fix hides uncertainty: write a proposed ADR, open or leave the PR draft, explain
the evidence and the decision needed, and stop the affected implementation. Already
approved decisions and ordinary engineering choices do not require this. Do not use xfail
or skip, loosen tolerances, or narrow scientific inputs to reach green. If the overall
plan looks unlikely to deliver accurate f within the runtime target, say so early and
explain why.
