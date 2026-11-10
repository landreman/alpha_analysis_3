# Lessons carried forward

Two earlier attempts exist: [landreman/alpha_analysis](https://github.com/landreman/alpha_analysis)
(a 3-D level-set mesh cut along companion curves, then a certified branch atlas) and
[landreman/alpha_analysis_2](https://github.com/landreman/alpha_analysis_2) (a numerical
first-return map, pursued in parallel). Their contracts do not govern this repository.
What their records teach:

- **Per-pitch geometry is the cost centre.** The legacy pipeline extracted the state space
  for every pitch as a level set of a 3-D mesh and cut it along curves found by backward
  tracing. Every reported failure traced to that representation: stepped-over contacts,
  zigzags that refinement relocates rather than removes, open endpoints, exhausted
  localization budgets, backward scans hitting a 1024-period cap. The final matrix resolved
  42 of 120 cases after 112 hours. Here all pitch-dependent structure is obtained by
  thresholding stored scalar fields on a pitch-independent surface.
- **Rigorous certificates do not scale.** The redesign bound itself to field-level
  enclosures; the L1 Fourier envelope made certified cells ~10⁻⁴ wide, so 0 of 240 coarse
  cells were certified, 0 of 17 queries classified, and after six milestones no Q(b)
  existed. Accuracy here comes from convergence and an explicit unknown-weight bracket.
- **End to end first.** Neither attempt computed f before its fifth milestone. The prototype
  behind this design computed f for all five files in its first day; every milestone here
  produces the five-file table.
- **The physics that was right stays.** b (not λₙ) as the conserved variable, undirected
  split/merge with A_W = A₁ + A₃, pitch-independent extremum heights to locate transitions,
  shared forward scans, an independent population ledger, the sign of G + ιI, the source at
  ρ, the b⁻² pitch weight. All of these are testable without certificates.
- **Reviewed corrections that remain binding.** Ordinary Reeb connectedness is not
  constant-action accessibility; small event measure does not bound global influence; a
  long-return cap is not passing; a surface-maximum trapping identity needs dense-line
  hypotheses.
- **Load-bearing mutations were missed.** Mutation testing found the easy mutations and
  a review bot found the ones that mattered. DESIGN §6 names the mutations each milestone
  must demonstrate.
- **Bound every data structure whose size is a data-dependent ratio.** The prototype's
  J-bin graph took its bin width from a per-sheet median that was once at rounding level;
  one cell asked for 5×10¹⁶ bins, eleven parallel runs filled 35 GB of swap and froze the
  researcher's 16 GB laptop. The graph is now dyadic and capped, every driver has a
  resident-memory watchdog, and whole-equilibrium runs are sequential.
- **Simple transition tests fail quietly.** A bounce-point jump threshold and an
  exit-identity change both produced false transitions that only the per-surface count
  check exposed. Keep the count check reported in every run.
