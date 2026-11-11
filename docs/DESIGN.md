# Design: ridge-surface computation of the alpha-loss surrogate f

Revision 2, 2026-11-11. This document is the scientific and numerical contract for
`alpha_analysis_3`. It supersedes the proposal that was iterated in chat on 2026-11-09
and carries its approved decisions. The original note defining the metric is
[docs/references/20261109-01 alpha loss fraction metric.tex](<references/20261109-01 alpha loss fraction metric.tex>).
Where this document and the note disagree about the physics, report the disagreement;
do not silently choose.

The algorithm was prototyped end to end before this repository existed
(`prototype/`, evidence in [validation/prototype](validation/prototype/README.md)).
The prototype is a reference for conventions and a source of measured numbers; the
production code is written fresh in `src/alpha_analysis_3/` against this document.

## 1. Objective and accepted decisions

Compute f, the fraction of alpha-particle births whose trapped state is connected to the
plasma boundary along contours of constant second adiabatic invariant J, including
discontinuous changes of J at transitions between trapping classes, for one stellarator
equilibrium given as a `boozmn` file.

The researcher has fixed the following. They are not reopened by implementation PRs.

| Decision | Value |
|---|---|
| Numerical accuracy target | δf ≲ max(0.02 f, 0.01), an error magnitude, not a full interval width |
| Error estimate | convergence (resolution N vs 2N) plus the unknown-weight bracket; no rigorous field enclosures |
| Runtime target | under 300 s per equilibrium on the researcher's laptop, including loading, setup and the N vs 2N check |
| Memory | under 2 GB resident per process; never a parallel batch of whole-equilibrium runs on the 16 GB laptop |
| Source profile | h(ρ) = (1 − ρ²)⁴ with s = ρ²; h = 1 kept as a second column in every report |
| Pitch variable | b = B_bounce conserved along a contour; λₙ only labels diagnostic levels |
| Diagnostic levels | λₙ ∈ {0.05, 0.1, 0.5, 0.8, 0.9, 0.95}, b = B_min + λₙ (B_max − B_min) with radially global extrema |
| Transitions | undirected; all branches allowed; J changes discontinuously at a split or merge |
| Boundary births | Θ = 1 for a state born at s = 1 |
| Passing orbits | never a bridge between trapped states (§2.4) |
| Inner cutoff | s_min = 10⁻³ with a graded grid; sensitivity in s_min is part of acceptance |
| Long wells | chains are capped at K_max periods; capped weight is reported as unknown, never as passing or confined |
| Dependencies | NumPy, SciPy, Numba (base), Matplotlib and contourpy (base, diagnostics and ridge extraction); multiple cores allowed |
| Reuse | `boozer_field.py` from `landreman/alpha_analysis` is the production field reader |

## 2. Physics and conventions

### 2.1 Coordinates and the field line

Boozer coordinates (s, θ, ζ), s ∈ (0, 1], field period L = 2π/N_fp, field-line label
α = θ − ιζ. Along a field line θ = θ₀ + ιζ. The along-line derivative is
D∥ = ι ∂_θ + ∂_ζ; the physical direction of motion along the line is sign(G + ιI), with
C = G + ιI. Arc length is dℓ = |C| / B · |dζ|. The field-period angle φ = N_fp ζ is used
for tables.

### 2.2 Pitch, wells, J and K

The bounce field b = B_bounce is conserved along a contour. A trapped well at pitch b is a
segment of a field line between two consecutive points where B = b with B < b in between.
Its half-bounce action length and time length are

J = ∫ √(1 − B/b) dℓ,   K = ∫ dℓ / √(1 − B/b),

both over the well. Both endpoint singularities are integrable.

### 2.3 The metric

With h the source profile, V_h the source-weighted volume and Θ ∈ {0, 1} the connectivity
indicator of a trapped state,

f = [ ∫ db/b² ∫_{Λ_b} h(√s) K Θ |ds dα| ] / [ 2 ∫ ds h(√s) |C| ∫∫ dθ dζ / B² ],

where Λ_b parametrizes the trapped states at pitch b (§3) and |ds dα| is the induced
measure. The denominator is the total birth weight; the numerator counts births whose
state is connected to s = 1. Equivalent forms: f = (1/2V_h) ∫ Q(b) b⁻² db with Q(b) the
connected trapped weight at pitch b.

### 2.4 Rules that are easy to get wrong

* **Transitions are undirected.** A contour may enter a split or merge from any branch
  and leave on any other, with the actions related by A_parent = A_child1 + A_child3.
* **Θ = 1 at s = 1.** Cells touching the outermost surface are seeds of connectivity.
* **Passing is never a bridge.** At fixed b a well can disappear by shrinking to zero length
  (J → 0, no finite-J contour reaches it), by an interior barrier reaching b (a transition,
  handled), or by every bounding maximum dropping below b (the orbit becomes passing). A
  passing particle may re-trap on a neighbouring surface with a different J; the note's rule
  Θ = 0 for never-bouncing particles excludes this, and in the constant-J formulation it is
  automatic because J diverges as the bounding maxima approach b. The only place the code
  meets that boundary is the chain cap, whose weight is reported as unknown.
* **Evaluate h at ρ = √s**, never at s.
* **The b⁻² weight** belongs to the outer pitch integral and must not be dropped or applied
  twice.
* **Global extrema** define λₙ; a radius-dependent normalization is never used inside a
  contour trace.

## 3. The state space: the ridge surface Λ

### 3.1 Ridge network

On each flux surface, Λ_s = {(θ, ζ) : D∥B = 0, D∥²B < 0} is the set of local maxima of B
along field lines: the zero contour of D∥B restricted to its concave arcs. Field lines cross
it transversally, from D∥B > 0 to D∥B < 0. Closed components are global cross-sections of
the field-line flow and carry a winding number (p, q) on the (θ, φ) torus; other components
are open arcs ending at folds (D∥²B = 0) where a maximum–minimum pair annihilates. In the
five reference equilibria the network is one to three closed curves plus a few arcs
(W7-X, DMercFail, TURBO) or a set of arcs that reconnect as s varies (PCA, n3are).
Λ = ∪_s Λ_s is a two-dimensional surface, **independent of pitch**.

### 3.2 Chains

From x ∈ Λ_s follow the field line forward: the successive crossings of Λ_s are the
along-line maxima x = x₀, x₁, x₂, … with heights H_k and intervening minima m_k.
[Accepted ADR 0001](adr/0001-long-well-cap-and-m1-ledger.md) replaces the prototype's
48–64-period assumption: production base chains use 64 periods; only capped nodes
are extended by doubling to a hard ceiling of 256. The measured PCA returns at
about 74 and 115 periods show why the prototype's longest observed 19-period well
was not a bound. The ceiling remains a cost guard. A node with no exit stays
unknown with NaN full J and K, and a separately reported partial-K lower bound.

### 3.3 Wells at pitch b

The trapped wells at pitch b correspond one-to-one to Λ_b = {x ∈ Λ : H₀(x) > b > m₀(x)}:
from such a maximum B falls through b at the incoming bounce point, and the well ends just
before x_n with n(x; b) = min{k ≥ 1 : H_k > b} the exit. Long wells, several wells per line
and the twisted α-periodicity are absorbed: a well is a point of Λ whatever its length.
Boundaries of Λ_b are H₀ = b (Γ_max, marginal maxima) and m₀ = b (Γ_min, zero-length wells).

### 3.4 Transitions

J jumps on T_b = {x ∈ Λ_b : H_k(x) = b for some 0 < k < n}: an interior barrier touches b
and the exit jumps from x_n back to x_k. Parent [x₋, x₊], child 1 [x₋, before x_k] and
child 3 [after x_k, x₊]: parent and child 1 are the *same point* x on the two banks of T_b;
child 3 is y = x_k ∈ Γ_max. Consequently

* the transition set at pitch b is the level curve {H = b} on Λ, and T_b is its preimage
  under the return map, found forward as a jump of the exit between neighbouring states;
* the gluing is pointwise and exact (same field line, same s), with no event matching;
* per surface, the number of T points equals the number of Γ_max crossings. This count is a
  built-in detector: a mismatch means two crossings fell inside one node spacing and the
  edge must be subdivided.

Below-b folds change the chain index but not the exit, so J is continuous and no transition
is created. Equal-height interior maxima on stellarator-symmetric lines give a crossing of
two T curves and show up as cells with one or three cut edges.

### 3.5 Measure

With u the parameter along a ridge curve and (θ_M, ζ_M)(u) the ridge point,
|ds ∧ dα| = |∂_u θ_M − ι ∂_u ζ_M| ds du, nonzero by transversality. For a closed (p, q)
curve the total |Δα| around the curve includes the seam term 2πp − ιLq.

## 4. Algorithm

Production refinements beyond the prototype are marked ►.

### 4.1 Pitch-independent stage (once per equilibrium)

1. **Surfaces.** s_i on a grid (base 24, production 48–96) with extra surfaces graded toward
   the axis down to s_min; ► adaptive insertion where the ridge census changes or where
   unpartnered weight is large.
2. **Tables.** Fourier coefficients → N_θ × N_φ tables (512² base) of B, D∥B, D∥²B by inverse
   FFT, accurate to roundoff at grid points. Values along field lines by periodic cubic
   interpolation with precomputed spline coefficients. ► Newton-polish extrema and bounce
   points on the exact series where the table error matters. Derivative tables are freed
   after ridge extraction.
3. **Ridge network.** Marching squares on D∥B over the wrap-padded grid, stitched across the
   periodic seams, split into concave arcs by the sign of D∥²B. Closed curves get nodes
   uniform in arc length; open arcs get the same linear density (base 128 nodes per 2π,
   production 256–512). Record θ_M, ζ_M, H, α and |Δα| per edge.
4. **Chains.** From every node, 128 samples per period to K_max periods, local extrema with
   derivative sign brackets and polishing. Vectorized or compiled over all nodes of a
   surface. Under ADR 0001, start base nodes at 64 periods and selectively double capped
   nodes to 128 and 256. Pitch-specific auxiliary quadrature points may stop as soon as
   their first exit is resolved; otherwise they double through the same finite ceiling.
5. **Radial partners.** ► Match nodes between adjacent surfaces in both directions (torus
   distance, threshold ~0.4 rad at base density, scaled with density). Nodes without a
   partner in either direction are frontier nodes and carry unknown weight; a top node not
   reached by any cell is an uncovered node and also carries unknown weight. Large
   unpartnered or uncovered weight triggers surface insertion.

### 4.2 Per-pitch stage

1. **Active states and exits.** H₀ > b > m₀; n = first chain maximum above b; a chain that
   reaches K_max without an exit is capped and its state is unknown.
2. **J and K.** Bounce points by bisection on the table (► polished); Gauss–Legendre with
   the sine substitution over the whole well (384 nodes base; ► composite per extremum
   segment), which removes both endpoint singularities.
   **Long-well population (ADR 0001).** Subdivide ridge edges adjacent to wells longer
   than 64 periods dyadically, retaining the original samples, up to 16× the base density.
   Compare successive strip-weight estimates against the slice tolerance and report the
   level and change. Activity-boundary edges are also refined for population quadrature;
   this is not a transition classification. A residual capped node contributes its
   observed partial-K lower bound, while its full weight remains unknown. Never fill that
   weight from the independent volume ledger.
   Add magnitudes of edge-weight changes before radial integration and source weighting;
   compare that full-slice change with its resolved-plus-partial weight for each source.
   Retain local ceiling failures separately. A stricter local quadrature control may be
   used without changing the slice acceptance tolerance. This estimates numerical change
   in known contributions, not an upper bound on the full weight of capped states.
   Composite trapezoid and composite Simpson rules use the same dyadic samples;
   both have positive weights and retain the full seam measure. Record the rule used.
3. **Transitions.** An edge between two active states is cut when the exit maximum of the
   shorter well appears as an interior maximum of the longer well (within a fraction of the
   local spacing between maxima). Two simpler tests were tried and rejected by the count
   check: a bounce-point jump threshold fires falsely where a shallow exit maximum makes the
   bounce point slide fast, and an exit-identity change fires falsely where a fold above b
   swallows the exit maximum after the well behind it has shrunk to zero. ► Bisect cut edges
   to locate T and the two bank actions; ► subdivide edges whose chain index changes by two
   or more without a transition. The partner y is the highest interior maximum of the
   long-side well, known exactly from the chain. ► Well-assignment check: adjacent nodes
   whose J differ by more than a few times the edge's expected variation without a detected
   transition are cut or unknown, never regular (the prototype found such a cell on n3are
   with J from 10⁻¹⁹ to 19).
4. **Cells.** For each edge (a_k, a_{k+1}) on surface i and its partners on surface i+1, a quad
   (or a triangle if the partners coincide). Weight = Δs · ½(|Δα|_bottom + |Δα|_top) · mean
   K with ► piecewise-linear clipping of J and interior K nodes. Cells with 0 cut edges are
   regular; with 2 cut edges they split into two sub-cells, one per bank; with 1 or 3 they
   are ► subdivided until each piece has 0 or 2 cuts; the residual is unknown. Unknown also:
   cells with a capped node, frontier edges, uncovered top nodes.
5. **Sheets, bins, reachability.** Sheets are connected components of cells sharing an uncut
   edge. Per sheet a bin width δJ from a robust scale (► median cell J-range weighted by
   measure, floored at 10⁻⁶ of the sheet's J span) divided by the bins-per-cell factor
   (4 base, 8 for the convergence check). Graph nodes are (cell, bin). Two adjacent cells are
   linked in bin k only if the *shared edge's* J range meets bin k; linking on the cells' own
   ranges over-connects across saddles and made f depend on δJ by 10–15 %. Transition links
   join the whole segment of T inside each cut cell: for t along the segment,
   (long bank, bin J_long(t)) ↔ (short bank, bin J_short(t)) ↔ (cells at y, bin J₃(t)) with
   J₃ = J_long − J_short. Seeds are all bins of cells touching s = 1. Axis rule: cells
   touching the innermost surface are linked where their inner-edge J ranges intersect. The
   reachable set is the connected component of the seeds. Motion along a Reeb edge is
   impossible by construction.
   **Bounded graph.** Bins are dyadic and multi-resolution: a cell whose range would span
   more than 128 sheet bins is binned at width 2^k δJ with the smallest k that fits. Coarse
   bins are exact unions of fine bins, so a link between cells of different levels mixes
   contours only inside the coarse bin; the J interval never diffuses along a path. The
   graph size is therefore bounded by 128 × (number of cells).
6. **Bracket.** Lower: unknown cells removed. Upper: each unknown cell is a wildcard joining
   every bin of every cell it touches, and its own weight counts as reachable. The point
   estimate is reported with the bracket; the bracket width is part of the error estimate.

### 4.3 Outer stage

f = (1/2V_h) ∫ Q(b) b⁻² db over b ∈ [B_min, B_max]. Q(b) is piecewise smooth with jumps
where connectivity flips; for the peaked source most of f comes from λₙ ≳ 0.8. Base:
midpoint rule with 24–48 slices; ► adaptive quadrature that locates the jumps. Slices are
independent and run in parallel. Convergence: repeat with doubled N_s, node density and
slice count and halved δJ; the reported error is |f_N − f_2N| plus the bracket width.

## 5. Numerical architecture and package layout

```
src/alpha_analysis_3/
  boozer_field.py   Boozer reader and exact Fourier evaluation (ported, see provenance)
  conventions.py    tolerances, λₙ levels, source profile, direction sign
  tables.py         per-surface FFT tables and periodic interpolation            (M1)
  ridge.py          ridge arcs, winding numbers, node placement, |Δα|            (M1)
  chains.py         extrema chains along field lines                            (M1)
  wells.py          active states, exits, bounce points, J and K                (M1)
  ledger.py         source-weighted volume, trapped-weight volume integral      (M1)
  transitions.py    interior-barrier test, bisection, counts, well assignment   (M2)
  cells.py          partners in both directions, cells, subdivision, coverage   (M2)
  graph.py          sheets, dyadic bins, links, seeds, reachability, bracket    (M3)
  weights.py        cell weights, clipping, interior K nodes, sources           (M3)
  pitch_slice.py    Q(b) with lower/upper/total/unknown for one pitch           (M3)
  outer.py          pitch quadrature, N vs 2N, error report                     (M4)
  kernels.py        Numba kernels with NumPy references                         (M5)
  diagnostics.py    ridge, chain and Θ-map plots                                (M1+)
  cli.py            compute f for one boozmn file                               (M5)
```

* The numerical core uses plain NumPy arrays and integer IDs. No graph, mesh or
  visualization framework objects inside it; SciPy sparse connected components is the
  one allowed graph primitive.
* Numba is a base dependency. Kernels use `cache=True`, are deterministic, and each one
  has a test against a NumPy reference or an analytic value.
* Every long-running driver checks resident memory and aborts above the cap rather than
  swapping. Whole-equilibrium runs on the laptop are sequential.
* Docstrings state units, coordinate conventions and the section of this document they
  implement.

## 6. Built-in checks

These run in every full computation and are reported, not only tested.

1. **Ledger.** The sum of all cell weights at pitch b, including unknown cells, equals the
   independent volume integral ∫ h |C| / (B √(1 − B/b)) [B < b] [B_max(s) > b] dV to within
   the stated tolerance (1 % at production resolution; §9 says what each milestone must
   reach). It catches the seam term, the measure, missing coverage and double counting.
   For residual capped nodes, ADR 0001 requires an incomplete ratio with the resolved
   sum, partial-weight lower bound and diagnostic gap, and checks
   resolved + lower bound ≤ reference × (1 + tolerance). The gap is never filled in.
   From M3 onward, capped birth weight has its own 0.5% of V_h budget, separate from
   frontier/uncovered weight; exceeding it requires a new ADR, not a higher ceiling.
   A partial lower bound alone cannot establish compliance with an upper weight budget.
2. **Transition count.** On every surface the number of T points equals the number of Γ_max
   crossings after subdivision; residual mismatches become unknown weight.
3. **Additivity.** At every located transition J_parent = J_child1 + J_child3 within the
   quadrature tolerance.
4. **Convergence.** N vs 2N in surfaces, node density, slices and δJ.
5. **Independent reachability check.** A birth-space quasi-Monte-Carlo sample classified
   through the same wells but independent weights (M4).

Load-bearing mutations every PR of the relevant milestone must demonstrate at least two of:
flip sign(G + ιI); drop the seam term; drop the b⁻² weight; evaluate h at s; link bins on
cell ranges instead of edge ranges; drop the transition links; treat a capped chain as
passing; skip the interior-barrier test.

## 7. Error accounting, acceptance and runtime

* **Point estimate** f at the production resolution, peaked source, with h = 1 alongside.
* **Error estimate** δf = |f_N − f_2N| + (f_upper − f_lower)/2 at resolution N. Acceptance
  requires δf ≤ max(0.02 f, 0.01) for each of the five files and the peaked source, and the
  same at N = 2N for the h = 1 source as a secondary check.
* **30-case matrix.** Five files × six λₙ levels are diagnostic probes: for each slice report
  the ledger ratio, T and Γ_max counts, unknown weight fraction, reachable fractions for both
  sources and the slice wall time. They are not six independent equilibria and do not replace
  the pitch integral.
* **Runtime.** Under 300 s per equilibrium including loading, tables, the N and 2N
  computations and JIT compilation on a cold cache; record 60-, 180- and 300-second snapshots
  for runs that are still active. A result whose bracket is wider than the tolerance when the
  budget expires is an accuracy failure, not a success.
* **Memory.** Under 2 GB resident.

## 8. Test budgets, evidence and prohibited shortcuts

* `make test` (fast tier) under 120 s of test-process wall time; a single fast test about
  20 s or less. `make test-full` under 300 s; a single slow test about 90 s or less. The
  runner enforces the totals.
* Every scientific claim keeps at least one fast test on the production path that fails under
  a mutation of the physics it checks. Reduce cost before marking a test slow: lower
  resolution, fewer surfaces, one file, shared fixtures with cached loaded fields.
* Full 30-case matrices and five-file convergence tables run outside unit CI with a committed
  driver and a concise report under `docs/validation/`.
* Never: replace a failed or capped well with zero weight, with passing, or with confined;
  drop unpartnered or uncovered nodes without accounting; treat a count mismatch as no
  transition; link cells by proximity in J; cap a long chain and assign Θ = 0; relax a
  tolerance, mark a test `xfail` or `skip`, or narrow the inputs to reach green. Each of these
  produces a plausible number, which is what makes it dangerous.

## 9. Milestones

Select the first incomplete milestone whose dependencies are complete in
[STATUS.md](STATUS.md). Every milestone from M1 on runs its machinery on all five files at
the six levels and ends with at least one diagnostic figure. Named tests are planned until
implemented; keep the names unless a better one is recorded in STATUS.

### M0 — Repository, scientific contract, data and workflow

Dependencies: none. This bootstrap: design, agent instructions, milestone skill, review
command and workflow, Tests workflow, original note, hashed five-file data, package skeleton
with conventions, isolated installation. Tests: `test_five_real_inputs_and_fourier_conventions`,
`test_reference_documents_unchanged`, `test_error_target_and_pitch_conventions`.

### M1 — Pitch-independent geometry and the per-pitch ledger

Dependencies: M0. Tables, ridge network with winding numbers, nodes and |Δα| with the seam
term, chains, one-directional partners (the full bidirectional matching is M2), active states
and exits, J and K, the source-weighted volume and the trapped-weight volume integral, the
ridge-network and chain plots, and a five-file ledger table at the six levels.

Named acceptance tests:

- `test_tables_match_exact_series`: B, D∥B, D∥²B tables agree with direct Fourier sums at
  grid points to roundoff and interpolated values agree to 10⁻⁷ relative at 512².
- `test_action_and_time_of_cosine_well`: J and K for a cosine well agree with the elliptic
  integral closed form to 10⁻⁸, for both signs of C.
- `test_axisymmetric_ridge_is_one_closed_curve`: an n = 0 field gives one closed ridge
  curve per surface with winding (0, 1) and constant H.
- `test_well_exit_and_cap_statuses`: multi-period wells, a hidden interior maximum above b,
  and a capped chain receive the right exit index or the capped status; no capped chain is
  reported passing.
- `test_adaptive_extension_only_retraces_capped_nodes` and
  `test_ceiling_keeps_unknown_with_partial_weight_lower_bound`: bounded extension and
  explicit lower-bound accounting without changing unknown status at the ceiling.
- `test_long_strip_refinement_preserves_measure_and_reports_level`: dyadic transverse
  refinement retains the closed-ridge seam measure and reports its convergence level.
- `test_ledger_closes` (fast, one file, two levels, reduced resolution, within 2 %) and
  `test_ledger_matrix` (slow, all 30 cases at base resolution, within 6 % for every case and
  within 1 % for DMercFail and d23p4). The 6 % allowance is the known coverage gap that M2
  closes; the table of ratios is committed.

Mutations: drop the seam term; flip sign(C); omit 1/B in dℓ.

### M2 — Transitions, cells and complete coverage

Dependencies: M1. Interior-barrier transition test; edge bisection for T and the bank
actions; subdivision of multi-jump edges and of 1-/3-cut cells; the per-surface count
check; the well-assignment check; bidirectional partner matching; repair of uncovered nodes
by subdivision and adaptive surface insertion; additivity at every T.

Named acceptance tests:

- `test_synthetic_split_merge`: a two-well model field with a controllable interior barrier;
  T is located within the bisection tolerance, A_parent = A_child1 + A_child3 to the
  quadrature tolerance, the count matches on every surface.
- `test_transition_count_matches_gamma_max`: all 30 cases at base resolution, counts match on
  every surface after subdivision; any residual is unknown weight and is listed.
- `test_coverage_closes_ledger`: all 30 cases within 1 % and unpartnered plus uncovered
  weight below 0.5 % on every file, including TURBO, PCA and n3are.
- `test_dmercfail_generic_event`: the DMercFail λₙ = 0.8 event of the old reports
  (`data/legacy_probes.json`) is found at the recorded pitch within a node spacing.
- `test_well_assignment_catches_inconsistent_edges`: an edge whose nodes carry wells of
  incompatible length without a detected transition becomes a cut or unknown cell.

Mutations: skip the interior-barrier test; skip bisection (additivity fails); match partners
in one direction only (coverage test fails on n3are).

### M3 — Reachability and the bracket: first complete coarse f

Dependencies: M2. Sheets, dyadic bins, edge-range links, transition segment links, axis
rule, seeds, connected components, lower and upper bounds, cell weights with clipping and
interior K nodes, both sources, Q(b) per slice, Θ maps, and a midpoint pitch integral giving
f for all five files at base resolution.

Named acceptance tests:

- `test_axisymmetric_field_has_no_losses`: f equals the weight of the outermost cell layer
  only and halves when the surface count doubles.
- `test_closed_contour_annulus_is_unreachable`: a model field whose contours inside an
  annulus are closed gives Θ = 0 there while the outside is lost.
- `test_unknown_connector_moves_only_the_upper_bound`: marking one small cell unknown leaves
  the lower bound unchanged and raises the upper bound by at most its reachable neighbourhood.
- `test_bin_halving_converges`: f at 4 and 8 bins per cell differ by less than 0.01 on
  DMercFail at reduced resolution.
- `test_coarse_f_five_files` (slow or offline evidence): Q(b) at the six levels and f for all
  five files with lower, upper, total and unknown weight. Feasibility gate: the bracket width
  for h = 1 is at most 0.02 on every file at base resolution. If a file cannot meet this, stop
  and write an ADR rather than relabel the width.

Mutations: link bins on cell ranges; drop the transition links; drop b⁻²; treat a capped
chain as passing.

### M4 — Outer integral, axis and convergence: accuracy acceptance

Dependencies: M3. Adaptive pitch quadrature that locates jumps of Q(b); graded axis grid to
s_min = 10⁻³ with a sensitivity study; N vs 2N tables for all five files and both sources;
source reweighting that reuses geometry; the birth-space quasi-Monte-Carlo check; a
restartable results record with provenance.

Named acceptance tests:

- `test_pitch_quadrature_handles_jumps`: a synthetic Q(b) with a jump integrates to 10⁻³.
- `test_axis_grading_sensitivity`: the peaked-source f of d23p4 changes by less than 0.005
  between s_min = 10⁻³ and 10⁻⁴ with the graded grid.
- `test_source_reweighting_reuses_geometry`: h = 1 and the peaked source from one geometry;
  the two denominators match independent quadrature.
- `test_birth_sampler_agrees_with_deterministic_f`: the quasi-Monte-Carlo estimate agrees
  with f within its statistical plus numerical scope on DMercFail.
- `test_five_file_accuracy_evidence` (offline evidence, committed table): for each file and
  the peaked source δf ≤ max(0.02 f, 0.01) as defined in §7, with every failure preserved.

Mutations: drop b⁻² in the outer integral; evaluate h at s.

### M5 — Speed and robustness within the runtime target

Dependencies: M4. Numba kernels for chains, wells and quadrature with NumPy references;
parallel slices; cold-start timing table for all five files at the resolution that met M4,
including the N vs 2N check, under 300 s and 2 GB; a fixed "optimization mode" resolution;
a command-line entry point; stress tests on any additional equilibria available.

Named acceptance tests:

- `test_kernels_match_numpy_reference`: every kernel against its reference on real data.
- `test_cli_end_to_end_small`: the command line computes f for one file at reduced
  resolution and writes the report.
- `test_runtime_budget_evidence` (offline evidence, committed table): five cold-start runs
  with 60/180/300 s snapshots, peak memory and the achieved δf.

## 10. Definition of done and decision process

A milestone is done when all of these hold:

1. Its dependencies are complete and `make check` is green inside the test budgets.
2. Each acceptance criterion in §9 is demonstrated by its named test or committed evidence,
   not by a plot alone and not by a residual that got small.
3. The one or two load-bearing mutations were applied, the suite went red, and they were
   reverted; the PR body names them.
4. Every new geometric object has a diagnostic plot.
5. GitHub `Tests` is green on the final revision; the Claude review ran or its absence is
   reported separately.
6. `STATUS.md` has the row marked and the notes the next milestone needs.
7. Every deviation from this document is an accepted ADR or is written up in the PR body.

STOP and write a proposed ADR in `docs/adr/`, leave the PR draft, and end the turn when:
this document is ambiguous or looks wrong about something that changes the physics or an
acceptance criterion; a criterion cannot be met without loosening a tolerance, widening the
bracket or deleting a check; a new base dependency is wanted; or a numerical failure will
not resolve and the only way to green is to hide it. Decisions already listed in §1 are not
reopened.

## 11. Prototype evidence and risks

The prototype (24 surfaces, 128 nodes per 2π, 24 slices, 4 bins per cell; doubled 48/256/48)
gave, single process on an Apple M1 laptop:

| file | f (h = 1) base → doubled | f (peaked) base → doubled | trapped fraction (h = 1) | time base / doubled [s] |
|---|---|---|---|---|
| DMercFail | 0.240..0.241 → 0.242..0.242 | 0.0381..0.0385 → 0.0402..0.0406 | 0.627 (old team 0.630) | 10 / 64 |
| d23p4 (W7-X) | 0.256..0.261 → 0.256..0.258 | 0.1431..0.1468 → 0.1415..0.1440 | 0.399 (0.396) | 8 / 53 |
| PCA | 0.153..0.218 → 0.143..0.205 | 0.0662..0.0968 → 0.0551..0.0877 | 0.388 (0.385) | 11 / 63 |
| TURBO | 0.094..0.098 → 0.083..0.085 | 0.0036..0.0090 → 0.0024..0.0056 | 0.368 (0.367) | 8 / 46 |
| n3are | 0.087..0.252 → 0.106..0.233 | 0.0065..0.2033 → 0.0173..0.1661 | 0.286 (0.278) | 10 / 52 |

DMercFail and d23p4 are converged at the target already; PCA and n3are have wide brackets
from uncovered nodes of reconnecting ridge arcs (the M2 coverage work), TURBO moves 12 %
between resolutions for the same reason. Full details, logs and Θ maps are in
[validation/prototype](validation/prototype/README.md).

| Risk | Containment |
|---|---|
| Ridge census changes between surfaces | bidirectional partners, uncovered nodes as unknown weight, adaptive surface insertion (M2) |
| Two transitions inside one edge | chain-index change ≥ 2 or a count mismatch triggers subdivision |
| Cells with 1 or 3 cut edges | subdivision until 0 or 2 cuts; residual unknown |
| Bin granularity near saddles of J | halving δJ is part of convergence |
| Logarithmic K near T and Γ_max | integrable; bisected cut cells refine it |
| Long wells, rational surfaces | cap at K_max, capped weight reported unknown |
| Near-axis | graded grid to s_min, axis rule, sensitivity study (M4) |
| Table interpolation error | polish extrema and bounce points on the exact series |
| Unbounded graph on a degenerate sheet | dyadic cap, robust δJ scale, memory watchdog |
| Over-pessimistic metric | by definition; per-λₙ reachable fractions are reported so the reader sees where f comes from |
