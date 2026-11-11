# Long wells and the M1 ledger gate

Status: accepted
Date: 2026-11-10, revised 2026-11-11.

## Evidence and problem

The M1 implementation in this PR reaches the contract's stop condition. It does not
complete M1. The [30-case report](../validation/m1/README.md) preserves all results,
including the failed first attempt and the subsequent numerical fixes.

The design anticipates 48–64 field periods per chain and says the longest observed
prototype well was 19 periods. Two actual PCA nodes at the required global pitches
have no exit in 64 periods. Diagnostic extensions find their first exits at
**73.80031579** and **114.84383316** periods. Their time lengths K are about 4065.88
and 5283.29 in the input length units. These are finite trapped wells, rather than
passing or confined states.

The extended scans agree at 128 and 512 samples per period. Direct Fourier-series
root finding verifies every sampled maximum through each exit: every preceding
maximum is below b and the exit maximum is above b. Composite quadrature at orders
32 and 64 agrees in K to 1.3e-7 relative or better. The small controlled experiment
takes 0.89 s and 0.14 GB on this cloud host; this is not an equilibrium timing.
See [long-wells.json](../validation/m1/long-wells.json) and the reproducible
[probe driver](../../scripts/probe_m1_long_wells.py).

At the unchanged 48-period setting, four nodes in three cases remain capped: PCA
lambda_n=0.8 and 0.95, and TURBO lambda_n=0.9. A capped node has unknown K. The
ridge ledger therefore reports its total and ratio as incomplete, rather than
inserting zero, copying a neighbouring K, or inferring its weight from the
independent ledger. A count of capped nodes is not their birth-weight fraction.
The current contract simultaneously permits capped unknown states and requires
every M1 ledger ratio to pass. How should a finite ledger total be obtained for
these states within the planned cap?

There is a separate, smaller unresolved gate: d23p4 at lambda_n=0.95 has the
uniform-source ridge/volume ratio 1.010127583, outside the required 1% tolerance.
This is a numerical refinement task; it is not grounds to round the ratio to a
pass or change the tolerance. All failed cases remain in the evidence.

### Why long wells occur and why one node matters

A well ends at the first along-line maximum above b. It is long when many successive
maxima lie just below b. Two mechanisms produce this:

1. **Marginal trapping.** Near Γ_max (H₀ ≈ b) at any pitch, most maxima on the
   surface are only slightly lower than the starting one.
2. **Near-rational ι/N_fp.** Near a low-order rational, for example 1/5, the line
   returns close to its start every few periods and meets nearly the same maxima again,
   drifting slowly in α. The well length scales roughly as the inverse of the α-drift
   per near-return. On the exact rational, a line whose maxima are all below b never
   exits; that is passing. No fixed cap therefore bounds every well.

K grows linearly with well length, and a well of 70 or more periods occupies a large
part of its surface. One node can therefore carry a large share of a surface's weight.
For PCA at s = 0.446, lambda_n = 0.8, the per-surface ridge sum is 203.0 at 48 periods
with one capped node, and 343.9 at 256 periods with none. The volume integral on that
surface is 347.4.

### How much weight is in long wells

`scripts/adr0001_compute.py` (diagnostic) traced every base ridge node to 256 periods on
29 surfaces for lambda_n = 0.40 to 0.99 in steps of 0.01.
`scripts/adr0001_tail_weight.py` integrates the per-slice weight ∑K|Δα| with ds and
db/b². The table gives the fraction of that integral carried by wells longer than each
length. The h=1 and peaked sources are separated by "/".

| File | > 48 periods | > 64 periods | > 128 periods | slices with no exit at 256 |
|---|---|---|---|---|
| PCA | 0.54 / 0.54 % | 0.37 / 0.31 % | 0.29 / 0.28 % | 2 of 1034 |
| TURBO | 0.24 / 0.01 % | 0.22 / 0.01 % | 0.11 / 0.00 % | 0 of 1119 |
| DMercFail | 0.12 / 0.00 % | 0.06 / 0.00 % | 0.06 / 0.00 % | 0 of 1502 |
| d23p4 | 0.11 / 0.18 % | 0.11 / 0.18 % | 0.00 / 0.00 % | 0 of 1037 |
| n3are | 0.57 / 0.92 % | 0.26 / 0.00 % | 0.05 / 0.00 % | 0 of 788 |

This is a coarse estimate, not a ledger, and it has three known biases:

- lambda_n < 0.4 is missing from the denominator, which makes the fractions too large;
- wells still capped at 256 have no K and are missing from the numerator, which makes
  them too small;
- one node carries a whole node spacing of measure, which is noisy where a long-well
  strip is narrower than that spacing.

The two no-exit PCA slices lie near ι/N_fp = 1/5 (s ≈ 0.36–0.40).

The conclusions are:

- Long wells are locally important but integrally small.
- The tail beyond 128 periods decays slowly on PCA because of the rational surface.
- Unknown births move the lower bound by at most their weight. They can move the upper
  bound by more, because an unknown cell is a wildcard connector and a long well spans
  much of its surface.

## Options and consequences

1. **Fixed cap, residual unknown (current contract).** Keep 48–64 periods and report
   capped weight as unknown.
   - Simple and honest.
   - Wells beyond 48 periods hold up to about 0.9% of births, which is a material part
     of the 0.01 absolute tolerance before any wildcard effect.
   - As written, M1 cannot complete: the four capped M1 nodes keep their ledgers
     incomplete.
2. **Bounded adaptive extension of capped nodes only.**
   - Trace every chain to 64 periods, then extend chains still capped up to a hard
     ceiling of 256 periods. Keep the unknown status for chains still capped, the 2 GB
     watchdog and the 300 s target.
   - All four M1 capped nodes exit by 115 periods, so every 30-case ledger becomes
     computable without changing a tolerance. The two TURBO lambda_n=0.9 nodes exit at
     99.69 periods with K ≈ 6019 and 5999. That is from the diagnostic trace in
     `scripts/adr0001_compute.py`; unlike the PCA pair, it is not yet verified on the
     exact series.
   - Extension is cheap: about 0.04 s per node at 256 periods in the probe, and only a
     few nodes need it.
   - It does not remove the near-rational tail: PCA has nodes with no exit at 256
     periods.
   - A large K at one node does not establish the weight of its strip. Long-well nodes
     need local transverse refinement, and that refinement, not the tracing, is the
     real cost.
3. **Keep the cap and account for the missing weight independently.** For example,
   estimate it as the volume integral minus the resolved weight.
   - No long tracing, and the totals always close.
   - Filling the gap from the reference integral makes the ledger agree by construction
     on exactly the slices it should test, which loses the bookkeeping check.
   - It gives no Θ for these states, so they remain in the bracket.
   - It changes M1's work and must specify how unknown connectivity propagates at M3.

## Decision

Choose option 2, with the stopping rule defined by a weight budget rather than by the
period count alone. All of this is subject to researcher approval.

1. **Base cap 64 periods. Extend capped nodes adaptively, doubling up to a ceiling of
   256 periods.** The ceiling is a cost guard, not an accuracy claim. Production
   defaults change only after the verification below.
2. **Transverse refinement.** Subdivide the ridge edges adjacent to any node whose well
   exceeds 64 periods (dyadically, up to 16× the base density) until the strip weight
   changes by less than the slice tolerance between successive levels. Record the level
   reached.
3. **Residual capped nodes stay unknown.**
   - Report their strip weight as a lower bound, using the partial K accumulated to the
     ceiling.
   - Never assign them zero, passing, confined, a neighbour's K, or weight inferred from
     the volume integral.
   - They enter the M3 bracket as wildcard cells.
4. **Ledger acceptance.**
   - On every slice with no residual capped node, the existing tolerances apply
     unchanged.
   - On a slice that still has residual capped nodes, report the ratio as incomplete
     with the resolved sum, the lower bound from item 3, and the gap to the volume
     integral. The check requires resolved sum + lower bound ≤ volume × (1 + tolerance).
     The gap is a diagnostic, never a filled-in value.
   - The current 30-case matrix has no such slice after extension. Every M1 ratio can
     therefore be judged against the unchanged tolerances. Whether the extended slices
     pass is still to be shown, and the separate d23p4 error remains.
5. **Unknown-weight budget (new acceptance item, M3 onward).** For each file and source,
   the total birth weight of residual capped cells must be at most 0.5% of V_h. Report it
   separately from frontier and uncovered weight. If a file exceeds the budget, stop and
   write an ADR; do not raise the ceiling silently.

The budget keeps the long-well share of the bracket half-width at about 0.0025 or less
by the cells' own weight. That share is a quarter of the absolute tolerance floor. The
coarse survey suggests every file meets the budget with the 256-period ceiling. The
wildcard effect on the upper bound must be measured at M3, not assumed.

Do not relax either ledger tolerance, classify a capped chain, or start M2 while M1
remains incomplete. The production code in this PR keeps the original 48-period default.
The 128- and 256-period runs are diagnostic experiments.

## Verification required after the decision

- Repeat the 30-case matrix with both sources. Retain the original failures, and use the
  actual global b values, chain-length and angular refinements, and explicit unresolved
  weight accounting.
- Show that the four previously capped nodes now exit and that their slices close
  within tolerance. Show that transverse refinement converges the strip weight.
- Rerun the long-well survey at production resolution with the adaptive rule. Report
  per file:
  - the residual capped weight against the 0.5% budget;
  - the ceiling actually reached;
  - the time and memory spent on extension and refinement.
- Close the separate d23p4 error within 1%.
- Compare the total ridge weight with the independent volume integral without using that
  integral to fill missing contributions.
- Run the normal unit, mutation, runtime and CI gates. Add a test that a chain with no
  exit at the ceiling is unknown with a lower-bound weight and is never passing.

This is not evidence that the ridge formulation should be abandoned. It shows that the
prototype's 19-period observation is not a safe bound. The present implementation
cannot yet claim accurate f or the full 300 s target; those claims require M3–M5
evidence and a usable treatment of these long wells.
