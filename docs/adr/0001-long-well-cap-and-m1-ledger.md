# Long wells and the M1 ledger gate

Status: **proposed**; no change to the accepted cap or tolerances is authorized here.
Date: 2026-11-10. Affected milestone: M1.

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

## Options and consequences

1. **Allow bounded adaptive chain extension**, initially up to 256 periods only
   for capped nodes. Keep the finite cap, unknown status when still capped, 2 GB
   watchdog and 300 s equilibrium target. Refine the transverse quadrature near
   long-well boundaries, since a large K at one node does not establish the
   population carried by its strip. This adds selective setup cost and needs
   five-file evidence; the small probe does not prove the overall budget.
2. **Keep 64 periods and add independent weight accounting for truncated wells.**
   This requires a separate method to measure the unresolved population without
   forcing the ledger to equal its reference integral. It changes M1's work and
   must specify how unknown connectivity propagates at M3.

## Proposed decision

Choose option 1, subject to researcher approval. Do not relax either ledger
tolerance, classify a capped chain, or start M2 while M1 remains incomplete.
The production code in this PR retains the original 48-period default. The
128/256-period runs are explicitly diagnostic experiments.

## Verification required after the decision

Repeat the 30-case matrix with both sources, original failures retained, actual
global b values, chain-length and angular refinements, and explicit unresolved
weight accounting. Close the separate d23p4 error within 1%. Compare the total
ridge weight to the independent volume integral without using that integral to
fill missing contributions. Run the normal unit, mutation, runtime and CI gates.

This is not evidence that the ridge formulation should be abandoned. It shows
that the prototype's 19-period observation is not a safe bound. The present
implementation cannot yet claim accurate f or the full 300 s target; those
claims require M3–M5 evidence and a usable treatment of these long wells.
