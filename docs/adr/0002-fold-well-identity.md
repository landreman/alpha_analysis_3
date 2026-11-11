# Well identity and continuation through generic folds

Status: proposed — researcher decision required; M2 is incomplete.
Date: 2026-11-11.

## Evidence and problem

The M2 implementation has reached the STOP condition in AGENTS. The current
candidate treats every polished ridge-fold endpoint as inactive with K=0.
That assumption is false. It affects `ledger._refine_layer` (including the
inherited M1 comment) and `transitions.surface_transitions`. The candidate
must not be accepted on its apparently passing population ratios.

Consider the positive analytic field, with N_fp=1 and iota=I=0:

    B(theta,phi) = 6 + (-2 + 0.5 cos(theta)) sin(phi) + sin(2 phi)
    C(s) = G(s) = 3 (1 + 0.2 s),    b = 5.

At theta=pi/2, phi=0, the field has a generic fold:
D_parallel B = D_parallel² B = 0, D_parallel³ B = -6, H=6>b.
The minimum which annihilates with the open maximum is **behind** it in
physical flow. The outgoing well survives. Its next minimum is at
phi=2pi/3, B=6-3sqrt(3)/2<5. Independent analytic-field bounce-root finding
and sine-substituted QUADPACK integration at s=0.5 give

- bounce pair (1.1124608386367665, 2.884599391840837);
- J=0.6601101743548781 and K=4.1614783045408625;
- estimated quadrature errors below 1.8e-12.

On the side with an open maximum, it represents this well; the preceding
closed-arc maximum is inactive because its first minimum remains above b.
On the other side the open maximum is absent and the closed-arc maximum
represents the same finite well. The physical bounce points and J/K have
matching finite limits. All maxima remain above b locally, so T=Gamma_max=0:
the existing transition-count check cannot detect a missing continuation.

Holding this J fixed, independent quadrature gives theta=1.6232166429393162
at s=0.4 (closed representative), theta=pi/2 at s=0.5 (fold boundary), and
theta=1.519739011305752 at s=0.6 (open representative). This is a concrete
constant-J curve crossing the change of representative without a split/merge
event. Its two representatives need not be adjacent on the ridge torus.

There is a separate extremum-resolution defect: at theta=pi/2-0.001, the
128-sample chain scan misses the preceding narrow max/min pair and reports
**both** representatives active, with identical J/K. A 512-sample scan resolves
the pair and reports only the open representative active. At pi/2-0.0001,
512 also misses it, while 4096 resolves it. A fixed scan increase cannot
certify arbitrarily close approaches to a generic fold. The false duplicate
can pass T=Gamma_max=0. This is a plausible source of excess ridge weight;
the proportion of the real-input excess caused by it has not been measured.

This phenomenon occurs in a required input: one polished n3are surface
(s=0.7922826086956523) has finite fold limits at lambda_n=0.1, 0.5 and 0.8,
with K approximately 35.47, 127.75 and 100.03 respectively. Those are diagnostic
one-sided evaluations, not strict maxima at the degenerate point or a repaired
population ledger. The separate n3are projection report also finds sampling
nodes projected onto minima. That is an engineering defect to repair after the
endpoint and incidence contract is settled; it is not itself an approval issue.

See the [reproducible JSON](../validation/m2/fold-identity.json),
[figure](../validation/m2/fold-identity.png), and
[driver](../../scripts/probe_m2_fold_identity.py). No tolerance or input subset
was changed. The diagnostic is below the 2 GB watchdog ceiling.

## Scientific decision needed

DESIGN sections 3.3–3.4 identify ordinary ridge states and interior-barrier
split/merge links, but do not specify the nonlocal identification needed when
an above-b fold removes a representative of a surviving well. Section 4's
shared-edge and T-link recipe therefore does not yet determine how this
constant-J curve continues. An endpoint's zero geometric area does not justify
deleting the connectivity of a finite family of curves through it.

Approve or revise this proposed clarification before implementing it:

1. Keep the accepted ridge/well formulation. Treat a disappearing maximum as
   a change of representative when its physical bounce interval survives.
   Link the limiting well to the appropriate surviving predecessor on the
   **same field line and surface**, modulo the physical periodicity. This is
   ordinary same-J continuation, distinct from a T split/merge. Require
   matching limiting bounce intervals and actions; do not infer identity
   from torus distance or a small endpoint weight.
2. Distinguish the two fold orientations using sign(C) D_parallel³ B.
   When the coalescing minimum is ahead, verify the appropriate inactive
   limit; when behind, compute the actual surviving outgoing-well limit.
   Check activity at the specified conserved b. Never assign K=0 to all folds.
3. Certify or bound unresolved first extrema near folds using local fold
   information and root isolation rather than trusting a uniform scan.
   Keep failed projection, unresolved identity and possible duplicates
   explicit. Propagate their unknown incidence to every possibly connected
   sheet in M3's upper graph; a population percentage alone cannot justify
   removing a connection. Any finite budget claim needs an actual bound.
4. Add named acceptance checks for finite above-b fold limits, both physical
   flow signs, matching representatives along a constant-J curve, and the
   missed narrow max/min pair. Extend the all-30-case fold inventory and
   retain failures. Re-run M1/M2 population and transverse checks after repair.

This proposal changes no metric, source, bounce-field normalization, input,
error tolerance, dependency, long-well ceiling or milestone order. It makes
the physical well quotient and its acceptance checks explicit. It is not a
request to relax a failing gate. Alternative formulations must demonstrate
the same physical continuation on this counterexample.

## Current disposition

[PR #2](https://github.com/landreman/alpha_analysis_3/pull/2) stays draft and
STATUS keeps M2 unchecked. The partial transition/cell implementation,
failed real-input evidence and diagnostic counterexample remain reviewable.
No M3 connectivity or f is claimed. No readiness-triggered Claude review is
requested until the scientific gates pass.

The PCA 0.95 refinement alone takes 434.5 s including 385.3 s setup, already
above the eventual 300 s full-equilibrium target on this host. Current code
has not demonstrated that target; adding a dense uniform chain scan is neither
a root certificate nor a credible speed solution. Local fold-aware geometry
may help, but no speed or full-f accuracy claim is supported yet.
