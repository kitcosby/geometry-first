# Rigorous Testing Plan — Symbolic Harmonics of Angular Deltas

Date: 2026-10-12. Status: PLAN (not yet executed).
Goal: Use Sage (exact arithmetic only) to determine the harmonic (Fourier)
structure of the 8 R-operators as functions of angle, and of their finite
differences ("deltas") under translation and doubling.

## 0. Discipline
- Exact fields only: Sage QQbar and CyclotomicField for values at rational
  multiples of pi; symbolic ring for identities.
- Floats permitted ONLY as a non-decisive sanity cross-check, never as a result.
- Every test ends PROVED (symbolic identity + human-readable argument),
  CHECKED (exact values at stated sample angles), FALSIFIED, or INCOMPLETE
  (with the timeout/failure disclosed).

## 1. Objects
Operators O in {srx, sxp, cxp, crx, urx, uxp, sawup, sawdown, FW} on
D = R \ {k pi/2}. On each open octant (k pi/4, (k+1) pi/4) all absolute-value
signs are constant, so each O is an honest ratio of trig polynomials.

Key variable: zeta = exp(i x). On a fixed octant each O is a rational
function R_O(zeta) (sin x = (zeta - zeta^-1)/2i, cos x = (zeta + zeta^-1)/2,
abs resolved by the octant's sign pattern alpha, beta).

## 2. Deltas under study
- Translation delta: D_h O(x) = O(x+h) - O(x), h in {pi/12, pi/6, pi/4, pi/2}.
  In zeta: x -> x+h is zeta -> zeta * exp(i h), an exact root-of-unity factor.
- Doubling delta: D_2 O(x) = O(2x) - O(x).  (zeta -> zeta^2.)
- Boundary delta (jump): J_O(k pi/2) = right limit - left limit at each
  singularity/cusp, computed exactly.

## 3. Sage procedures
1. Octant normal forms: for each octant, derive R_O(zeta) symbolically in
   Sage's SR with the octant sign assumptions; verify against QQbar values
   at the octant midpoint (exact equality test).
2. Fourier coefficients (global, piecewise): a_k, b_k by exact symbolic
   integration over the 8 octants (integrate with assumptions), summed.
   Report coefficients as exact surds / rational multiples of 1/pi.
3. Laurent content (per octant): expand R_O(zeta) as a Laurent series in
   zeta about a sample point; record which powers occur (the "harmonics"
   present) and whether the rational form is finite (Möbius) or infinite.
4. Delta identities: compute D_h O, D_2 O symbolically per octant; test
   candidate closed forms by exact simplification to zero.

## 4. Tests with predictions
- H1 (FW square wave). FW = sgn(sin 2x) is a square wave. Prediction:
  FW = (4/pi) * sum_{k>=0} sin(2(2k+1) x)/(2k+1); only odd harmonics of
  sin 2x, coefficients 4/(pi(2k+1)). Verify k=0,1,2 coefficients exactly.
  Falsifier: any even harmonic or cos term with nonzero coefficient.
- H2 (Primitive normal forms). Prediction: on QI, srx = cot(x/2) is Möbius
  in zeta^{1/2} (finite rational, degree 1), so its Laurent content is a
  geometric tail, not new harmonics. Verify R_srx in closed form.
- H3 (Doubling delta). On the FW=+1 branch, srx(2x) = (srx - sxp)/2 PROVED,
  so D_2 srx = -(srx + sxp)/2. Prediction: harmonic content of D_2 srx is
  exactly the union of the contents of srx and sxp on that branch (no new
  frequencies). On FW=-1 use (crx - cxp)/2 — test the family swap exactly.
- H4 (Quarter-turn delta). Prediction: translation by pi/2 swaps the
  families (srx <-> cxp up to sign/reciprocal per octant). Compute
  D_{pi/2} srx symbolically per octant and identify the closed form.
- H5 (Boundary deltas). Compute J_O at 0, pi/2, pi, 3pi/2 exactly.
  Known from 2026-10-12 numeric work: cxp jumps at pi/2 (left +inf,
  right ~0); srx, sxp continuous (limit 1). Promote the jump table to
  exact one-sided limits (Laurent expansions at the boundary, not Taylor —
  several operators have poles/cusps there, so Taylor does not exist).
- H6 (Reciprocal pairs). For sxp = 1/srx, crx = 1/cxp: test whether the
  delta of a reciprocal is determined by the delta of the base operator:
  D_h(1/O) = -D_h O / (O(x) O(x+h)) — exact algebra, verify symbolically.

## 5. What counts as "discerning the harmonics"
Success = for each operator (or at minimum FW, srx, cxp), an exact,
octant-resolved statement of (a) the rational normal form in zeta,
(b) the global Fourier coefficients in closed form or a proved pattern,
(c) the delta identities D_h, D_2 in closed form.
Partial success with honest scope labels is acceptable; a numeric FFT is NOT.

## 6. Known risks
- sawup = 1/urx, sawdown = 1/uxp have poles where urx/uxp vanish; their
  Fourier coefficients may fail to exist in L^2. Test FW and the 4
  primitives first; treat saw pair separately (principal value / distribution
  sense only if needed, labeled as such).
- Sage `integrate` can time out on piecewise abs expressions; mitigation:
  integrate per octant with signs already resolved (Procedure 1 first).
- zeta^{1/2} appears via half-angles; work in CyclotomicField(2m) when the
  angle grid is pi/m, to stay exact.

## 7. Order of execution
1. Procedure 1 (normal forms) -> enables everything else.
2. H5 (boundary deltas; settles the pi/2 equilibrium question rigorously).
3. H1 (FW harmonics). 4. H3, H4, H6 (delta identities). 5. H2 write-up.
