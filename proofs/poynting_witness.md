# Poynting-Flux Dichotomy Inside the Maxwell Exercise's Mode System

Date: 2026-10-26. Status: PROVED as mode algebra (24 counted checks,
exit 0, exact at three points + one general symbolic identity);
the physics reading stays ASSERTED-WITNESS under the exercise's own
dictionary (maxwell_exercise.md Step 5). Script:
proofs/poynting_witness.py (sympy, this machine; ~7 s). Kit's
Poynting tangent, 2026-10-26: test whether the exercise's families
reproduce the classical flux dichotomy — constant flux for circular
polarization, pulsing flux for linear.

## Setup (the exercise's objects only)

Placement vector omega = (sxp, cxp, uxp); its cross matrix Y = [omega]x
(charpoly lambda(lambda^2 + theta^2), theta^2 = |omega|^2, re-verified
as a control). Step 5's mode equations k x E = w B, k x B = -w E
(c = 1) make B a positive scalar multiple of khat x E, so the
instantaneous Poynting vector along k is S(phi) ~ |E(phi)|^2 khat,
E(phi) = Re{E0 e^{i phi}}. Only the phi-dependence is tested;
positive scalar factors are irrelevant to it.

## Verified

- General identity (symbolic, any complex amplitude E0 = a + ib):
    |E(phi)|^2 = |E0|^2/2 + pulse,  pulse at 2 phi with amplitude
    |E0.E0|/2 (exactly: ((|a|^2-|b|^2)/2) cos 2phi - (a.b) sin 2phi).
  So the flux is constant iff E0.E0 = 0 — the circularity condition.
- At three exact points (tie pi/4 in Q(sqrt2, i, sqrt10); pi/6 in
  Q(sqrt3, i); pi/3), with v = a + (i/theta) b built from an
  unnormalized frame (a perp omega, b = Y a):
  Y v = -i theta v exactly (helicity eigenvector), v.v = 0 exactly,
  the circular flux is constant in phi, and a real transverse mode's
  flux is cos^2(phi) — pulsing, with zeros.
- Model-specific closed form: theta^2 = |omega|^2
    = 2 (sxp^2 + cxp^2 - uxp + 1),
  equivalent to the census lock sxp.cxp = uxp - 1 (re-verified
  exactly at all three points; the closed form also checked against
  the definition at three rational x by 50-digit evaluation,
  residual < 1e-40). Values: theta^2(tie) = 10,
  theta^2(pi/6) = 26 - 12 sqrt3.

## Scope

The dichotomy itself is general cross-matrix algebra — true for any
placement vector — so what is verified is witness CONSISTENCY: the
exercise's helicity sectors are exactly its constant-flux modes,
and the linear modes pulse, inside the model's own mode system,
with the model-specific theta^2(x) above. This does not derive a
Poynting vector, an energy, or a spacetime flux; BR-FIELD's status
is unchanged.

## Disclosures (three checker defects, all in the checker)

1. First run's charpoly control compared factored polynomials in
   mismatched symbols (false FAIL); fixed by comparing coefficients.
2. First run's helicity check used a NORMALIZED frame; the nested
   radicals defeated sympy's simplify and returned false FAILs at
   both points, while v.v = 0 passed. Rebuilt with an unnormalized
   frame (v = a + (i/theta) Y a); the eigenvector identity then
   verifies exactly.
3. The x-space numeric control used mpmath, which in this
   environment ignores mp.dps and computes tan at float53
   (diagnosed: mpmath's c agreed with a 50-digit evaluation only to
   1.8e-15, constant across dps 40/60/80; the sympy-50 residual of
   the identity is 8.6e-50). Controls moved to sympy N at 50 digits.
