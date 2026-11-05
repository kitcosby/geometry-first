# O_prop next step — correction and reframed choice (DRAFT)

Date: 2026-10-12. Prepared after the Algebra-first chat message of
19:20 UTC, which offered Kit a choice between (a) running the pending
R(theta)=0 leakage test on the corrected 2026-10-09 O_prop candidate and
(b) writing the self-propagation involution T's covariance law up as a
formal selection constraint. Option (a) as offered is NOT executable,
and one status in that message was understated. Correction below, with
the audit record.

## What the record actually says

Sources: io-word sweep RUN229_THETA_TYPING.md and
RUN233_OPROP_AUDIT.md (2026-10-09, exact sympy, 9/9 and 13/13 PASS).

1. R(theta) has NO definition in the extracted transfer-matrix
   formalization (run 229). The day-log's offered step "run the
   corrected O_prop through R(theta)=0" was already audited in run 233
   and returned BLOCKED: "cannot be executed as stated... No definition
   invented." So there is no pending runnable leakage test; there is a
   pending DEFINITION, and it is Kit's to supply (or the sector residual
   R denotes).

2. R(theta)=0 is also the wrong frame even once R is defined: for
   scalar theta it is 12 equations in 1 unknown (11 excess, generically
   overdetermined, run 229 C1). The formalization's own argmin frame —
   theta_phys = argmin ||R(theta)||^2 — is the consistent reading. The
   correct scalar-theta uniqueness object is the scalar Hessian
   f''(theta*) = 376 (J^T J + R^T H) > 0, not a "non-singular Jacobian"
   (the Jacobian is 12x1, non-square — ill-typed as stated, run 229 J1).

3. Minimal consistent repair (run 229, exact typing): per-sector triple
   (theta_E, theta_B, theta_O), W_O using theta_O only. Every stated
   formula stays typed; the selection becomes a genuine 3-dim argmin
   with a square 3x3 Hessian; the G-orthogonal points theta_O in {0,2}
   survive.

4. Gram-compatibility is not "likely false" — it is FALSIFIED exactly
   for both readings audited (run 233 E1/B1): diagonal-O_E reading has
   defect -3/4 at x=pi/4; fixed complex-structure J reading has
   M^T M - I = [[I/4, J], [-J, I/4]] != 0 at pi/4. And a proved
   no-constant-J theorem (run 233 B2): no x-independent J makes
   O_B(x) = [[D,J],[-J,D]] orthogonal on an interval.

5. Constructive repair already proved (run 233 B3):
   J(x) = diag(sqrt(1-1/urx^2), sqrt(1-1/uxp^2)) gives M^T M = I_4
   exactly wherever |urx|, |uxp| >= 1. An orthogonal O_B exists, but
   only with x-varying J.

6. One new positive already on record (run 233 T1, sweep proposal, not
   Kit's equation): identifying theta := omega(x) = sin(2x)/4 gives
   det(T) = (k_10/k_0)^2 (1-omega)^{10r} with 1-omega >= 3/4 > 0 for all
   real x — the run-224 theta=1 singularity is avoided on the
   omega-locus.

## The reframed choice for Kit

Option A (leakage/argmin test) — executable only after Kit supplies:
  (i) a definition of R (or the sector residual it denotes);
  (ii) the exact O_E entries and J for the candidate (gaps flagged in
      run 233, not invented);
  (iii) a decision on the theta reading: scalar-argmin (run 229 repair)
      vs per-sector triple;
  (iv) adopt or reject the B3 J(x) repair and the theta := omega(x)
      identification.
  With (i)-(iv) the argmin computation is a concrete exact-sympy run.

Option B (T as a selection constraint) — writable now, from proved
material only. The verified form of T (quadrant_flip_polar.md, PROVED):
the pi-shift zeta -> -zeta, w -> i w, i.e. the Mobius involution
t -> -1/t, maps EVERY operator's QI formula to the SAME operator's QIII
formula (the involution preserves each operator; it does not swap
urx <-> uxp — that first-pass reading was corrected in the polar file).
At block level: M1_QIII = diag(1,-1) M2_QI diag(1,-1). FW =
sgn(Im zeta^2); alpha = sgn(Im zeta) selects the perfect square.
Stating this as a constraint ON candidates in the 12-dim b still needs
one declared input: how the quadrant involution acts on b's 12
dimensions. That action is not in the record; it would be Kit's
declaration or a separate derivation, not an invention.

## Status

CP-OPROP stays OPEN. The 2026-10-09 candidate is exactly audited, not
verified. Nothing above is a new candidate matrix.
