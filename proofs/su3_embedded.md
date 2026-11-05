# The su(3) Embedded in the SU(5) Stabilizer of the 16, Verified

Date: 2026-10-26. Status: PROVED (explicit exact matrices over
Q(i, sqrt(3)), 16 counted checks, exit 0). Script:
proofs/su3_embedded.py (sympy, this machine; ~2 s). Stage 1 of the
gluon-index route approved by Kit 2026-10-26 (the "Map physics
notations to model" chat): the target is the algebra of the gluon
octet — su(3) with its standard structure constants, compact form —
sought inside proved material, with every choice flagged as a choice.

## Construction

Work in the stabilizer build's own model (su5_stabilizer.md): the
16 is Lambda^even(W'), W' the creation 5-plane, and for X in the
stabilizer algebra with W-block A(X),

    rho(X)|_piece = dLambda^m(-A(X)^T) + mu(X) I,   mu = (1/2) tr A,

on the pieces m = |S| = 0, 2, 4 (dims 1, 10, 5), subset basis.
On su(5), A is anti-Hermitian and traceless, so mu = 0.

THE SPLIT IS A CHOICE: W' = span{e0,e1,e2} (+) span{e3,e4}
(C^3 (+) C^2). Nothing in the record selects it (Stage 2's problem,
left open here). With E_a = i lambda_a/2 (a = 1..8, lambda the
Gell-Mann matrices), set A_a = diag(E_a, 0_2) in su(5) and
R_a = dLambda^even(-A_a^T).

## Verified exactly

- The f^{abc} table is computed in-script from the lambda matrices
  ([T^a, T^b] = i f^{abc} T^c, T = lambda/2): f^123 = 1,
  f^458 = sqrt(3)/2, totally antisymmetric; and the E-basis law
  [E_a, E_b] = -f^{abc} E_c holds in su(3) itself (sample pairs).
  Controls: the Cartan rho's commute (known-true); a scaled
  generator FAILS the bracket test (known-false) — the checker can
  say no.
- S1: the 8 R_a are linearly independent.
- S2: every R_a is anti-Hermitian and preserves the pieces
  1 (+) 10 (+) 5.
- S3 (the point): [R_a, R_b] = -f^{abc} R_c for ALL 28 pairs, with
  the standard su(3) constants. The contragredient twist in -A^T
  is absorbed by the R_a themselves (each R_a is the image of the
  diagonal automorphism E_a |-> -E_a^T of su(3)), so the constants
  land on the nose, not up to signs.
- S4: tr_16(R_a R_b) = -2 delta_ab (Dynkin index of the 16 under
  this su(3): four triplet copies x -1/2); the abstract identity
  sum_cd f^{acd} f^{bcd} = 3 delta_ab holds, so the Killing form
  in the E-basis is -3 delta_ab — negative-definite: the COMPACT
  real form, not sl(3, R).
- S5: the Casimir C = sum_a R_a^2 on the 16 satisfies the exact
  identity C (C + 4/3) = 0, with tr C = -16 and rank C = 12:
  spectrum {0 (x4), -4/3 (x12)}. So as an su(3) module,
      16 = 3 (+) 3 (+) 3bar (+) 3bar (+) 1 (+) 1 (+) 1 (+) 1,
  four triplets with C_2 = 4/3 and four singlets.
- S6 (the (8,1) placement): the R_a commute with the SU(2) block
  (A = diag(0_3, i sigma_k/2); 24 commutators) and with the
  commutant U(1)' (A = i diag(2,2,2,-3,-3); 8 commutators). The
  octet is the (8,1) of su(5) under SU(3) x SU(2) x U(1)'.

## Scope

PROVED: an su(3) subalgebra with the standard structure constants,
in the compact real form, exists inside the proved SU(5)
stabilizer algebra of the 16, occupies the (8,1) slot, and the 16
branches under it as 3+3+3bar+3bar+1^4 — computed in the
stabilizer build's own spinor model, exactly.

NOT established: (i) any native selection of the 3 (+) 2 split —
it is a choice, one of a conjugate family (Stage 2); (ii) any
physical identification — this is algebra inside the Spin(10)
record, not a field, not a gauge action on spacetime, and it does
not touch BR-FIELD or the ledger's gluon row status (a CANDIDATE
note is filed there, status unchanged); (iii) nativeness in the
core-operator sense — the construction lives in the standard
Cl(10) spinor record, whose operator content is the path weights
only (cf. the so8_triality verdict).

Disclosure: the S4 trace-form expectation was first written
-4 delta (wrong index guess); the computed form is -2 delta,
matching the branching above. The build was right; the
expectation was corrected, not the build.
