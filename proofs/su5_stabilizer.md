# The SU(5) x U(1) Stabilizer of the 16, Built

Date: 2026-10-23. Status: PROVED (explicit exact matrices over Q
and Q(i), 34 counted checks, exit 0). Script:
proofs/su5_stabilizer.py (sympy, this machine; ~9 s). Companion to
pure_spinor_orbit.md: the orbit build produced the COMPLEX line
stabilizer (parabolic, Levi GL(5, C)); this build produces the
COMPACT real form of that stabilizer and the branching under it.

## Claim

Let e_0 be the Fock vacuum of the 16 (Gamma = -1 Weyl block of
spin10_spinor.md), W its annihilator 5-plane, W' = W-bar the
creation 5-plane (pure_spinor_orbit.md). Then:

1. The stabilizer of the line C e_0 in the compact real form
   so(10) (real span of the Sigma_ab) has dimension 25 and is
   exactly u(5): every element preserves both W and W' (blocks
   B = C = 0 in the decomposition C^10 = W (+) W'), the W-block
   A satisfies A* + A = 0, and the A-blocks span all of u(5).
   It is the compact real form of the orbit build's Levi GL(5, C).
2. Centre dim 1; derived algebra dim 24, equal to ker mu (mu the
   orbit build's character, mu = (1/2) tr A); derived A-blocks are
   traceless and span su(5). So the algebra is su(5) (+) u(1).
   The stabilizer of the VECTOR e_0 (not the line) inside the
   compact form is the su(5) part: X e_0 = 0 iff mu(X) = 0.
3. Branching of the 16, with the centre generator Z normalised by
   A(Z) = -i I on W (so rho(Z)/i is the charge operator):
   the three Z-eigenspaces are the occupation pieces |S| = 0, 2, 4,
   of dims 1, 10, 5, with charges -5/2, -1/2, +3/2 (doubled:
   -5, -1, +3), and as SU(5) representations they are
   1, 10bar, 5 -- NOT 1, 10, 5bar. Precisely: the pieces are
   Lambda^even of the creation space W', which carries the
   contragredient action -A^T (a 5bar): Lambda^0 = 1,
   Lambda^2(5bar) = 10bar, Lambda^4(5bar) = 5 (x) det with det
   trivial on SU(5). All three pieces are irreducible under the
   derived algebra (Burnside span = End of each piece).
4. The conjugate sector 16' under the SAME stabilizer and the
   SAME Z carries the conjugate branching: odd pieces |S| = 1, 3, 5
   of dims 5, 10, 1 with charges -3/2, +1/2, +5/2.

## Verified exactly (script sections)

- The stabilizer is SOLVED, not assumed: real nullspace of the
  W- and W'-preservation equations in the adapted basis gives
  dim 25. Controls: Y0 = sum_k E_{2k,2k+1} is in (true control);
  E_{0,2} is not (false control); comparator controls as usual.
- Block form on all 25 basis elements: B = C = 0, D = -A^T,
  A* + A = 0; spin action preserves C e_0; mu = (1/2) tr A
  re-verified in this basis; A-projection has real rank 25.
- Brackets: bracket span dim 24; centre nullspace dim 1;
  ker mu dim 24; derived = ker mu (rank of the union 24);
  derived A-blocks traceless, real rank 24 inside the
  24-dimensional su(5).
- Centre: A-block a nonzero scalar lam * I; normalised to
  A(Z) = -i I. rho(Z) diagonal on the 32; [rho(Y_j), rho(Z)] = 0
  on the 16 for all j; charges and piece dims as in item 3;
  16' charges as in item 4.
- Branching operator identity, all 25 basis elements, all pieces:
    rho(X)|_piece = dLambda^m(-A(X)^T) + mu(X) I,   m = 0, 2, 4,
  with dLambda^m computed by epsilon-minors (themselves
  cross-checked against the explicit wedge formula for m = 2).
- Epsilon intertwining: the complement-pairing matrix intertwines
  the Lambda^4 action with A + (mu - tr A) I -- the defining
  action with a scalar twist; on su(5) (tr A = mu = 0) it is the
  defining 5 itself.
- Weights (decisive, convention-free): for the Cartan element
  with A = diag(i, 2i, 3i, 4i, -10i) (verified real, i.e. in the
  compact algebra), the spin eigenvalues on piece 4 are exactly
  {a_k} and on piece 2 exactly {-(a_i + a_j)}.
- Irreducibility: for each piece, {I} u {rho(D_t)} u
  {rho(D_t) rho(D_s)} over the derived basis spans End(piece)
  (ranks 25, 100, 25); no invariant vector in the 10 piece.

## Disclosures (expectations corrected en route)

- First draft of the operator identity used -A in place of -A^T:
  it failed on 10-15 of 25 basis elements per piece; the -A^T form
  passes all 25 x 3. The contragredient is the content: the spinor
  is Lambda^even of the CREATION space W', whose action is -A^T.
- Consequently the draft's branching guess (1 (+) 10 (+) 5bar,
  the form usually quoted) is wrong FOR THIS SECTOR AND THIS
  EMBEDDING: the weights computation (no convention freedom --
  distinct eigenvalues) gives 1 (+) 10bar (+) 5 with the charges
  above. The usually-quoted form is the conjugate branching,
  carried here by the 16' (item 4) up to the charge normalisation;
  which sector is called "the 16" and which space (W or W') is
  called "the 5" is exactly the labelling choice that flips it.
  Both sectors and both labellings are pinned in this note, so
  nothing is left to convention.
- An earlier commutant test (simultaneous-commutant nullspace)
  was mis-built (a real-split system that conflated the unknown
  vector's real and imaginary parts) and slow; it was discarded
  before producing any verdict and replaced by the Burnside span
  test, which is decisive and runs in seconds. No result in this
  note rests on the discarded test.

## Scope

Algebra level throughout: the stabilizer GROUP is
(SU(5) x U(1))/Z_5 (the diagonal Z_5 quotient is the standard
group-level statement; what is computed here is its Lie algebra
u(5) = su(5) (+) u(1) and the representation branching, which is
where the content is). "GUT breaking" is the physics name for a
spinor vev along the singlet direction breaking Spin(10) to this
stabilizer; the branching itself is the mathematics proved here.
No FW/saw operator content: this is the standard Spin(10)
geometry of the 16, built under spin10_spinor.md and
pure_spinor_orbit.md. Local file; not pushed.
