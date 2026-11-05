# E6 Seed: the Algebra of the Pair (Spin(10) x U(1), 16 + 16')

Date: 2026-10-23. Status: CONSTRUCTED and PROVED exact
(17 counted checks, exit 0, ~33 s; script proofs/e6_seed.py),
with the identification as e6 by the classification step
cited below. This is the seed of a new ladder, not a rung of
the DIII ladder (diii_ladder.md): the DIII n = 5 rung supplies
the so(10) and both spinors; everything else is built here.

## The construction (complexified)

e = k + p+ + p-, dim 46 + 16 + 16 = 78, where
  k = so(10) (E-coords) + C Z6,   p+ = the vacuum-sector 16,
  p- = the conjugate sector 16',
brackets: [k, k] standard, Z6 central in k acting as +-1 on
p+-; [p+, p+] = [p-, p-] = 0 (Z6 charge); and
  [psi, chi] = A(psi, chi) + lam B(psi, chi) Z6
for psi in p+, chi in p-, where:
- B is the invariant pairing 16 x 16' -> C. Found
  constructively: B(psi, chi) = psi^T M chi with
  M = g1 g3 g5 g7 g9 (charge conjugation by the odd gammas),
  verified sector-swapping (M Gamma = -Gamma M), invariant
  (Sig_j^T M + M Sig_j = 0 for all 45 j), and nondegenerate
  (block det != 0). Uniqueness up to scale is Schur's lemma
  (cited); existence is the computation.
- A is the moment map onto so(10): kappa(A(psi, chi), Y) =
  B(rho(Y) psi, chi), kappa = -tr_32(rho rho) (Gram = 8 I,
  verified). The full 45 x 256 table is computed.
- lam = -3/32 in these normalisations, SOLVED from the
  Jacobi identity (below), not fitted: the identity is linear
  in lam, one nonzero component fixes it, the rest verify.

## Jacobi, by cases

(k,k,k): k is a Lie algebra. (k,k,p): p+- are k-modules
(corpus representations; Z6 by scalars). (k,p+,p-):
equivariance of B (verified in full) and of A. A's
equivariance is DERIVED (B-invariance + kappa-invariance by
trace cyclicity + rho a representation) and additionally
spot-verified exactly for the 5 torus generators x all
(a, b) -- the assembly conventions are what the spot-check
guards. (p,p,p) same-sector: zero. Mixed: the Fierz identity
  rho([x, y]) z + rho([y, z]) x = 0   (x, z in p+, y in p-)
and its mirror -- the one non-formal input -- verified on
all 4096 basis triples each, with the solved lam.

## Identification data (all computed)

- [p+, p-] surjects onto k: extended A-table rank 46.
- Centre = 0, derived from checked inputs: bracketing a
  central element with p+- kills its p-parts (A injective in
  each slot, checked) and forces rho(Y) = -+ t I with
  tr rho(Y) = 0, so t = 0 and Y = 0 by rho injective
  (checked: the 45 spin matrices are independent).
- Perfect: [p+, p-] = k (above) and [k, p+-] spans p+-
  (checked ranks 16, 16).
- Simple, by the standard case argument on the checked
  inputs: an ideal is Z6-graded and a k-submodule in p;
  16 and 16' are irreducible and non-isomorphic (their
  u(5)-branchings in the DIII ladder are conjugate and
  unequal, and Z6 separates the charges), so an ideal
  meeting p in either summand contains [p+, p-] = k and
  then the other summand; an ideal inside k is an ideal of
  so(10) + C normalised by p, and rho != 0 in each case
  kills all four candidates.
- Killing form: Killing_e6 restricted to so(10) is exactly
  -3 kappa (= 3 tr_32(rho . rho)), a single rational
  multiple on all entries (checked); tr(ad_{Z6}^2) = 32 on p.
- A 78-dimensional complex simple Lie algebra is e6 --
  the classification step, cited (78 occurs for no other
  simple type). The pair (k, p) is accordingly the EIII
  symmetric pair, and Z6 is its grading centre.

## Scope

Built: the complexified algebra, its bracket in full
(structure constants implicit in the tables: so(10)
structure constants, the two spin actions, the A-table, B,
lam), Jacobi, simplicity, centre, Killing data. NOT built:
the 27-dimensional representation and its cubic invariant
(the charge pattern 16 + 10 + 1 under this k is the next
layer); real forms (the compact form is a further reality
check on these tables, not made here); the EIII homogeneous
geometry and any Kostant-type operator on it. No FW/saw
content. Local file; not pushed.
