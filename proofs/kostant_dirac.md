# Kostant's Cubic Dirac Operator for Spin(10)/U(5), Built

Date: 2026-10-23. Status: the operator is ASSEMBLED from computed
parts and every finite algebraic input is PROVED exactly (18
counted checks, exit 0); the D^2 Casimir formula itself is
Kostant's theorem, CITED (see Scope). Script:
proofs/kostant_dirac.py (sympy + an exact sparse evaluation
engine over the 1024 monomials of Lambda*(n), coefficients in
Q(i) carried as Fraction pairs; ~35 s). Sixth and last layer of
the 16's stack, on spinor_bundle.md (the pair (g, h), m, the
tangent isotropy) and su5_stabilizer.md (h and its action).

## Conventions (fixed once)

kappa(X, Y) = -tr(XY) (fundamental trace) on so(10); note
B_Killing = 8 tr(XY) for so(10), so kappa = -B_Killing / 8 --
the rescaling is recorded because the cited D^2 constants are
normalisation-sensitive. gamma is the Clifford action for
kappa|_m on S = Lambda*(n), n the B-block radical of the
parabolic (dim_C n = 10, dim S = 1024 = 512 even + 512 odd;
gamma is odd, swapping the halves). With {M_j} the real basis
of m from spinor_bundle.py and {M^{j*}} its kappa-dual basis,
Kostant's operator is the explicit element of U(g) (x) Cl(m):
    D = sum_j M_j (x) gamma(M^{j*})  +  cubic term.

## Verified exactly (script sections)

1. Metric data. kappa Gram on m nondegenerate; dual basis
   computed and checked (Gram x inverse = I); kappa(h, m) = 0
   (25 x 20); kappa invariant under h on m (the trace-cyclicity
   identity, checked on the sparse h basis x all m pairs).
2. Clifford model. The m_C basis {e_t, e-bar_t} (radical
   matrices; e-bar_t normalised as the kappa-dual of e_t) has
   Gram [[0, I], [I, 0]] computed from the matrices. gamma(e_t)
   = wedge, gamma(e-bar_t) = contraction on Lambda*(n) then
   satisfies {gamma_a, gamma_b} = kappa_{ab} I for all a <= b,
   verified by evaluation on all 1024 monomials. Diagonal
   Clifford relations for the real basis gamma(M_j) (built from
   its B/C blocks) verified the same way; off-diagonal real
   relations follow by C-linearity from the m_C ones.
3. Spin lift. sigma(X) = sum_{pq} M_{pq} a+_p a_q - (1/2) tr(M),
   M = ad(X)|_n = dLambda^2(A), for X in a sparse A-unit basis
   of h (verified real, in the solved h, A recovered). The
   intertwining [sigma(X), gamma(v)] = gamma([X, v]) holds for
   all 25 X and all 20 basis v, by evaluation on all monomials
   (the e-bar side uses the contragredient action, checked
   consistent: N = -M^T by the block formulas). tr sigma(X) = 0
   on S for all X. A lift satisfying the intertwining is a Lie
   homomorphism: the defect commutes with every gamma, hence
   with Cl(m_C) = End(S) (standard, even dimension), hence is
   scalar, and trace 0 kills it.
4. gamma-contraction. sum_j gamma(M_j) gamma(M^{j*}) = 10 I on
   all 1024 monomials -- the canonical metric contraction
   (dim_C n = 10 pairs each contributing {gamma, gamma} = I).
5. The cubic term is ZERO. phi(M_j, M_k, M_l) =
   kappa([M_j, M_k], M_l) = 0 for all j < k < l (the rest by
   antisymmetry). Visible reason: [m, m] subset h
   (spinor_bundle.md) and kappa(h, m) = 0 (item 1) -- for a
   Hermitian symmetric pair the Cartan 3-form of m vanishes,
   so Kostant's cubic operator here is the plain Dirac sum
   above, with no cubic correction.

## Disclosures

First run: 13 PASS / 5 FAIL. Two were my errors, both fixed and
re-run to 18/18. (i) The dual radical basis was normalised with
C-block U/4, giving kappa(e, e-bar) = 1/2; the correct dual is
U/2 (tr(U * U/2) = -1). That single factor accounted for the
Gram, Clifford, diagonal-Clifford, and contraction failures.
(ii) The kappa-invariance check was twice written as the wrong
trace expression (it tested tr(X [M_j, M_k]), a different,
generally nonzero quantity); the correct cyclicity form passes.
The mathematics was not adjusted to fit either error -- the
checks were.

## Scope

What is proved: the operator's finite algebraic content --
metric and dual bases, the Clifford module and its relations,
the spin lift of the isotropy action, the metric contraction
scalar, and the vanishing of the cubic term. What is cited:
Kostant's theorem (with Parthasarathy's formula) that for this
D, D^2 = C_g - C_{h, diagonal} + scalar (rho-norm constants,
in the conventions where the Casimirs and gamma use the same
invariant form). That identity lives on sections of the bundle
(the induced picture); it is not an identity on a finite
V (x) S module, and no such finite evaluation is claimed here.
Reproducing the scalar constant in the kappa = -tr
normalisation from the Killing-form statement is a rescaling
computation not performed. No FW/saw operator content: this is
the standard Dirac geometry of OG+(5,10), completing the layer
stack of the 16. Local file; not pushed.
