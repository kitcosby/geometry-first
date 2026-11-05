# The Pure Spinor Orbit: Spin(10) on the 16 is OG(5,10)

Date: 2026-10-22. Status: PROVED (explicit exact matrices over Q(i),
34 counted checks, exit 0; the two standard inputs used are named in
the Scope section). Script: proofs/pure_spinor_orbit.py (sympy, this
machine; runs in ~3 s).

## Claim

Let 16 be the Spin(10) Weyl spinor containing the Fock vacuum e_0
(the Gamma = -1 block of the corpus Cl(10) construction,
spin10_spinor.md). Then:

1. e_0 is a pure spinor: its annihilator
   W = { v in C^10 : (v . gamma) e_0 = 0 }
   is a maximal isotropic 5-plane, W = span of the five annihilation
   directions a_k = (gamma_{2k} + i gamma_{2k+1})/2.
2. The stabilizer of the line C e_0 in so(10, C) is the maximal
   parabolic of the isotropic 5-plane: dim 35 = Levi gl(5) (25) +
   unipotent radical Lambda^2 (10), in the block form
   [[A, B], [0, -A^T]], B skew, relative to C^10 = W (+) W'.
   Hence the orbit has dim 45 - 35 = 10.
3. The orbit of [e_0] in P^15 is OG^+(5,10), the component of the
   orthogonal Grassmannian of maximal isotropic 5-planes in the
   family of W.
4. The other Weyl spinor's vacuum e_31 (Gamma = +1 block, the 16')
   is also pure, with annihilator W' spanned by the creation
   directions a_k^+. Since dim(W n W') = 0 is even while n = 5 is
   odd, W' lies in the opposite family: the 16' orbit is OG^-(5,10).
5. Line vs vector stabilizer: the parabolic acts on the vacuum line
   by the character mu = +(1/2) tr(A) (A = the Levi block on the
   annihilator space W; the often-quoted -1/2 is the same statement
   parametrized on the creation space). The stabilizer of the vector
   e_0 itself is ker mu, dim 34 (Levi sl(5), 24, + radical 10).

## Verified exactly (script sections)

- Checker controls first: {gamma_1, gamma_2} = 0 accepted (true
  control); gamma_1 gamma_2 - gamma_2 gamma_1 rejected as 0 (false
  control). Clifford relations on all 100 pairs; Gamma^2 = I,
  tr Gamma = 0, sectors 16 + 16, vacuum in Gamma = -1, top vector
  in Gamma = +1.
- Casimir on BOTH Weyl blocks = -45/4 I_16 (the layer-1 value of
  2026-10-22, now computed on the nose for so(10), not only by
  the n(n-1)/8 formula of so8_triality.md).
- Purity, two independent routes each: dim ker(v |-> (v.gamma) e_0)
  = 5 by nullspace, and the five explicit a_k coefficient vectors
  lie in that kernel and annihilate e_0 directly; Gram matrix of W
  is exactly 0. Same for e_31 with the a_k^+ vectors. Pairing
  W x W' has Gram (1/2) I in the explicit a-bases; W (+) W' = C^10
  (rank 10, intersection 0).
- False control for purity: psi_c = e_0 + e_{11110} has annihilator
  dim 1, not 5 -- the test discriminates.
- Family parity sanity: for the graph W_S of the skew map
  S = E_01 - E_10 (W -> W'), dim(W n W_S) = 3, odd, i.e. same
  family as W -- consistent with the criterion used in item 4.
- Stabilizer, solved twice: spin side (X e_0 in C e_0 inside the
  16) gives nullity 35; vector side (Y W subset W) gives nullity
  35; rho carries the vector-side algebra into the line stabilizer
  (all 35 basis elements), so the two coincide under rho.
- Levi form on all 35 vector-side basis elements in the adapted
  basis [W | W']: lower-left block 0, D = -A^T, B = -B^T;
  A-projection rank 25 (onto gl(5)), B-projection rank 10
  (onto Lambda^2).
- Character: mu = +(1/2) tr(A) on all 35 basis elements;
  dim ker mu = 34 on the spin-side basis.

## Argument (orbit = OG^+)

W_psi is isotropic for every psi (v, w in W_psi implies
2 B(v,w) psi = {v.gamma, w.gamma} psi = 0), so dim W_psi <= 5 with
equality iff psi is pure. The Clifford action intertwines the
vector representation (rho of spin10_spinor.md), so
W_{g psi} = g W_psi: the orbit of [e_0] maps Spin(10, C)-
equivariantly into the isotropic Grassmannian. Its stabilizer is
the parabolic computed above, so the orbit is Spin(10, C)/P, a
projective variety, hence closed in P^15, of dim 10. By the pure
spinor correspondence (Scope) the pure spinor lines of the 16 are
exactly OG^+(5,10), irreducible of dim 10; a closed irreducible
10-dimensional subvariety of it is all of it.

## Disclosures (checker expectations corrected en route)

- First run, one FAIL: the W x W' Gram check demanded (1/2) I from
  sympy's nullspace bases, which are not the a_k bases; the matrix
  is basis-dependent. Split into det != 0 (nullspace bases) and
  (1/2) I (explicit a-bases): both pass. The claim was never at
  issue; the check was.
- First character run, one FAIL: expected mu = -(1/2) tr(A) from
  the standard Fock formula; the computation gives +(1/2) tr(A)
  uniformly (whenever tr A = i, mu = i/2; whenever tr A = 0,
  mu = 0), because A here acts on the annihilator space W.
  Functional identity verified on a basis, hence as functionals.
  The script and this note record the corrected sign.

## Scope

Computed here, exactly: every matrix statement above (purity,
stabilizers, block forms, character, Casimir, parity numbers).
Cited as standard inputs, not re-proved: (a) Chevalley's
correspondence between pure spinor lines of one Weyl spinor and
one family of maximal isotropic subspaces, with the family
criterion dim(W n W') = n mod 2; (b) a quotient G/P by a parabolic
is projective, hence its orbit in P(V) is closed. No operator
(FW/saw) content: this is the standard Spin(10) geometry of the 16,
built to sit under spin10_spinor.md. Local file; not pushed.
