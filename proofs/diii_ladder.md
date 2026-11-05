# The DIII Ladder: Spin(2n)/U(n) for General n, Derived

Date: 2026-10-23. Status: the n = 5 stack (orbit, stabilizer,
bundle, Kostant) is the n-th rung of a uniform construction,
DERIVED below for all n >= 2 and VERIFIED exactly rung by rung
for n = 2..6: proofs/diii_ladder.py, 150 counted checks, 0
failures, exit 0 (~65 s). Per-rung scope in the table at the
end. The n = 5 rung reproduces the three filed builds
(pure_spinor_orbit, su5_stabilizer, spinor_bundle,
kostant_dirac) number for number.

Notation: nu = n(n-1)/2 throughout. The spinor is the Weyl
sector CONTAINING the Fock vacuum e_0 -- its Gamma sign is not
fixed across the ladder: Gamma e_0 = i^{n+1} e_0 (n odd,
Gamma = i prod) resp. i^n e_0 (n even, Gamma = prod), giving
signs -, +, +, -, - for n = 2, 3, 4, 5, 6 (period 4). The
script reads the sector off the matrix; the n = 5 builds'
"Gamma = -1 sector" is this same sector at n = 5.

## The uniform construction (the derivation)

1. Spinor. The Pauli recursion gives Cl(2n) on C^{2^n}; the
   a-basis polarises C^{2n} = W + W' with W, W' maximal
   isotropic; the vacuum e_0 is pure with annihilator W
   (dim n, Gram 0 -- the same two-route computation as n = 5,
   uniform in the recursion). Sector dim 2^{n-1}; spinor
   Casimir (kappa = -tr convention, Sigma_ab generators)
   -n(2n-1)/4.
2. Orbit / parabolic. The line stabiliser is the
   W-preserving algebra: block form [[A, B], [0, -A^T]], B
   skew; dim_C n^2 + nu; the spinor character is
   mu = (1/2) tr A (the vacuum is the top exterior power of
   the creation space, half-density corrected -- per-rung
   verified on the sparse basis). Orbit = OG+(n, 2n),
   dim_C nu.
3. Compact stabilizer. Reality <=> preserve W and W' <=>
   B = C = 0 and A anti-Hermitian: u(n), dim n^2. Centre:
   Z = -sum_k E_{2k,2k+1}, A(Z) = -i I_n. Derived algebra
   = ker mu = su(n), dim n^2 - 1.
4. Branching. The sector is Lambda^even(W'), W' carrying
   -A^T: rho(X)|_m = dLambda^m(-A^T) + mu I. Piece m:
   dim C(n, m), rho(Z)/i charge m - n/2. As SU(n) modules the
   pieces are Lambda^m(n-bar) (m even); Lambda^n = 1.
   At n = 5: 1 + 10-bar + 5, as filed.
5. Bundle data. m = {A-block 0}: real dim 2 nu;
   [h, m] subset m, [m, m] subset h (block multiplication;
   Hermitian symmetric, type DIII); kappa(h, m) = 0;
   tangent isotropy B |-> A B + B A^T.
6. Kostant layer. Cubic form phi = kappa([., .], .) = 0 on m
   for every n (from 5: [m, m] subset h, kappa(h, m) = 0), so
   D = sum_j M_j (x) gamma(M^{j*}) with no cubic correction at
   every rung. The Clifford/spin-lift statements
   ({gamma, gamma} = kappa, [sigma, gamma] = gamma o ad,
   tr sigma = 0, contraction = nu I) reduce, via the
   canonical anticommutation relations of the monomial model
   (pure combinatorics of wedge/contraction, n-independent),
   to two per-rung matrix inputs: ad|_n = dLambda^2(A) in the
   form B |-> A B + B A^T, and the contragredient identity
   N = -M^T. Both inputs are verified at EVERY rung,
   including n = 6.

## Per-rung verification (script sections)

Every rung, exact over Q(i): gamma controls; sector dim;
Casimir; sparse u(n) basis real and line-preserving with
mu = (1/2) tr A; parabolic dim by complex nullity; compact
stabilizer by real solve (dim, anti-Hermiticity, A-span rank
n^2); derived = ker mu (dim n^2 - 1, bracket rank and mu
annihilation); centre dim 1 and Z identified; charges and
piece dims; the branching operator identity on the sparse
basis, all pieces; the distinct-weight Cartan weights test
(w_k = 2k - n + 1, piece-m eigenvalues exactly
{-i sum_{k in S} w_k}); m by solve (dim 2 nu); both bracket
inclusions on full bases; kappa nondegenerate on m;
kappa(h, m) = 0; cubic phi = 0 on all triples; radical
commutator identity.

Rung table (all checks pass; counts exclude the RESULT line):
  n   so(2n)  spinor  nu   checks  engine (Lambda* n)
  2   so(4)   2       1    31      2 monomials, full
  3   so(6)   4       3    31      8 monomials, full
  4   so(8)   8       6    31      64 monomials, full
  5   so(10)  16      10   31      1024 monomials, full
  6   so(12)  32      15   26      not evaluated (32768)

## Scope and disclosures

- n = 6 engine: the monomial evaluation (gamma relations,
  sigma intertwining, contraction on 2^15 monomials) was NOT
  run -- the sparse engine's cost scales with the monomial
  count and the operator count together. Those three
  statements at n = 6 rest on the uniform CAR derivation in
  item 6, whose per-rung inputs ARE verified at n = 6. Status:
  PROVED by the derivation with inputs CHECKED; at n <= 5
  additionally verified by direct evaluation. This is the
  proof-over-sampling split, stated rather than smoothed
  over.
- First script run crashed at n = 2 before any check beyond
  the controls: Gamma = i prod is the wrong chirality
  operator for even n (prod^2 = (-1)^n). Fixed to the
  standard Gamma (i prod for n odd, prod for n even) and the
  vacuum sector read from the matrix. No n = 5 number was
  affected (n odd there).
- Group-level statements ((SU(n) x U(1))/Z_n as the
  stabilizer group, OG+ projective hence orbit closed) are
  the standard ones, cited as in the n = 5 notes; the
  computed content is Lie-algebra level throughout.
- The ladder is the DIII family only. The other pure-spinor
  families (e.g. the E_6/(Spin(10) x U(1)) stratum at the
  27, or quaternionic analogues) are different constructions,
  not rungs of this ladder. No FW/saw content. Local files;
  not pushed.
