# The n = 6 fermion miniature — SU(4) = Spin(6) assembled and checked

Date: 2026-10-18. Status: exact finite-dimensional verification, COMPLETE.
Script: `proofs/fermion_miniature_n6.py` (local only, not uploaded).

Question (Kit, 2026-10-18): does SU(4) tie into the fermion problem?
This note assembles the standard Clifford miniature at dimension 6 into
one object and checks every piece exactly. It does not connect the
miniature to the FW/saw operator algebra; see the scope section.

## Construction

Cl(6) gamma matrices on C^8 by the Jordan–Wigner Pauli-monomial
construction (the so-ladder spinor engine's construction):
gamma_{2k-1}, gamma_{2k} carry sigma_x / sigma_y on tensor factor k
with a sigma_z string on the earlier factors. Spin generators
Sigma_ab = (1/2) gamma_a gamma_b (a < b). Chirality
Gamma = i * gamma_1 gamma_2 gamma_3 gamma_4 gamma_5 gamma_6
      = sigma_z (x) sigma_z (x) sigma_z in this basis.

Engines: desktop sympy, exact over Q(i) — 654/654 checks passed.
Cross-run on the laptop (WSL lampy, Sage Python route): identical
output, 654/654. All checks are exhaustive matrix identities over all
indices, not samples.

## Results (all exact)

A. Clifford and chirality.
{gamma_a, gamma_b} = 2 delta_ab I for all pairs; gamma_a Hermitian.
Gamma Hermitian, Gamma^2 = I, tr Gamma = 0, {Gamma, gamma_a} = 0,
[Sigma_ab, Gamma] = 0 for all 15 pairs. Gamma sectors are 4 + 4:
Gamma = +1 on basis indices {0,3,5,6}, -1 on {1,2,4,7} (occupation
parity in this basis).

B. The two blocks are su(4) and its conjugate.
Each Gamma block of the 15 Sigma_ab: anti-Hermitian, traceless,
real rank 15, closed under the bracket (all 225 ordered pairs per
block), Killing form exactly -8 * delta_ab in the Sigma basis
(signature (0,15), the compact form), Casimir sum_a Sigma_a^2
= -15/4 * I_4 per block. Controls against the 2026-10-13 record
(`proofs/so6_so5.md`): Dirac Casimir -15/4 * I_8 and vector Casimir
-5 * I_6 both reproduced in this normalization.
Conjugacy: the symmetric triple trace (d-symbol) on the triple
(Sigma_12, Sigma_34, Sigma_56) is -i for the Gamma = +1 block and
+i for the Gamma = -1 block. Opposite signs: the two blocks are
complex-conjugate, non-equivalent representations — the 4 and the
4bar. (Which block carries which name is fixed by the orientation
choice — gamma ordering plus the sign of Gamma; reversing the
orientation exchanges the names. The invariant content is the
conjugacy and the sign split.)

C. Lambda^2 of the 4 is the vector 6, as a real representation.
The induced action on Lambda^2 C^4 has rank 15, is anti-Hermitian,
lies in so(6, C) w.r.t. the wedge form Q (alpha wedge beta =
Q(alpha, beta) vol), and commutes with the Hodge real structure J
(J^2 = +I): the 6 is a real SO(6) representation, as required.

D. Fock form / statistics.
With the declared pairing (gamma_1,gamma_2), (gamma_3,gamma_4),
(gamma_5,gamma_6), set a_k = (gamma_{2k-1} + i gamma_{2k})/2 and
a_k^dagger its adjoint. Then exactly:
{a_k, a_l^dagger} = delta_kl I, {a_k, a_l} = {a_k^dagger, a_l^dagger} = 0
(so (a_k^dagger)^2 = 0 — Pauli nilpotency as an identity);
N_k = a_k^dagger a_k are commuting projectors; and
Gamma = prod_k (1 - 2 N_k) = (-1)^(N_1+N_2+N_3):
chirality IS occupation parity in this polarization. The common
kernel of a_1, a_2, a_3 is 1-dimensional (unique vacuum); the 8 Fock
states (a^dagger)^S |vac> are an orthonormal basis of C^8, with the
even-occupation states spanning the Gamma = +1 sector and the odd
the Gamma = -1 sector.

E. Double cover signs.
For the 2pi rotation in the 12-plane (via the exact identities
exp(tM) = cos(t/2) I + 2 sin(t/2) M when M^2 = -I/4, and
exp(tG) = I + sin t G + (1 - cos t) G^2 when G^3 + G = 0):
exp(2 pi Sigma_12) = -I_8 on the spinor, -I_4 on the 4-block,
+I_6 on Lambda^2(4), +I_6 on the vector representation.
The nontrivial kernel element of the covering acts as -1 on the
spinor and +1 on everything built from tensor powers of the vector
side checked here.

## Scope — what this does and does not establish

Established (PROVED, exact): at dimension 6 there is one proved
object in which the double cover, the chirality split 4 (+) 4bar,
and Fock statistics with Pauli nilpotency coexist, and its algebra
is su(4) = so(6), already on record since 2026-10-13. Fermion gap 2
(statistics) and the mechanism behind gap 3's sign now have a
complete miniature; the chirality grading of gap 4 has a realized
standard-Clifford instance.

Not established: any connection to the FW/saw operator algebra.
No action of the operator algebra on this C^8 (or on V_32) is
constructed — OQ1a stands. The differential Dirac operator (SC-5,
gap 1), masses (gap 5), and the seam twist prescription (gap 3 as
posed in the program) are untouched: this miniature is internal
algebra and representation theory, with one declared input, the
gamma pairing (polarization) of section D. WALL-FERMION stays OPEN.
The identification of the miniature's spinors with physical
fermions remains ASSERTED, per algebra-first discipline.
