# The Spin(10) Spinor Representation: 32 = 16 + 16'

Date: 2026-10-13. Status: PROVED (explicit exact matrices, Sage 9.5
over QQbar). Standard Clifford construction; operator tie-in stated
at the end.

## Construction
Gamma matrices for Cl(10) by the Pauli recursion: start
gamma_1 = sigma_1, gamma_2 = sigma_2 (2x2); step from Cl(2k) to
Cl(2k+2) by g -> g (x) sigma_3 for old gammas, plus I (x) sigma_1,
I (x) sigma_2. After 5 steps: 10 gamma matrices, 32x32, entries in
{0, +/-1, +/-i}.
Spin generators: Sigma_ab = (1/4)[gamma_a, gamma_b], a < b.
Chirality: Gamma = i * gamma_1 gamma_2 ... gamma_10.

## Verified exactly
- {gamma_a, gamma_b} = 2 delta_ab I_32, all 100 pairs.
- The 45 Sigma_ab are linearly independent (rank 45).
- rho: so(10) -> End(C^32), e_ab |-> Sigma_ab, is a Lie algebra
  homomorphism: [rho(u), rho(v)] = rho([u, v]) verified on all 990
  basis pairs.
- Gamma^2 = I, Gamma Hermitian, tr Gamma = 0; chirality eigenspaces
  have dims 16 and 16.
- In the chirality eigenbasis every Sigma_ab is block-diagonal:
  the Dirac spinor 32 splits as 16 (+) 16' under so(10) — the two
  Weyl spinor representations of Spin(10).

## Operator tie-in
Kit's path generators Y_i = f_i (E_{i,i+1} - E_{i+1,i}) map to
rho(Y_i) = (f_i/2) gamma_i gamma_{i+1}; all 9 (weights at 30 deg)
verified block-diagonal in the chirality basis with nonzero
16-blocks. The operators weight the spinor action through the
vector-rep generators; the spinor representation itself is the
standard Cl(10) one.
Script: /tmp/spin10.py (this run; ~15 min wall, dominated by the
990 homomorphism pair checks over QQbar).
