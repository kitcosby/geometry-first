# so(9) Spinor: the Ladder's Gap Closed

Date: 2026-10-13. Status: PROVED (Sage 9.5 exact, QQbar).
Closes the gap flagged in so7_so4_so3.md ("so(9) spinor not
built"). The orthogonal ladder so(3)..so(10) is now complete,
algebras and spinors alike.

## Construction and verification
- Cl(9) from the first 9 gammas of the Cl(10) build (32x32):
  anticommutation exact.
- Sigma_ab (36 of them): rank 36; rho a Lie homomorphism on all
  630 basis pairs.
- Volume w9 = g1...g9: w9^2 = I, commutes with every Sigma,
  tr w9 = 0, eigenspaces 16 + 16; Sigma block-diagonal; EACH block
  has rank 36 — each is the spinor 16 of so(9) (the two blocks are
  the two Cl(9) modules, same so(9) representation).
- Real type: invariant symmetric bilinear form on the 16, solved
  on the distance<=2 generators: solution space dim exactly 1,
  invertible, verified preserved by all 36 Sigma. (Mod-8 pattern:
  n = 9 = 1 mod 8, real — matching so(7)'s symmetric form, against
  so(5)'s symplectic one.)
- Casimir control: vector -8 I9, spinor -9 I32. Not equal — and
  note the order has swapped past the crossing: n=7 gives
  -6 vs -21/4 (spinor smaller), n=8 gives -7 = -7, n=9 gives
  -8 vs -9 (spinor larger). The two scalar sequences meet only
  at n = 8.
Script: /tmp/so9spin.py (this run; ~12 min wall, dominated by the
630 homomorphism pairs).
