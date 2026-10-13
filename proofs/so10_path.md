# so(10) from the Antisymmetric Path (Spin(10) algebra)

Date: 2026-10-13. Status: PROVED (Kit's length-generic path induction;
Sage 9.5 exact cross-check over QQbar).

## Construction
Y_i = f_i (E_{i,i+1} - E_{i+1,i}), i = 1..9, in 10x10 matrices, with
weights the 8 operators plus FlatWave as the ninth edge:
f = (srx, sxp, cxp, crx, urx, sawup, uxp, sawdown, FW).
All f_i are nonzero on the admissible domain (FW = +/-1).

## Proof
Kit's so(9) path induction is length-generic: induction on the
distance k = j - i produces every E_{ij} - E_{ji} (i < j) as a
commutator word, for any path length and any nonzero weights.
Hence the Y_i generate the full space of antisymmetric 10x10
matrices: so(10), dim C(10,2) = 45, type D5, rank 5.

## Sage cross-check (x = 30 deg, exact surd weights)
- Commutator closure dim = 45.
- Every basis element antisymmetric (closure subset of so(10);
  dimension equality gives the whole algebra).
- Closed under all pairwise commutators.
- First 8 edges alone close at dim 36 = so(9) on the first 9
  coordinates: so(9) subset so(10) by construction.

## Note on Spin(10) vs SO(10)
The Lie algebra is the same for Spin(10) and SO(10); the algebra
result above is the complete algebra-level statement. The spinor
representation (16 of Spin(10)) is a separate construction via the
Clifford algebra Cl(10) — not claimed here.
Script: /tmp/so10.py (this run).
