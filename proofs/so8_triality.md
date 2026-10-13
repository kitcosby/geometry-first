# so(8): Path, Spinors, and the Triality Test

Date: 2026-10-13. Status: PROVED (Sage 9.5 exact, QQbar). Kit's aim
was falsification of so(8)'s specialness; the attempt failed —
see the verdict at the end.

## Baseline (forced by the length-generic induction)
- Antisymmetric path, 8x8, 7 edges weighted by the first 7 operators
  at 30 deg: closure dim 28, closed. (= so(8), D4, rank 4.)
- Cl(8) by the Pauli recursion: 8 gamma matrices 16x16,
  anticommutation exact; Sigma_ab rank 28; rho a Lie homomorphism on
  all 378 basis pairs.
- Chirality Gamma = gamma_1...gamma_8 (no factor i needed at n=8):
  Gamma^2 = I, Hermitian, tr 0, eigenspaces 8 + 8; all Sigma_ab
  block-diagonal: Dirac 16 = 8_s (+) 8_c.

## The triality test (the falsifiable edge)
Quadratic Casimir scalar C(rho) = sum_{a<b} rho(e_ab)^2:
  vector rep:  -(n-1) I        (each index in n-1 pairs)
  spinor reps: -n(n-1)/8 I     (C(n,2) pairs x Sigma_ab^2 = -1/4 I)
Equal iff n(n-1)/8 = n-1 iff n = 8. Verified exactly at n=8:
  8_v: -7 I8;  8_s block: -7 I8;  8_c block: -7 I8.
All three 8-dimensional representations carry the SAME scalar —
the Casimir fingerprint of triality (equal Dynkin indices). For
every other n the spinor scalar differs (e.g. so(10): vector -9,
spinor -45/4). The test could have failed; it did not.

## Verdict (two-sided, scoped)
- so(8)'s specialness: CONFIRMED as algebra. Triality's numerical
  fingerprint holds exactly.
- Operator content in that specialness: NONE. The operators entered
  only as nonzero path weights; any nonzero weights (integers
  included) give the same result. Triality is a property of D4
  itself, not of the trig system. What is falsified is any claim
  that the operators MAKE so(8) special.
Script: /tmp/so8.py (this run).
