# so(7), so(4), so(3): The Orthogonal Ladder Completed

Date: 2026-10-13. Status: PROVED (Sage 9.5 exact, QQbar).

## so(7) (B3)
- Path 7x7 (6 edges): closure dim 21, antisymmetric, closed.
- Spinors from Cl(8) gammas restricted to 7: Sigma rank 21,
  homomorphism; volume w7 commutes with all Sigma, (i w7)^2 = I,
  eigenspaces 8 + 8, Sigma block-diagonal, each block rank 21:
  the spinor 8 (two copies from the two Cl(7) modules).
- Type test: the spinor 8 carries an invariant SYMMETRIC bilinear
  form, unique (solution dim 1), invertible, preserved by all 21
  generators — the spinor of B3 is of real (orthogonal) type.
- Casimir: vector -6 I, spinor -21/4 I. Not equal.

## so(4) (D2 = A1 + A1, not simple)
- Path 4x4 (3 edges): closure dim 6, antisymmetric, closed.
- Hodge star on so(4) (via the Levi-Civita tensor): star^2 = I,
  tr 0, eigenspaces 3 + 3 (self-dual / anti-self-dual).
- Split verified: [+,+] in +, [-,-] in -, [+,-] = 0 identically,
  both factors with nonzero brackets: so(4) = so(3) (+) so(3).
- Spinors (Cl(4), chirality Gamma = g1g2g3g4 — note omega^2 = +1
  at n=4, no factor i): eigenspaces 2 + 2; every self-dual
  generator acts ONLY on the second chiral block, every
  anti-self-dual ONLY on the first: spinor = (2,1) (+) (1,2).

## so(3) (A1)
- Path 3x3 (2 edges): closure dim 3, antisymmetric, closed.
- Spinors from the Pauli matrices: Sigma rank 3, traceless
  anti-Hermitian, homomorphism: image = su(2), the spinor 2 is the
  defining representation.

## Disclosure
First so(4) spinor attempt used Gamma = i*omega (correct at
n = 6, 10 where omega^2 = -1; wrong at n = 4 where omega^2 = +1):
chirality eigenspaces came out empty (0, 0) and the run crashed
downstream. Corrected to Gamma = omega; factorization then
verified. The algebra-level split was unaffected.

## The completed ladder (all PROVED in this campaign)
n:    3   4     5     6     7   8        9   10
dim:  3   6     10    15    21  28       36  45
spin: 2   2+2   4     4+4'  8   8s+8c    -   16+16'
note: su2 split sp4   su4   sym triality path Spin(10) 16
(so(9) spinor not built; path only. Casimir coincidence of the
three fundamental reps occurs only at n = 8.)
Scripts: /tmp/so743.py, /tmp/so4spin.py (this run).
