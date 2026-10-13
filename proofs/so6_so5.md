# so(6) and so(5): Paths, Spinors, Exceptional Isomorphisms

Date: 2026-10-13. Status: PROVED (Sage 9.5 exact, QQbar). Same
battery as so(10)/so(8).

## so(6) (D3 = A3)
- Path 6x6 (5 edges, operator weights at 30 deg): closure dim 15,
  basis antisymmetric, closed => so(6).
- Cl(6): Sigma rank 15, rho a homomorphism; chirality Gamma = i*w
  (w = g1...g6 has w^2 = -1 at n=6), Gamma^2 = I, tr 0, eigenspaces
  4 + 4; Sigma block-diagonal.
- The 4-block: rank 15, all traceless, all anti-Hermitian => image
  inside su(4) with dim 15 = dim su(4): the spinor 4 IS the defining
  representation of su(4). so(6) = su(4) realized explicitly.
- Casimir control: vector -5 I, spinor -15/4 I. NOT equal — the
  triality coincidence of n=8 is absent, as the formula predicts.

## so(5) (B2 = C2)
- Path 5x5 (4 edges): closure dim 10, antisymmetric, closed => so(5).
- Cl(6) gammas restricted to the first 5: Sigma rank 10, rho a
  homomorphism; volume w5 = g1...g5 has w5^2 = +1 and commutes with
  every Sigma; eigenspaces 4 + 4; Sigma block-diagonal; each block
  rank 10.
- Symplectic test (so(5) = sp(4)): solving Sig^T J + J Sig = 0 for
  antisymmetric J on each 4-block: solution space dim exactly 1,
  J invertible (det != 0), preserved by all 10 generators. The
  spinor 4 lands in sp(4) (dim 10) and fills it: so(5) = sp(4)
  realized explicitly.
- Casimir control: vector -4 I, spinor -5/2 I. NOT equal.

## Disclosures
- The symplectic solve took three attempts: (1) equation system
  built transposed (meaningless 55-dim kernel, discarded);
  (2) silent no-run — Sage reading a script from stdin silently
  skips a trailing compound block at EOF; fixed by ending the
  script with a top-level print. Both disclosed in chat.
Scripts: /tmp/so65.py, /tmp/so5b.py, /tmp/so5d.py (this run).
