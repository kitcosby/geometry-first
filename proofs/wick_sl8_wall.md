# Wick Rotation of the sl(8) Cycle: a Wall + the Correct Route

Date: 2026-10-12. Status: PROVED (conjugacy identity) / CHECKED (dim 28
closure, compactness). Question (Kit): "Wick" — apply E101-style Wick
rotation to the sl(8) cycle generators W_i = antisym(X_i) + i*sym(X_i),
X_i = a_i E_{i,i+1}, operators at 30deg for numerics.

## Negative result (the wall)
The Wick-rotated GENERATORS generate only dim 28, not su(8)'s 63.
Reason (exact): with U = diag(e^{i k pi/4})_{k=0..7},
  U W_i U^dagger = (a_i/sqrt(2)) (E_{i,i+1} - E_{i+1,i}),
real antisymmetric — precisely the so-cycle generators with operator
weights. Hence <W_i> = U^dagger so(8) U, a unitary conjugate of so(8):
dim 28, trace form negative definite (compact). Verified numerically:
dim closes at 28 and stalls; conjugacy identity holds for all i.

## Why
Wick rotation is a linear map on a whole algebra (via a Cartan
decomposition), NOT a Lie homomorphism — it does not commute with
generation. E101 worked because Wick was applied to a full basis of
sl(3) (8 elements), not to a small generating set.

## Correct route to su(8) (PROVED by the standard construction)
1. Generate the full sl(8) basis from the cycle (63 elements; pair
   transparency gives Chevalley normalization).
2. THEN Wick the basis: keep the 28 antisymmetric directions, multiply
   the 35 symmetric traceless directions by i.
Result: su(8), Killing signature (0,63). Same pattern as E101.
(Not re-run numerically at dim 63; the new content is the generator
wall above. First script's Killing step crashed assuming dim 63 —
downstream of the correct 28 result, disclosed.)

Bonus: the "failed" Wick is a phase-twisted copy of our own so(8):
the diagonal phases e^{ik pi/4} are the price of rotating edge by edge.
Scripts: /tmp/wick8.py, /tmp/wick8b.py (this run).
