# The FW-seeded ladder: an induction attempt (all-n)

**Date:** 2026-10-16 · **Status: INCOMPLETE — reduced to one named lemma.**
Local note, written at Kit's invitation ("try your hand at induction").
Companion to `proofs/fw_ladder.md` (rungs n = 3..10 PROVED Sage-exact,
both variants; n = 11 CHECKED in GF(p) at two primes, both variants,
2026-10-16). Target: the all-n theorem for the two specific diagonal
families. The universal E104 statement is untouched and stays OPEN.

## 1. The reformulation (today's synthesis applied)

WALL-ARG/VAL (proved 2026-10-16, scoped) says no argument-side
rescaling argument can prove generation; the proof must live on the
value side. The value-side normal form of the ladder family is a
**coefficient decomposition**: the family matrix is a linear
combination of seven *constant* matrices with operator-function
coefficients. Both variants use the same seven functions
`{srx, crx, cxp, sxp, FW, cos 2x, sin 2x}`:

- Ladder FW: `M(x) = srx·A₁ + crx·A₂ + cxp·A₃ + sxp·A₄ + FW·A₅ +
  cos 2x·A₆ + sin 2x·A₇`, where A₁..A₄ carry the alternating path
  (odd/even super- and subdiagonal) plus the diagonal parts of
  `(urx, uxp, FW, 0, …, 0)` and the `−(tr/n) I` correction split by
  coefficient; A₅ is the FW diagonal part; A₆ = E_{1n}; A₇ = E_{n1}.
- Ladder DA: same path matrices A₁..A₄; the sin 2x / cos 2x diagonal
  coefficient matrices and the FW corner matrix in place of
  A₅, A₆, A₇. (Exact tables reproducible from `fw_ladder.py`'s
  `build()` by collecting coefficients.)

## 2. Lemma A (PROVED) — the coefficient space is the sample space

The seven functions are linearly independent over ℚ, and the ten
samples span the same space as the family. Proof: the 7×7 evaluation
determinant at degrees {15, 30, 45, 60, 75, 105, 120} (a subset of
the sample set, both FW signs included), computed in GF(p) with
ℚ(√2, √3) embedded, is **205676 ≢ 0 (mod 1000033)** and **797465 ≢ 0
(mod 999983)**; a rational determinant divisible by neither prime is
nonzero. Hence the evaluation matrix has rank 7 over ℚ, the functions
are independent, and span(samples) = span(family) = span{A₁..A₇}.

Verified alongside (exact rational arithmetic, desktop):
rank{A₁..A₇} = 7 and Lie closure{A₁..A₇} = n² − 1 at n = 4, 5, 6
(Ladder FW) and closure = 24 at n = 5 (Ladder DA) — agreeing with the
Sage rungs, through an independent formulation.

**Consequence.** The all-n theorem is equivalent to: *for every
n ≥ 3, the seven explicit constant matrices A₁(n)..A₇(n) Lie-generate
sl(n, ℝ).* No functions remain in the statement.

## 3. Why block induction fails (obstruction, identified exactly)

A(n) does not extend A(n−1) as a block: the traceless correction
`−(tr/n) I` puts `1/n` in every diagonal coefficient, and Ladder DA's
last diagonal `d_n = −Σ_{i<n} d_i` likewise depends on n. The top-left
(n−1)-block of A(n) is not A(n−1). Naive "embed the previous rung and
add connectors" induction is therefore unavailable — the ladder is a
**cycle** (path closed by the two corners), not a path with a tail.
The induction must be on *distance around the cycle*, in the style of
Kit's sl(8) cycle theorem, not on the rung index by blocks.

## 4. The corner walk — first steps (VERIFIED exactly, n = 6, rationals)

The two corners close the cycle, and bracketing them against the
coefficient matrices walks roots in from both ends. Exact brackets at
n = 6 (1-based positions; the identities are position-local):

- `[A₆, A₇] = E₁₁ − E₆₆` — the long diagonal H_{1n} (identity valid
  for all n by inspection: A₆ = E_{1n}, A₇ = E_{n1}).
- `[A₆, A₃] = E_{1,5} − E_{2,6}` — exactly the two roots adjacent to
  the corner (1, n): E_{1,n−1} (through the last subdiagonal edge)
  and −E_{2,n} (through the first). No other support: the diagonal
  part of A₃ has equal entries at positions 1 and n, so its bracket
  with the corner vanishes identically.
- `[A₇, A₁] = −E_{5,1} + E_{6,1} + E_{6,2}` — the mirror walk from
  (n, 1) along the superdiagonal: E_{n,2} produced, −E_{n−1,1}
  produced, plus a multiple of A₇ itself (absorbable).
- `[A₆, A₁] = −E_{1,6}`, `[A₆, A₂] = +E_{1,6}`, `[A₇, A₂] = −E_{6,1}`
  — corner multiples only (diagonal-part brackets); `[A₆, A₄] =
  [A₆, A₅] = [A₇, A₃] = [A₇, A₄] = [A₇, A₅] = 0` at n = 6.

Parity note: the (1, n) corner engages the odd subdiagonal chain
first; which of A₃/A₄ carries the last edge flips with the parity of
n. The first-step pattern — one bracket producing the two
corner-adjacent roots and nothing else — is n-uniform in form.

## 4b. Addendum (2026-10-16, later same day) — the first step is parity-split

Exact-rational check of the §4 brackets at every n = 5..10 (script:
`fw_ladder_corner_walk.py`, A₁..A₇ read off `fw_ladder.py`'s
`build()`; Fractions throughout, no Sage). This **corrects** §4's
closing sentence: the first-step pattern is NOT n-uniform in form —
it splits by parity, and the split is good news:

- `[A₆, A₇] = E₁₁ − Eₙₙ` at every n checked (as §4 said).
- **Even n** (the §4 case, n = 6): `[A₆, A₃] = E_{1,n−1} − E_{2,n}`
  — the pair, in one bracket; `[A₆, A₄] = 0`. Mirror corner:
  `[A₇, A₁] = −E_{n−1,1} + E_{n,1} + E_{n,2}`, `[A₇, A₂] = −E_{n,1}`.
- **Odd n**: the pair never forms. `[A₆, A₃] = −E_{2,n}` and
  `[A₆, A₄] = +E_{1,n−1}` — each corner-adjacent root arrives
  **individually**, because for odd n the last path edge is even, so
  its two sides live in different coefficient matrices. Mirror
  corner: `[A₇, A₁] = E_{n,1} + E_{n,2}` and
  `[A₇, A₂] = −E_{n−1,1} − E_{n,1}` — the roots separate modulo the
  corner multiple `E_{n,1} = A₇` itself (absorbable, as in §4).

Consequence for §5: the uniform triangular-separation gap is an
**even-n-only** gap at the first step. For odd n the separation the
gap asks for is already done by the A₃/A₄ (resp. A₁/A₂) split at
step one; what remains for odd n is only the front propagation.
Scope: VERIFIED exactly at n = 5..10, both parities represented
three times each; the parity mechanism (which coefficient matrix
carries the last edge) is visible in the construction, but a
written all-n parity argument is still part of the unwritten lemma
— this addendum does not upgrade the all-n statement, which stays
INCOMPLETE.

## 5. The remaining gap (named)

**Gap (uniform triangular separation).** The first step produces the
*sum/difference pair* E_{1,n−1} − E_{2,n}, not the individual roots:
ad(H_{1n}) gives both roots the same weight (+1) and does not
separate them; the separating diagonals (e.g. the diagonal part of
A₁ − A₂ ≈ 2E₁₁ modulo center) act *triangularly* — their brackets
also move roots along the walk (e.g. produce E_{2,n−1} terms). What
remains to be written is the triangularity argument: an ordering of
the cycle's root positions under which iterated brackets with
{A₁ − A₂, A₃ − A₄-type differences, H_{1n}} peel the pair apart and
propagate the two fronts (from (1, n) and from (n, 1)) until they
meet, for all n and both parities, followed by the classical
path/cycle lemma (root vectors on all edges + H's ⇒ all E_{ij} by
distance induction ⇒ sl(n)). The classical lemma is standard; the
separation bookkeeping specific to these seven matrices is the open
work. Estimated: a few pages of careful bracket accounting, no new
ideas required beyond the walk of §4 — but it is NOT written, and
until it is, the all-n statement stays INCOMPLETE, with rungs
n = 3..10 PROVED, n = 11 CHECKED (GF(p)), both variants.

## 6. Scope

- Nothing here touches E104's universal diagonal claim, other
  skeletons, or any physical identification.
- The determinant of Lemma A and the rational closures of §2 are
  exact computations, reproducible from `fw_ladder.py`'s definitions
  (computation scripts were run on the desktop, 2026-10-16; method:
  Tonelli–Shanks roots of 2, 3 in GF(p); Fraction-based RREF closure,
  the Pascal TSpanBasis algorithm).
