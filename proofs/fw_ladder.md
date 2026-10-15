# The FW-Seeded Ladder: urx/uxp branches joined by FlatWave, rungs 3–8

**Date:** 2026-10-14 · **Status: PROVED (Sage-exact) at every rung, both variants.**
**Script:** `proofs/fw_ladder.py` (byte-identical to the executed file; see
Reproduction). **Engine:** SageMath 9.5, exact arithmetic in QQbar.

## 1. The ask

Kit's construction (main chat, 2026-10-14): *"Can we use the urx matrix and the
uxp matrix as the branches, joining with the flatwave matrix to oscillate, and
climb the ladder six more times to the full 8"* — refined by Kit to *"Perhaps
with the double angles in the diagonals."* Staged run, confirmed by Kit:
rungs n = 3 … 8, targets sl(n, R), the FlatWave matrix family (E96/E98) as the
seed pattern, the urx and uxp branches, FlatWave as the joining entry, and —
in the second variant — double-angle functions on the diagonals.

What is at stake is **generation by this specific pattern** (the
self-propagation statement for these families), nothing more: the existence
of sl(8) by Kit's operator cycle is already PROVED (cycle theorem), and the
universal "any traceless diagonal" statement remains OPEN in its E104 scope.

## 2. The two constructions (fixed before the run)

Both ladders share the E96/E106 path skeleton at size n:

- superdiagonal (i, i+1): `srx` if i odd, `crx` if i even;
- subdiagonal (i+1, i): `cxp` if i odd, `sxp` if i even;
- all other off-diagonal entries 0 except the two corners below.

**Ladder FW** (E96 generalized). Diagonal `(urx, uxp, FW, 0, …, 0)`;
corners `(1,n) = cos 2x`, `(n,1) = sin 2x`; then the traceless correction
`−(tr/n) I` (mandatory — the E95 raw matrix is not traceless and fails
closure). At n = 3 this is exactly E96. FlatWave appears only as a matrix
entry, never in auxiliary slots (E80 discipline).

**Ladder DA** (the double-angle variant). Diagonal
`d_i = sin 2x` (i odd) / `cos 2x` (i even) for `i < n`,
`d_n = −Σ_{i<n} d_i` — traceless by construction, no correction needed
(at n = 3 the diagonal is exactly `(sin 2x, cos 2x, −sin 2x − cos 2x)`);
corners `(1,n) = (n,1) = FW`.

Operators from signed sin/cos (valid in QI and QII, all samples have
sin x > 0): `srx = 1/|sin x| + cos x/sin x`, `cxp = 1/|cos x| + sin x/cos x`,
`sxp = 1/srx`, `crx = 1/cxp`, `urx = srx − crx`, `uxp = cxp − sxp`,
`sin 2x = 2 sin x cos x`, `cos 2x = cos²x − sin²x`, `FW = sgn(sin 2x)`.

## 3. Samples, domain, method

- **Samples:** x ∈ {15°, 30°, 45°, 60°, 75°, 105°, 120°, 135°, 150°, 165°} —
  quadrant I **and** quadrant II, so the FlatWave join is sampled at both
  signs (+1 below 90°, −1 above): the oscillation is part of the test.
  Exact values in Q(√2, √3).
- **Self-checks (asserted exactly at every sample, every rung):**
  `FW = sgn(sin 2x)` and `urx · uxp = 4/|sin 2x|` (the PROVED identities).
- **Domain:** x ∈ D = R ∖ {kπ/2}; sin 2x ≠ 0 for all operator formulas.
  Degeneracy angles named: `sin 2x = cos 2x` at x = 22.5° + k·45° (DA
  diagonal entries collide; not sampled); `cos 2x = 0` at 45° + k·45°
  (45° *is* sampled deliberately, as in E96, as a boundary probe);
  the boundary points sin 2x = 0 are excluded.
- **Closure:** exact incremental RREF span basis over flattened matrices
  (the Pascal TSpanBasis algorithm); repeated passes over basis pairs until
  a full pass adds nothing — that fixed point is the Lie algebra generated
  by the samples.
- **Identification:** every family matrix is traceless at *every* x in D
  (FW ladder by the correction, DA ladder by construction) and real, so
  ⟪family⟫ ⊆ sl(n, R); the samples' closure is contained in ⟪family⟫; hence
  closure(samples) = n² − 1 forces ⟪family⟫ = sl(n, R) (E91 dimension
  precedent).
- **Killing signature:** Gram `2n·tr(B_i B_j)` on the closure basis, exact
  LDL pivot signs (no zero pivot occurred at any rung), compared against
  the split signature of sl(n, R), `((n(n+1)/2 − 1), n(n−1)/2, 0)`.

## 4. Results — all 12 rungs PROVED

At every rung: E96-subset span (15°,30°,45°,60°) = 4; 10-sample span = 7;
all closure basis matrices traceless; Killing signature = the split
signature. Zero stalls, zero failures.

| n | FW closure | FW Killing | DA closure | DA Killing | target | elapsed FW / DA |
|---|---|---|---|---|---|---|
| 3 | 8  | (5, 3, 0)   | 8  | (5, 3, 0)   | 8  | 0.1s / 0.1s |
| 4 | 15 | (9, 6, 0)   | 15 | (9, 6, 0)   | 15 | 0.3s / 0.3s |
| 5 | 24 | (14, 10, 0) | 24 | (14, 10, 0) | 24 | 1.1s / 1.0s |
| 6 | 35 | (20, 15, 0) | 35 | (20, 15, 0) | 35 | 3.2s / 4.1s |
| 7 | 48 | (27, 21, 0) | 48 | (27, 21, 0) | 48 | 8.5s / 11.1s |
| 8 | 63 | (35, 28, 0) | 63 | (35, 28, 0) | 63 | 28.3s / 24.5s |

**Theorem (PROVED, Sage-exact).** For each n = 3, …, 8, the samples of
Ladder FW generate sl(n, R), and the samples of Ladder DA generate
sl(n, R); by the containment argument of §3, each *family* (all x ∈ D)
generates sl(n, R).

## 5. Cross-checks against the existing record

- **FW rung 3 = E96**, reproduced from scratch: subset span 4, closure 8,
  Killing (5, 3) — the machinery's sanity gate, passed before climbing.
- **Rung 8 Killing (35, 28, 0)** equals the sl(8) Killing signature of
  Kit's cycle theorem, `(35, 28)`: two different operator constructions
  of sl(8) — the cornered cycle and the urx/uxp-branch ladder — agree.
- Rungs 4–6 Killing signatures match the block-ladder record exactly
  (E88 (9,6), E90 (14,10), E91 (20,15)).
- The 10-sample span is **7 at every rung**: the sampled matrices always
  live in a 7-dimensional subspace, and the brackets do the rest of the
  climbing. (Computed datum, stated as such.)

## 6. Scope — what is NOT proved here

- The **universal** statement "any traceless diagonal works with this
  skeleton" is NOT proved and stays in its E104 conditional scope. This
  run proves two *specific* diagonal choices per rung (the urx/uxp/FW
  diagonal, and the double-angle diagonal), at n = 3…8.
- No claim for n > 8, other skeletons, other corner choices, or other
  diagonals.
- The two ladders give the *same* algebra at every rung; as with the
  cycle/path constructions, the algebra does not distinguish the
  operator content — what is proved is that this operator pattern
  suffices at every rung, not that it is distinguished by the result.
- No physical identification of any kind is made (algebra-first
  discipline). The E8 wall is untouched.

## 7. Disclosed tooling failures (run log)

Two Sage paths on this install are broken and were routed around; no
result depends on either:

1. `number_field_elements_from_algebraics` (QQbar → number-field
   conversion) fails inside its own homomorphism construction
   (`ValueError: relations do not all (canonically) map to 0 …`,
   Sage 9.5). First version of the script died there, after the sample
   self-checks had passed.
2. `NumberField(polynomial, embedding=…)` imports
   `sage.rings.polynomial.real_roots`, which hits the install's broken
   numpy (`ValueError: numpy.dtype size changed …`). Second version
   died there.

Resolution: all arithmetic performed in QQbar directly (exact; real
elements compare exactly, which covers the LDL pivot signs). The
occasional numpy ValueError also appears as non-fatal noise in the Sage
startup of the run logs; it comes from the same broken module import
chain and does not touch the QQbar/matrix code paths used.

## 8. Reproduction

Each rung was executed as the two-line prefix `KIND = '<FW|DA>'`,
`N = <n>` followed by `proofs/fw_ladder.py`, piped to Sage 9.5
(stdin mode) on the laptop WSL `lampy` distro. The file ends with a
top-level `main()` call (required: a piped Sage script whose last
statement is a compound block may silently never execute it).

## Addendum 2026-10-15 — extension to n = 9 and n = 10

At Kit's order, the ladder was extended by two rungs. Same
script (`proofs/fw_ladder.py`), same 10 samples, same method and
identification criteria as §§1–6; the run configuration is the
two-line prefix of §8 with `N = 9` and `N = 10`. All four rungs
PROVED (exact closure over QQbar):

| n | Variant | Sample span | Closure | Killing (pos, neg, zero) |
|---|---|---|---|---|
| 9 | FW | 7 | **80** | (44, 36, 0) |
| 9 | DA | 7 | **80** | (44, 36, 0) |
| 10 | FW | 7 | **99** | (54, 45, 0) |
| 10 | DA | 7 | **99** | (54, 45, 0) |

Every closure equals $n^2 - 1$, every closure basis is entirely
traceless, and every Killing signature equals the split form
$\big(\tfrac{n^2+n-2}{2}, \tfrac{n^2-n}{2}, 0\big)$ exactly.
Run times 103–221 s per rung.

The ladder now stands at **n = 3..10 for both variants (16
rungs), all PROVED**: both diagonal families generate
$\mathfrak{sl}(n,\mathbb{R})$ at every rung through
$\mathfrak{sl}(10,\mathbb{R})$. Scope unchanged from §1: the two
specific diagonals only — the universal E104 statement stays
OPEN, and no claim is made for $n > 10$.
