# Blending intertwiner hunt — odd-candidate solve (steps a/b + c): RESULT = MISS (proved), with a unique partial HIT underneath

**Date:** 2026-10-17 · Run at Kit's order (queued 2026-10-17). Spec:
`proofs/blending_intertwiner_hunt_spec_DRAFT.md`, steps (a)–(c) under the
2026-10-16 addendum (impose {C, Ω} = 0 from the start). Step (d) was the
even screen: clean MISS, `proofs/blending_step_d_screen.md`.
**Outcome: no element of M(16, ℚ) satisfies conditions (1)–(4)
simultaneously. The obstruction is identified exactly and proved:
condition (2) (sector swap) and condition (3) (chirality flip) are
incompatible for any involution — the failing combination is (2)∧(3)
given (1). Underneath the wall sits a unique partial pattern: (1)∧(2)
is satisfiable by exactly one involution up to sign (§3).**

## 1. Setup, gates (PROVED — reproduction before search)

Realization rebuilt exactly as in `proofs/d2_structure.py`: e = (σx⊗I,
Y⊗I, σz⊗σx, σz⊗Y) in Cl(2,2); E_μ = e_μ⊗I₄ (μ = 0..3, R sector),
E_{4+ν} = Ω₄⊗e_ν (J sector), Ω₄ = e₀e₁e₂e₃; signature
(+,−,+,−,+,−,+,−); Ω = E₀E₁···E₇. Gates, all exact: gamma squares and
anticommutation; Ω² = I; Ω anticommutes with every gamma; the E6 record
([S₀₄,S₁₅] = [S₀₅,S₁₄] = 0, other mixed pairs nonzero, T² = 4I,
tr T = 0, T = −2E₀E₁E₄E₅); dim V = 16 recomputed (d2's value).
Parity gate: every Clifford monomial E_S anticommutes with Ω iff |S| is
odd and commutes iff |S| is even (all 256 monomials checked) — so the
odd candidate space is exactly the odd-grade subspace.

Method: coefficient space over the 256 monomials; left/right
multiplication by a gamma is a signed permutation of coefficients, so
each condition below is an exact rational linear system until the
single quadratic step (C² = I) on a 1-dimensional space.

## 2. The odd candidate space and the linear solve (PROVED)

Candidate spaces (dimensions over ℚ; V = one-direction intertwiners
{X : X E_μ = E_{4+μ} X}, W = the reverse {X : X E_{4+μ} = E_μ X} — with
C² = I imposed, membership in V∩W is exactly condition (2) in
intertwining form; O = odd part):

| space | dim | even part | odd part |
|---|---|---|---|
| O (odd part of Cl(4,4)) | 128 | — | 128 (monomials E_S, |S| odd) |
| V | 16 | 8 | 8 |
| W | 16 | 8 | 8 |
| V∩W (the (2)-space) | 1 | 1 | **0** |

The (2)-space V∩W was computed twice (W-conditions restricted to a V
basis, and independently V-conditions restricted to a W basis): dim 1
both ways. **Its odd part is 0: no odd element satisfies condition (2)
in both directions at all** — the odd solve's solution set is empty
already at the linear stage, before C² = I is even imposed. Odd
intertwiners exist one-directionally (odd parts of V and W are
8-dimensional each) but, by the lemma in §4, none of them can square
to +I.

## 3. The unique swap involution (partial pattern — satisfies (1), (2), (4)-data; fails (3) necessarily)

V∩W is spanned by a single element B₀, a sum of 16 grade-4 monomials
with coefficients ±1:

B₀ = E₀₁₂₃ − E₀₁₂₇ + E₀₁₃₆ + E₀₁₆₇ − E₀₂₃₅ − E₀₂₅₇ + E₀₃₅₆ − E₀₅₆₇
   + E₁₂₃₄ + E₁₂₄₇ − E₁₃₄₆ + E₁₄₆₇ + E₂₃₄₅ − E₂₄₅₇ + E₃₄₅₆ + E₄₅₆₇

(subscripts list the monomial's gamma indices). Since B₀ intertwines
both ways, B₀² commutes with all eight gammas, hence is scalar:
**B₀² = 16·I** exactly. Therefore (1)∧(2) has exactly two solutions,
**C = ±B₀/4**. For C = B₀/4 (verified identically on both engines):

- (1) C² = +I. ✓
- (2) Ad_C(E_a) = E_{a⊕4} for all eight a — the strict, unsigned
  sector swap. ✓
- (3) C Ω C⁻¹ = **+Ω**. ✗ — [C, Ω] = 0, {C, Ω} ≠ 0: C is even, as
  V∩W forces. tr C = 0.
- (4) {C, C} = 2I ✓. [C, E_a] is a sum of 8 grade-3 monomials with
  coefficients ±1/2; {C, E_a} a sum of 8 grade-5 monomials with
  coefficients ±1/2 (full support tables are printed by
  `blending_odd_solve2.py`). Derived relations, matching the tables:
  [C, E_{4+μ}] = −[C, E_μ] and {C, E_{4+μ}} = {C, E_μ}. So the brackets
  and anticommutators land on odd (mixed-grade) monomials in the
  ambient algebra, and C paired with itself lands on the
  sector-preserving scalar 2I — the pairing-back *pattern* is present;
  only the chirality flip is absent.

Logged as a partial pattern per the spec, not rounded up: the
sector-blending involution exists and is unique, and it provably
cannot flip chirality.

## 4. Why (3) must fail — the wall (PROVED)

**Theorem.** Let C ∈ M(16, ℚ) with C² = I, and suppose Ad_C maps the
R-sector gamma span onto the J-sector gamma span (hence also conversely,
since Ad_C is an involution). This covers condition (2) in its stated
pairing form and every signed, permuted, or linearly mixed involutive
span swap. Then C Ω C⁻¹ = +Ω.

*Proof.* Ad_C acts on the 8-dimensional gamma span as a linear map L
with L² = I of block anti-diagonal form [[0, N], [M, 0]] (R↔J blocks).
L² = I gives NM = I₄, so det M · det N = 1, and
det L = (−1)^{4·4} det M det N = +1 (the two 4-blocks contribute an even
reordering sign). Conjugation sends the volume element Ω = E₀···E₇ to
det(L)·Ω = +Ω. ∎

**Lemma (one-direction form).** If X ∈ V (or X ∈ W) and X² = I, then
Ad_X is the exact sector swap, so X Ω X⁻¹ = +Ω by the same product
computation (E₄E₅E₆E₇E₀E₁E₂E₃ = +Ω, 16 transpositions). Hence no odd
element of V or W is an involution: the 8-dimensional odd parts of V
and W contain no C with C² = +I at all.

So the step-(d) obstruction — "every candidate was even" — was not an
artifact of the screened directions. It is structural at two levels:
the intertwiner space implementing the swap is purely even in the
relevant direction (V∩W), and any span-swapping involution whatsoever
acts on Ω with determinant +1. A chirality flip requires determinant
−1 on the gamma span (a reflection-type action); a 4+4 sector swap has
determinant +1. One involution cannot be both.

## 5. Step (c) — Lie closure of the partial pattern (PROVED, exact)

For the unique swap involution C of §3, the Lie closure of
{C} ∪ {E₀..E₇} in M(16, ℚ) is **sl(16, ℚ), dimension 255** — the full
traceless algebra (parity split: odd part 128, even part 127).
Certificate: iterative bracket closure saturated, then an explicit
final pass verified all 2295 brackets [basis, generator] lie in the
span. For contrast, the gammas alone close at dimension 36
(8 gammas + 28 bivectors: γ-span ⊕ so(4,4)). So C's brackets do NOT
collapse back to sector-preserving pieces: adjoining the blending
involution generates the entire traceless ambient algebra. The
blending is Lie-algebraically maximal even though the chirality flip
fails.

A mod-p (p = 10007) desktop computation of both closures reproduced
36 and 255 independently (CHECKED corroboration).

## 6. Scope, labels, guardrails

- PROVED (exact rational arithmetic, Sage on toe; second engine:
  desktop exact integer/Fraction verification of B₀, C, the conjugation
  table and the Ω test, plus mod-p rank cross-checks at p = 10007 and
  p = 10009 reproducing dim V = 16, dim V∩W = 1, odd part of V∩W = 0):
  §§2–5 — the empty odd solution set, the unique (1)∧(2) involution,
  its condition table, the determinant wall, the closure dimension.
- The wall is proved for condition (2) as the spec states it:
  conjugation pairing of the sector generators, equivalently an
  involutive swap of the two sector gamma spans. A looser condition
  (e.g. swapping the sector *algebras* while moving gamma spans to
  non-gamma generating sets) is a different condition and was not
  searched; it would first need its own spec.
- Nothing here is INCOMPLETE: steps (a), (b), (c) all ran to
  completion. Combined with step (d), the hunt is closed: conditions
  (1)–(4) are jointly unsatisfiable, and the exact failing combination
  is (2)∧(3) under (1).
- "Antineutrino" remains Kit's heuristic search label only; no physics
  identification is made or implied.

## 7. Artifacts and reproduction

- `proofs/blending_odd_solve.py` — part 1: gates, monomial machinery,
  the linear solves and parity splits, the B₀ basis.
- `proofs/blending_odd_solve2.py` — part 2: C's condition table,
  condition-(4) support tables, Lie closures with certificates.
- `proofs/blending_odd_check.py` — desktop cross-check (exact integer
  arithmetic; B₀ transcribed from part 1's output and re-verified
  against both intertwining systems; mod-p rank cross-checks; mod-p
  Lie closure scout).
- Toe run note: scripts were run on toe as `.sage` copies — plain
  `sage file.py` fails there on a rational_field import-order error and
  the interactive stdin pipe mangles compound blocks; the `.sage` path
  initializes Sage correctly. No GitHub push.
