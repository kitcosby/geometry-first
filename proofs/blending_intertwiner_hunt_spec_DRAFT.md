# DRAFT — Hunt spec: a sector-blending intertwiner C ("antineutrino-like bracket")

Status: SPEC ONLY — nothing here is a result. Assembled 2026-10-14 from the
main-chat discussion of the same day. "Antineutrino-like" is Kit's heuristic
search label only (algebra-first: no physical identification without a
derived Γ map).

## Background already in the record (scope as labeled there)
- A0.3: the imaginary sector carries the same structure as the real sector.
  Consequence drawn in-thread 2026-10-14: if no internal test distinguishes
  "real volume" from "imaginary volume," the real/imaginary designation is
  ONE global orientation choice, and chirality bookkeeping downstream of it
  is globally dependent, not volume-by-volume local. (Consistent with T1:
  global chirality is not derivable from the math; it is fixed by
  stipulation or a Γ-side witness.)
- T2: the real↔imaginary intertwiner space V = Hom(R, J) is 16-dimensional,
  M₄(R) as a vector space. T3: an involution with X² = +I exists there
  (mirror-type blending). T4: no real X with X² = −I (no rotation-type).
- D2: V is NOT closed under ordinary matrix multiplication (B₀B₀ ∉ V,
  I₁₆ ∉ V); V carries an exact inner product (Gram 4I₁₆, orthogonal basis
  of norm 2) — a corridor, not a bracket space.
- Ambient bracket: Cl(4,4) ≅ M(16, R), so [A, B] = AB − BA exists for every
  pair in the ambient algebra; the whole Lie ladder is commutator closure.
- E6 (in record): mixed bivectors S₀₄, S₀₅, S₁₄, S₁₅ map the two (1,1)
  planes into one another under commutator; [S₀₄,S₁₅] = [S₀₅,S₁₄] = 0, the
  other mixed pairs nonzero; T = S₀₄S₁₅ − S₀₅S₁₄ satisfies T² = 4I,
  tr T = 0. (Script: proofs/e6_interplay.py per 2026-10-08 log.)

## The object to search for
A single conjugation-like intertwiner C ∈ V (equivalently in the ambient
M(16, R) action) satisfying ALL of:
  1. C² = +I  (mirror-type; T3 says this kind exists; T4 rules out −I).
  2. Sector swap by conjugation: C E_μ C⁻¹ = E_{4+μ}-type pairing for the
     sector generators — C maps each sector's generators to the OTHER
     sector's, not merely coexisting with them.
  3. Chirality flip: C Ω C⁻¹ = −Ω for the global chirality element Ω, so C
     reverses the global handedness set by the A0.3 designation choice.
  4. Pairing-back: brackets/anticommutators of C with sector elements land
     in the mixed space; C paired with itself lands on sector-preserving
     (diagonal) terms. Note [C, C] = 0 trivially — the content is in the
     anticommutator {C, C} = 2I and in [C, sector] / {C, sector} landing
     spots, computed in M(16, R) since V itself has no internal bracket (D2).

## Search plan (bounded, exact-first per standing rules)
  a. Fix the T2 orthogonal basis of V (Gram 4I₁₆) and parametrize
     C = Σ cᵢBᵢ with symbolic/rational cᵢ.
  b. Impose C² = +I and the conjugation-pairing (2) as exact polynomial
     constraints; solve over ℚ first (Sage on toe, targeted imports).
  c. For any solution, test (3) exactly, then compute the Lie closure of
     {C} ∪ {sector generators} in M(16, R): does it close, at what
     dimension, and does it collapse back to sector-preserving pieces?
  d. Screen first against the known mixed directions: the E6 bivectors and
     T = S₀₄S₁₅ − S₀₅S₁₄ (T² = 4I, so T/2 is an involution candidate
     up to normalization) — check whether any E6 element or combination
     already satisfies (2)–(3) before a blind 16-parameter solve.

## Outcomes — both are clean
- HIT: an explicit C + the closed blending subalgebra it generates, all
  scope-labeled PROVED/CHECKED by the usual two-engine/exact standards.
- MISS: a wall result — no intertwiner satisfies (1)–(4) simultaneously —
  stated with the exact constraint set that fails, published with the wall
  per the external-review standing advice.
- Partial patterns (satisfies some conditions, fails a named one) are
  logged as such, never rounded up.

## Addendum 2026-10-16 — step (d) is RUN; constraint for steps (a)/(b)

Step (d) was run and is a clean MISS, filed at
`proofs/blending_step_d_screen.md` (two engines: desktop exact
rationals + Sage on toe). None of the screened candidates (S_04, S_05,
S_14, S_15, both products, T, T/2, the four pairwise sums) lies in V,
none swaps sectors by conjugation, and none flips chirality. T/2
satisfies (1) exactly as the W-plane volume involution
(T = −2·E_0E_1E_4E_5) — a partial pattern, logged, not rounded up.

Structural consequence (PROVED in the screen, §3): every screened
candidate is EVEN, and Ω commutes with every even element, so
C Ω C⁻¹ = +Ω for any invertible even C — condition (3) is unreachable
by even elements at all. Therefore, before the step (a)/(b)
16-parameter solve is run:

- add {C, Ω} = 0 (C anticommutes with Ω) as an explicit constraint
  from the start — with C invertible it is sufficient for the flip in
  (3), and it is necessary for any C satisfying (3);
- record the parity of each T2 basis element B_i of V first, and do
  not spend the solve on an even-only parametrization;
- screening further even elements cannot succeed and is closed.

Steps (a)–(c) remain unrun. This addendum changes the search plan
only; it is not a result about existence — HIT/MISS for (1)–(4)
simultaneously is still open until (a)–(c) run under the added
constraint.

## Not claimed / guardrails
- No physics identification of C (the "antineutrino" name stays a search
  heuristic until a Γ map is derived).
- No claim that the blending subalgebra exists; existence is the question.
- The saw/octant orientation datum (companion draft) does NOT feed this
  hunt: it provably cannot see the spinor-chirality sign C would flip.
