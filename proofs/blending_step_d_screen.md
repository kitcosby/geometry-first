# Blending intertwiner hunt — step (d) screen: RESULT = MISS (clean)

**Date:** 2026-10-16 · Run at Kit's order (main chat). Spec:
`proofs/blending_intertwiner_hunt_spec_DRAFT.md`, step (d): screen the
E6 mixed bivectors S_04, S_05, S_14, S_15 and T = S_04·S_15 − S_05·S_14
against conditions (1) C² = +I, (2) sector-swap conjugation,
(3) C Ω C⁻¹ = −Ω, (4) pairing-back; V membership also tested.
**Outcome: no candidate satisfies the conditions. The miss is
structural and is proved, not merely observed (see §3).**

## 1. Setup and reproduction gate (PROVED)

Realization rebuilt exactly as in `proofs/d2_structure.py`:
e = (σx⊗I, Y⊗I, σz⊗σx, σz⊗Y) in Cl(2,2); E_μ = e_μ⊗I₄ (μ = 0..3, the
R sector), E_{4+ν} = Ω₄⊗e_ν (the J sector); signature
(+,−,+,−,+,−,+,−); S_ab = E_aE_b. Gate against the E6 record, exact:
[S_04, S_15] = [S_05, S_14] = 0, the other four mixed pairs nonzero,
T² = 4I, tr T = 0 — all reproduced. **New identification: T =
−2·E_0E_1E_4E_5** (i.e. T/2 is minus the W-plane volume element for
W = {0,1,4,5}); also S_04·S_15 = −E_0E_1E_4E_5, S_05·S_14 =
+E_0E_1E_4E_5. Two engines: desktop exact integer/rational
arithmetic; Sage 9.5 on toe with targeted imports — all key items
cross-checked in both (gate, T identification, (T/2)² = I, the T/2
conjugation table, Ω tests, V tests).

## 2. Screen results (exact)

Conjugation is diagonal for every candidate: each E_μ maps to ±E_μ
(itself), never to the sector partner E_{μ±4}.

| candidate | C² | conj. signs on E_0..E_7 | Ω | in V | (2) swap |
|---|---|---|---|---|---|
| S_04 | −I | (−,+,+,+,−,+,+,+) | +Ω | no | no |
| S_05 | +I | (−,+,+,+,+,−,+,+) | +Ω | no | no |
| S_14 | +I | (+,−,+,+,−,+,+,+) | +Ω | no | no |
| S_15 | −I | (+,−,+,+,+,−,+,+) | +Ω | no | no |
| S_04·S_15 = −vol_W | +I | (−,−,+,+,−,−,+,+) | +Ω | no | no |
| S_05·S_14 = +vol_W | +I | (−,−,+,+,−,−,+,+) | +Ω | no | no |
| T = −2·vol_W | 4I | (−,−,+,+,−,−,+,+) | +Ω | no | no |
| T/2 | +I | (−,−,+,+,−,−,+,+) | +Ω | no | no |
| S_04 ± S_15, S_05 ± S_14 | nonscalar | — (singular: no inverse, no conjugation) | — | no | no |

Pairing-back (4), from the tables: {C, C} = 2C² (so 2I exactly for
the +I rows; −2I for S_04, S_15). With s_μ the conjugation sign,
[C, E_μ] = (1 − s_μ)·C·E_μ and {C, E_μ} = (1 + s_μ)·C·E_μ: brackets
with negated gammas land on single monomials — e.g. [S_04, E_0] =
−2E_4 (the E6 plane-mapping phenomenon) — and anticommutators with
fixed gammas land on C·E_μ monomials. No bracket produces a sector
exchange by conjugation; the swap condition (2) is about conjugation
and fails in every row.

## 3. Why the miss is structural (PROVED)

Every screened candidate is an EVEN element of Cl(4,4). The volume
element Ω commutes with every even element, so for any invertible
even C, C Ω C⁻¹ = +Ω identically: **condition (3) is unreachable by
any even element**, not just by these. Condition (2) also fails
maximally here: monomial conjugations are sign patterns (plane
π-rotations and plane-volume actions), which fix every gamma's line.
And no candidate lies in V = Hom(R, J) to begin with (intertwiner
test X·R_μ = J_μ·X fails for all).

**Consequence for the hunt (steps a/b):** a C satisfying (1)–(3)
must anticommute with Ω (C Ω + Ω C = 0 with C invertible suffices for
the flip) — in particular it cannot be even. The 16-parameter solve
of step (b) should be run with that constraint built in; screening
further even elements cannot succeed. Partial patterns logged per
spec: T/2 satisfies (1) exactly and is the W-plane volume involution;
S_05 and S_14 satisfy (1) as mixed-plane π-rotations. None is
rounded up.

## 4. Scope

Exact rational/integer arithmetic throughout; no floating point.
"Antineutrino" remains Kit's heuristic search label only; no physics
identification is made or implied. Steps (a)–(c) of the spec are not
run here.
