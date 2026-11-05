# T4 admission proposal — geometry-first theorems/INDEX.md (DRAFT, nothing applied)

Prepared 2026-10-20 (PT) after the family-bridge directed thread closed.
Status: awaiting Kit's word (rides with the next geometry-first sync,
which is separately Kit-gated). No index line has been edited.

## Why the withheld status is now stale

`theorems/INDEX.md` still carries, in "Withheld — strong evidence only":

> | T4* | X²=−I: 2 solutions over ℚ̄, both non-real; 0 real found —
> strong evidence, NOT a complete real-algebraic proof of emptiness |
> `proofs/q6b_groebner.py`, `proofs/q6b_solutions.py` |

and the line "**Total: 27 PROVED.** T4* withheld — strong evidence
only, not admitted under the proved-only rule."

The gap in the original q6b evidence was its reality screen: numerical
(1e-10 thresholds, first 10 solutions). That gap is now closed twice
over, both exact:

1. **Run 148** (2026-10-08, `goals/r-theory-volume-i-audit/hidden_files/
   io_word_sweep/RUN148_T4_NONEXISTENCE.md`): complete proof of
   nonexistence — {R_μ, J_ν} = 0 exactly (16/16); S_μ, D_μ satisfy the
   Cl(4,4) relations over ℚ(√2) (8/8 squares, 28/28 anticommutators);
   any X ∈ A satisfies [X, S_μ] = 0; X² = −I would embed
   M₁₆(ℂ) ↪ M₈(ℂ), contradiction. Labeled PROVED in its record.
2. **Family-bridge thread** (2026-10-20, `io_word_sweep/
   FAMILY_BRIDGE.md` §3, case C2): independent exact enumeration —
   parametrizing X over the 16-dim basis of A, the system X² = −I has
   184 unique equations, a 0-dimensional ideal, and its variety over
   ℚ̄ consists of exactly 2 solutions, both non-real under the exact
   conjugation test. A real solution would be a ℚ̄-solution, so the
   real solution set is empty. The thread also re-verified every
   computational premise of run 148 exactly in Sage.

Nothing in either proof rests at CHECKED level.

## Proposed index edits (exact)

1. Move T4 out of the "Withheld" section into the Clifford-structure
   table as a PROVED entry, statement: "No real X in the
   sector-intertwiner space A with X² = −I (nonexistence proved:
   run 148 commutant argument; independent exact enumeration —
   variety over ℚ̄ has exactly 2 points, both non-real)".
   Proof column: `io_word_sweep/RUN148_T4_NONEXISTENCE.md`;
   `io_word_sweep/FAMILY_BRIDGE.md` §3 (scripts fb_t4variety.py,
   fb_native.py — working copies were in /tmp on toe; re-derive from
   the FAMILY_BRIDGE record if the sync wants them filed in-repo).
2. Total line: "**Total: 27 PROVED.** T4* withheld …" →
   "**Total: 28 PROVED.**" (T4 admitted on the two proofs above.)
3. Delete the now-empty "Withheld — strong evidence only" section, or
   leave the heading with "None currently withheld." — Kit's
   preference at sync time.
4. Downstream note (no edit proposed here): the physics notation
   mapping ledger's X-row records "X² = −I withheld (T4*)" — it should
   be updated to match in the same sync pass so the two ledgers do not
   disagree.

## Scope kept exact

- This admits only the nonexistence statement for the intertwiner
  space A. It does not admit any bridge to the Spin(10) Weyl 16: the
  family-bridge thread proved the bridge obstructed in shapes
  (C1)–(C3), and (C4) remains fenced — it needs a Kit-declared so(10)
  action on a new native module (real 32 with complexification
  16 ⊕ 16′, or a genuinely complex 16), or Kit closes C4 as
  type-obstructed for the native ℝ¹⁶. That choice is separate from
  this index admission.
