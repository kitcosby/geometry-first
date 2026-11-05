# Consolidation — 2026-10-12 evening results (DRAFT, for Kit's approval)

Prepared 2026-10-12 ~15:35 PDT. Extends `dependency_graph_v2_post_cycle_DRAFT.md`
(16:38 UTC) with the results recorded later the same day (tracking updates
effective 17:07–20:22 UTC). It does NOT replace `theorems/INDEX.md` (still
"27 PROVED", last updated 2026-10-09, no Lie-campaign entries) — the v2 draft's
suggested merge still needs Kit's approval; this draft only makes the evening
batch auditable in one place, with each result's exact scope as written in its
own proof file.

## Evening batch — statement, scope, proof file

1. **H3 / H4 / H6 — PROVED** (exact symbolic, all 4 quadrants).
   H4: quarter-turn swaps families, srx(x+π/2)=crx(x), cxp(x+π/2)=sxp(x);
   hence urx and uxp ARE the negated quarter-turn deltas of srx and cxp.
   H3: doubling delta, branch-dependent ((srx−sxp)/2 on FW=+1,
   (crx−cxp)/2 on FW=−1). H6: reciprocal delta.
   File: `harmonics_h3_h4_h6_deltas.md`.
2. **Cross-family diagonal — PROVED** (elementary algebra from the
   Pythagorean Foundation definitions): (srx+crx)−(sxp+cxp) = 4 cot(2x),
   the E77 numerator. File: `cross_family_sum_difference_DRAFT.md`
   (filename says DRAFT; the file's own status line says PROVED — Kit may
   want the filename cleaned up at merge time).
3. **Jacobi pair transparency — PROVED** (symbolic): every complete
   reciprocal pair contributes exactly 1 to bracket/Jacobi coefficients;
   yields Chevalley normalization for sl(8). File: `jacobi_pairs.md`.
4. **Cross-family couplings — PROVED** (symbolic in s,c,u,v and in t):
   connector brackets carry srx·cxp as a Cartan-scale coupling (no
   cancellation), while the connector cross-products telescope to 1
   around the full cycle. File: `jacobi_cross_family.md`.
5. **Cycle product C^8 = I — PROVED** items 1–3, 6 / **CHECKED** items 4–5
   of `expmap_operators.md`: ∏aᵢ=1 (four reciprocal pairs), so C is an
   exact order-8 Coxeter element, charpoly λ⁸−1, spectrum independent
   of x; under exp the operators act as rotation angles.
6. **M1–M2 scalar invariant — PROVED** (trace identity, all quadrants) /
   **CHECKED** (entry forms): tr(M1·M2) = 3 identically, equivalently
   urx·uxp − srx·crx − sxp·cxp = 2. Recorded near-miss: det M2 at 30° is
   1.6188, NOT the golden ratio. File: `cross_m1_m2.md`.
7. **Quadrant flip = perfect-square swap — PROVED** (symbolic in w and t):
   the π-shift is w → iw, i.e. t ↦ −1/t; α = sgn(Im ζ) selects which
   Euler square (w+w⁻¹)² / (w−w⁻¹)² sits upstairs; at block level
   M1_QIII = diag(1,−1) M2_QI diag(1,−1). This is the mechanism Kit
   named "self propagation" (2026-10-12 19:18 UTC).
   File: `quadrant_flip_polar.md`.
8. **SO(3) exponential — PROVED** (structure, θ² identity) / **CHECKED**
   (numeric matrix at 30°): Rodrigues closes exp(Y) for
   ω=(sxp,cxp,uxp), and θ² = 2(sxp² − sxp·cxp + cxp²), the A2 Eisenstein
   norm of the cores, because uxp = cxp − sxp.
   File: `expmap_uxp_so3.md`.
9. **Wick-on-generators wall — PROVED** (conjugacy identity) / **CHECKED**
   (dim-28 closure, compactness): Wick-rotating the 8 cycle generators
   gives a unitary conjugate of so(8), dim 28, NOT su(8). Wick does not
   commute with generation; the correct route Wicks the full sl(8)
   basis (28 keep, 35 × i). File: `wick_sl8_wall.md`.
10. **M_refined generates sl(4,R) — PROVED** (Sage 9.5 exact: QQ(t) for
    det, QQbar for closure/Killing): raw M_refined closes at dim 16
    (gl(4,R)); its traceless correction closes at dim 15, Killing
    (9,6,0) ⇒ sl(4,R); det factors exactly in t=tan(x/2), vanishing at
    t²=1/5 in QI. Scripts archived 2026-10-12 (see pending M_refined
    rescue note). File: `m_refined_sl4.md`.

## Explicitly NOT theorems in this batch

- **Maxwell exercise** (`maxwell_exercise.md`): operator identities
  PROVED; the E/B identification is a declared dictionary (ASSERTED);
  Maxwell's equations as spacetime PDEs are NOT derived. A
  representation witness only, per Kit's certified spine.
- **O_prop**: still OPEN. No new candidate emerged today; the proved
  transport laws model its shape as selection criteria only, and no
  12×12 candidate passed a leakage test. See also
  `oprop_next_step_correction_DRAFT.md`: the R(θ)=0 leakage test as
  previously offered is NOT executable (R has no definition in the
  record — run 233 returned BLOCKED); the pending item is a definition
  from Kit, plus the reframed choices in that draft.

## Script archive status (this run)

The proof files for items 6, 8, 9, 7 and the Maxwell exercise cited
scripts living only in /tmp (ephemeral). Rescued 2026-10-12 by this
research run, byte-identical, sha256-verified against the /tmp originals:
`cross_m1_m2_a/b/c.py` (was /tmp/cross2*.py), `expmap_uxp_so3_a/b.py`
(/tmp/uxpexp*.py), `wick_sl8_a/b.py` (/tmp/wick8*.py),
`quadrant_flip_polar_a/b.py` (/tmp/polar*.py), `maxwell_exercise.py`
(/tmp/maxwell_ex.py) — all now in `geometry-first/proofs/`.
Items 1, 3, 4, 5 in the list above cite no separate script files in
their proof notes (symbolic work recorded inline); item 2 is hand
algebra from the Foundation. /tmp/mref.py (19:41 UTC, an earlier
M_refined attempt) was NOT rescued — the four archived m_refined_*.py
scripts are the ones the proof file cites.

## Remaining decision for Kit (one, unchanged from the v2 draft)

Approve the INDEX merge: add a Lie-campaign section built on the v2
minimal core plus this evening batch (items 1–10 above), keep the 27
pre-Lie theorems counted separately, and update `sl_n_exact.tex` so the
cycle theorem is the main theorem with the block ladder as
cross-verification. Nothing here pushes to GitHub.
