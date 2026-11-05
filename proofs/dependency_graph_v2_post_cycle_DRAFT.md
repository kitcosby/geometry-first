# Theorem Dependency Graph v2 — Post-Cycle Consolidation (DRAFT, 2026-10-12)

Status of this draft: prepared as a consolidation note after the sl(8) cycle theorem
(archived 2026-10-12 as `proofs/sl8_cycle_theorem.tex`, Kit's human proof).
It does NOT replace `dependency_graph_e70_e110.md` or `theorems/INDEX.md` yet —
it records what changed after both were written, so Kit can approve a merge.

## What changed after the v1 graph (built 2026-10-12 ~15:29 UTC)

1. **Cycle theorem (Kit, human proof, PROVED).** For any n, and any nonzero
   a_i on an interval where the operators are defined and nonzero,
   X_i = a_i E_{i,i+1} (indices mod n) generate sl(n). Proved by induction
   around the cycle: [X_i, X_{i+1}] gives E_{i,i+2} up to a unit, induction
   gives all E_{i,j} (i != j), and [E_{i,j}, E_{j,i}] gives the Cartan.
   Instantiated at n = 8 with the eight operators
   (srx, sxp, cxp, crx, urx, sawup, uxp, sawdown): sl(8), dim 63, rank 7, A_7.
   Killing signature for the split form sl(8,R): (35, 28, 0)
   by the general formula ((n^2+n-2)/2, (n^2-n)/2, 0); B(X,Y) = 16 tr(XY).

2. **Mobius / anharmonic cross-check of sl(6) (PROVED symbolically).**
   The six anharmonic transforms of z = srx are nonzero units in F = Q(z);
   placed in a 6-cycle they generate sl(6, F) by the same cycle argument.
   Specializing z -> sqrt(3) (x = 60 deg) gives sl(6, R). This is an
   independent route to sl(6) that does not use the block ladder.

3. **E100 / E103 / E104 promoted CHECKED -> PROVED** (see
   `e100_e103_e104_proved.md`): E103 by the fixed 4-dim span
   {H, E12, E21, K}; E100 by the dimension argument (an 8-dim real
   subalgebra of the 8-dim sl(3,R) is sl(3,R), so (0,8) is unreachable
   with real matrices); E104 conditional on the rank-4 hypothesis.

4. **Pythagorean Foundation archived** to `proofs/pythagorean_foundation.md`,
   and Kit's direct-algebra proof gives E77 a second independent proof.

## Revised minimal core (v2)

1. **Pythagorean Foundation** — operator definitions and core identities.
2. **Cycle theorem (general n)** — the main route upward: sl(n) for any n
   from any nonzero cycle. Subsumes the *existence* claims of the block
   ladder (E88 -> E90 -> E91/E94) and of E93/E96/E97/E98 as corollaries
   in spirit; the ladder remains valuable as independent exact
   cross-verification (Sage, QQbar), not as the primary proof route.
3. **E107 -> E106** — the *minimality* result (2 independent functions,
   4 positions, overlapping sl(2)s lift to sl(3)). Different virtue from
   the cycle theorem: fewest ingredients, not greatest generality.
4. **E101 (Wick rotation)** — the only route to the compact form su(3),
   (0,8); explicitly imports i, dagger, and the target real form.
5. **Mobius extension** — ONEGEN in action: one generator z = srx, its
   Mobius orbit supplies arbitrarily many nonzero cycle entries.

## Wall clarification (important for any paper revision)

E108 / E109 / E110 are NOT walls for sl(n) in general. They are walls for
one specific construction: the *sampled minimal chain* that reuses the
same 2 functions in overlapping sl(2) positions. That construction stalls
at rank 10/15 for sl(4) because the sampled matrices live in a 4-dim
function space (E103, PROVED) and neither more samples (E109) nor Re/Im
complex variants (E110) add real independence. The cycle construction
bypasses the wall entirely: its generators are *placed* at distinct cycle
positions with nonzero units, not sampled from the T_a(x) curve.
Any write-up that lists E108-E110 as "the ladder terminates" would now
be wrong; the correct statement is "the sampled minimal chain terminates
at n = 3; the cycle construction does not terminate."

## Index gap

`theorems/INDEX.md` (last updated 2026-10-09) still totals "27 PROVED"
and contains no Lie-campaign entries at all (no E70-E110, no cycle
theorem, no Mobius sl(6)). The external review advice Kit asked to have
remembered on 2026-10-12 (points 3 and 5: build the dependency graph,
reduce to the minimal generating structure) cannot be acted on from that
index as it stands. Suggested merge, pending Kit's approval:
add a "Lie campaign" section to the index using the v2 minimal core
above, keep the 27 pre-Lie theorems counted separately, and update
`sl_n_exact.tex` (2026-10-12 05:46 version, block ladder to sl(6) only)
to make the cycle theorem the main theorem with the ladder as
cross-verification.

## Still open (unchanged, do not re-push)

- Z/W complex-structure moduli: Kit approved the research order
  "parameter theorem first, then Z/W" on 2026-10-12. The cycle theorem
  and E104 cover much of the parameter-theorem ground for the cycle and
  diagonal cases, but an explicit statement for general traceless 3x3
  operator families, and the Z/W moduli study itself, remain open.
- Doubling map z(2x): branch-dependent ((srx-sxp)/2 on FW = +1,
  (crx-cxp)/2 on FW = -1); the absolute-value form is disproved.
  Piecewise Mobius covariance is the correct framing.
