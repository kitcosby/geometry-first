# Theorem Dependency Graph: Lie Results E70–E110

## Root: Pythagorean Foundation
**Source:** `proofs/pythagorean_foundation.md` (Kit Cosby, 2026-10-12)
**Status:** PROVED
**Provides:** Operator definitions (srx, sxp, cxp, crx, urx, uxp), reciprocal conjugates, UNA forms, double-angle, FlatWave=sgn(sin2x), octant cycle.

---

## Level 1: Core Identities (independent of Lie)

### E77 — Double-Angle Identities
**Status:** PROVED (t=tan(x/2), reviewer sketch + analytic verification; also Kit's direct algebraic proof 2026-10-12)
**Depends on:** Pythagorean Foundation (operator definitions)
**Provides:** sin(2x)/4 = 1/(urx+uxp), cos(2x) = (srx+crx−sxp−cxp)/(urx+uxp)
**Used by:** (none directly; foundational identity)

---

## Level 2: Base Lie Algebras

### E70–E72 — sl(3,R) via T_a
**Status:** PROVED (exact rank 8, 28 commutators close, Killing (5,3), C2=4/9, Cartan 2-dim, A2 roots, all over QQbar)
**Depends on:** Pythagorean Foundation
**Provides:** First exact sl(3,R) from operators; reference for all sl(3) variants
**Used by:** E96, E98, E101, E103, E104

### E107 — sl(2,R) via 2×2
**Status:** PROVED (rank 2→3, closure YES, Killing (2,1))
**Depends on:** Pythagorean Foundation (srx, cxp only)
**Provides:** Minimal building block; 2 functions suffice for sl(2)
**Used by:** E106 (lift via overlapping)

### E88 — sl(4,R) via Blocks
**Status:** PROVED (exact rank 15, 105 commutators, Killing (9,6,0))
**Depends on:** Pythagorean Foundation (block construction pattern)
**Provides:** Second rung; block method template
**Used by:** E93, E99, E90b/c

---

## Level 3: Higher Rungs (block method)

### E90b–E90c — sl(5,R)
**Status:** PROVED (rank 24, 276 commutators, Killing (14,10,0))
**Depends on:** E88 (block pattern)
**Used by:** E91

### E91 / E94 — sl(6,R)
**Status:** PROVED (E91: rank 35, 595 commutators; E94: independent rank-35 re-verification; 35-dim subalgebra of 35-dim sl(6,R) ⇒ spans)
**Depends on:** E90b/c
**Note:** E92c numeric Killing (18,17,0) discarded as unstable; identification via dimension count, not Killing.

---

## Level 3: sl(3) Variants (all PROVED, independent constructions)

### E96 — FlatWave Matrix sl(3)
**Status:** PROVED (traceless-corrected M−tr/3·I; rank 8, 28 commutators)
**Depends on:** E70 (reference; independent construction)
**Provides:** Second sl(3); first to use FlatWave explicitly
**Used by:** E97

### E97 — Corner Variants
**Status:** PROVED (corners (0,0), (sin x,cos x), (sin2x,cos2x) all give sl(3,R))
**Depends on:** E96 (FW skeleton)
**Provides:** Corners are free (for tested cases)

### E98 — Sawtooth Matrix sl(3)
**Status:** PROVED (diagonal sawup,sawdown,FW; rank 8, closure)
**Depends on:** E70 (reference; independent)
**Provides:** Third sl(3); uses quasi-linear sawtooths natively

### E105 — Zero-Diagonal sl(3)
**Status:** PROVED (6 off-diagonals, zero diagonal; rank 4→8, closure YES, Killing (5,3,0))
**Depends on:** E104 (diagonal freedom)
**Provides:** Fourth sl(3); sparsest (no diagonal)

### E106 — Minimal 4-Operator sl(3)
**Status:** PROVED ([[0,srx,0],[cxp,0,crx],[0,sxp,0]]; rank 4→8, closure YES, Killing (5,3,0))
**Depends on:** E105 (remove UNA) + E107 (sl2 lift)
**Provides:** Fifth sl(3); minimal (2 independent functions, 4 positions); parameter theorem base case

---

## Level 3: Parameter Exploration

### E93 — Nonzero Extension (sl4)
**Status:** PROVED (specific nonzero v-entry still gives sl(4,R); rank 15, Killing (9,6))
**Depends on:** E88
**Provides:** Robustness evidence (not universal proof)

### E99 — w=−v (sl4)
**Status:** PROVED (sign flip preserves (9,6) sl(4,R))
**Depends on:** E88
**Provides:** Real form stability

### E103 — 5° Steps
**Status:** CHECKED (numeric; 17 angles → rank 4)
**Depends on:** E70
**Provides:** Finer sampling doesn't increase rank; curve is 4-dim

### E104 — Diagonal Parameters
**Status:** CHECKED (5 patterns tested: urx, constant, varying, zero, all give rank 8)
**Depends on:** E70
**Provides:** Evidence that any traceless diagonal works (including zero)
**Used by:** E105

---

## Level 4: Real Forms and Complexification

### E100 — Random Permutations
**Status:** CHECKED (100 samples; all closed gave (5,3,0); zero gave (0,8))
**Depends on:** E70 (concept)
**Provides:** Evidence that (5,3) is structural to real operators, not arrangement
**Note:** Not a proof; structural argument is dimension-based (8-dim subalgebra of 8-dim space)

### E101 — Wick Rotation → su(3)
**Status:** PROVED (3 antisym + 5 i·sym; rank 8, closure YES, Killing (0,8), eigenvalues [−6,−18,−12×6])
**Depends on:** E70 (Chevalley basis)
**Provides:** (0,8) reachable via Wick; compact form from split via i
**Note:** Imports i, †, and target (not native)

---

## Failed / Stalled (Negative Results)

### E108 — Minimal sl(4) Chain (4 samples)
**Status:** FAILED (rank 10/15, stalled)
**Depends on:** E106, E107 (extrapolation attempt)
**Provides:** Wall: minimal chain insufficient for n≥4 with 4 samples

### E109 — Minimal sl(4) Chain (8 samples)
**Status:** FAILED (rank 10/15, stalled; 8 samples → still rank 4 initially)
**Depends on:** E108
**Provides:** Wall: more samples don't help; function space is 4-dim

### E110 — Complexified sl(4)
**Status:** FAILED (rank 10/15; 8 Re/Im variants → still rank 4)
**Depends on:** E108
**Provides:** Wall: complexification doesn't add real independence

---

## Summary Counts
- **PROVED:** E70-72, E77, E88, E90b/c, E91/E94, E93, E96, E97, E98, E99, E101, E105, E106, E107 (14)
- **CHECKED:** E100, E103, E104 (3)
- **FAILED (wall):** E108, E109, E110 (3)

## Minimal Core
The dependency graph collapses to:
1. **Pythagorean Foundation** (definitions + core identities)
2. **E107** (sl2 from 2 functions) → **E106** (sl3 from overlapping sl2s)
3. **E88** (sl4 via blocks) → **E90** (sl5) → **E91** (sl6)
4. **E101** (Wick → su3)

Everything else is variant, robustness check, or wall.
