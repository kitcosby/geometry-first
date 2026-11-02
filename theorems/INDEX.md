# Theorems Index

Each entry: statement, status, proof script in `../proofs/`. The pre-Lie
sections below are all PROVED (27 total). The Lie campaign section that
follows is counted separately and carries per-entry scope labels, including
CHECKED parts and OPEN items, exactly as established.

## Clifford structure

| ID | Statement | Proof |
|----|-----------|-------|
| T1 | Cl(4,4) 16-spinor is irreducible (commutant dim 1, exact over ℚ) | `proofs/q4_exact.py` |
| T2 | Sector intertwiner space dim = 16 (exact over ℚ) | `proofs/q5_exact.py` |
| T3 | ∃ X in intertwiner space with X²=+I (explicit, over ℚ̄, verified) | `proofs/q6a_groebner.py` |

## Operator algebra

| ID | Statement | Proof |
|----|-----------|-------|
| RECIP | sawup + sawdown = FlatWave (exact) | algebra |
| ONEGEN | Every operator is a function of z=srx alone | algebra |
| FLIP1 | crx = \|sec x\| − tan x (exact) | sec²−tan²=1 |
| FLIP2 | sxp = \|csc x\| − cot x (exact) | csc²−cot²=1 |
| INVPAIR | c₁·c₄=1, c₂·c₃=1 (exact) | definitions |

## Symmetry

| ID | Statement | Proof |
|----|-----------|-------|
| V4ACT | V₄ acts on (srx,sxp,cxp,crx) via ι=(0 1)(2 3), τ=(0 3)(1 2), ι_d=(0 2)(1 3) | trig + reciprocal |
| V4GRP | {id,ι,τ,ι_d} is Klein V₄ (closure verified) | permutation computation |
| CHIF | FlatWave character χ_F=(+1,−1,−1,+1) | FW definition |

## Möbius and double-angle

| ID | Statement | Proof |
|----|-----------|-------|
| MOBIUS | cxp = (srx+ε)/(ε·srx−1), ε=FlatWave | equivalent to RECIP |
| KRED | u=(urx,uxp)=K·P, K=[[1,0,0,−1],[0,−1,1,0]] | definitions |
| AB | A=srx−sxp=2cot x, B=cxp−crx=2tan x | FLIP1, FLIP2 |
| COSRAT | K₂/D₂=cos(2x), K₂=A−B, D₂=A+B | AB + trig |
| SINPROD | urx·uxp=4·\|csc(2x)\| globally (+4csc(2x) on Q1/Q3, −4csc(2x) on Q2/Q4); the unsigned form holds only where sin(2x)>0 | `proofs/run150b_doubleangle_product.py` (sympy trigsimp→0 on both \|\.\| branch forms; numeric rel ≤1.2e-11) |
| HYPER | D₂²−K₂²=16 | COSRAT |
| CIRCLE | (sin2x,cos2x)=(4/D₂,K₂/D₂) | HYPER |

## Conic geometry

| ID | Statement | Proof |
|----|-----------|-------|
| QCONIC | Q_conic=diag(16,4,−1); xᵀQx=0 ⟺ 4H_R²+V_R²=1/4 | algebra |
| QHYP | Q_hyp=diag(1,−1,−16); yᵀQy=0 ⟺ Σ²−Δ²=16 | algebra |
| HYP2ELL | (Σ,Δ)↦(1/Σ,Δ/2Σ)=(H_R,V_R) maps hyperbola→ellipse; seam Σ→∞ ↦ (0,1/2) | algebra |
| VERON | ν₂(cos x,sin x) satisfies X₀X₂−X₁²=0 | algebra |
| SPIN1 | ⟨J_z⟩=cos2x, ⟨J_x⟩=sin2x; ⟨J_z⟩/2=V_R, ⟨J_x⟩/4=H_R | spin-1 representation |

## Dynamics

| ID | Statement | Proof |
|----|-----------|-------|
| Q7 | dS/dx=(srx·\|csc x\|−crx·\|sec x\|)·S², S=sawup (Bernoulli, exact per-octant) | `proofs/q7_dSdx.py` |
| DURX | d(urx)/dx=crx·\|sec x\|−srx·\|csc x\| (exact) | Q7 + S=1/urx |
| RUN144 | 8 (M,P,G) types, uniquely identified; G flips for c₁,c₄ at octant boundaries | exact lemma + census |
| HARM | F=(sin2x/4,cos2x/2), dF/dx=[[0,1],[−4,0]]·F | trig |

## Intertwiner geometry (D2/D3)

| ID | Statement | Proof |
|----|-----------|-------|
| D2CLOSE | V not closed under 16×16 multiplication (B₀·B₀∉V; I₁₆∉V) | `proofs/d2_structure.py` |
| D2GRAM | Gram G[i,j]=tr(B_iᵀB_j)=4·I₁₆ (orthogonal, uniform norm 2) | `proofs/d2_structure.py` |
| D3FORM | Exactly 2-dim Spin(4,4)-invariant symmetric forms, both rank 8; 0 skew | `proofs/d3b_forms.py` |
| D3VOL | vol=E₀···E₇: vol²=+I, not scalar, anticommutes with E_μ | `proofs/d3b_forms.py` |

**Total: 27 PROVED.** T4* withheld — strong evidence only, not admitted under the proved-only rule.

## Lie campaign (2026-10-12) — counted separately from the 27 above

Full dependency graph: `proofs/dependency_graph_e70_e110.md`, updated by
`proofs/dependency_graph_v2_post_cycle_DRAFT.md` and
`proofs/consolidation_2026-10-12_evening_DRAFT.md` (both merged into this
index 2026-10-13 with Kit's approval).

### Foundation and main generation theorems

| ID | Statement | Status | Proof |
|----|-----------|--------|-------|
| FOUNDATION | Pythagorean Foundation: operator definitions, reciprocal conjugates, UNA forms, urx·uxp=4·\|csc(2x)\|, FlatWave=sgn(sin 2x), \|saw_r\|+\|saw_x\|=1 partition of unity, half-angle forms, octant dominance cycle | PROVED | `proofs/pythagorean_foundation.md` |
| E77 | Double-angle identities: sin(2x)/4 = 1/(urx+uxp), cos(2x) = (srx+crx−sxp−cxp)/(urx+uxp) | PROVED | analytic (t=tan(x/2)); also Kit's direct algebraic proof in the Foundation |
| CYCLE | For any n and any nonzero a_i, X_i = a_i·E_{i,i+1} (indices mod n) generate sl(n); at n=8 with the eight operators (srx, sxp, cxp, crx, urx, sawup, uxp, sawdown): sl(8), dim 63, rank 7, A7, Killing (35,28,0) | PROVED (human proof, Kit) | `proofs/sl8_cycle_theorem.tex` |
| SO9PATH | The same eight operators as Y_i = f_i·(E_{i,i+1} − E_{i+1,i}) in 9×9 generate so(9), dim 36, rank 4, B4 — the algebra is fixed by quiver/symmetry type, not by the values (nonzeroness only) | PROVED (human proof, Kit) | `proofs/so9_path_theorem.tex` |
| MOBSL6 | The six anharmonic transforms of z = srx are units in ℚ(z) and generate sl(6, ℚ(z)) by the cycle argument; specializing z = √3 (x = 60°) gives sl(6,ℝ) — an independent route to sl(6) that does not use the block ladder | PROVED (symbolic) | recorded 2026-10-12; see dependency graph v2 |

### Block ladder — exact cross-verification (Sage, exact over ℚ̄)

The ladder independently cross-verifies the A-series existence claims for
n = 3–6 by a different construction; it is not the primary proof route
(the cycle theorem is).

| ID | Statement | Status | Proof |
|----|-----------|--------|-------|
| E70–E72 | sl(3,ℝ) via T_a: rank 8, all 28 commutators close, Killing (5,3,0), C2 = 4/9, Cartan 2-dim, A2 roots | PROVED | `proofs/e70_exact.py`, `proofs/e71_exact_casimir.py`, `proofs/e72_exact_rest.py` |
| E88 | sl(4,ℝ) via blocks: rank 15, 105 commutators close, Killing (9,6,0) | PROVED | E88 (see `sl_n_exact.tex`) |
| E90b–E90c | sl(5,ℝ): rank 24, 276 commutators close, Killing (14,10,0) | PROVED | E90b, E90c (see `sl_n_exact.tex`) |
| E91 / E94 | sl(6,ℝ): rank 35, 595 commutators close (E91), independent rank-35 re-verification (E94); identification by dimension count (35-dim subalgebra of the 35-dim sl(6,ℝ)). The E92c numeric Killing (18,17) is discarded as unstable; the exact signature is (20,15,0) | PROVED | E91, E94 (see `sl_n_exact.tex`) |
| E107 | sl(2,ℝ) building block: [[0,srx],[cxp,0]], rank 2→3, closure, Killing (2,1) | PROVED | E107 script archived on GitHub (`proofs/e107_sl2.py`, sha b682e95a8bd0) |
| E106 | Minimal sl(3,ℝ): [[0,srx,0],[cxp,0,crx],[0,sxp,0]], rank 4→8, closure, Killing (5,3,0); 2 independent functions in 4 positions, via two overlapping sl(2)s (E107) | PROVED | E106 script archived on GitHub (`proofs/e106_minimal.py`, sha 191e285ecdc4) |

### sl(3) variants and parameter results

| ID | Statement | Status | Proof |
|----|-----------|--------|-------|
| E96 | FlatWave matrix (traceless-corrected) generates sl(3,ℝ): rank 8, 28 commutators | PROVED | dependency graph E96 |
| E97 | Corner variants (0,0), (sin x, cos x), (sin 2x, cos 2x) all give sl(3,ℝ) — corners free for the tested cases | PROVED | dependency graph E97 |
| E98 | Sawtooth matrix (diagonal sawup, sawdown, FlatWave) generates sl(3,ℝ): rank 8, closure | PROVED | dependency graph E98 |
| E93 | A specific nonzero v-entry still gives sl(4,ℝ): rank 15, Killing (9,6) — robustness evidence for that case, not a universal proof | PROVED | dependency graph E93 |
| E99 | The w = −v sign flip preserves sl(4,ℝ), Killing (9,6) | PROVED | dependency graph E99 |
| E105 | Zero-diagonal sl(3,ℝ): six off-diagonal operators, zero diagonal; rank 4→8, closure, Killing (5,3,0). The specific construction is one of E104's five tested diagonal patterns | PROVED | dependency graph E105 |
| E100 | Structural no-go: an 8-dimensional real Lie subalgebra of sl(3,ℝ) equals sl(3,ℝ), so Killing (0,8,0) is unreachable by real traceless matrices. The 100-permutation search illustrates this; it is not itself the proof | PROVED (dimension argument) | `proofs/e100_e103_e104_proved.md` |
| E103 | T_a(x) lies in the fixed 4-dimensional span {H, E12, E21, E23−E32} for all x; no sampling of the T_a family exceeds rank 4 | PROVED (algebraic) | `proofs/e100_e103_e104_proved.md` |
| E104 | Five specific traceless diagonal patterns — (urx,−urx,0), (1,−1,0), (srx,cxp,−srx−cxp), (0,0,0), (urx,uxp,−urx−uxp) — each generate rank 8 | PROVED for the five tested patterns only; the universal "any traceless diagonal" claim is **FALSE** — falsified 2026-10-26 by verified counterexample (every sample shares eigenvector (1,1,1); rank ≤ 6 < 8, sharp); the constant-pattern subcase is PROVED for every traceless constant D on affinely spanning sample sets — see `proofs/e104_universal_falsified.md` | `proofs/e100_e103_e104_proved.md` |
| E101 | Wick rotation to su(3): 3 antisymmetric + 5 i·symmetric generators, rank 8, closure, Killing (0,8). Imports i, †, and the target real form — not native to the operators | PROVED | dependency graph E101 |

### Evening batch (2026-10-12)

| ID | Statement | Status | Proof |
|----|-----------|--------|-------|
| H3/H4/H6 | Deltas: quarter-turn swaps families — srx(x+π/2)=crx(x), cxp(x+π/2)=sxp(x) — so urx and uxp are the negated quarter-turn deltas of srx and cxp; doubling delta exact per FlatWave branch; reciprocal delta | PROVED (exact symbolic, all 4 quadrants) | `proofs/harmonics_h3_h4_h6_deltas.md` |
| XFAM | Cross-family diagonal: (srx+crx) − (sxp+cxp) = 4·cot(2x), the E77 numerator | PROVED (elementary algebra from the Foundation definitions) | `proofs/cross_family_sum_difference_DRAFT.md` (status line in file: PROVED) |
| JACPAIR | Jacobi pair transparency: every complete reciprocal pair contributes exactly 1 to bracket/Jacobi coefficients; yields the Chevalley normalization for sl(8) | PROVED (symbolic) | `proofs/jacobi_pairs.md` |
| JACXFAM | Cross-family couplings: connector brackets carry srx·cxp as a Cartan-scale coupling (no cancellation), while the connector cross-products telescope to 1 around the full cycle | PROVED (symbolic in s,c,u,v and in t) | `proofs/jacobi_cross_family.md` |
| CYCPROD | Cycle product C^8 = I: ∏aᵢ = 1 (four reciprocal pairs), so C is an exact order-8 Coxeter element, characteristic polynomial λ⁸−1, spectrum independent of x | PROVED (items 1–3, 6 of the source note) / CHECKED (items 4–5) | `proofs/expmap_operators.md` |
| M1M2 | M1–M2 scalar invariant: tr(M1·M2) = 3 identically, equivalently urx·uxp − srx·crx − sxp·cxp = 2. (Recorded near-miss: det M2 at 30° is 1.6188, NOT the golden ratio) | PROVED (trace identity, all quadrants) / CHECKED (entry forms) | `proofs/cross_m1_m2.md` |
| QFLIP | Quadrant flip = perfect-square swap: the π-shift is w → iw, i.e. t ↦ −1/t; α = sgn(Im ζ) selects which Euler square (w+w⁻¹)² / (w−w⁻¹)² sits upstairs; at block level M1_QIII = diag(1,−1)·M2_QI·diag(1,−1). This is the mechanism Kit named "self propagation" | PROVED (symbolic in w and t) | `proofs/quadrant_flip_polar.md` |
| SO3EXP | SO(3) exponential: Rodrigues closes exp(Y) for ω = (sxp, cxp, uxp), and θ² = 2(sxp² − sxp·cxp + cxp²), the A2 Eisenstein norm of the cores (since uxp = cxp − sxp) | PROVED (structure, θ² identity) / CHECKED (numeric matrix at 30°) | `proofs/expmap_uxp_so3.md` |
| WICKWALL | Wick-on-generators wall: Wick-rotating the 8 cycle generators gives a unitary conjugate of so(8), dim 28, NOT su(8). Wick does not commute with generation; the correct su(8) route Wicks the full sl(8) basis (28 keep, 35 × i) | PROVED (conjugacy identity) / CHECKED (dim-28 closure, compactness) | `proofs/wick_sl8_wall.md` |
| MREF4 | M_refined generates sl(4,ℝ): raw M_refined closes at dim 16 (gl(4,ℝ)); its traceless correction closes at dim 15, Killing (9,6,0); det factors exactly in t = tan(x/2), vanishing at t² = 1/5 in Quadrant I | PROVED (Sage 9.5 exact: ℚ(t) for det, ℚ̄ for closure/Killing) | `proofs/m_refined_sl4.md` |

### Walls — scoped negative results

E108 / E109 / E110 are NOT walls for sl(n) in general. They are walls for
one specific construction: the sampled minimal chain that reuses the same
2 functions in overlapping sl(2) positions. That chain stalls at rank 10
of 15 for sl(4) — with 4 samples (E108), with 8 samples (E109), and with
Re/Im complex variants (E110) — because the sampled matrices live in a
4-dimensional function space (E103, PROVED), so more samples and complex
packaging add no real independence. The cycle construction bypasses the
wall entirely: its generators are placed at distinct cycle positions with
nonzero units, not sampled from the T_a(x) curve. The sampled minimal
chain terminates at n = 3; the cycle construction does not terminate.

### Explicitly not theorems (Lie campaign)

- **Maxwell exercise** (`proofs/maxwell_exercise.md`): operator identities
  PROVED; the E/B identification is a declared dictionary (ASSERTED);
  Maxwell's equations as spacetime PDEs are NOT derived. A representation
  witness only, per Kit's certified spine.
- **O_prop**: OPEN. No 12×12 candidate passed a leakage test; the proved
  transport laws model its shape as selection criteria only. The R(θ) = 0
  leakage test as previously offered is NOT executable (R has no
  definition in the record — run 233 returned BLOCKED); the pending item
  is a definition of R from Kit, plus the reframed choices in
  `proofs/oprop_next_step_correction_DRAFT.md`.
- **Universal parameter theorem** for general traceless 3×3 operator
  families, and the **Z/W complex-structure moduli** study: OPEN
  (research order approved 2026-10-12: parameter theorem first, then Z/W).

## Withheld — strong evidence only (not proved)

Excluded from the proved count. Admit only if a complete real-algebraic proof is established.

| ID | Statement | Proof |
|----|-----------|-------|
| T4* | X²=−I: 2 solutions over ℚ̄, both non-real; 0 real found — strong evidence, NOT a complete real-algebraic proof of emptiness | `proofs/q6b_groebner.py`, `proofs/q6b_solutions.py` |
