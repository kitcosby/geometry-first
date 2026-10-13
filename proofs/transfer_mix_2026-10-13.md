# TRANSFER — The 2×2 Mix, the FlatWave Cartan Extension, and the Saw Forms

**Date:** 2026-10-13 · **Source:** R Theory trigonometry side chat (Kit Cosby, directed session)
**Status key:** PROVED = exact derivation or exact rational/integer computation ·
CHECKED = verified numerically at 30-digit precision (mpmath) or by finite scan ·
All computations this session: desktop SymPy (exact) + mpmath (30-digit).

## 0. Definitions (all previously established)

Core operators, on the domain where sin x · cos x ≠ 0:

- srx = |csc x| + cot x, sxp = 1/srx = |csc x| − cot x
- cxp = |sec x| + tan x, crx = 1/cxp = |sec x| − tan x
- urx = srx − crx, uxp = cxp − sxp (the double-angle pair; urx + uxp = 4/sin(2x))
- FW = FlatWave = sgn(sin 2x) ∈ {+1, −1} on each octant (pattern over octants 0–7: +,+,−,−,+,+,−,−)
- sawup p = 1/urx, sawdown q = 1/uxp (Foundation D.4); p + q = FW, pq = |sin 2x|/4

Previously proved facts used: run 245 (urx·uxp = 4|csc 2x|, hence urx·uxp ≥ 4 and neither
factor ever vanishes); run 366 (T(x) generates sl(3,ℝ)); run 368 (Cartan split of sl(3,ℝ):
𝔨 = so(3) dim 3, 𝔭 dim 5, Killing signature (5,3), θ(X) = −Xᵀ, C ∈ 𝔨, H, D, K ∈ 𝔭);
H4/H6 (quarter-turn family swap; uxp = urx/(srx·crx) in all quadrants).

Notation: H = diag(1,−1,0), K = diag(0,1,−1), C = E₂₃−E₃₂, D = E₂₃+E₃₂ (run-366 basis of sl(3,ℝ)).

## 1. The two single pieces and their pair algebra — PROVED (exact)

Inside T(x) = urx·H + srx·E₁₂ + cxp·E₂₁ + uxp·C:

- The urx piece (direction H) alone and the uxp piece (direction C) alone each have
  Lie closure of dimension 1 (abelian).
- Together: closure{H, C} = span{H, C, D, K}, dimension 4, with [H,C] = −D, [C,D] = 2K,
  [K,C] = 2D, [K,D] = 2C. The triple {K, C, D} is sl(2,ℝ) (the 2–3 block, so(2,1)), and
  2H + K = diag(2,−1,−1) is central. Hence closure{H, C} ≅ sl(2,ℝ) ⊕ ℝ.
- Position in the run-368 Cartan split: C ∈ 𝔨 (the 3), H, D, K ∈ 𝔭 (the 5) — the pair's
  4 dimensions land 1 in the compact part and 3 in the noncompact part.
- The pair alone does not generate sl(3,ℝ); the srx/cxp pieces are also required (run 366).

## 2. Kit's 2×2 — M = [[srx, cxp],[sxp, crx]] = [[s, c],[1/s, 1/c]], s = srx, c = cxp

- **det M = s/c − c/s — PROVED** (definitions only). det M = 0 iff s = c.
- **Rank-1 collapse at handoffs — PROVED at π/4, extended by the proved quarter-turn
  symmetry:** at x = π/4, s = c = √2+1, the columns coincide, M = [[√2+1, √2+1],[√2−1, √2−1]].
  Same collapse at every x ≡ π/4 (mod π/2); M is invertible strictly inside octants.
  |det M| at octant midpoints = 3.061467459, sign alternating by octant, π-periodic (CHECKED).
- **span{M(x)} = gl(2,ℝ) and span{traceless parts} = sl(2,ℝ) — CHECKED** (sampled rank 4
  resp. 3, capped by ambient dimension). The traceless part is
  N(x) = [[urx/2, cxp],[sxp, −urx/2]].
- **det N(x) = −(urx²/4 + sxp·cxp) < 0 at every x — PROVED:** all four core operators are
  positive (|csc| > |cot|, |sec| > |tan|), so sxp·cxp > 0. N(x) is therefore always
  hyperbolic (real eigenvalues of opposite sign), never elliptic.

## 3. The FlatWave Cartan extension — A(x) = [[urx/2, cxp, 0],[sxp, −urx/2, FW],[0, −FW, 0]]

- The direction set {H, E₁₂, E₂₁, C} closes at dimension **8 = sl(3,ℝ)** — PROVED (exact
  integer closure; same bracket tree as run 366; generation is weight-blind, and FW = ±1 ≠ 0
  inside each octant is the only property of the edge used).
- Family-span step — CHECKED: coefficient functions {urx, cxp, sxp} are independent on an
  octant (sample determinant 0.21816368 ≠ 0).
- **Cartan rank 2:** {H, K}, K = [C,D]/2 in the closure, [H,K] = 0 — exactly one rank more
  than the 2×2's rank 1.
- **Oscillation:** the FW-weighted direction C is the compact so(3) generator (run 368),
  eigenvalues 0, ±i — a pure rotation, of a kind the 2×2 alone cannot contain (§2, N always
  hyperbolic). FW reverses the rotation's sign quadrant to quadrant; sign-flipping C is
  conjugation by diag(1,1,−1) (the quadrant-flip / self-propagation mechanism).
  Limit, disclosed: A(x)'s own eigenvalues at all 8 octant midpoints are real (CHECKED);
  the oscillation lives in the algebra element FW·C, not in A's midpoint spectrum.
- **Bonus — PROVED (exact closure):** a second FW edge into a 4th coordinate closes at
  dimension 15 = sl(4,ℝ), rank 3. Each added FW edge buys one rank.
- **det A = urx/2 — PROVED** (cofactor expansion; FW² = 1 kills all other terms).

## 4. The diagonal double-angle variant — B(x) = [[urx, srx, 0],[cxp, uxp, FW],[0, −FW, −(urx+uxp)]]

- Closure of directions {D₁ = diag(1,0,−1), D₂ = diag(0,1,−1), E₁₂, E₂₁, C} = sl(3,ℝ),
  dimension 8 — PROVED (exact). Coefficient independence {urx, uxp, srx, cxp} CHECKED
  (sample determinant 0.010998671 ≠ 0).
- **span{D₁, D₂} = span{H, K} — PROVED (exact):** the double angles re-coordinatize the
  *same* Cartan subalgebra; no new algebra, no new rank. The gain is visibility: the full
  rank-2 Cartan sits on the diagonal, and (3,3) = −(urx+uxp) = −4/sin(2x).
- Cost: det B has no single-operator closed form (values 14.9210225, 11.8595551,
  26.8889923, 29.9504598 at the quadrant representatives — CHECKED), vs det A = urx/2.
  B's midpoint spectra are also all real (CHECKED). B dissolves the 2×2 block that A keeps.

## 5. The uxp/2 variant — Ã(x) = [[uxp/2, srx, 0],[crx, −uxp/2, FW],[0, −FW, 0]]

- Full twin of A: coefficient functions {uxp, srx, crx} independent (CHECKED, sample
  determinant 1.1778068 ≠ 0); spans {H, E₁₂, E₂₁, C}; closure sl(3,ℝ); Cartan {H, K}.
- **det Ã = uxp/2 — PROVED** (same cofactor derivation as det A).
- Spectra swap with A's between neighboring octants (CHECKED) — the s↔c family exchange.
- **Pairing identity — PROVED** (exact algebra + run 245): tr(A·Ã) = (urx·uxp)/2 = 2|csc 2x|;
  value 2√2 at octant midpoints.
- **DISCLOSED TRAP:** the naive placement with off-diagonals (sxp, cxp) has coefficient
  functions {uxp, sxp, cxp} that are linearly dependent *by definition* (uxp = cxp − sxp;
  sample determinant exactly 0). Its family span collapses to 3 dimensions
  ({E₁₂ − H/2, E₂₁ + H/2, C}). Bracket closure still recovers sl(3,ℝ) (PROVED, exact),
  but the clean span argument is lost. Use the srx/crx placement (§5 as written).

## 6. Inverses — closed forms PROVED by cofactor algebra, verified CHECKED (~1e-30 residuals, 4 sample points)

- **M⁻¹** = 1/(s² − c²) · [[s, −s·c²],[−c, s²·c]]. Exists iff s ≠ c — fails exactly at the
  handoffs (§2).
- **A-shape inverse:** for [[a, b, 0],[c, −a, FW],[0, −FW, 0]] (determinant a):
  inverse = [[1/a, 0, b·FW/a],[0, 0, −FW],[−c·FW/a, FW, (−a² − bc)/a]].
  (−a² − bc is the determinant of the 2×2 traceless block.)
  - A⁻¹: a = urx/2, b = cxp, c = sxp.
  - Ã⁻¹: a = uxp/2, b = srx, c = crx.
  - **A⁻¹ and Ã⁻¹ exist everywhere on the domain — PROVED:** determinants urx/2, uxp/2,
    and urx·uxp = 4|csc 2x| ≥ 4 (run 245) forbids either factor vanishing.
- **B⁻¹** = adjugate(B)/det B, cofactors:
  C₁₁ = 1 − uxp(urx+uxp), C₁₂ = cxp(urx+uxp), C₁₃ = −cxp·FW,
  C₂₁ = srx(urx+uxp), C₂₂ = −urx(urx+uxp), C₂₃ = urx·FW,
  C₃₁ = srx·FW, C₃₂ = −urx·FW, C₃₃ = urx·uxp − srx·cxp.
  Existence: det B ≠ 0 on a 1,600-point scan of the full domain (min |det B| ≈ 9.21) —
  CHECKED only; no analytic proof yet.

## 7. The mix in saw terms — p = sawup = 1/urx, q = sawdown = 1/uxp, ε = p + q = FW, R = √(1 + 4pq) = √(1 + |sin 2x|)

**Dictionary (principal chart, ε = +1) — PROVED by exact derivation**
(srx − crx = 1/p and cxp − sxp = 1/q plus the reciprocal pairs force):
srx = (1+R)/(2p), cxp = (1+R)/(2q), sxp = (R−1)/(2q), crx = (R−1)/(2p).
On ε = −1 charts the (1+R)/(R−1) assignments swap with signs (CHECKED at quadrant
representatives) — chart-local, the same status λ has. The FW edge is ε = p + q.

**Matrices (principal chart; verified CHECKED on QI samples, residuals ~1e-30):**
- M = [[(1+R)/(2p), (1+R)/(2q)],[(R−1)/(2q), (R−1)/(2p)]]
- A = [[1/(2p), (1+R)/(2q), 0],[(R−1)/(2q), −1/(2p), ε],[0, −ε, 0]]
- Ã = A with p↔q; B = [[1/p, (1+R)/(2p), 0],[(1+R)/(2q), 1/q, ε],[0, −ε, −ε/(pq)]]

**Global scalar identities in saw form — PROVED** (definitions + H6 + INVPAIR + run 245;
verified at all 8 octant midpoints, residuals ~1e-31):
- det M = q/p − p/q · det A = 1/(2p) · det Ã = 1/(2q)
- sxp·cxp = p/q · srx·crx = q/p · tr(A·Ã) = 1/(2pq) · pq = |sin 2x|/4

**Inverses in saw form (principal chart):**
- A⁻¹ = [[2p, 0, εp(1+R)/q],[0, 0, −ε],[−εp(R−1)/q, ε, −(1/(2p) + 2p²/q)]] (verified CHECKED)
- Ã⁻¹ by the p↔q mirror.
- M⁻¹ = [pq/(ε(q−p))] · [[(R−1)/(2p), −(1+R)/(2q)], [−(R−1)/(2q), (1+R)/(2p)]]

## 8. What is NOT claimed

- No new Lie algebra appears anywhere in this session: every closure is sl(2,ℝ), sl(3,ℝ),
  or sl(4,ℝ) already in the proved record. The content is coordinatization, invariants,
  and the saw parametrization.
- det B ≠ 0 globally is CHECKED (scan), not proved.
- Family-span steps resting on sampled coefficient determinants are CHECKED; the bracket
  closures of the resulting direction sets are exact.
- No physics identification is asserted for any matrix here (algebra-first discipline).

*End of transfer.*
