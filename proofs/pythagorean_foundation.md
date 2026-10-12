# R-Theory: Pythagorean Foundation of the Bounded Phase Operators

## IDENTIFICATION DIVISION

**PROGRAM-ID:** R-THEORY-PYTHAGOREAN-CORE

**PURPOSE:** To derive the exact algebraic properties of the R-operators from the Pythagorean identity sin²x + cos²x = 1, together with standard definitions of the reciprocal and quotient trigonometric functions.

**NON-PURPOSE:** This division does not establish:
* physical conservation of energy, mass, or information;
* uniqueness of the R-operator family;
* empirical applicability;
* replacement of exponential or logarithmic laws;
* a law of temporal evolution.

Those require additional arguments beyond the Pythagorean foundation.

---

## ENVIRONMENT DIVISION

**Definition E.1 — Phase variable:** Let x be a real phase variable measured in radians, with periodic identification x ≡ x + 2π.

**Definition E.2 — Admissible domain:** D = {x ∈ R : sin x ≠ 0 and cos x ≠ 0} = R \ {kπ/2 : k ∈ Z}.

**Definition E.3 — Unit-circle coordinates:** s = sin x, c = cos x. P(x) = (c,s) lies on the unit circle.

---

## AXIOM DIVISION

**Inherited Theorem A.1 — Pythagorean theorem:** For a right triangle with legs a,b and hypotenuse h, a² + b² = h².

**Corollary A.2 — Unit-circle Pythagorean identity:** s² + c² = 1, i.e. sin²x + cos²x = 1.

---

## DATA DIVISION

**Definition D.1 — Absolute reciprocal functions:** M = |sec x| = 1/|c|, N = |csc x| = 1/|s|. Both > 0.

**Definition D.2 — Primitive displayed operators:**
* srx(x) = N + cot x (Sininverex)
* sxp(x) = 1/srx(x) (Sinexpon)
* cxp(x) = M + tan x (Cosinexpon)
* crx(x) = 1/cxp(x) (Cosininverex)

**Definition D.3 — UNA operators:**
* urx(x) = srx(x) − crx(x) (Unainverex)
* uxp(x) = cxp(x) − sxp(x) (Unaexpon)

**Definition D.4 — Saw operators:** saw_r(x) = 1/urx(x), saw_x(x) = 1/uxp(x).

**Definition D.5 — Flatwave:** Flatwave(x) = saw_r(x) + saw_x(x) = 1/urx(x) + 1/uxp(x).

---

## PROCEDURE DIVISION

### SECTION 1 — Reciprocal Trigonometric Identities

**Lemma 1.1:** M² − tan²x = 1. (Proof: M² = sec²x, sec²x − tan²x = (1−s²)/c² = c²/c² = 1.)

**Lemma 1.2:** N² − cot²x = 1. (Proof: analogous.)

**Corollary 1.3:** M > |tan x|, N > |cot x|, so M±tan x > 0, N±cot x > 0. No denominator vanishes on D.

### SECTION 2 — Exact Reciprocal Conjugates

**Theorem 2.1:** 1/(M+tan x) = M − tan x. (From (M+tan x)(M−tan x) = 1.)

**Theorem 2.2:** 1/(N+cot x) = N − cot x.

**Corollary 2.3:** crx(x) = M − tan x, sxp(x) = N − cot x.

**Corollary 2.4:** cxp(x)·crx(x) = 1, srx(x)·sxp(x) = 1.

**Corollary 2.5:** srx, sxp, cxp, crx > 0 on D. Signed direction first appears in UNA differences.

### SECTION 3 — Sign Variables and Branch Forms

**Definition 3.1:** α = sgn(s), β = sgn(c). α,β ∈ {−1,+1}. N = α/s, M = β/c.

**Lemma 3.2:** srx = (α+c)/s, sxp = (α−c)/s, cxp = (β+s)/c, crx = (β−s)/c.

### SECTION 4 — UNA Structure

**Definition 4.1:** p = αc − βs.

**Lemma 4.2:** urx(x) = (1+p)/sc, uxp(x) = (1−p)/sc.
(Proof: urx = (α+c)/s − (β−s)/c = [c(α+c) − s(β−s)]/sc = (1+αc−βs)/sc.)

**Lemma 4.3:** |p| < 1. (Proof: p² = 1 − 2|sc| < 1 since |sc| > 0 on D.)

**Corollary 4.4:** 1+p > 0, 1−p > 0.

**Theorem 4.5:** sgn(urx) = sgn(uxp) = sgn(sc) = sgn(sin 2x).

### SECTION 5 — Double-Angle Closure

**Theorem 5.1:** urx(x) + uxp(x) = 4/sin(2x).
(Proof: (1+p)/sc + (1−p)/sc = 2/sc = 4/sin(2x).)

**Corollary 5.2:** 1/(urx+uxp) = sin(2x)/4.

### SECTION 6 — UNA Product

**Theorem 6.1:** urx(x)·uxp(x) = 2/|sc| = 4/|sin 2x|.
(Proof: (1−p²)/s²c², with 1−p² = 2|sc|.)

### SECTION 7 — Flatwave

**Theorem 7.1:** Flatwave(x) = sgn(sin 2x).
(Proof: (urx+uxp)/(urx·uxp) = (2/sc)/(2/|sc|) = |sc|/sc.)

**Corollary 7.2:** |Flatwave(x)| = 1.

**Corollary 7.3:** Flatwave = +1 in Quadrants I/III, −1 in Quadrants II/IV.

**Corollary 7.4:** d/dx Flatwave(x) = 0 on every open quadrant.

### SECTION 8 — Saw Partition

**Theorem 8.1:** saw_r(x) + saw_x(x) = sgn(sin 2x).

**Corollary 8.2:** |saw_r(x)| + |saw_x(x)| = 1.
(Proof: urx, uxp same sign ⇒ |saw_r + saw_x| = |saw_r| + |saw_x| = |sgn| = 1.)

### SECTION 9 — Half-Angle Representation

**Theorem 9.1 (Quadrant I):** srx(x) = cot(x/2), sxp(x) = tan(x/2), cxp(x) = tan(π/4 + x/2), crx(x) = tan(π/4 − x/2).

### SECTION 10 — Octant Dominance

**Theorem 10.1:** For 0 < x < π/4, srx is the unique dominant operator.

**Theorem 10.2:** For π/4 < x < π/2, cxp is dominant.

**Lemma 10.3 — Quarter-turn permutation:** srx(x+π/2) = crx(x), cxp(x+π/2) = sxp(x), etc.

**Theorem 10.4 — Dominance cycle:** srx → cxp → crx → sxp → srx over successive octants.

| Octant | Dominant |
|--------|----------|
| (0,π/4) | srx |
| (π/4,π/2) | cxp |
| (π/2,3π/4) | crx |
| (3π/4,π) | sxp |
| (π,5π/4) | srx |
| (5π/4,3π/2) | cxp |
| (3π/2,7π/4) | crx |
| (7π/4,2π) | sxp |

---

## CLOSURE DIVISION

From the Pythagorean identity alone, rigorously established:

1. Reciprocal conjugates: crx = M−tan x, sxp = N−cot x.
2. Positivity of srx, sxp, cxp, crx.
3. Reciprocal pairing: cxp·crx = 1, srx·sxp = 1.
4. UNA forms: urx = (1+p)/sc, uxp = (1−p)/sc.
5. Common sign: sgn(urx) = sgn(uxp) = sgn(sin 2x).
6. Double-angle: urx + uxp = 4/sin(2x).
7. UNA product: urx·uxp = 4/|sin 2x|.
8. Flatwave: Flatwave = sgn(sin 2x).
9. Saw partitions: saw_r + saw_x = sgn(sin 2x), |saw_r| + |saw_x| = 1.
10. Octant dominance cycle: srx → cxp → crx → sxp.

These establish the operator geometry. They do not yet establish physical interpretation.

---
*Authored by Kit Cosby. Archived 2026-10-12.*
