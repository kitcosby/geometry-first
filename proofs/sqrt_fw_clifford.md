# √FlatWave: the Double Cover, Its Clifford Pair, and the z′ Junction

**Date:** 2026-10-14 · **Status: constructions VERIFIED exactly (finite
octant model + the PROVED record); the z′ identification is REFUTED as
stated and replaced by the precise relation in §4.**
Sources for the retrieved definitions: `SYMBOLIC_CALIBRATION.md` Results
11–13; `RUN143_RESULT11_EPS_BRANCH.md`; `RUN144_RESULT13_BIT_TABLE.md`.

## 1. The square root of the FlatWave

$FW = \operatorname{sgn}(\sin 2x)$ is constant per quadrant:
$(+1, -1, +1, -1)$. Its consistent square root around the circle is

$$g = (1,\; i,\; -1,\; -i) \quad \text{per quadrant} \qquad
(\text{or the conjugate branch } \bar g),$$

with $g^2 = FW$ pointwise (verified at all 8 octants). But
$g(x+\pi) = -g(x)$: **$g$ is antiperiodic over the operator period
$\pi$** — it lives on the double cover, period $2\pi$. This is an
explicit function realizing the antiperiodic sector that SC-7 PROVED
the fermionic twist requires ($\pi$-periodic gluing forces holonomy
$+I$; the twist needs $\pi$-periodicity broken).

Transport relations (verified exactly; $T$ = quarter-turn, $S = T^2$,
$G$ = multiplication by $g$):

$$T G = \mp i\, G T \quad (\text{sign by transport convention}), \qquad
S G = -\, G S .$$

The quarter-turn commutes with $\sqrt{FW}$-multiplication only up to
the phase $\pm i$ (cycling with period 4 quadrants); the half-turn
anticommutes with it, commutator the central $-1$. The lifted symmetry
group is dicyclic (binary dihedral), with the quaternion $Q_8$ pattern
over the Klein subgroup.

## 2. The Clifford pair

With $S$ = $\pi$-transport and $G$ as above (all verified):

$$S^2 = 1,\qquad G^2 = M_{FW},\qquad \{S, G\} = 0,\qquad
(SG)^2 = -M_{FW},$$

$M_{FW}$ central. This is a 2-generator Clifford algebra
$\{1, S, G, SG\}$ canonically attached to the quadrant sign. The
function space splits into $FW = \varepsilon$ sectors ($S$ preserves
them); the generator squares there are $(+1, \varepsilon)$ — the
$Cl(2,0)$ pattern on $FW = +1$ quadrants, the $Cl(1,1)$ pattern on
$FW = -1$ quadrants (over $\mathbb C$ both are $M(2,\mathbb C)$).
**The FlatWave sign is the signature of this Clifford algebra.**
The dicyclic lift is a finite Spin-type group — the double cover
arrives inside a Clifford algebra, its canonical home.

Reflection law (verified; corrects a first guess): $FW$ is
reflection-odd, so reflection sends $g \mapsto -i\,\bar g$, a square
root of $-FW$. The branch choice is reflection-sensitive with this
exact law; the choice itself remains a stipulation.

## 3. What z′ is (retrieved, not assumed)

Result 11 (PROVED): on each octant all operators are rational
functions of $z = srx(x)$ alone, with branch parameter
$\varepsilon = FW$ entering through

$$cxp = \varepsilon \cdot M_\varepsilon(z),\qquad
M_{\pm 1}(z) = \frac{z \pm 1}{z \mp 1},$$

and $\varepsilon$ is not independent data: $\varepsilon =
\operatorname{sgn}(z - 1)$ (the $|z| \gtreqless 1$ encoding, PROVED).
The Cl(2,2) 4-spinor is the graded tensor product $\psi_\ell \otimes
\psi_\eta$ with components (RUN144, PROVED table)

$$c_1 = \sqrt{srx \cdot cxp} = e^{(\ell+\eta)/2},\quad
c_2 = \sqrt{srx/cxp} = e^{(\ell-\eta)/2},\quad
c_3 = \sqrt{cxp/srx},\quad c_4 = 1/\sqrt{srx \cdot cxp},$$

principal (positive) roots — the components are the half-rapidity
exponentials. $z'$ is the **imaginary sector's own generator**: under
Kit's A0.3 (ASSUMED isomorphism), the partner sector carries the same
definitions in its own variable ("same definitions, applied to $z'$",
Result 13; its instance of the closure is conditional, not probed,
Result 11 doc). The 16-spinor components are functions of $(z, z')$;
Result 12: sector holonomies multiply, $h_{16} = h_R \cdot h_I$, and
the 16 is fermionic iff exactly one sector is twisted. The twist
itself ($\Psi(x+2\pi) = -\Psi(x)$ in a twisted sector) is Kit's
ASSERTED seam prescription, not a derived formula.

## 4. The junction — verdict

**$g$ is not $z'$, and $g$ is not the spinor twist factor.** Three
distinct objects, with distinct period structures:

| object | period | nature |
|---|---|---|
| operators, $FW$, principal components $c_i$ | $\pi$ | values / branch labels |
| $g = \sqrt{FW}$ | $2\pi$ (antiperiod $\pi$) | $\mu_4$-valued, canonical up to conjugation |
| twisted spinor $\Psi_{tw}$ (Result 12) | $4\pi$ (antiperiod $2\pi$) | ASSERTED holonomy choice |

The precise relations that DO hold:

1. **$g^2 = \varepsilon$ is the branch label inside the component
   formulas.** The $(z, z')$ formulas carry $\varepsilon$ through
   $cxp = \varepsilon M_\varepsilon(z)$ (and its primed twin). $g$ is
   the square root of that branch parameter — the components see
   $\varepsilon$; $g$ sees it at half weight. But $\varepsilon$ is
   itself $\operatorname{sgn}(z-1)$, so no independent $g$-factor
   enters the components: they are functions of $z$ (resp. $z'$)
   alone per sector.
2. **$g$ is the canonical middle rung of a cover ladder** the record
   previously had only the endpoints of: period-$\pi$ tensor level
   (operators, principal spinor) → $g$'s $2\pi$ level → Result 12's
   $4\pi$ twisted-spinor level. A function-level twist factor for
   Result 12 would have to be a fourth-root-of-sign object
   ($\mu_8$-valued, antiperiod $2\pi$); no canonical one is on
   record, and Result 12 treats the twist as a discrete choice.
3. **$z'$ remains an independent generator** (of an ASSUMED-isomorphic
   sector). Nothing in the record defines $z'$ from $z$, $\varepsilon$,
   or $g$; each sector would carry its OWN $g$ (root of its own
   $\varepsilon$, resp. $\varepsilon'$), constructed identically.

## 5. Scope

- VERIFIED here (finite exact computation on proved data): §1, §2,
  the reflection law.
- PROVED on record (cited, not re-derived): Result 11 closure and
  $\varepsilon = \operatorname{sgn}(z-1)$; RUN144 component formulas
  and $(P,G)$ table; Result 12 holonomy product rule.
- ASSERTED on record: A0.3 sector isomorphism; the twisted-seam
  prescription.
- NOT claimed: any identification of $g$ with $z'$ or with the
  Result-12 twist; any bearing on the Ansatz-4 obstruction, the
  $Cl(4,4)$ 16's origin, or WALL-FERMION (all unchanged).

## Addendum 2026-10-15 — saw antiperiodicity and the cube root

### 6. The saws are antiperiodic (PROVED)

From the RUN144 branch lemma (PROVED), the quarter-turn swaps the
primitive pairs:

$$srx(x+\tfrac\pi2) = crx(x), \qquad cxp(x+\tfrac\pi2) = sxp(x)$$

(verified exactly at the octant midpoints), and reciprocation gives
the other two swaps. Hence both UNA differences negate:

$$urx(x+\tfrac\pi2) = -urx(x), \qquad uxp(x+\tfrac\pi2) = -uxp(x),$$

$$\boxed{\;sawup(x+\tfrac\pi2) = -sawup(x), \qquad
sawdown(x+\tfrac\pi2) = -sawdown(x).\;}$$

The saws are real-valued antiperiodic objects with antiperiod $\pi/2$
— they do at the quarter-turn what $g$ does at the half-turn, with no
$i$ imported. With $T$ = quarter-turn transport, $U$ = multiplication
by $sawup$, $D$ = multiplication by $sawdown$:

$$TU = -UT, \qquad TD = -DT, \qquad [T,U] = 2TU, \quad [T,D] = 2TD,$$

both commutators nonzero. Factor bookkeeping: an anticommuting pair's
commutator carries $|1 - (-1)| = 2$; and $U + D = M_{FW}$, whose
values are exactly $\pm 1$. (Both saw commutators carry the 2; the
$\pm 1$ is the FlatWave sum, not a per-saw split.)

### 7. The cube root of the FlatWave (VERIFIED)

Construction (Kit): adjust the sign's period if desired, then take
the cube root. Let $\rho = e^{i\pi/3}$, so $\rho^3 = -1$ — the
per-sector sign flip, cube-rooted. Set $h(\text{step } k) = \rho^k$.
Verified exactly:

- $h^3 = FW$ on the quadrant grid, and $h^3 = \operatorname{sgn}(\sin
  8x)$ on the 16-sector grid of $FW(4x)$. The $4x$ adjustment rescales
  the grid only; the root's structure is unchanged.
- $h$ has period 6 steps: $3\pi$ — a **triple cover** (on the $4x$
  grid, $3\pi/4$ over the $\pi/4$ base).
- Transport is a Weyl clock-shift pair, $T H = \rho^{\pm 1} H T$, and
  the commutator factor is $|1 - \rho| = 1$ **exactly**
  ($\cos\frac\pi3 = \frac12$): the cube-root channel's commutator has
  modulus 1, beside the square root's antipodal factor 2.
- General pattern: the $n$-th root has per-step phase $e^{i\pi/n}
  \neq 1$, so its transport commutator is automatically nonzero at
  every order.
- **Non-closure on the $4\pi$ lap (proved arithmetic):**
  $\rho^8 = e^{2i\pi/3} \neq 1$ — the phase has order 6 and the $4\pi$
  lap has 8 quadrant steps. The cube root's home is the triple cover,
  a cover distinct from both the double cover ($2\pi$, §1) and the
  Result-12 twist ($4\pi$).

Scope of the addendum: §6 PROVED (derivation from the PROVED RUN144
lemma; swaps confirmed by exact evaluation). §7 VERIFIED exactly on
the finite sector models; the non-closure is exact arithmetic. No
identification with $z'$, the Result-12 twist, or any $Cl(4,4)$
structure is claimed; WALL-FERMION unchanged.

## Addendum 2026-10-15 (II) — the 8x level

### 8. The fourth root, the order-8 pseudoscalar, and the binary tower

Let $F_8 = FW(8x) = \operatorname{sgn}(\sin 16x)$: 32 sectors of
width $\pi/16$ around the circle, $F_8 = (-1)^j$ on sector $j$. All
statements below VERIFIED exactly on the 32-sector model.

**Root tower at a single argument.** With $\zeta = e^{i\pi/4}$,
the accumulated roots are

$$G_8(j) = i^j, \qquad Q_8(j) = \zeta^j, \qquad
Q_8^2 = G_8, \quad G_8^2 = F_8, \quad Q_8^4 = F_8.$$

So at the $8x$ level the **fourth root is canonical** (same
accumulation rule as $g$, one level higher): $Q_8$ has period
8 sectors $= \pi/2$ and antiperiod $\pi/4$, with values in $\mu_8$.
(The value group $\mu_8$ is also the spectrum of the sl(8) cycle
matrix, PROVED $x$-independent — recorded here as a shared
order-8 cyclic structure, not an identification.) The sector-step
transport is a Weyl pair: $U Q_8 U^{-1} = \zeta^{-1} Q_8$, and the
antiperiod transport anticommutes with $Q_8$:
$\{\text{shift}_4, Q_8\} = 0$.

**The pseudoscalar has order 8.** $J_8 = \text{shift}_4 \cdot Q_8$
satisfies

$$J_8^2 = -\,(\pi/2\text{-transport}) \cdot G_8, \qquad
J_8^4 = S \cdot M_{F_8}, \qquad J_8^8 = I,$$

with $S$ the $\pi$-transport. Against §2 ($J^2 = -M_{FW}$, order 4)
and the intermediate $2x$ case (computed 2026-10-15, session record:
$J_2^2 = -S \cdot M_{FW(2x)}$, order 4), the pattern is that each
level's $J^2$ is dressed by the level below: bare sign, then
transport × sign, then transport × phase.

**The bits are binary coordinates.** At the 32 sector centers,
$(FW(x), FW(2x), FW(4x), FW(8x))$ are exactly the binary digits of
the sector index: the FlatWave doubling tower is a complete binary
coordinate system on the circle (the Rademacher system), one bit of
position per doubling.

**The 1820 ε-pattern is a Walsh product of the first two bits.**
The 1820 lift's octant ε-pattern
$\{-1, +1, +1, -1, -1, +1, +1, -1\}$ (RUN143) equals, per octant,

$$\varepsilon_{1820}(k) = -\,FW(x)\cdot FW(2x)\Big|_{\text{octant } k}
\quad \text{exactly, all 8 octants.}$$

Pattern-level identity only: the chart-level branch Möbius and the
1820's $M_\varepsilon$ lift remain distinct structures (RUN143's
anti-conflation datum stands); their sign patterns are related by
this product.

**Refinement of §4.2's wording.** The $\mu_8$ accumulation root of
$FW$ itself (period $4\pi$) is likewise canonical up to conjugation;
what fails is only its identification with Result 12's twist — same
$2\pi$ holonomy $-1$, but phase-valued within the lap, where the
record's twist keeps components real and negates across laps.

*Disclosure:* the raw check $\{\text{shift}_4, Q_8\} = 0$ first
printed False under sympy structural equality on $e^{i\pi/4}$
powers; rerun with $\zeta$ in algebraic form confirms exactly 0
(presentation-layer artifact, the math was re-verified, not
assumed).

Scope: finite-model VERIFIED as stated; the binary-digit and
Walsh-product identities are exact finite arithmetic on PROVED sign
patterns. No identification of $Q_8$ with the sl(8) cycle, the 1820
lift, or any spinor structure is claimed.

## Addendum 2026-10-15 (III) — the 2x square root, and FW^{1/4}

### 9. The 2x square root (VERIFIED)

$F_2 := FW(2x) = \operatorname{sgn}(\sin 4x) = (-1)^k$ per octant.
$F_2$ is the **saw-dominance sign**: $F_2 = +1$ exactly where
$|sawdown| > |sawup|$ (exact at all 8 octant midpoints; the ties at
the $45^\circ$ points are its zeros). Hence the complete octant
label mod $\pi$ is $(FW(x), FW(2x))$ — both entries FlatWaves.

The square root interlocks with doubling:
$\sqrt{FW(2x)} = g(2x)$ (substitution: $(g \circ 2)^2 = FW \circ 2$;
confirmed at 8 exact interior points). Per octant $H_2$ has values
$i^k$; period $\pi$, antiperiod $\pi/2$. The octant-step transport
is a Weyl pair, $T_8 H_2 T_8^{-1} = i^{-1} H_2$, and with $A$ the
$\pi/2$-transport, $\{A, H_2\} = 0$, $[A, H_2] = 2 A H_2 \neq 0$.

Contrast with §2: the pseudoscalar $J_2 = A H_2$ has

$$J_2^2 = -\,S \cdot M_{F_2}, \qquad J_2^4 = I,$$

$S$ the $\pi$-transport. The clean sector complex structure
$J^2 = -M_{FW}$ is special to the $1x$ root; at $2x$ the square is
the doubled sign dressed with the $\pi$-transport. (This is the
intermediate case cited in §8.)

### 10. The fourth root FW^{1/4} (VERIFIED)

On the $4\pi$ lap (8 quadrant cells), with $\zeta = e^{i\pi/4}$,
the accumulated fourth root is $q(k) = \zeta^k$ (conjugate branch
$\bar q$). Verified exactly on the lap model:

- **Tower:** $M_q^2 = M_g$, $M_q^4 = M_{FW}$ — $q$ squares to the
  square root and fourth-powers to the FlatWave, globally on the
  lap. $q(x + 2\pi) = -q(x)$: period $4\pi$, antiperiod $2\pi$.
- **Transport:** quarter-turn Weyl pair with the finest phase yet,
  $T M_q T^{-1} = \zeta^{-1} M_q$; and with $\Lambda$ the
  $2\pi$-transport, $\{\Lambda, M_q\} = 0$,
  $[\Lambda, M_q] = 2\Lambda M_q \neq 0$.
- **Pseudoscalar:** $J_q = \Lambda M_q$ satisfies
  $J_q^2 = -M_g$ — the fourth root's pseudoscalar squares to
  (minus) the *previous root*, where §2's squared to the sign
  itself. $J_q^4 = M_{FW}$, $J_q^8 = I$: order 8, the same order
  §8 reached by refining the argument instead of the root.
- **Group:** $\langle \Lambda, M_q \rangle$ has order **32**
  (exact enumeration on the $(a,b,c)$ normal forms
  $\Lambda^a \cdot \mathrm{diag}(\zeta^{bk+c})$): a dicyclic-type
  central extension, one level above the $Q_8$ pattern of §1.
- **Reflection:** $q \mapsto \zeta^7 \bar q$, a fourth root of
  $-FW$ — consistent with $FW$ being reflection-odd (§2's law for
  $g$ is the level below).

**Not the Result-12 twist.** $q \cdot \Psi_{\mathrm{prin}}$ obeys
the twisted lap law $\Psi(x+2\pi) = -\Psi(x)$, but $q$ is
phase-valued *within* the lap (factor $\zeta$ on QII), where the
record's twist keeps first-lap components real and negates only
across laps (the real lap-step $\tau = \pm 1$, $\tau^2 = 1$, while
$q^2 = g \neq \tau$). Same holonomy, different object; the twist
identification remains declined, as in §4.

Remark: at fixed argument the root tower has periods
$\pi \cdot 2^k$ and value groups $\mu_{2^{k+1}}$ — the 2-adic
tower of covers of the operator circle.

Scope: §§9–10 VERIFIED exactly on the finite sector/lap models
from the PROVED sign patterns. No identification of $q$ with $z'$,
the Result-12 twist, or any $Cl(4,4)$ structure is claimed.
