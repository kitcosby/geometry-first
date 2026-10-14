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
