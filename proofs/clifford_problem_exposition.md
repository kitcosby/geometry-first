# The Clifford Problem for the FlatWave Operator System: An Exposition

**Kit Cosby — working exposition for peer review, 2026-10-16**

Status labels are used throughout in the sense of this repository:
**PROVED** (complete proof, analytic or exact-computational with a
dimension/containment argument), **CHECKED** (exact finite verification
without a general proof), **ASSERTED** (stated, not established),
**OPEN** (well-posed, unresolved), **FALSIFIED** (disproved by proof or
counterexample). No physical identifications are claimed anywhere in
this document; the subject is pure algebra and analysis.

---

## 1. The operator system

On the domain $\mathcal D = \mathbb R \setminus \{k\pi/2 : k \in \mathbb Z\}$
define
$$
srx = |\csc x| + \cot x, \qquad cxp = |\sec x| + \tan x,
$$
$$
sxp = 1/srx, \qquad crx = 1/cxp,
$$
$$
urx = srx - crx, \qquad uxp = cxp - sxp \quad \text{(the UNA pair)},
$$
$$
sawup = 1/urx, \qquad sawdown = 1/uxp, \qquad
FW = sawup + sawdown.
$$

**Facts (PROVED; Pythagorean Foundation document, this repository).**
Each of $srx, sxp, cxp, crx$ is strictly positive on $\mathcal D$ and
$\pi$-periodic. One has
$$
urx + uxp = \frac{4}{\sin 2x}, \qquad
urx \cdot uxp = \frac{4}{|\sin 2x|}, \qquad
FW(x) = \operatorname{sgn}(\sin 2x) \in \{+1, -1\},
$$
$$
|sawup| + |sawdown| = 1,
$$
and, with $\cos 2x = \big[(srx + crx) - (sxp + cxp)\big]/(urx + uxp)$,
the quadrature pair $(\sin 2x, \cos 2x)$ is carried by UNA composites:
$\sin 2x = 4/(urx + uxp)$. The UNA operators and the saws are
antiperiodic with antiperiod $\pi/2$: $f(x + \pi/2) = -f(x)$ for
$f \in \{urx, uxp, sawup, sawdown\}$. All operators are $\pi$-periodic,
so no operator-internal quantity distinguishes $x$ from $x + \pi$
(PROVED; the system factors through the quotient circle
$\mathbb R / \pi \mathbb Z$).

The primitives are rapidity exponentials (PROVED): with
$\eta = \ln cxp$ and $\ell = \ln srx$ one has $\sinh \eta = \tan x$,
$\cosh \eta = |\sec x|$, $\tanh \eta = \sin x$ on QI with the signed
global form $\tanh \eta = \operatorname{sgn}(\cos x)\sin x$;
$d\eta/dx = \cosh \eta = |\sec x|$ (SQ1).

## 2. The square root of the FlatWave and its Clifford pair

**Construction (VERIFIED exactly on the quadrant/octant models;
`proofs/sqrt_fw_clifford.md` §§1–2).** There is a function $g$ on the
$2\pi$-circle, $g = (1, i, -1, -i)$ per quadrant (and a conjugate branch
$\bar g$), with
$$
g^2 = FW \ \text{pointwise}, \qquad g(x + \pi) = -g(x).
$$
So $\sqrt{FW}$ exists as an explicit function, but only on the double
cover: the operators live on the $\pi$-circle, their square root lives
on the $2\pi$-circle and is antiperiodic over the operator period.

Let $S$ be the half-turn transport on functions and $G$ multiplication
by $g$. Then (VERIFIED exactly):
$$
S^2 = 1, \qquad G^2 = M_{FW}, \qquad \{S, G\} = 0,
$$
where $M_{FW}$ is multiplication by $FW$, central in this algebra. The
pseudoscalar $J = SG$ satisfies
$$
J^2 = -M_{FW}, \qquad
\exp(\theta J) = A(\theta) + J B(\theta),
$$
with $(A, B) = (\cos\theta, \sin\theta)$ on $FW = +1$ sectors and
$(\cosh\theta, \sinh\theta)$ on $FW = -1$ sectors (PROVED in closed
form by centrality of $M_{FW}$; §14). The function space splits by the
FlatWave sign, and the generator squares $(+1, \varepsilon)$ give the
**Clifford signature**
$$
FW = +1: \mathrm{Cl}(2,0)\text{-pattern}, \qquad
FW = -1: \mathrm{Cl}(1,1)\text{-pattern}.
$$
The FlatWave sign is, in this precise sense, a Clifford signature.

**The tower (§§6–13).** The same construction iterates. $FW^{1/4}$
exists on the $4\pi$-lap with pseudoscalar square $J_q^2 = -M_g$ (the
*previous root* becomes the center); the saws are real antiperiodic
objects anticommuting with the quarter-turn ($[T, U] = 2TU \neq 0$);
the $n$th root carries Weyl phase $e^{i\pi/n}$ with commutator size
$2\sin(\pi/2n) \to 0$, so the tower's classical (commutative) top is
the dyadic solenoid, the inverse limit of the laps (§12). Clifford is
the $n = 2$ slice of a Weyl family. The fourth root
$FW(8x)^{1/4}$ and the diagonal root $FW(8x)^{1/8}$ (which is again
$\pi$-periodic) are VERIFIED on the 32-sector model.

**What the pair $J$ does not do (FALSIFIED).** $J$ is a complex
structure on $FW = +1$ function sectors ($J^2 = -1$ there) and a
product structure on $FW = -1$ sectors ($J^2 = +1$): there is no
canonical complex structure on the full operator system from this
construction (WALL-CPLX, partial).

**WALL-ARG/VAL (J-route FALSIFIED by exact test; general form
OPEN).** The Wick-via-$J$ test (`proofs/wick_via_j.md`) establishes
more than the failure of one route. The argument-space grading and
the value-space grading are distinct $\mathbb Z/2$-structures on the
FlatWave system, and the test separates them:

1. *Argument-space grading.* The half-turn $S$ and the root-phase
$G$ of §2 act on the domain of the operators. Their anticommutator
$\{S, G\} = 0$ lives in the algebra of endomorphisms of function
space $\mathcal F(\mathcal D)$. The pseudoscalar $J = SG$ satisfies
$J^2 = -M_{FW}$ — a relation whose sign is controlled by the value
of $FW$.
2. *Value-space grading.* The Wick rotation acts on the target of
the generators as matrices. It produces the compact real form
$\mathfrak{su}(n) \subset \mathfrak{sl}(n, \mathbb C)$ exactly when
the basis is split into symmetric and antisymmetric parts under the
transpose (the E101 route, PROVED, Killing form $(0,8)$) — at the
cost of importing $i$ and the split.
3. *The identification fails for the $J$-route.* In the Sage-exact
test over QQbar, adjoining the $J$-partners of the
$\mathfrak{sl}(3, \mathbb R)$ generators produces a real algebra
of dimension 16 — the full $\mathfrak{sl}(3, \mathbb C)$ — not
dimension 8 ($\mathfrak{su}(3)$). The argument-space grading closes
over the complexification rather than selecting a real form. The
control arm (transpose split applied to the value space) correctly
yields dimension 8 with Killing signature $(0,8)$.

Consequence, with its scope: the wall is not "we lack a complex
structure." For the $J$-route it is proved that the argument-side
grading does not perform the value-side selection. The general
statement — that *no* operator-native construction identifies the
two gradings — is **OPEN**; it is OQ4 of §5.

## 3. Internal orientation: the saw bit (PROVED 2026-10-16)

Orientation questions cannot be answered by values at a tie: at the
four equilibrium points $x \equiv \pi/4 \pmod{\pi/2}$ the two saws
coincide, $sawup = sawdown = FW/2$. They can be answered by trend.
From $d(sawup)/dx = (srx + crx)/(2\,urx)$ and
$d(sawdown)/dx = -(cxp + sxp)/(2\,uxp)$, together with primitive
positivity and $\operatorname{sgn}(urx) = \operatorname{sgn}(uxp) = FW$:
$$
\operatorname{sgn}(sawup') = FW, \qquad
\operatorname{sgn}(sawdown') = -FW \quad \text{on all of } \mathcal D,
$$
neither slope ever vanishing, with tie-point slopes exactly
$+FW\sqrt2/2$ (sawup) and $-FW\sqrt2/2$ (sawdown). Consequently the
pair (value of either saw, whose sign is $FW$; trend through the tie)
identifies each saw uniquely in every sector: **the operator system
carries a derived internal orientation bit** — derived, not declared,
and not importable from any scalar anchor (a scalar cannot supply
handedness; this is a derivative-sign datum). Circle-wide propagation
is by the system's own transport (the quarter-turn swaps and
$\pi$-periodicity).

## 4. The tie-gluing dichotomy (OPEN; identified 2026-10-16)

Section 3 identifies the *functions*. A *trajectory* riding a saw
faces, at each equilibrium, two continuous continuations:
**pass-through** (exit on the same function; the X-crossing) and
**bounce** (enter as $sawup$, leave as $sawdown$; a $\Lambda$- or
V-path). Value data at the tie cannot prefer one — the values are
equal there — and trend data exhibits the exit arm without settling
the rider's identity under the other convention. The choice between
them is a gluing rule at the seam, of the same species as the
antiperiodic seam prescription $\Psi(x + 2\pi) = -\Psi(x)$
(bounce-type) familiar from spin structures. No gluing class for the
ties is currently declared; the dichotomy is recorded as an open
structural decision, not a missing datum.

**Corollary (resonance at the ties; PROVED, 2026-10-17).** The
symmetric power sums of the six operator pairs (expon, inver, the
two sum pairs, UNA, the saws, the plus-saws) are governed by two
sequences $Q_n, Sp_n$ depending only on $a = |\sin 2x|$: after the
pair-dependent rescaling of the all-orders theorem (audit workspace,
runs 963–964), the power sums $P_n = 2Q_n$ and the difference
quotients $Sp_n$ both obey $P_n = 2P_{n-1} - aP_{n-2}$, with
Chebyshev closed forms $Q_n = a^{n/2} T_n(1/\sqrt a)$,
$Sp_n = a^{(n-1)/2} U_{n-1}(1/\sqrt a)$. The characteristic
equation $r^2 - 2r + a = 0$ has discriminant $4(1 - a)$, which
vanishes exactly at $a = 1$ — that is, exactly at the equilibria
$x \equiv \pi/4 \pmod{\pi/2}$ of §3, where the characteristic root
is double ($r = 1$) and the whole multiple-angle family degenerates
to $Q_n = 1$, $Sp_n = n$ for every $n$. The power-sum family is
therefore resonant precisely at the ties: the same points where the
two saws coincide and where this section's gluing choice must be
made. The corollary identifies the locus of the degeneracy; it
does not prefer pass-through or bounce.

## 5. The bridge problem (OPEN — the Clifford problem proper)

On the Clifford side, the orthogonal ladder built from the same
operators closes $\mathfrak{so}(10)$ (path construction, PROVED;
`proofs/so10_path.md`), and its spinor representation splits under
the chirality operator $\Gamma$ as $32 = 16 \oplus 16'$
(`proofs/spin10_spinor.md`, PROVED). The two halves are a value-tie:
same dimension, same Casimir data, exchanged by the mirror. No
scalar computed inside the Clifford algebra separates them; only an
orientation datum can — $\Gamma$ itself is the volume element, the
product of all gamma matrices at once.

The operator side now possesses a derived orientation datum (§3),
and a global one: the phasor $(\cos 2x, \sin 2x)$ of §1 has winding
number $+2$ about its circle ($-2$ under the mirror) — a property of
the whole closed wave, not of any point. The bridge problem is
therefore a **global-to-global transport**, sharpened as follows.

Two corrections of fact frame the sharpening. First, the punctures
of the quotient circle $\mathcal D_\pi = \mathcal D / \pi\mathbb Z$
are the images of the quadrant *boundaries* $x \equiv 0, \pi/2$ —
two punctures, not four; the equilibria $x \equiv \pi/4$ are
interior points of $\mathcal D$ where every operator is finite (the
saws tie there; nothing is removed). Second, the monodromy that
genuinely exists on the operator side belongs to the root tower,
not to the saw bit: antiperiodic functions on the $2\pi$-circle are
sections of the nontrivial $\mathbb Z/2$ local system over the
$\pi$-circle, the generator of $H^1(S^1; \mathbb Z/2) \cong
\mathbb Z/2$; $g$ of §2 is such a section. The saw-trend bit of §3
is a derivative-sign *field* on the sectors, not a monodromy class,
and the phasor winding is an integer class — three distinct objects
at three distinct levels, not to be conflated.

> **OQ1 (Bridge, sharpened).** The spinor module $V_{32}$ of the
> path construction carries the chirality grading $V_{32} = V_{16}^+
> \oplus V_{16}^-$. In the construction as it stands, the operator
> algebra $\mathcal A$ (generated by $srx, sxp, cxp, crx$ and their
> composites) has **no action** on $V_{32}$: the operators enter the
> path construction as scalar weights evaluated at sample points,
> not as operators acting on the module. The bridge therefore has
> two ordered sub-questions. **OQ1a:** construct an action of
> $\mathcal A$ (or of the function-space algebra generated by $S$
> and $G$ of §2) on $V_{32}$, or prove that the weight-based
> construction admits none. **OQ1b:** given such an action,
> determine whether the saw-trend bit (§3) or the root-tower
> monodromy class above is carried by it to the chirality grading —
> i.e., whether a grading-preserving intertwining selects $V_{16}^+$
> vs $V_{16}^-$ from operator data. Whether the phasor winding
> factors into the picture is a separate sub-question, not part of
> OQ1. A functorial formulation (monodromy data $\to$ graded
> $\mathcal A$-modules) is well-posed only after OQ1a is answered,
> since "$\mathcal A$-linear" is undefined for $V_{32}$ until an
> action exists.

> **OQ2 (Miniature).** In $\mathfrak{so}(4)$ the operator-weighted
> path construction splits $6 = 3 + 3$ under the Hodge star
> (`proofs/so7_so4_so3.md`, PROVED). The star is the volume element
> acting on 2-forms. Does any operator datum — the weights, their
> arrangement, the field's circulation — have grip on the star's
> sign, or is the split provably weight- and circulation-blind?

> **OQ3 (Gluing).** Declare (§4) the trajectory gluing class at the
> four equilibria — pass-through, bounce, or a stated mixture — and
> determine whether the chirality bit of §3 is invariant under the
> choice or is itself the gluing in disguise.

> **OQ4 (Complex structure).** WALL-CPLX, restated: $J$ of §2 is a
> complex structure exactly on $FW = +1$ sectors. Does any operator-
> native construction give a complex structure, or a real structure
> pair, on the full system? (The 1820 four-form analysis in the audit
> workspace shows two incompatible real structures are the linear-
> algebraic content of the analogous obstruction there.)

**On the wall's location (ASSERTED synthesis, 2026-10-16).** The
negative results of §2 and the questions of §5 suggest a single
joint: the argument side of the operator algebra (where the root
tower, the half-turn $S$, and the phase function $g$ live) is not
canonically identified with the value side (where the Clifford
generators, the chirality volume element $\Gamma$, and the spinor
module live). The campaign's negative or open items all sit at this
joint: Wick-via-$J$ (§2; FALSIFIED — the argument-side $J$ does not
act on the value-side matrix algebra); WALL-CPLX (§2; the
argument-side $J$ is a complex structure only on $FW = +1$
sectors); the $16$ vs $16'$ tie (§5; no value-side scalar selects
the chirality split); OQ2 (§5; OPEN, not yet a negative result —
whether any argument-side datum grips the star's sign is untested).
This is an organizing synthesis, not a proved reduction: the four
items are not derived from a single non-identification theorem, and
the general non-identification statement is OPEN (see §2). If the
wall is a single wall, this is where it stands.

## 6. What is not claimed

No identification of any of these objects with physical spinors,
particles, spacetime, or forces is made or implied; such
identifications in the wider manuscript are ASSERTED and lie outside
this exposition. The E8 host hypothesis is excluded by the
repository's clean-break policy. The negative results of §2
(Wick-via-$J$) and §5 are considered part of the result set: the wall,
precisely located, is a contribution.

## References (this repository / audit workspace)

- `proofs/pythagorean_foundation.md` — the operator foundation (§1).
- `proofs/sqrt_fw_clifford.md` §§1–14 — the root tower (§2).
- `proofs/wick_via_j.md` + `proofs/wick_via_j.py` — §2, negative test.
- `proofs/so10_path.md`, `proofs/spin10_spinor.md`,
  `proofs/so7_so4_so3.md` — the orthogonal ladder and spinors (§5).
- `proofs/rapidity_dynamics.md` — §1 rapidity dictionary.
- Audit workspace, IO-word sweep runs 963–964
  (`RUN963_ALLORDERS.md`, `RUN964_CLOSEDFORMS.md`) — the all-orders
  power-sum theorem and its closed forms (§4 corollary).
- Session record 2026-10-16 — §§3–4 (saw bit; tie gluing), this chat's
  log in `~/memory/2026-10-16.md`.
