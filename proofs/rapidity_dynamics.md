# Rapidity Dynamics: the Operator Families in the η Parameter

**Date:** 2026-10-14 · **Status: velocity PROVED symbolically in QI;
transport classified (2×2); QII dictionary CHECKED exactly at points;
speed-of-light screen: NO MATCH (N = 11).**
**Scripts:** `proofs/rapidity_dynamics.py` (28/29 checks passed — the one
FAIL is a malformed test in the script itself, disclosed in §6),
`proofs/rapidity_dynamics_fixup.py` (the corrected checks, all passed).
**Engine:** sympy 1.12, exact symbolic (identities verified by rewriting
to rational functions of $e^\eta$).

## 1. Question

Following the FW-ladder run and the sampling-vs-dynamics question: rewrite
the operator families in the rapidity parameter $\eta = \ln cxp$ (the
earlier rapidity work: SQ1/SQ4), compute the velocity $dM/d\eta$ in closed
operator form, determine whether the transport
$U = \mathcal P \exp \int M\,d\eta$ closes in form, and screen the results
for a rational factor of the speed of light.

## 2. Setup (QI; all identities on record or verified here)

On QI, $\eta > 0$ and the dictionary is
$\sinh\eta = \tan x$, $\cosh\eta = \sec x$, $\tanh\eta = \sin x$, with

$$cxp = e^{\eta},\quad crx = e^{-\eta},\quad
srx = \coth\tfrac{\eta}{2} = \frac{1+\operatorname{sech}\eta}{\tanh\eta},
\quad sxp = \tanh\tfrac{\eta}{2},\quad
urx = srx - crx,\quad uxp = cxp - sxp,$$
$$\sin 2x = 2\tanh\eta\,\operatorname{sech}\eta,\qquad
\cos 2x = \operatorname{sech}^2\!\eta - \tanh^2\!\eta .$$

Each line verified symbolically in the script (Part A, 6/6).

## 3. The velocity closes in form — PROVED (QI, symbolic)

Every entry of $dM/d\eta$ is again an operator expression:

$$\frac{d\,srx}{d\eta} = -\frac{srx^2-1}{2},\qquad
\frac{d\,sxp}{d\eta} = \frac{1-sxp^2}{2},\qquad
\frac{d\,cxp}{d\eta} = cxp,\qquad \frac{d\,crx}{d\eta} = -crx,$$
$$\frac{d\,urx}{d\eta} = -\frac{srx^2-1}{2} + crx,\qquad
\frac{d\,uxp}{d\eta} = cxp - \frac{1-sxp^2}{2},$$

and the double-angle block rotates against itself with an operator rate:

$$\frac{d\sin 2x}{d\eta} = \frac{4\cos 2x}{cxp+crx} = 2\,\operatorname{sech}\eta\,\cos 2x,
\qquad
\frac{d\cos 2x}{d\eta} = -\frac{4\sin 2x}{cxp+crx}.$$

All 8 identities verified symbolically (Part B), plus the chain-rule
consistency check $d\sin 2x/d\eta = (d\sin 2x/dx)/(d\eta/dx)$.
Consequently the full E96 velocity $dM_{96}/d\eta$ (traceless-corrected
matrix, $FW = +1$ on QI) is an explicit operator matrix — all 9 entries
verified, and $\operatorname{tr}(dM_{96}/d\eta) = 0$: the velocity stays
inside $\mathfrak{sl}(3,\mathbb R)$, as it must. **The curve's velocity in
rapidity is a closed operator form, not a series.**

## 4. The transport does not close elementarily — classified (2×2)

For the $2\times2$ family $M = \begin{pmatrix} 0 & a \\ b & 0\end{pmatrix}$,
$a = \coth\frac\eta2$, $b = e^\eta$, the transport reduces channel-wise to

$$u'' + P u' + Q u = 0,\qquad
P = \frac{1}{\sinh\eta},\qquad Q = -e^{\eta}\coth\tfrac{\eta}{2}$$
($P = -a'/a$ verified symbolically). In the variable $w = e^\eta$:

$$u_{ww} + \frac{w^2 + 2w - 1}{w(w-1)(w+1)}\,u_w
- \frac{w+1}{w(w-1)}\,u = 0,$$

coefficients verified against the $\eta$-form at $w = e^\eta$. Singular
points, by exact denominator factorization and Fuchs limits:

- $w = 0$, $w = 1$, $w = -1$: **regular singular** (both Fuchs limits
  finite at each).
- $w = \infty$: **irregular singular** — $Q(1/z)/z^2 \to -\infty$; the
  $\eta \to +\infty$ end is essential, with WKB growth $\sim e^{\eta}$
  in $Q$.

On the real $\eta$-line this says: the transport is regular at
$\eta = 0$ (the $x = 0$ end, $\beta = \tanh\eta = 0$) and has its
essential singularity at $\eta = +\infty$ (the $x = \pi/2$ boundary,
$\beta = 1$). The 3×3 (E96) ordered transport was **not** attempted in
closed form — OPEN; §3 is the closed-form content at that size.

*Observation (exact, significance OPEN):* the friction coefficient's
numerator $w^2 + 2w - 1$ vanishes at $w = \sqrt2 - 1$, i.e. at
$\eta = -\operatorname{arsinh} 1$, where the $u_w$ term drops out of the
ODE. $\sqrt2 - 1$ is the value $cxp$ takes at the $135^\circ$ tie point
(record: $srx = cxp = \sqrt2 - 1$ at $135^\circ$). This is a value-level
coincidence between the ODE's analytic form and the tie value; no
dynamical significance is claimed.

## 5. QII: the signed dictionary (CHECKED exactly)

The SQ1 parenthetical $\tanh\eta = \sin x$ holds only where $\cos x > 0$.
The global signed form is $\tanh\eta = \operatorname{sgn}(\cos x)\sin x$
($\sinh\eta = \tan x$ and $\cosh\eta = |\sec x|$ are global). Verified by
exact rational arithmetic at $x = 120^\circ$ ($w_0 = e^{\eta_0} = 2-\sqrt3$):
$\sinh\eta_0 = -\sqrt3$, $\cosh\eta_0 = 2$, $\tanh\eta_0 = -\sqrt3/2$,
and the QII branch form $srx = -\frac{1}{\tanh\eta} + \frac{1}{\sinh\eta}$
reproduces $srx(120^\circ) = \sqrt3/3$. The $2\times2$ velocity was
checked against the $x$-side at $120^\circ$ (agreement to $10^{-165}$ in
exact-algebraic evaluation), and at $30^\circ$, $60^\circ$ on the QI side.

## 6. Disclosed script defect

`rapidity_dynamics.py` reports 28/29: the failing check
(`sinh eta = -1 at w = sqrt2-1`) tested $(w - 1/w)/2 + 1 \equiv 0$ as an
identity in the symbol $w$ instead of substituting $w = \sqrt2 - 1$.
The corrected substitution is in `rapidity_dynamics_fixup.py` and gives
exactly $-1$. No other result depends on that check.

## 7. Speed-of-light screen (H-EIG discipline)

Pre-registered form: $c = 299792458$ m/s is dimensionful, so a pure-algebra
result can be a "rational factor of $c$" only as a rational multiple of
that bare number, with zero free parameters (no unit choice may be
introduced). Screened quantities produced by this computation (N = 11):
rational factors $\tfrac12$ (half-angle/Riccati velocities), $2$
(double-angle rate numerator over $cxp + crx$), $4$ (same rate,
unnormalized), $\tfrac13$ (E96 trace correction); singular points
$w \in \{0, 1, -1\}$; friction zero $w = \sqrt2 - 1$; tie-point data at
$\eta^* = \ln(1+\sqrt2)$: eigenvalues of $M_2$ equal $\pm(1+\sqrt2)$,
eigenvalues of $dM_2/d\eta$ equal $\pm i(1+\sqrt2)$; and $\eta^*$ itself.

**Result: NO MATCH.** No screened quantity is a rational multiple of
$299792458$, and the only rationals present are the structural factors
$\{\tfrac12, \tfrac13, 2, 4\}$. The remaining quantities are irrational
($1+\sqrt2$, $\sqrt2-1$) or logarithmic ($\ln(1+\sqrt2)$). Near-misses:
none reported as leads, per discipline.

The one place $c$-type structure genuinely lives in this formalism is
not a factor at all: $\beta = \tanh\eta \to 1$ at the irregular boundary
$\eta = \infty$ (§4). In any relativistic reading, $c$ is the boundary
value of the dictionary, not a coefficient produced by it. BR-TIME and
all physical identifications remain OPEN/undeclared, untouched by this
computation.

## 8. Scope

- PROVED (symbolic, QI): the η-dictionary as used, all velocity formulas,
  the E96 velocity and its tracelessness, the 2×2 ODE and its
  $w$-form, the singular-point classification.
- CHECKED (exact point evaluations): QI points 30°, 60°; QII signed
  dictionary and 2×2 velocity at 120°.
- OPEN: E96 (and ladder) ordered transport in any form; the ℓ-mirror of
  §3–§4 (expected by the complementary-angle symmetry, not computed);
  any significance of the §4 friction zero; BR-TIME.
