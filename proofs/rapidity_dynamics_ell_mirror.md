# Rapidity Dynamics — the ℓ-Mirror

**Date:** 2026-10-14 · **Status: PROVED symbolically in QI (29/29 checks);
point checks at 30°/60° exact to ~10⁻¹⁶⁵. Companion to
`proofs/rapidity_dynamics.md` §8, which listed this mirror as OPEN
("expected by the complementary-angle symmetry, not computed").**
**Script:** `proofs/rapidity_dynamics_ell_mirror.py` (sympy 1.12, exact;
identities verified by rewriting to rational functions of $e^\ell$).
**Not pushed — local file, Kit's call on sync like the other geometry-first work.**

## 1. Question

`rapidity_dynamics.md` computed the velocity $dM/d\eta$ and the 2×2
transport classification for $\eta = \ln cxp$. Compute the same objects
for the second rapidity variable $\ell = \ln srx$ (SQ4, PROVED:
$d\ell/dx = -|\csc x| = -(srx+sxp)/2$).

## 2. The ℓ-dictionary (QI) — PROVED

The mirror is exact, with $\sin x$ and $\cos x$ swapped against the
$\eta$-dictionary:

$$srx = e^{\ell},\quad sxp = e^{-\ell},\quad
cxp = \coth\tfrac{\ell}{2} = \frac{1+\operatorname{sech}\ell}{\tanh\ell},
\quad crx = \tanh\tfrac{\ell}{2},$$
$$\sinh\ell = \cot x,\quad \cosh\ell = |\csc x| = \tfrac{srx+sxp}{2},
\quad \sin x = \operatorname{sech}\ell,\quad \cos x = \tanh\ell,$$
$$\sin 2x = 2\tanh\ell\,\operatorname{sech}\ell,\qquad
\cos 2x = \tanh^2\!\ell - \operatorname{sech}^2\!\ell .$$

(The $\eta$ forms were $\sin x = \tanh\eta$, $\cos x = \operatorname{sech}\eta$,
$\cos 2x = \operatorname{sech}^2\!\eta - \tanh^2\!\eta$.) This is the
complementary-angle symmetry made literal: $\ell(x) = \eta(\pi/2 - x)$,
and under $x \to \pi/2 - x$ the pairs $(srx, sxp)$ and $(cxp, crx)$ swap.

## 3. The velocity closes in form — PROVED (QI, symbolic)

$$\frac{d\,srx}{d\ell} = srx,\qquad \frac{d\,sxp}{d\ell} = -sxp,\qquad
\frac{d\,cxp}{d\ell} = -\frac{cxp^2-1}{2},\qquad
\frac{d\,crx}{d\ell} = \frac{1-crx^2}{2},$$
$$\frac{d\,urx}{d\ell} = srx - \frac{1-crx^2}{2},\qquad
\frac{d\,uxp}{d\ell} = -\frac{cxp^2-1}{2} + sxp,$$

i.e. the Riccati pair and the exponential pair exchange roles against §3
of the $\eta$ note. The double-angle block rotates against itself with
the **opposite sign** (because $d\ell/dx < 0$ — $\ell$ runs backwards
along the quadrant):

$$\frac{d\sin 2x}{d\ell} = -\frac{4\cos 2x}{srx+sxp}
= -2\,\operatorname{sech}\ell\,\cos 2x,\qquad
\frac{d\cos 2x}{d\ell} = +\frac{4\sin 2x}{srx+sxp}.$$

Consequently the full E96 velocity $dM_{96}/d\ell$ is an explicit
operator matrix — all 9 entries verified, and
$\operatorname{tr}(dM_{96}/d\ell) = 0$: the velocity in $\ell$ also
stays inside $\mathfrak{sl}(3,\mathbb R)$.

## 4. The 2×2 ℓ-transport — classified, and simpler than the η form

For $M = \begin{pmatrix} 0 & a \\ b & 0\end{pmatrix}$ with
$a = srx = e^\ell$, $b = cxp = \coth\frac\ell2$, the same channel
reduction as the $\eta$ note ($u'' + Pu' + Qu = 0$, $P = -a'/a$,
$Q = -ab$) gives

$$P \equiv -1,\qquad Q = -e^\ell \coth\tfrac{\ell}{2}.$$

In the variable $w = e^\ell$ the friction coefficient is
$(1+P)/w \equiv 0$ — **it vanishes identically**, unlike the $\eta$ form
whose friction numerator $w^2 + 2w - 1$ vanished only at the single
point $w = \sqrt2 - 1$. The ℓ-transport channel is exactly

$$u_{ww} - \frac{w+1}{w(w-1)}\,u = 0.$$

Singular points (Fuchs limits, exact): $w = 0$ and $w = 1$ are
**regular singular** (both limits finite, in fact 0); $w = \infty$ is
**irregular singular** ($Q_w(1/z)/z^2 \to -\infty$). Endpoint reading
*in the ℓ variable*: $w = 1$ is $\ell = 0$, the $x = \pi/2$ ($\beta = 1$)
end — regular here; the irregular end is $w = \infty$, i.e.
$\ell = +\infty$, the $x = 0$ ($\beta = 0$) end. This does not contradict
the $\eta$ classification: $\mathcal P\exp\int M\,d\ell$ and
$\mathcal P\exp\int M\,d\eta$ are different transports (different
parameters, different generators after reparametrization), and the
$\ell$ form is, via §2, the $\eta$-transport of the complementary-angle
(swapped-entry) family.

## 5. Scope

- PROVED (symbolic, QI): the ℓ-dictionary, all velocity formulas, the
  E96 ℓ-velocity and its tracelessness, the 2×2 ℓ-ODE, its $w$-form
  and singular-point classification.
- CHECKED (exact point evaluations): the 2×2 ℓ-velocity against the
  $x$-side at 30° and 60° (agreement ~10⁻¹⁶⁵).
- NOT computed: QII signed ℓ-dictionary (the $\eta$ note's §5 mirror);
  E96/ladder ordered transport in either parameter (still OPEN);
  any significance of the friction-free $w$-form.
- BR-TIME and all physical identifications remain OPEN/undeclared,
  untouched by this computation. No H-EIG screen: this computation
  produced no new spectrum (the structural rationals are the same
  $\{\tfrac12, \tfrac13, 2, 4\}$ class already screened).
