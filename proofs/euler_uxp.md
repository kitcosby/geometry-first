# Euler mapped to the UNA pair: the theta parametrization, by quadrant

**Kit Cosby — 2026-10-16. Status: PROVED (derivation from the
Pythagorean Foundation's branch forms; tripwire numerics in §6).
§5 is a proved remark. No physical identifications.**

Notation. On $\mathcal D = \mathbb R \setminus \{k\pi/2\}$: the core
operators $srx, sxp = 1/srx, cxp, crx = 1/cxp$; the UNA pair
$urx = srx - crx$, $uxp = cxp - sxp$; their reciprocals, the saws,
$sawup = 1/urx$, $sawdown = 1/uxp$, satisfying (PROVED, Pythagorean
Foundation) $sawup + sawdown = FW$ and $|sawup| + |sawdown| = 1$,
where $FW = \operatorname{sgn}(\sin 2x)$. The Foundation's branch
forms give, with $sc = \sin x \cos x$ (signed) and the branch
parameter $p = \sigma_s \cos x - \sigma_c \sin x$,
$$
urx = \frac{1 + p}{sc}, \qquad uxp = \frac{1 - p}{sc},
\qquad |p| < 1 \ \text{on all of } \mathcal D \ \text{(PROVED).}
$$

## 1. The direct form on QI. PROVED.

On QI ($0 < x < \pi/2$), $uxp = \sec x + \tan x - \csc x + \cot x$,
and clearing denominators,
$$
uxp \;=\; \frac{1 + \sin x - \cos x}{\sin x \cos x}
\;=\; \frac{2}{1 + \cos x - \sin x},
$$
the last step by $(1 + s - c)(1 + c - s) = 1 - (s - c)^2 = 2sc$.
Since $1 + \cos x - \sin x = 1 + \sqrt2 \cos(x + \pi/4)
= 1 + \sqrt2\,\mathrm{Re}\, e^{i(x + \pi/4)}$,
$$
uxp(x) \;=\; \frac{2}{\,1 + \sqrt2\,\mathrm{Re}\, e^{i(x+\pi/4)}\,}
\;=\; \frac{2}{1 + \mathrm{Re}\big[(1+i)e^{ix}\big]}.
$$
This is a true restatement; its Euler content is the rewriting of
$\cos x - \sin x$ as $\sqrt2\cos(x + \pi/4)$. The structural content
is §2.

## 2. The structural form on QI: the FlatWave identity is Pythagoras. PROVED.

Define $\theta(x) \in (0, \pi/2)$ on QI by
$$
\cos 2\theta \;=\; \cos x - \sin x .
$$
Then $\cos^2\theta = \frac{1 + \cos x - \sin x}{2}$, so by §1
$uxp = \sec^2\theta$; and $1/urx = \frac{\sin x \cos x}{1 + \cos x
- \sin x} = \sin^2\theta$ by the same $2sc = 1 - (c - s)^2$
identity, so $urx = \csc^2\theta$. Hence on QI
$$
\frac{1}{uxp} + \frac{1}{urx} \;=\; \cos^2\theta + \sin^2\theta \;=\; 1:
$$
the FlatWave identity is the Pythagorean identity in operator
clothing, and Euler's formula reads
$$
e^{i\theta} \;=\; \frac{1}{\sqrt{uxp}} + i\,\frac{1}{\sqrt{urx}}
\;=\; \sqrt{sawdown} + i\,\sqrt{sawup},
$$
with $\sin 2\theta = \sqrt{\sin 2x}$. The components are square
roots of the saws. $\theta$ is a bijection of QI onto $(0, \pi/2)$,
increasing, fixing $\pi/4$, with the mirror symmetry
$\theta(\pi/2 - x) = \pi/2 - \theta(x)$ inherited from
$uxp(\pi/2 - x) = urx(x)$.

Values (the first and last rows are **limits**; $0, \pi/2 \notin
\mathcal D$):

| $x$ | $uxp$ | $urx$ | $1/uxp = \cos^2\theta$ | $1/urx = \sin^2\theta$ | $\theta$ |
|---|---|---|---|---|---|
| $0$ (limit) | $1$ | $\infty$ | $1$ | $0$ | $0$ |
| $\pi/6$ | $2(\sqrt3 - 1) \approx 1.464$ | $2 + \frac{2\sqrt3}{3} \approx 3.155$ | $\approx 0.683$ | $\approx 0.317$ | $\approx 0.599$ |
| $\pi/4$ | $2$ | $2$ | $1/2$ | $1/2$ | $\pi/4$ |
| $\pi/3$ | $2 + \frac{2\sqrt3}{3}$ | $2(\sqrt3 - 1)$ | $\approx 0.317$ | $\approx 0.683$ | $\approx 0.972$ |
| $\pi/2$ (limit) | $\infty$ | $1$ | $0$ | $1$ | $\pi/2$ |

## 3. Theorem (proof by quadrant). PROVED.

The parametrization is global. Since $p^2 = 1 - 2\sigma_s\sigma_c
\, sc$ and $\sigma_s\sigma_c = FW$, we have $1 - p^2 = 2 \cdot FW
\cdot sc$ in every quadrant, and therefore
$$
sawup = \frac{FW\,(1 - p)}{2}, \qquad
sawdown = \frac{FW\,(1 + p)}{2}.
$$
Defining $\theta = \frac12 \arccos p \in (0, \pi/2)$ — legitimate
everywhere, since $|p| < 1$ globally — this reads
$sawup = FW \sin^2\theta$, $sawdown = FW \cos^2\theta$, and
inverting,
$$
\boxed{\;uxp = FW \cdot \sec^2\theta, \qquad
urx = FW \cdot \csc^2\theta\;}
$$
in all four quadrants, with $p$ spelled per quadrant:

| quadrant | $p$ | $FW$ |
|---|---|---|
| QI | $\cos x - \sin x$ | $+1$ |
| QII | $\cos x + \sin x$ | $-1$ |
| QIII | $\sin x - \cos x$ | $+1$ |
| QIV | $-(\cos x + \sin x)$ | $-1$ |

Per sector: QI is §2. QIII is *identical* to QI — $p$ and all
operator values are $\pi$-periodic. QII and QIV: $uxp =
-\sec^2\theta$, $urx = -\csc^2\theta$; the saws are negative,
$sawdown = -\cos^2\theta$, $sawup = -\sin^2\theta$, and Euler reads
$e^{i\theta} = \sqrt{|sawdown|} + i\sqrt{|sawup|}$. In every
quadrant $\theta$ is determined by the saws alone:
$p = \frac{sawdown - sawup}{sawdown + sawup}$.

## 4. Corrections forced by the quadrant proof.

Two earlier statements of 2026-10-16 do not survive §3 and are
retracted here so the file carries the corrected record:

1. **No QIII swap.** An audit remark suggested the construction
   swaps the $uxp/urx$ roles in QIII. It does not: the remark came
   from extending the QI-only variable $\cos x - \sin x$ into QIII,
   where $p$ is actually $\sin x - \cos x$. With $p$ as the
   variable there is no swap anywhere.
2. **No hyperbolic sectors.** The claim that $FW = -1$ sectors
   become hyperbolic (via $J^2 = +1$) does not occur for this
   parametrization: $\theta$ stays real and trigonometric in QII
   and QIV because $|p| < 1$ is global. The sector signature
   appears only as the overall $FW$ factor in §3. The hyperbolic
   appearance is an artifact of insisting on $\cos x - \sin x$
   outside QI, where it leaves $[-1, 1]$ and $\arccos$ fails.

## 5. Remark (polar traversal). PROVED.

$p$ is $\pi/2$-periodic: under $x \mapsto x + \pi/2$,
$(\sigma_s, \sigma_c) \mapsto (\sigma_c, -\sigma_s)$ and
$(\sin x, \cos x) \mapsto (\cos x, -\sin x)$, so $p \mapsto
-\sigma_c \sin x + \sigma_s \cos x = p$; hence so is $\theta$.
And $\theta$ sweeps $(0, \pi/2)$ monotonically across each
quadrant: $\theta \to 0$ at each quadrant's opening seam, $\theta
= \pi/4$ at the tie, $\theta \to \pi/2$ at the closing seam. So
the polar point $\Phi(x) = e^{i\theta(x)}$ sweeps the quarter
circle once per quadrant, identically — four sweeps per lap — and
the only data distinguishing one sweep from another live outside
the $\theta$ map: $FW$ itself (on record as
$FW = \operatorname{sgn}(\operatorname{Im} \zeta^2)$ for
$\zeta = e^{ix}$) and the lap bit. $FW$ is the circle's sign
readout, not the circle itself.

## 6. Tripwire numerics.

Thread verification (2026-10-16): the §3 formulas checked at
twelve points — all four quadrants plus the four ties, where
$p = 0$ and $\theta = \pi/4$ — to $10^{-15}$.

Independent recheck (2026-10-16, 40-digit mpmath, operators built
directly from $|\sec x| \pm \tan x$, $|\csc x| \pm \cot x$): at
the twelve points $15°, 30°, 45°$ offset into each quadrant,
$|uxp - FW\sec^2\theta|$, $|urx - FW\csc^2\theta|$,
$|1/urx - sawup|$ and $|1/uxp - sawdown|$ all $\le 2.8 \times
10^{-40}$.
