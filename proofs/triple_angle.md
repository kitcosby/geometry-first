# Multiple-angle identities in the operator system: the 1x dictionary, 3x, 4x, and the ladder pattern

**Kit Cosby — 2026-10-16. Status: PROVED (assembly from proved
ingredients; derivations given). Tripwire numerics in §5.**

Notation. On $\mathcal D = \mathbb R \setminus \{k\pi/2\}$: the core
operators $srx, sxp = 1/srx, cxp, crx = 1/cxp$; the UNA pair
$urx = srx - crx$, $uxp = cxp - sxp$; the sign bits
$\sigma_s = \operatorname{sgn}(\sin x)$, $\sigma_c =
\operatorname{sgn}(\cos x)$; and the three denominators
$$
D = urx + uxp, \qquad
N = (srx + crx) - (sxp + cxp), \qquad
K = srx \cdot cxp - 1.
$$
Recall (PROVED, E77 and the Pythagorean Foundation): $\sin 2x =
4/D$ and $\cos 2x = N/D$; equivalently $D = 4/\sin 2x$ and
$N = 4\cot 2x$.

## 1. Lemma (the 1x dictionary). PROVED.

$$
\sin x = \frac{\sigma_s\, cxp + \sigma_c}{K}, \qquad
\cos x = \frac{\sigma_c\, srx + \sigma_s}{K}.
$$

*Proof.* The branch lemma (RUN144, PROVED) gives, with signed
$\sin x, \cos x$,
$$
srx = \frac{\cos x + \sigma_s}{\sin x}, \qquad
cxp = \frac{\sin x + \sigma_c}{\cos x}
$$
on every open octant. Cleared of denominators this is the linear
system $srx \sin x - \cos x = \sigma_s$,
$-\sin x + cxp \cos x = \sigma_c$, whose determinant is
$srx \cdot cxp - 1 = K$. Solving gives the two formulas. ∎

So $K$ is the *1x denominator*, in the same sense that $D$ is the
2x denominator. Note $K$ already appears in the FlatWave identity
$\frac{srx + cxp}{srx \cdot cxp - 1} = sawup + sawdown$.

*Bit accounting.* $FW = \sigma_s \sigma_c$ is operator-internal.
The individual signs are not: $\sin x$ and $\cos x$ are
$2\pi$-periodic, not $\pi$-periodic, so no formula in the
$\pi$-periodic operators alone can produce them (cf. the
octant-label result of 2026-10-14: operator data determine the
octant modulo $\pi$; one external bit, e.g. $\sigma_s$, completes
it). The dictionary consumes exactly that one external **lap bit**;
given $\sigma_s$, the other sign is $\sigma_c = FW \cdot \sigma_s$.
The 1x objects live one cover up from the operators — the same
double cover on which $g = \sqrt{FW}$ lives
(`proofs/sqrt_fw_clifford.md` §1).

## 2. Theorem (3x). PROVED.

Expanded form:
$$
\sin 3x = \frac{4\,(\sigma_c\, srx + \sigma_s) + N\,(\sigma_s\, cxp + \sigma_c)}{D \cdot K},
\qquad
\cos 3x = \frac{N\,(\sigma_c\, srx + \sigma_s) - 4\,(\sigma_s\, cxp + \sigma_c)}{D \cdot K}.
$$
Dressed form:
$$
\sin 3x = \sin x \cdot \frac{2N + D}{D}, \qquad
\cos 3x = \cos x \cdot \frac{2N - D}{D},
$$
with $\sin x, \cos x$ from §1.

*Proof.* Angle addition $3x = 2x + x$ with $\sin 2x = 4/D$,
$\cos 2x = N/D$ gives $\sin 3x = (4\cos x + N \sin x)/D$ and
$\cos 3x = (N \cos x - 4 \sin x)/D$; substituting §1 gives the
expanded form. For the dressed form, use
$\cos x \sin 2x = \sin x (1 + \cos 2x)$, whence $\sin 3x = \sin x
(2\cos 2x + 1)$, and $\cos 3x = \cos x (2 \cos 2x - 1)$ similarly;
substituting $\cos 2x = N/D$ gives the claim. ∎

*Remark.* The dressed form is the classical triple-angle identity
in operator costume: with $\cos 2x = 1 - 2\sin^2 x$,
$\sin x (2\cos 2x + 1) = 3 \sin x - 4 \sin^3 x$, and
$\cos x (2\cos 2x - 1) = 4\cos^3 x - 3\cos x$. Note the dressed
factor is $\frac{2N + D}{D} = 2\cos 2x + 1$ — **not**
$\frac{2N + D}{2D}$ (which equals $\cos 2x + \tfrac12$ and would
halve the result; this factor-of-2 slip appeared in an audit draft
of 2026-10-16 and is corrected here).

## 3. Theorem (4x) and the ladder pattern. PROVED.

$$
\sin 4x = \frac{8N}{D^2}, \qquad \cos 4x = \frac{N^2 - 16}{D^2}.
$$

*Proof.* $\sin 4x = 2 \sin 2x \cos 2x = 2 \cdot \frac4D \cdot
\frac ND$; $\cos 4x = \cos^2 2x - \sin^2 2x = \frac{N^2 - 16}{D^2}$.
∎

**Ladder pattern.** For every $n \ge 1$, $\sin nx$ and $\cos nx$
have operator forms with denominator $D^{\lfloor n/2 \rfloor}$ if
$n$ is even, and $D^{\lfloor n/2 \rfloor} \cdot K$ if $n$ is odd.

*Proof sketch.* Even case $n = 2m$: $\cos 2mx = T_m(\cos 2x)$ and
$\sin 2mx = \sin 2x \cdot U_{m-1}(\cos 2x)$ (Chebyshev), so
substituting $\sin 2x = 4/D$, $\cos 2x = N/D$ gives denominator
$D^m$. Odd case $n = 2m + 1$: angle addition against the even case
introduces exactly one factor of $\sin x$ or $\cos x$, hence exactly
one factor $K$ by §1. ∎

In particular $K$ never accumulates: it enters once, through the 1x
recovery, exactly when the multiple is odd. (The guess that 4x
would carry $D \cdot K \cdot L$ for some further operator $L$ was
raised in the same 2026-10-16 audit draft and is FALSIFIED by §3:
$4x$ stays in the doubling world, over $D^2$.)

## 4. Orientation accounting (corrected corollary)

The multiple-angle formulas consume: the operator data themselves;
the internal bit $FW = \sigma_s \sigma_c$; and, for odd multiples
only, the single external lap bit of §1. They do **not** consume the
saw-trend orientation bit of the Clifford exposition §3 — that bit
is a different object (a derivative-sign field, operator-internal
and $\pi$-periodic). Even multiples need no external bit at all,
consistent with $\sin 2mx$ and $\cos 2mx$ being $\pi$-periodic —
functions on the quotient circle, as they must be.

## 5. Verification (tripwire)

Exact hand checks (from the 2026-10-16 audit, independently
recomputed): at $x = \pi/6$: $srx = 2 + \sqrt3$, $cxp = \sqrt3$,
$D = 8/\sqrt3$, $N = 4/\sqrt3$, $K = 2(1 + \sqrt3)$; both forms give
$\sin 3x = 1$, $\cos 3x = 0$. At $x = 2\pi/3$ ($\sigma_s = +1$,
$\sigma_c = -1$): $srx = 1/\sqrt3$, $cxp = 2 - \sqrt3$,
$D = -8\sqrt3/3$, $N = 4/\sqrt3$, $K = 2(\sqrt3 - 3)/3$; both forms
give $\sin 3x = 0$, $\cos 3x = 1$.

Numerical tripwire (mpmath, 10 points across all four quadrants,
degrees 15–345): the §1 dictionary, the §2 pair in both forms, and
the §3 pair all agree with the direct values to $\le 2 \times
10^{-15}$. Per the corpus rule, the derivations are the
establishment; the numerics are a tripwire only.

## Cross-references

- `proofs/pythagorean_foundation.md` — E77 double-angle pair; the
  branch lemma's home is RUN144 (audit workspace).
- `proofs/sqrt_fw_clifford.md` §1 — the double cover and $g$.
- `proofs/clifford_problem_exposition.md` §3 — the saw-trend bit
  (a distinct orientation datum; see §4 above).
