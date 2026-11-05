# The CTN–FW continued fraction: a conjugate-pair theorem, corrected

**Kit Cosby — 2026-10-17. Status: DRAFT, not yet filed on Kit's word.
Algebraic identities PROVED (Vieta, verified symbolically in the main
chat 2026-10-17); convergence statement CORRECTED as in §3 — the
draft Kit brought in claimed convergence to the positive root on
both sectors, and that claim is FALSE for FW = −1. No physical
identifications.**

Notation. $CTN$ is the curtain function,
$CTN = 1/\sqrt{1+|\sin 2x|}$ (PROVED, sweep runs 930–935), so on the
operator domain ($\sin 2x \ne 0$) $CTN \in [1/\sqrt2, 1)$, with
$CTN \to 1$ only as a limit at the excluded axes. $FW =
\operatorname{sgn}(\sin 2x) \in \{+1, -1\}$. The theorem itself is
stated for a general $CTN \in (0,1]$.

## 1. The positive roots. PROVED.

For $FW \in \{+1,-1\}$ let $X_{FW}$ be the positive root of
$$
X^2 + FW\cdot X - CTN = 0,
\qquad
X_\pm := X_{\pm 1} = \frac{\mp 1 + \sqrt{1+4\,CTN}}{2}.
$$
Then, by Vieta / direct substitution:
$$
\boxed{\;X_- - X_+ = 1, \qquad X_+ \cdot X_- = CTN.\;}
$$
Further $X_+ + X_- = \sqrt{1+4\,CTN}$, $X_+^2 + X_+ = CTN$, and
$X_-^2 - X_- = CTN$. All five identities were verified symbolically
in CTN in the main chat, 2026-10-17.

Reading: the two sector values differ by exactly one unit — the
discrete (FW) bit, independent of $x$ — and their product is exactly
the curtain — the continuous bit. On positive roots, the sector
flip is a **unit shift**.

Ranges over $CTN \in (0,1]$: $X_+ \in (0, 1/\varphi]$,
$X_- \in (1, \varphi]$ (corrected endpoints: $X_-$ does not attain 1).

## 2. Golden limits at the boundary. PROVED (as limits).

At $CTN \to 1$ (the excluded-axis limit $|\sin 2x| \to 0^+$):
$$
X_+ \to \frac{\sqrt5-1}{2} = \frac1\varphi \approx 0.6180,
\qquad
X_- \to \frac{\sqrt5+1}{2} = \varphi \approx 1.6180.
$$
At $CTN = 1$ the two sector quadratics *are* the golden equations
$X^2 + X - 1 = 0$ and $X^2 - X - 1 = 0$, and $\varphi - 1/\varphi = 1$
is the unit-shift identity of §1 at that point. $\varphi$ is a
**limit at the excluded axis, never an attained value** on the
operator domain. This is exact algebra, of a different kind from
the logged near-miss det $M_2$ at $30^\circ = 1.6188$, which remains
NOT a match (H-EIG registry).

Interior point ($|\sin 2x| = 1$, $CTN = 1/\sqrt2$), corrected
decimals (40-digit check, 2026-10-17):
$$
X_+ = \frac{\sqrt{1+2\sqrt2}-1}{2} \approx 0.4783,
\qquad
X_- \approx 1.4783.
$$
(The decimals 0.277 / 1.277 in the incoming draft came from
evaluating $\sqrt{1+\sqrt2}$ in place of $\sqrt{1+2\sqrt2}$; the
symbolic forms were already right.)

## 3. Convergence — the correction. PROVED.

Consider the periodic continued fraction, i.e. iteration of
$f(t) = CTN/(FW + t)$:
$$
X \;=\; \frac{CTN}{FW + \dfrac{CTN}{FW + \dfrac{CTN}{FW + \cdots}}}
$$

- **FW = +1:** $f$ attracts at $X_+$: $|f'(X_+)| = X_+/(1+X_+) < 1$.
  The CF converges to $X_+$, as drafted.
- **FW = −1:** $X_-$ is a **repelling** fixed point. Indeed
  $FW + X_- = X_- - 1 = X_+$ by §1, so
  $|f'(X_-)| = X_-/X_+ > 1$ for every $CTN > 0$. The attracting
  fixed point is the other root of $X^2 - X - CTN = 0$, namely
  $-X_+$, with $|f'| = X_+/X_- < 1$. Numerically (40-digit, several
  seeds, $CTN = 1$ and $1/\sqrt2$, main chat 2026-10-17) the
  FW = −1 iterates go to $-X_+$ every time; at $CTN = 1$ the
  convergents are $-1, -1/2, -2/3, -3/5, \dots \to -1/\varphi$.
  (One seed hits the pole exactly and is excluded.)

Corrected reading: on **CF limits** the sector flip is a **sign**
(the limits are $X_+$ and $-X_+$); on **positive roots** it is a
**unit shift** ($X_- = X_+ + 1$); and the product identity of §1 is
exactly the converter between the two readings, since the roots of
$X^2 - X - CTN$ are $X_-$ and $-X_+$, with product $-CTN$, so
$-X_+ = -CTN/X_-$.

## 4. Structural reading.

The CF packages the skeleton and the branch in one object: the
skeleton ($CTN$, a function of $a = |\sin 2x|$ alone) is recovered
as the product of the two positive roots; the branch ($FW$) is
recovered as their difference (unit shift) or, at the level of the
CF's own limits, as a sign. This is the skeleton/dressing relation
of the ledger stated in a single map — with the §3 correction that
the packaging lives in the pair (positive root, CF limit), not in
the two positive roots alone.

## 5. Remark: the doubly-infinite Fibonacci skeleton at CTN = 1.

Kit's sequence (2026-10-17): $\dots 13, -8, 5, -3, 2, -1, 1, 0, 1,
1, 2, 3, 5, 8, 13 \dots$ is $F_n$ for all integers $n$, with
$F_{-n} = (-1)^{n+1} F_n$. Checked in the main chat 2026-10-17:
the recurrence $F_n = F_{n-1} + F_{n-2}$ holds across 0 with no seam
($n = -25 \dots 27$) and Cassini's identity
$F_{n-1}F_{n+1} - F_n^2 = (-1)^n$ survives for negative $n$.
Right-tail ratios tend to $\varphi$; left-tail ratios tend to
$-1/\varphi$ — the conjugate root $\psi$, and exactly the FW = −1
CF limit of §3 at $CTN = 1$. At the boundary the discrete skeleton
of the sector pair is thus a single doubly-infinite integer
sequence; Binet's formula $F_n = (\varphi^n - \psi^n)/\sqrt5$
uses both sector roots at once.

## 6. Ledger row (corrected, ready to file on Kit's approval).

"CTN–FW continued fraction. The periodic CF $X = CTN/(FW + X)$
converges for FW = +1 to $X_+ = (\sqrt{1+4CTN}-1)/2$; for FW = −1
it converges to $-X_+$, the negative root of $X^2 - X - CTN$
($X_-$ is repelling, $|f'| = X_-/X_+ > 1$). The positive roots obey
$X_- - X_+ = 1$ and $X_+ X_- = CTN$: the sector flip is a unit
shift on positive roots, a sign on CF limits, and the product
identity converts between them. Golden limits $1/\varphi$ and
$\varphi$ (roots), $-1/\varphi$ (FW = −1 CF) at $CTN \to 1$."

**Scope:** all sectors where $\sin 2x \ne 0$; golden values are
excluded-axis limits, stated as limits, not values.
