# Cross-Family Sum-Difference (srx+crx) − (sxp+cxp)

Date: 2026-10-12. Status: PROVED (elementary algebra from the definitions
in pythagorean_foundation.md). Question (Kit, Algebra-first chat,
2026-10-12 17:30 UTC): "And (srx+crx)-(sxp+cxp) cross families?"

## Result

On the whole admissible domain D = {x : sin x != 0, cos x != 0}:

    (srx + crx) - (sxp + cxp) = 4 cot(2x) = cos(2x) * (urx + uxp)

This is exactly the numerator of Kit's already-PROVED E77 identity
cos(2x) = (srx + crx - sxp - cxp) / (urx + uxp) — so it is not a new
independent quantity; it is the cos(2x) channel, in the same way
urx + uxp = 4/sin(2x) is the sin(2x) channel.

## Proof

From definitions D.2/D.3 and Cor. 2.3 of the Pythagorean Foundation,
with M = |sec x|, N = |csc x|:

    srx = N + cot x,   sxp = N - cot x
    cxp = M + tan x,   crx = M - tan x

Regroup the expression by family:

    (srx + crx) - (sxp + cxp)
      = (srx - sxp) + (crx - cxp)
      = 2 cot x - 2 tan x            (M and N cancel exactly)

    2 cot x - 2 tan x
      = 2 (cos^2 x - sin^2 x)/(sin x cos x)
      = 2 cos(2x) / (sin(2x)/2)
      = 4 cot(2x).

And since urx + uxp = 4/sin(2x) (Thm 5.1), dividing gives cos(2x),
recovering E77. No branch/sign case split is needed: the absolute
values M, N cancel before any sign enters. ∎

## Consequences

1. Zero exactly at the octant boundaries x = pi/4 + k*pi/2 — the
   dominance-switch points: there srx + crx = sxp + cxp.
2. Sign: sign of 4 cot(2x) = sign(cos 2x / sin 2x); poles are simple,
   at x = k*pi/2 (the domain boundaries), inherited from cot/tan.
3. Cross-family reading: within-family differences (srx-sxp = 2 cot x,
   crx-cxp = -2 tan x) ADD here, whereas the UNA operators are the
   cross-family differences (urx = srx-crx, uxp = cxp-sxp).
   The sum-difference and the UNA sum are the two channels:
   E = 4 cot(2x), urx+uxp = 4/sin(2x), E/(urx+uxp) = cos(2x).
4. The complementary total (srx+crx)+(sxp+cxp) = 2(M+N) =
   2(|sec x| + |csc x|) — the only place M, N survive.

Sanity check (CHECKED, not the proof): numeric evaluation at 11 points
spanning all four quadrants agrees with 4 cot(2x) and with
cos(2x)*(urx+uxp) to machine precision.
