# The Quadrant Flip in Polar/Euler Notation

Date: 2026-10-12. Status: PROVED (symbolic in w = e^{ix/2} and in t).
Question (Kit): represent the quadrant flip (QI<->QIII block swap;
srx form flip 1/t <-> -t; FW sign) in polar notation — is there a
complex reason the operators flip in the negative quadrants?

## Setup
zeta = e^{ix} (Euler), w = zeta^{1/2} = e^{ix/2}. The three sign data
are argument data of zeta: alpha = sgn(Im zeta), beta = sgn(Re zeta),
FW = sgn(Im zeta^2). Every flip is zeta (or zeta^2, or w) crossing an
axis of the complex plane.

## The alpha-flip: two perfect squares
srx = i (zeta + zeta^{-1} + 2 alpha)/(zeta - zeta^{-1}).
The numerator is a perfect square in w, in exactly two ways:
  alpha=+1 (upper half-plane): zeta+zeta^{-1}+2 = (w + w^{-1})^2,
    srx = i (w+w^{-1})/(w-w^{-1}) = cot(x/2) = 1/t.
  alpha=-1 (lower half-plane): zeta+zeta^{-1}-2 = (w - w^{-1})^2,
    srx = i (w-w^{-1})/(w+w^{-1}) = -tan(x/2) = -t.
So the flip 1/t <-> -t is the swap of the Euler factors
(w + w^{-1}) = 2cos(x/2) and (w - w^{-1}) = 2i sin(x/2). The absolute
value |csc x| in the definition IS the half-plane datum sgn(Im zeta),
and it selects which square — hence which factor sits upstairs.

## The pi-flip: multiplication by i on the double cover
x -> x+pi sends zeta -> -zeta and w -> i w. On t = tan(x/2) this is
the Möbius involution t -> -1/t, and it maps EVERY operator's QI
formula to the same operator's QIII formula (all 8 verified).
At block level the flip is a sign conjugation:
  M1_QIII = diag(1,-1) * M2_QI * diag(1,-1)   (verified),
so traces/dets swap, as observed in the crossing computation.

## FW
FW = sgn(sin 2x) = sgn(Im zeta^2): the negative quadrants are exactly
where zeta^2 crosses into the lower half-plane.

## Conclusion
The operators live on the half-angle double cover w = e^{ix/2} and
are ratios of the two Euler factors (w +/- w^{-1}). Quadrant flips
are: alpha = choice of square (half-plane of zeta), pi-shift = factor
i on w (Möbius t -> -1/t), block swap = residual sign conjugation,
FW = half-plane of zeta^2. That is the complex reason.
Scripts: /tmp/polar.py, /tmp/polar2.py (this run; three first-pass
test targets were mis-stated and corrected — the involution preserves
each operator, it does not swap urx<->uxp).
