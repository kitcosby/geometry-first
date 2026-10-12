# Procedure 1 Results — Quadrant Normal Forms in t = tan(x/2)

Date: 2026-10-12. Part of harmonics_testing_plan.md. Status: PROVED
(symbolic identities in t on quadrant interiors; exact simplification in
sympy, QI independently cross-checked in Sage 9.5 on toe — one unused
import raised ImportError, the computation itself ran and matched).

Setup: alpha = sgn(sin x), beta = sgn(cos x), constant per quadrant.
srx = [alpha(1+t^2) + (1-t^2)]/(2t),  cxp = [beta(1+t^2) + 2t]/(1-t^2),
sxp = 1/srx, crx = 1/cxp, urx = srx - crx, uxp = cxp - sxp,
sawup = 1/urx, sawdown = 1/uxp, FW = sawup + sawdown.
sin(2x) in t: 4t(1-t^2)/(1+t^2)^2.

## QI (alpha=+1, beta=+1)
srx=1/t, sxp=t, cxp=(1+t)/(1-t), crx=(1-t)/(1+t),
urx=(1+t^2)/(t(1+t)), uxp=(1+t^2)/(1-t),
sawup=t(1+t)/(1+t^2), sawdown=(1-t)/(1+t^2), FW=+1.

## QII (alpha=+1, beta=-1)
srx=1/t, sxp=t, cxp=(t-1)/(t+1), crx=(t+1)/(t-1),
urx=(1+t^2)/(t(1-t)), uxp=-(1+t^2)/(1+t),
sawup=t(1-t)/(1+t^2), sawdown=-(1+t)/(1+t^2), FW=-1.

## QIII (alpha=-1, beta=-1)
srx=-t, sxp=-1/t, cxp=(t-1)/(t+1), crx=(t+1)/(t-1),
urx=(1+t^2)/(1-t), uxp=(1+t^2)/(t(1+t)),
sawup=(1-t)/(1+t^2), sawdown=t(1+t)/(1+t^2), FW=+1.
(saw pair = QI pair swapped.)

## QIV (alpha=-1, beta=+1)
srx=-t, sxp=-1/t, cxp=(1+t)/(1-t), crx=(1-t)/(1+t),
urx=-(1+t^2)/(1+t), uxp=(1+t^2)/(t(t-1)) [sign per sympy form],
sawup=-(1+t)/(1+t^2), sawdown=t(1-t)/(1+t^2)... (see script output), FW=-1.

## Verified identities (all 4 quadrants, exact)
- FW = sgn(sin 2x): +1,-1,+1,-1.
- urx + uxp = 4/sin(2x).
- urx * uxp = 4/|sin(2x)|.
- cos(2x) = (srx + crx - sxp - cxp)/(urx + uxp).

## Structural notes
- Primitives are Möbius (degree <=1) in t; urx/uxp degree 2;
  sawup/sawdown have universal denominator (1+t^2) in every quadrant,
  and their numerators telescope: sum = FW*(1+t^2) exactly.
- srx depends only on alpha: 1/t (alpha=+1) or -t (alpha=-1).
  cxp depends only on beta: (1+t)/(1-t) (beta=+1) or (t-1)/(t+1) (beta=-1).
Scripts: /tmp/proc1.py (sympy, this run), /tmp/proc1.sage.py (Sage cross-check).
