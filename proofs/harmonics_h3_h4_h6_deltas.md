# H3/H4/H6 Results — Delta Identities (exact, in t = tan(x/2))

Date: 2026-10-12. Part of harmonics_testing_plan.md. Status: PROVED
(exact symbolic identities, all 4 quadrants, sympy simplification to zero).
H1 (Fourier) dropped by Kit's ruling 2026-10-12: series = approximation
(Gibbs); the limit method is exact.

## H4 (headline): the quarter-turn swaps the families
srx(x + pi/2) = crx(x),   cxp(x + pi/2) = sxp(x)   for all x in D.
Consequences (all verified):
- Delta_{pi/2} srx = -urx,  Delta_{pi/2} cxp = -uxp.
  => The UNA functions ARE the negated quarter-turn deltas of the
  primitives: urx = srx(x) - srx(x+pi/2), uxp = cxp(x) - cxp(x+pi/2).
- sxp(x+pi/2) = cxp(x), crx(x+pi/2) = srx(x): the quarter-turn is the
  involution (srx sxp cxp crx) -> (crx cxp sxp srx); period pi follows.
- urx(x+pi/2) = -urx(x), uxp(x+pi/2) = -uxp(x) (anti-period pi/2);
  FW(x+pi/2) = -FW(x), consistent with sgn(sin(2x+pi)).

## H3: doubling delta (branching confirmed exactly)
srx(2x) equals (srx-sxp)/2 on FW=+1 quadrants (QI, QIII) and
(crx-cxp)/2 on FW=-1 quadrants (QII, QIV) — symbolically, via
t2 = 2t/(1-t^2) and alpha' = FW.
Delta_2 srx closed forms in t:
QI: -(1+t^2)/(2t);  QIII: +(1+t^2)/(2t)  [both = -(srx+sxp)/2]
QII: (1+t^2)/(t(t-1)(t+1));  QIV: t(1+t^2)/((t-1)(t+1)).

## H6: reciprocal delta
Delta_h(1/O) = -Delta_h O / (O(x) O(x+h)) verified for O = srx, cxp
at h = pi/2, all quadrants (pure algebra, as predicted).
Bonus identity (verified all quadrants): uxp = urx/(srx * crx).
Derived partner (follows from verified H4+H6): urx = uxp/(cxp * sxp).
Script: /tmp/h346.py (this run).
