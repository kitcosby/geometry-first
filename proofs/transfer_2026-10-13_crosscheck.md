# Cross-Check of the 2026-10-13 Transfer (2x2 Mix / Cartan Extension / Saw Forms)

Date: 2026-10-13, Algebra-first chat. The transfer from the
trigonometry side chat was independently re-verified with this
chat's exact methods (SymPy in t = tan(x/2), Proc-1 normal forms).
No conflicts with this chat's proved record were found.

## Re-verified exactly (all PROVED-consistent)
- S1 brackets: [H,C] = -D, [C,D] = 2K, [K,C] = 2D, [K,D] = 2C,
  [H,K] = 0; 2H+K central in the span; closure{H,C} dim = 4.
- S2/S3/S5 formulas in ALL FOUR quadrants (t normal forms):
  det M = srx/cxp - cxp/srx; det N = -(urx^2/4 + sxp*cxp);
  det A = urx/2; det A-tilde = uxp/2; tr(A*A-tilde) = (urx*uxp)/2.
- S6 inverses: the A-shape inverse multiplies to I for both
  FW = +1 and FW = -1 (det = a); the M inverse formula verified
  symbolically in s, c.
- S7 saw dictionary: all four assignments verified in QI.
  Bonus: 1 + 4pq = (t^2 - 2t - 1)^2/(1 + t^2)^2 in QI — a perfect
  square, so R = (1 + 2t - t^2)/(1 + t^2) is rational in t there.
  Saw-form identities det M = q/p - p/q, sxp*cxp = p/q,
  srx*crx = q/p verified.

## PROMOTION: det B != 0 is PROVED (was CHECKED by 1600-point scan)
det B(t) per quadrant is rational with numerator:
  QI:   -(t^2+1)(t^4 + t^3 + 2t^2 + t + 3)
  QII:   (t^2+1)(t^4 - t^3 + 2t^2 - t + 3)
  QIII: -(t^2+1)(3t^4 - t^3 + 2t^2 - t + 1)
  QIV:  -(t^2+1)(3t^4 + t^3 + 2t^2 + t + 1)
Exact Sturm counts: 0 roots in each quadrant's t-range
(0,1), (1,oo), (-1,0), (-oo,-1). Each quartic is moreover
manifestly positive on its range (QI: all coefficients > 0;
QII: t^3(t-1) + t(2t-1) + 3; QIII/QIV by u = -t substitution).
Signs of det B by quadrant (one exact sample each, no roots =>
constant sign): QI +, QII +, QIII -, QIV +
(samples: 335/36, 85/6, -185/18, 235/12). Hence B^{-1} exists on
the whole domain — S6's existence statement is now PROVED.

## Accepted as labeled (not re-run)
Weight-blind bracket closures (S3, S4, S5) — same bracket trees
as this chat's E96/E88 record. Sampled items remain CHECKED as
the transfer labels them. S8's disclaimer is accurate: no new
Lie algebra; the content is coordinatization and invariants.
Scripts: /tmp/xfer_check.py, /tmp/detB.py (this run).
