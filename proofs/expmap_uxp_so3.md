# Exponential Map: Orthogonal Placement of uxp and Its Cores

Date: 2026-10-12. Status: PROVED (structure, theta^2 identity) /
CHECKED (numeric matrix at 30deg).
Placement: so(3) cross-product matrix Y for omega = (sxp, cxp, uxp),
the cores of uxp = cxp - sxp. exp(Y) in SO(3) via Rodrigues.

## Results
1. Y^3 = -theta^2 Y with theta^2 = sxp^2 + cxp^2 + uxp^2, so
   exp(Y) = I + (sin theta/theta) Y + ((1-cos theta)/theta^2) Y^2
   exactly (Rodrigues). Verified against scipy expm at 30deg.
2. theta^2 identity (symbolic, QI t-form):
   theta^2 = 2 (sxp^2 - sxp*cxp + cxp^2)
           = 2 (t^4 - t^3 + 2 t^2 + t + 1)/(t-1)^2.
   The bracket is the A2 / Eisenstein norm form of the two cores —
   the uxp term contributes exactly the -sxp*cxp cross term doubled,
   because uxp = cxp - sxp.
3. At x=30deg: omega = (2-sqrt3, sqrt3, 2sqrt3-2),
   theta = sqrt(26 - 12 sqrt3) = 2.2837 rad, axis (0.117, 0.758, 0.641).
   exp(Y) orthogonal, det 1, tr = 1 + 2 cos theta. (Matrix entries
   numeric: sin/cos of theta are not algebraic.)
4. Companion fact (standard Pfaffian): the so(4) path placement with
   edges (sxp, uxp, cxp) has pf = sxp*cxp — the cross product of the
   cores is the Pfaffian.
Scripts: /tmp/uxpexp.py, /tmp/uxpexp2.py (this run).
