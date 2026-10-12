# Cross-Family Jacobi Structure

Date: 2026-10-12. Status: PROVED (symbolic in s,c,u,v and in t).
Question (Kit): pairs cancel by inversion — what about cross-family?

## Findings
1. Connector brackets (only cross-family links in the sl(8) cycle):
   [X2,X3]=(c/s) E24, [X4,X5]=(u/c) E46, [X6,X7]=(v/u) E68,
   [X8,X1]=(s/v) E82 — coefficients are the adjacent cross products
   (sxp*cxp, crx*urx, sawup*uxp, sawdown*srx). NO cancellation.
   BUT their product telescopes: (c/s)(u/c)(v/u)(s/v) = 1 — cross-family
   holonomy around the full cycle is also trivial.
2. Cross-family sl(2) (E107 block, X=srx E12, Y=cxp E21, s,c independent):
   [X,Y] = (srx*cxp)(E11-E22); with H=[X,Y]: [H,X]=2(srx*cxp) X,
   [H,Y]=-2(srx*cxp) Y; Jacobi(X,Y,H)=0 verified.
   => The cross product srx*cxp is the COUPLING CONSTANT: it sets the
   scale of the Cartan element and does not drop out (contrast pairs).
3. Cross products pair among themselves (two splittings of prod=1):
   (srx*cxp)(sxp*crx) = 1 and (srx*crx)(sxp*cxp) = 1.
4. Mediators: urx = srx - crx, uxp = cxp - sxp are themselves
   cross-family differences; urx+uxp = 4/sin2x, urx*uxp = 4/|sin2x|,
   uxp = urx/(srx*crx) (proved earlier). H4 quarter-turn is the
   cross-family symmetry: srx(x+pi/2)=crx, cxp(x+pi/2)=sxp.

Summary: intra-family pairs are Jacobi-transparent (factor 1);
cross-family brackets carry genuine coupling constants (cross products),
whose cycle-product is nevertheless 1.
Script: /tmp/cross.py (this run).
