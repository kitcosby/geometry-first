# Crossing M1 = [[urx,srx],[-cxp,FW]] with M2 = [[uxp,sxp],[-crx,FW]]

Date: 2026-10-12. Status: PROVED (trace identity, all quadrants) /
CHECKED (entry forms). Kit's two 2x2 cross-family blocks, crossed by
commutator, product, anticommutator, Kronecker.

## Findings
1. No entry-wise collapse: [M1,M2] entries are genuine rational
   functions of t (contrast the Jacobi pair transparency). tr[M1,M2]=0.
2. PROVED scalar invariant, all 4 quadrants:
     tr(M1 * M2) = 3   identically.
   Equivalent operator identity (using FW^2 = 1):
     urx*uxp - srx*crx - sxp*cxp = 2.
   (Hand-check in t: numerator collapses to 2t(1-t^2) over t(1-t^2).)
   Hence tr{M1,M2} = 6 as well.
3. [M1,M2]_(1,2) in QI = 4t/(t^2-1) = -(urx+uxp) sin^2(x).
4. Quadrant symmetry: QIII is QI with M1 <-> M2 (traces/dets swap);
   QIV is QII swapped — the pi-shift exchanges the two blocks.
5. Exact at 30deg: det M1 = 5 + 8 sqrt(3)/3 = 9.6188,
   det M2 = -3 + 8 sqrt(3)/3 = 1.6188. NOTE: det M2 is NOT the golden
   ratio (1.61803) — a near-miss at the 3rd decimal, recorded to avoid
   a false pattern. M1 is elliptic at 30deg (complex eigenvalues).
6. Kronecker cross: det/tr factor standardly (det^2 det^2, trM1 trM2);
   nothing special.
Scripts: /tmp/cross2.py, /tmp/cross2b.py, /tmp/cross2c.py (this run).
