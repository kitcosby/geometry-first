# H5 Results — Exact Boundary Deltas (one-sided limits + residues)

Date: 2026-10-12. Part of harmonics_testing_plan.md. Status: PROVED
(exact symbolic one-sided limits from the Procedure-1 normal forms in
t = tan(x/2); residues = lim (x-x0) O(x), exact).

## Jump table at sec-boundaries x = pi/2 and 3pi/2 (left | right)
srx: 1|1 continuous. sxp: 1|1 continuous.
cxp: +inf|0. crx: 0|+inf. urx: 1|-inf. uxp: +inf|-1.
sawup: 1|0. sawdown: 0|-1. FW: +1|-1.

## Jump table at csc-boundaries x = 0 and pi (left | right)
srx: 0|+inf. sxp: +inf|0.
cxp: 1|1 continuous. crx: 1|1 continuous.
urx: -1|+inf. uxp: -inf|1.
sawup: -1|0. sawdown: 0|1. FW: -1|+1.

## Residues (every pole is simple, residue exactly +/-2)
x=0: srx right +2, urx right +2, sxp left -2, uxp left +2.
x=pi/2: cxp left -2, uxp left -2, crx right +2, urx right -2.
x=pi: sxp left -2, uxp left +2, srx right +2, urx right +2.
x=3pi/2: cxp left -2, uxp left -2, crx right +2, urx right -2.

## Structural conclusions
1. At every boundary exactly one family is continuous, and it sits at
   value exactly 1 (the equilibrium): csc-boundaries keep cxp=crx=1;
   sec-boundaries keep srx=sxp=1.
2. The singular family swaps 0 <-> +inf across the boundary.
3. FW flips sign at every boundary; the saw pair hands off its
   partition ((1,0) -> (0,-1) at pi/2, etc.); saws never blow up —
   they only jump (numerators finite, urx/uxp nonzero inside quadrants
   since their numerators are +/-(1+t^2)).
4. All finite one-sided boundary values lie in {-1, 0, 1}.
5. Uniform pole strength: residue +/-2 everywhere (the 2 is the
   half-angle Jacobian: t ~ x/2 at 0/pi, t-1 ~ (x-pi/2) at pi/2).
Scripts: /tmp/h5.py (limits), /tmp/h5b.py (residues), this run.
