# Jacobi Transparency of the Four Pairs

Date: 2026-10-12. Status: PROVED (symbolic in s,c,u,v; pair relations
a_{2k} = 1/a_{2k-1} are definitional for the operators).
Setting: sl(8) cycle X_i = a_i E_{i,i+1}, a = (s,1/s, c,1/c, u,1/u, v,1/v)
= (srx,sxp, cxp,crx, urx,sawup, uxp,sawdown).

## Results (all verified symbolically)
1. Intra-pair brackets are BARE:
   [X1,X2]=E13, [X3,X4]=E35, [X5,X6]=E57, [X7,X8]=E71 — coefficient 1.
2. Jacobi (X1,X2,X3): [X1,[X2,X3]] = +cxp E14, [X3,[X1,X2]] = -cxp E14,
   middle term 0. Pair 1 (srx*sxp=1) cancels completely; the surviving
   coefficient is the unpaired third operator alone.
3. Jacobi (X2,X3,X4): nonzero terms are +/-(1/s) E25; pair 2 (cxp*crx=1)
   cancels; the connector operator sxp=1/s survives.
4. Cartan is integrally normalized: E31 via the long path (pairs 2,3,4
   all collapse) has coefficient 1, so [E13,E31] = E11 - E33 exactly.
   Integer control 1..8 would give Cartan coefficient 2 * 20160 = 40320.

## Statement
Every complete reciprocal pair is JACOBI-TRANSPARENT: it contributes a
factor of exactly 1 to every bracket and Jacobi coefficient containing
it. Only unpaired (connector-adjacent) operators survive in coefficients.
Equivalently: the operator generators yield the standard Chevalley
normalization of sl(8) (the diagonal conjugacy D of the exp-map note has
trivial holonomy on each pair, and total holonomy prod a_i = 1).
Non-adjacent pairs commute pairwise; pairs interact only through the
four connector brackets [X2,X3], [X4,X5], [X6,X7], [X8,X1].
Script: /tmp/jacobi.py (this run).
