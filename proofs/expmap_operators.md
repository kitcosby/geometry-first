# Exponential Map with the 8 Operators Restored

Date: 2026-10-12. Status: PROVED (items 1-3, 6) / CHECKED numeric (4-5).

Setup: C = sum_i a_i E_{i,i+1} (8x8 cycle, sl(8) construction) with a_i the
8 operators at angle x; Y_i = f_i (E_{i,i+1} - E_{i+1,i}) (so(9) path).

1. PROVED: prod a_i = 1 for all x. The 8 operators are 4 reciprocal pairs
   (srx*sxp = cxp*crx = urx*sawup = uxp*sawdown = 1). Integers 1..8 give
   8! = 40320 instead.
2. PROVED: C^8 = (prod a_i) I = I for all x. So the operator-weighted
   cycle is an exact ORDER-8 element. (Integer control: C^8 = 40320 I.)
   Verified exactly at x=30deg in sympy (C^8 == I, operators in sqrt(3)).
3. PROVED: C is diagonally conjugate to the bare cyclic permutation
   matrix (consistency around the cycle is exactly prod a_i = 1), hence
   charpoly(C) = lambda^8 - 1 and spec(C) = {8th roots of unity},
   INDEPENDENT of x. C is a Coxeter element of the Weyl group S_8.
4. CHECKED: det exp(C) = 1 (tr C = 0); spec(exp C) = {exp(omega^j)}.
5. CHECKED: exp(C) = sum_{r=0}^7 s_r C^r with s_r = sum_{m>=0} 1/(8m+r)!
   s_r ~= (1.000025, 1.000003, 0.5, 1/6, 1/24, 1/120, 1/720, 1/5040)
   — i.e. 1/r! plus tail corrections from C^8 = I.
6. PROVED (standard formula, verified at 30deg): exp(Y_i) is the rotation
   of the (i,i+1)-plane by angle exactly f_i(x) — under exp, the operators
   ARE angles. E.g. exp(Y_1) at 30deg rotates by srx = 2+sqrt(3) rad.

Conclusion: Lie generation only saw nonzeroness, but the exponential map
sees the values: the operators are the canonical normalization making
the cycle a root-of-unity (Coxeter) element, and they become literal
rotation angles on the so(9) side.
Scripts: /tmp/expmap.py (exact), /tmp/expmap2.py, /tmp/expmap3.py (numeric).
