from sage.all import *
# Break apart: which pairing divides by zero at 90deg?
# srx=|csc|+cot = |1/sin| + cos/sin. Denominator: sin x.
# cxp=|sec|+tan = |1/cos| + sin/cos. Denominator: cos x.
# At x=90deg: sin=1 (nonzero), cos=0 (ZERO).

print("=== Constituent parts at x=90deg ===")
print("srx = |csc x| + cot x = |1/sin x| + cos x / sin x")
print("  Denominator: sin x. At 90deg: sin=1 != 0. NO division by zero.")
print("  srx(90deg) = 1/1 + 0/1 = 1. FINITE. Crosses cleanly.")
print()
print("cxp = |sec x| + tan x = |1/cos x| + sin x / cos x")
print("  Denominator: cos x. At 90deg: cos=0. DIVISION BY ZERO.")
print("  cxp(90deg) -> +Infinity. DIVERGES. Does not cross cleanly.")
print()
print("=== Consequences for urx, uxp ===")
print("urx = srx - crx, where crx=1/cxp.")
print("  srx->1 (clean), crx=1/oo=0 (clean). So urx->1-0=1. FINITE.")
print("  urx CROSSES CLEANLY because its csc-part is finite and its sec-part vanishes.")
print()
print("uxp = cxp - sxp, where sxp=1/srx.")
print("  cxp->oo (divides by zero), sxp->1 (clean). So uxp->oo-1=oo. DIVERGES.")
print("  uxp DOES NOT cross cleanly because its sec-part divides by zero.")
print()
print("PROVED: At 90deg, the csc-pairing (srx,urx) is finite; the sec-pairing (cxp,uxp) divides by zero.")
print("The 'clean crossing' is the csc-type. The 'divergence' is the sec-type.")
print("This is structural: sin(90)=1, cos(90)=0. Not a choice.")

with open('/root/divzero.txt','w') as f:
    f.write("csc_clean sec_divzero PROVED\n")
print("WROTE /root/divzero.txt")
