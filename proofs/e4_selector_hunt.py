from sage.all import *
# Natural selector hunt: octant pairing, reversal, singularity obstruction.
# Operators: srx=|sec|+tan, crx=|sec|-tan, cxp=|csc|+cot, sxp=|csc|-cot.

print("=== Operator definitions and singularities ===")
# Singularities: sec/tan at pi/2, 3pi/2; csc/cot at 0, pi.
# srx, crx (sec-type): singular at pi/2, 3pi/2.
# cxp, sxp (csc-type): singular at 0, pi.
print("srx,crx singular at pi/2, 3pi/2 (sec-type)")
print("cxp,sxp singular at 0, pi (csc-type)")

print()
print("=== Quadrant pairs ===")
print("P+ = {Q_I (oct 0,1), Q_III (oct 4,5)}: operators (srx,cxp) [sec,csc], sign +")
print("P- = {Q_II (oct 2,3), Q_IV (oct 6,7)}: operators (crx,sxp) [sec,csc], sign -")

print()
print("=== Test x -> x+pi (Q_I <-> Q_III, Q_II <-> Q_IV) ===")
# srx(x+pi) = |sec(x+pi)|+tan(x+pi) = |-sec x|+tan x = |sec x|+tan x = srx(x). PI-PERIODIC.
# Similarly all four are pi-periodic (since |sec|,|csc| have period pi, tan,cot have period pi).
print("srx, crx, cxp, sxp are ALL pi-periodic (|sec|,|csc|,tan,cot all pi-periodic).")
print("So x->x+pi preserves ALL operators. Both pairs 'allow' this map.")
print("This is NOT the distinguisher.")

print()
print("=== Test x -> -x (reversal) ===")
# srx(-x) = |sec(-x)|+tan(-x) = |sec x|-tan x = crx(x).
# cxp(-x) = |csc(-x)|+cot(-x) = |csc x|-cot x = sxp(x).
# So x->-x swaps: srx<->crx, cxp<->sxp. I.e., P+ <-> P-.
print("x->-x swaps srx<->crx, cxp<->sxp. So it swaps P+ <-> P-.")
print("Fixed points of x->-x: x=0, x=pi (mod 2pi).")
print("At x=0,pi: csc-type (cxp,sxp) SINGULAR. sec-type (srx,crx) finite.")
print("So the P+<->P- swap via x->-x is BLOCKED at x=0,pi by csc singularities.")
print("The 'point singularities' are x=0 and x=pi.")

print()
print("=== Which pair is 'allowed' vs 'blocked'? ===")
print("Consider staying WITHIN a pair via x->x+pi:")
print("  P+: (srx,cxp) -> (srx,cxp). Both pi-periodic. No singularity in the MAP.")
print("  P-: (crx,sxp) -> (crx,sxp). Both pi-periodic. No singularity in the MAP.")
print("Consider the REVERSAL x->-x (which swaps pairs):")
print("  Blocked at x=0 (affects cxp in P+, sxp in P-).")
print("  Blocked at x=pi (affects cxp in P+, sxp in P-).")
print()
print("Alternative: 'reversal' = operator order swap within octant?")
print("  Oct 0: (srx,cxp). Oct 1: (cxp,srx). The pair SWAPS across the octant boundary (45 deg).")
print("  This is smooth (no singularity at 45 deg).")
print("  Similarly for all quadrants. So intra-quadrant reversal is always allowed.")
print()
print("=== KEY INSIGHT: sign of tan/cot ===")
print("srx = |sec|+tan, crx = |sec|-tan. They differ by sign of tan.")
print("cxp = |csc|+cot, sxp = |csc|-cot. They differ by sign of cot.")
print("P+ uses (+tan, +cot). P- uses (-tan, -cot).")
print("The 'reversal' tan->-tan maps P+ to P-.")
print("tan = 0 at x=0, pi. At these POINTS, srx=crx and the reversal does nothing (fixed).")
print("But csc blows up at x=0,pi, so cxp,sxp are singular there.")
print("The obstruction is specifically in the csc-type operators at the tan-fixed points.")
print()
print("=== Natural selector candidate ===")
print("P+ = {oct 0,1,4,5} <-> vector directions {E_0,E_1,E_4,E_5}.")
print("P- = {oct 2,3,6,7} <-> vector directions {E_2,E_3,E_6,E_7}.")
print("If P+ is the 'allowed' pair (reversal works), it selects span{E_0,E_1,E_4,E_5}.")
print("Signature: E_0^2=+1,E_1^2=-1,E_4^2=+1,E_5^2=-1 => (+,-,+,-) = (2,2).")
print("NOT Lorentz (3,1), but it IS a natural 4-dim (selected by dynamics, not by choice).")
print()
print("=== Is (2,2) a problem? ===")
print("A (2,2) 4-plane is a neutral signature. In 4D, (2,2) admits:")
print("  - Two independent null directions.")
print("  - A natural complex structure? (2,2) ~ C^2 as real?")
print("  - Or: (2,2) could be (1,1)+(1,1), and we need to break further.")
print("The octant structure gives 4 quadrants, not 4 spacetime dims directly.")
print("Maybe: 4 quadrants -> 4 dims, and the P+/P- pairing gives the signature?")

with open('/root/selector_hunt.txt','w') as f:
    f.write("P+ allowed (pi-periodic, no obstruction in x->x+pi)\n")
    f.write("P- <-> P+ swap via x->-x blocked at x=0,pi (csc singularities)\n")
    f.write("Natural 4-dim: span{E_0,E_1,E_4,E_5}, signature (2,2)\n")
print("WROTE /root/selector_hunt.txt")
