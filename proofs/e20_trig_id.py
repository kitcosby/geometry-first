from sage.all import *
# Verify:
# (1) sin(2x)/4 = 1/(urx+uxp)
# (2) cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp)   [assumed parsing]
# Definitions: srx=|csc|+cot, sxp=1/srx, cxp=|sec|+tan, crx=1/cxp,
#              urx=srx-crx, uxp=cxp-sxp.
# Work in (0,pi/2) where |.| drops, then note quadrant generality.

x = var('x')
assume(x>0, x<pi/2)
srx = 1/sin(x) + cos(x)/sin(x)
cxp = 1/cos(x) + sin(x)/cos(x)
sxp = 1/srx
crx = 1/cxp
urx = srx - crx
uxp = cxp - sxp

print("=== Identity 1: sin(2x)/4 = 1/(urx+uxp) ===")
lhs1 = sin(2*x)/4
rhs1 = 1/(urx+uxp)
diff1 = (lhs1 - rhs1).simplify_trig()
print(f"LHS - RHS simplified: {diff1}")
print(f"Identity 1: {'PROVED' if diff1==0 else 'NOT proved (diff != 0)'}")
# If not zero, try numeric.
if diff1 != 0:
    for deg in [15,30,45,60]:
        r = RDF(deg*pi/180)
        l = float(lhs1.subs(x=r)); rr = float(rhs1.subs(x=r))
        print(f"  x={deg}: LHS={l:.6f}, RHS={rr:.6f}, match={abs(l-rr)<1e-9}")

print()
print("=== Identity 2: cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp) ===")
lhs2 = cos(2*x)
rhs2 = (srx+crx-sxp-cxp)/(urx+uxp)
diff2 = (lhs2 - rhs2).simplify_trig()
print(f"LHS - RHS simplified: {diff2}")
print(f"Identity 2: {'PROVED' if diff2==0 else 'NOT proved (diff != 0)'}")
if diff2 != 0:
    for deg in [15,30,45,60]:
        r = RDF(deg*pi/180)
        l = float(lhs2.subs(x=r)); rr = float(rhs2.subs(x=r))
        print(f"  x={deg}: LHS={l:.6f}, RHS={rr:.6f}, match={abs(l-rr)<1e-9}")

with open('/root/trig_id.txt','w') as f:
    f.write(f"id1={diff1==0} id2={diff2==0}\n")
print("WROTE /root/trig_id.txt")
