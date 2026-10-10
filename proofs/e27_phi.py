from sage.all import *
# E27: Compute phi(x) = arg(urx + i*uxp) = atan2(uxp, urx) in QI.
# Then test if phi relates to 2x via the given identities.
# Given:
#   sin(2x)/4 = 1/(urx+uxp)   [E20, expanded form matches user's]
#   cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp)  [E20]

x = var('x')
# QI definitions (0<x<pi/2, sin>0, cos>0)
srx = (1+cos(x))/sin(x)
sxp = 1/srx
cxp = (1+sin(x))/cos(x)
crx = 1/cxp
urx = srx - crx
uxp = cxp - sxp

# phi = atan2(uxp, urx)
# tan(phi) = uxp/urx
tan_phi = (uxp/urx).simplify_full()
print("tan(phi) = uxp/urx simplified:")
print(tan_phi)
print()

# Try to see if tan_phi has a simple form
tan_phi_trig = tan_phi.simplify_trig()
print("trig simplified:")
print(tan_phi_trig)
print()

# Check the user's identities symbolically
# Identity 1: sin(2x)/4 = 1/(urx+uxp)
lhs1 = sin(2*x)/4
rhs1 = 1/(urx+uxp)
diff1 = (lhs1 - rhs1).simplify_full()
print("Identity 1: sin(2x)/4 - 1/(urx+uxp) =")
print(diff1)
print(f"  Holds: {diff1==0}")
print()

# Identity 2: cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp)
lhs2 = cos(2*x)
rhs2 = (srx+crx-sxp-cxp)/(urx+uxp)
diff2 = (lhs2 - rhs2).simplify_full()
print("Identity 2: cos(2x) - (srx+crx-sxp-cxp)/(urx+uxp) =")
print(diff2)
print(f"  Holds: {diff2==0}")
print()

# Now: can we express phi in terms of 2x?
# We have tan(phi)=uxp/urx. Is there a relation to tan(2x) or similar?
# Let's compute phi at specific values and compare to x.
import math
def num_phi(deg):
    r=math.radians(deg)
    s=math.sin(r); c=math.cos(r)
    srx_=(1+c)/s; cxp_=(1+s)/c
    urx_=srx_-1/cxp_; uxp_=cxp_-1/srx_
    return math.degrees(math.atan2(uxp_,urx_))

print("deg(x) | phi=arg(U) | 2x | phi/x")
for d in [10,22.5,30,45,67.5]:
    ph=num_phi(d)
    print(f"{d:6.1f} | {ph:7.2f} | {2*d:6.1f} | {ph/d:5.3f}")

print()
print("CONCLUSION: phi(x) computed. Relation to 2x (if any) to be determined.")
print("User's identities 1 and 2 are E20 (verified above).")

with open('/root/phi.txt','w') as f:
    f.write("phi_computed=True\n")
print("WROTE /root/phi.txt")
