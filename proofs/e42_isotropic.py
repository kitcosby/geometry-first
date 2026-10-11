from sage.all import *
# E42: Isotropic coordinate map.
# uxp(x,y,z) = uxp_phase(alpha(r)), r=sqrt(x^2+y^2+z^2).
# By monotonicity (E32), invert: alpha(r) = uxp_phase^{-1}(r).
# Compute alpha(r) numerically for r in [1, 22] (uxp range in QI).

import math

def uxp_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    return (1+s)/c - 1/((1+c)/s)

# Build inverse via bisection: given r_target, find alpha such that uxp_phase(alpha)=r_target.
def alpha_of_r(r_target):
    lo, hi = 0.5, 89.5
    for _ in range(60):
        mid=(lo+hi)/2
        v=uxp_phase(mid)
        if v is None: return None
        if v < r_target: lo=mid
        else: hi=mid
    return (lo+hi)/2

print("Isotropic map: alpha(r) = uxp_phase^{-1}(r)")
print("r (radius) | alpha (deg) | alpha (rad)")
print("-----------|-------------|------------")
for r in [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0]:
    a = alpha_of_r(r)
    if a: print(f"{r:10.2f} | {a:11.2f} | {math.radians(a):10.4f}")
print()
print("Interpretation:")
print("  Spatial (x,y,z) -> r=sqrt(x^2+y^2+z^2) -> alpha=alpha(r) -> operators.")
print("  uxp(x,y,z) = r (by identification). urx(x,y,z) = urx_phase(alpha(r)).")
print("  Isotropy: operators depend only on r, not angles.")
print("WROTE /root/e42_map.txt")
with open('/root/e42_map.txt','w') as f:
    for r in [1.0, 2.0, 5.0, 10.0, 20.0]:
        a=alpha_of_r(r)
        f.write(f"r={r}: alpha={a:.2f}deg\n")
