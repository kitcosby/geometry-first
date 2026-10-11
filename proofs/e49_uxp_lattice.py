from sage.all import *
# E49: uxp-adapted finite difference.
# Define Delta_alpha such that uxp(alpha+Delta) = 2*uxp(alpha) (radius doubling).
# D_uxp f(alpha) = (f(alpha+Delta)-f(alpha))/Delta.
# This uses uxp to set the lattice spacing, not as a multiplicative q.

import math

def uxp_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    return (1+s)/c - 1/((1+c)/s)

def srx_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9: return None
    return abs(1/s)+c/s

def find_delta(alpha, factor=2.0):
    # Find Delta>0 such that uxp(alpha+Delta)=factor*uxp(alpha).
    u0=uxp_phase(alpha)
    if u0 is None: return None
    target=factor*u0
    lo, hi=0.01, 89.9-alpha
    if uxp_phase(alpha+hi) is None or uxp_phase(alpha+hi)<target:
        return None  # cannot double within domain
    for _ in range(50):
        mid=(lo+hi)/2
        v=uxp_phase(alpha+mid)
        if v is None: return None
        if v<target: lo=mid
        else: hi=mid
    return (lo+hi)/2

def D_uxp(f, alpha, factor=2.0):
    d=find_delta(alpha, factor)
    if d is None: return None
    f0=f(alpha); f1=f(alpha+d)
    if f0 is None or f1 is None: return None
    return (f1-f0)/d

print("uxp-adapted finite difference: step = radius doubling.")
print("alpha | uxp(alpha) | Delta (deg) | D_uxp(srx) | d/dx(srx)")
print("------|------------|-------------|-----------|----------")
for a in [10, 20, 30, 45, 60, 75]:
    u=uxp_phase(a)
    d=find_delta(a, 2.0)
    if d is None:
        print(f"{a:5.1f} | {u:10.4f} | {'cannot double':>11s} |")
        continue
    dq=D_uxp(srx_phase, a, 2.0)
    h=0.01
    d_ord=(srx_phase(a+h)-srx_phase(a-h))/(2*h)
    print(f"{a:5.1f} | {u:10.4f} | {d:11.4f} | {dq:9.4f} | {d_ord:8.4f}")
print()
print("Interpretation:")
print("  The lattice spacing Delta is set by uxp doubling (physical radius).")
print("  At small alpha (early), uxp~1, Delta is large (coarse).")
print("  At large alpha (late), uxp large, Delta is small (fine).")
print("  This is INVERSE to q-calculus: discreteness decreases with time.")
print("WROTE /root/e49_uxp_lattice.txt")
