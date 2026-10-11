from sage.all import *
# E73: Exponentiate T_a(x). What does exp(T_a(x)) do to R^3?
# T_a traceless => exp(T_a) in SL(3,R) (det=1, volume-preserving).

import math
import numpy as np
from scipy.linalg import expm

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta(deg):
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)

print("exp(T_a(x)) for x = 15, 30, 45, 60:")
print()
for deg in [15,30,45,60]:
    T=Ta(deg)
    E=expm(T)
    det=np.linalg.det(E)
    ev=np.linalg.eigvals(E)
    # Eigenvalues of T.
    evT=np.linalg.eigvals(T)
    print(f"x={deg}deg:")
    print(f"  T_a eigenvalues: {np.sort(evT)}")
    print(f"  exp(T_a) det: {det:.6f} (should be 1)")
    print(f"  exp(T_a) eigenvalues: {np.sort(ev)}")
    # Classify.
    if np.all(np.abs(ev.imag)<1e-8):
        print(f"  Type: HYPERBOLIC (all real positive eigenvalues). Pure stretch.")
        print(f"  Stretch factors: {sorted(ev.real, reverse=True)}")
    else:
        # Complex pair => rotation.
        print(f"  Type: MIXED (rotation + stretch).")
    print()

# 1-parameter subgroup: t -> exp(t*T_a(15)).
print("Flow t -> exp(t*T_a(15deg)) for t=0, 0.5, 1.0, 2.0:")
print("  (tracks how unit sphere deforms)")
T15=Ta(15)
for t in [0, 0.5, 1.0, 2.0]:
    E=expm(t*T15)
    ev=np.linalg.eigvals(E)
    # Singular values = stretch factors along principal axes.
    sv=np.linalg.svd(E, compute_uv=False)
    print(f"  t={t}: singular values (principal stretches) = {sorted(sv, reverse=True)}")
    print(f"         det = {np.linalg.det(E):.6f}")
print()
print("Interpretation:")
print("  exp(t*T_a) is a volume-preserving linear flow on R^3.")
print("  For x=15deg (real T_a eigenvalues): pure anisotropic stretch.")
print("  The three singular values are the stretch factors along orthogonal axes.")
print("  Their product = 1 (volume preserved).")
print()
print("WROTE /root/e73_exp.txt")
with open('/root/e73_exp.txt','w') as f:
    f.write("done\n")
