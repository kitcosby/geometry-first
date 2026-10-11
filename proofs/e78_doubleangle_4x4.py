from sage.all import *
# E78: 4x4 with double-angle in v,w.
# X(x) = [[T_a(x), v(x)], [w(x)^T, 0]]
# v = w = [sin(2x)/4, cos(2x), 0]^T
# Does this generate a Lie algebra? What dimension?

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Xmat(deg):
    r=math.radians(deg)
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    # T_a (3x3)
    A=np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
    # Double-angle via operators (CHECKED identities).
    s2 = 1/(urx+uxp)  # sin(2x)/4
    c2 = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    v=np.array([s2, c2, 0], dtype=float)
    w=np.array([s2, c2, 0], dtype=float)
    # 4x4 block.
    X=np.zeros((4,4))
    X[0:3,0:3]=A
    X[0:3,3]=v
    X[3,0:3]=w
    # X[3,3]=0 (trace zero: tr(A)=0, so tr(X)=0).
    return X

print("X(x) = [[T_a, v],[w^T, 0]], v=w=[sin(2x)/4, cos(2x), 0]^T")
print()
# Build generators: 4 X(x_i) + commutators.
G=[]
for p in [15,30,45,60]:
    G.append(Xmat(p))
# Add commutators.
G.append(G[0]@G[1]-G[1]@G[0])
G.append(G[0]@G[2]-G[2]@G[0])
G.append(G[0]@G[3]-G[3]@G[0])
A=G[0]@G[1]-G[1]@G[0]
G.append(A@G[2]-G[2]@A)
# That's 8. But we're in 4x4, so max dim is 15.
# Let's see what dimension we actually get.
M=np.array([g.flatten() for g in G]).T  # 16x8
r=np.linalg.matrix_rank(M, tol=1e-8)
print(f"Rank of 8 generators: {r}")
print()
# Try adding more commutators to see if it grows.
# Compute all pairwise commutators of the 8, add independent ones.
Gfull=list(G)
for i in range(len(G)):
    for j in range(i+1, len(G)):
        comm=G[i]@G[j]-G[j]@G[i]
        # Check if independent.
        Mall=np.array([g.flatten() for g in Gfull]).T
        r0=np.linalg.matrix_rank(Mall, tol=1e-8)
        Mall2=np.array([g.flatten() for g in Gfull+[comm]]).T
        r1=np.linalg.matrix_rank(Mall2, tol=1e-8)
        if r1>r0:
            Gfull.append(comm)
            # print(f"  Added [{i},{j}], rank now {r1}")
print(f"After adding independent commutators: {len(Gfull)} generators, rank {np.linalg.matrix_rank(np.array([g.flatten() for g in Gfull]).T, tol=1e-8)}")
print()
# Check closure of Gfull.
n=len(Gfull)
Mfull=np.array([g.flatten() for g in Gfull]).T
r_full=np.linalg.matrix_rank(Mfull, tol=1e-8)
closed=True
for i in range(n):
    for j in range(i+1, n):
        comm=(Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]).flatten()
        coeffs=np.linalg.lstsq(Mfull, comm, rcond=None)[0]
        recon=Mfull@coeffs
        if np.linalg.norm(recon-comm)>1e-6:
            closed=False
            break
    if not closed:
        break
print(f"Closed: {closed}")
print(f"Dimension: {r_full}")
if closed:
    print(f"This is a {r_full}-dimensional Lie algebra.")
    if r_full==15:
        print("  -> sl(4) ?")
    elif r_full==8:
        print("  -> sl(3) (embedded) ?")
    else:
        print(f"  -> {r_full}-dim, not sl(4) or sl(3).")
