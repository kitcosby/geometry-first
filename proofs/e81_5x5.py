from sage.all import *
# E81: 5x5 via 4x (2^2). Pattern: 2^n.
# X4 from E78 (4x4 with 2x in v,w). Now X5 = [[X4, v4],[w4^T, 0]].
# v4 = w4 = [sin(4x), cos(4x), 0, 0]^T, via operator double-angle twice.

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def X4mat(deg):
    # 4x4 from E78.
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    A=np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
    s2 = 1/(urx+uxp)  # sin(2x)/4
    c2 = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    v=np.array([s2, c2, 0], dtype=float)
    w=np.array([s2, c2, 0], dtype=float)
    X=np.zeros((4,4))
    X[0:3,0:3]=A; X[0:3,3]=v; X[3,0:3]=w
    return X

def X5mat(deg):
    X4=X4mat(deg)
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    # 2x via operators.
    S2 = 1/(urx+uxp)  # sin(2x)/4
    C2 = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    # 4x = 2*(2x).
    # sin(4x) = 2*sin(2x)*cos(2x) = 2*(4*S2)*C2 = 8*S2*C2.
    S4 = 8*S2*C2
    # cos(4x) = cos^2(2x) - sin^2(2x) = C2^2 - (4*S2)^2.
    C4 = C2**2 - 16*S2**2
    v4=np.array([S4, C4, 0, 0], dtype=float)
    w4=np.array([S4, C4, 0, 0], dtype=float)
    X5=np.zeros((5,5))
    X5[0:4,0:4]=X4
    X5[0:4,4]=v4
    X5[4,0:4]=w4
    # X5[4,4]=0 (trace: tr(X4)=0).
    return X5

print("5x5 via 4x (2^2). X5=[[X4,v4],[w4^T,0]], v4=w4=[sin4x,cos4x,0,0]^T")
print()
# Verify 4x identities at a test point.
deg=30
srx,sxp,cxp,crx,urx,uxp=ops(deg)
S2=1/(urx+uxp); C2=(srx+crx-sxp-cxp)/(urx+uxp)
S4=8*S2*C2; C4=C2**2-16*S2**2
r=math.radians(deg)
print(f"At x=30deg:")
print(f"  sin(4x) via ops: {S4:.6f}, actual: {math.sin(4*r):.6f}")
print(f"  cos(4x) via ops: {C4:.6f}, actual: {math.cos(4*r):.6f}")
print()

# Build generators.
G=[]
for p in [15,30,45,60]:
    G.append(X5mat(p))
# Add commutators.
G.append(G[0]@G[1]-G[1]@G[0])
G.append(G[0]@G[2]-G[2]@G[0])
G.append(G[0]@G[3]-G[3]@G[0])
A=G[0]@G[1]-G[1]@G[0]
G.append(A@G[2]-G[2]@A)

M=np.array([g.flatten() for g in G]).T  # 25x8
r=np.linalg.matrix_rank(M, tol=1e-8)
print(f"Rank of 8 generators: {r}")

# Grow.
Gfull=list(G)
for i in range(len(G)):
    for j in range(i+1, len(G)):
        comm=G[i]@G[j]-G[j]@G[i]
        Mall=np.array([g.flatten() for g in Gfull]).T
        r0=np.linalg.matrix_rank(Mall, tol=1e-8)
        Mall2=np.array([g.flatten() for g in Gfull+[comm]]).T
        r1=np.linalg.matrix_rank(Mall2, tol=1e-8)
        if r1>r0:
            Gfull.append(comm)

Mfull=np.array([g.flatten() for g in Gfull]).T
r_full=np.linalg.matrix_rank(Mfull, tol=1e-8)
print(f"After commutators: {len(Gfull)} generators, rank {r_full} (expect 24 for sl(5))")

# Closure.
n=len(Gfull)
closed=True
for i in range(n):
    for j in range(i+1, n):
        comm=(Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]).flatten()
        coeffs=np.linalg.lstsq(Mfull, comm, rcond=None)[0]
        recon=Mfull@coeffs
        if np.linalg.norm(recon-comm)>1e-6:
            closed=False; break
    if not closed: break
print(f"Closed: {closed}, Dimension: {r_full}")
if closed and r_full==24:
    print("  -> sl(5) (rank 4). The 2^n pattern holds!")
elif not closed:
    print("  -> Not closed. The 4x choice may need adjustment.")
