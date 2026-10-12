from sage.all import *
# E84: Perfect fifth. Try 3x (and 3/2 x) in 4x4 v,w.
# Does angle-tripling also give sl(4)? Or is 2 special?

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Xmat(deg, mult):
    # 4x4 with v=w=[sin(mult*x), cos(mult*x), 0]^T.
    r=math.radians(deg)
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    A=np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
    s_m = math.sin(mult*r)
    c_m = math.cos(mult*r)
    v=np.array([s_m, c_m, 0], dtype=float)
    w=np.array([s_m, c_m, 0], dtype=float)
    X=np.zeros((4,4))
    X[0:3,0:3]=A; X[0:3,3]=v; X[3,0:3]=w
    return X

def test_mult(mult, label):
    G=[]
    for p in [15,30,45,60]:
        G.append(Xmat(p, mult))
    G.append(G[0]@G[1]-G[1]@G[0])
    G.append(G[0]@G[2]-G[2]@G[0])
    G.append(G[0]@G[3]-G[3]@G[0])
    A=G[0]@G[1]-G[1]@G[0]
    G.append(A@G[2]-G[2]@A)
    # Grow.
    Gfull=list(G)
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            comm=G[i]@G[j]-G[j]@G[i]
            Mall=np.array([g.flatten() for g in Gfull]).T
            r0=np.linalg.matrix_rank(Mall, tol=1e-8)
            if r0>=15: break
            Mall2=np.array([g.flatten() for g in Gfull+[comm]]).T
            r1=np.linalg.matrix_rank(Mall2, tol=1e-8)
            if r1>r0:
                Gfull.append(comm)
    Mfull=np.array([g.flatten() for g in Gfull]).T
    r_full=np.linalg.matrix_rank(Mfull, tol=1e-8)
    # Closure.
    n=len(Gfull); closed=True
    for i in range(n):
        for j in range(i+1, n):
            comm=(Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]).flatten()
            coeffs=np.linalg.lstsq(Mfull, comm, rcond=None)[0]
            recon=Mfull@coeffs
            if np.linalg.norm(recon-comm)>1e-6:
                closed=False; break
        if not closed: break
    print(f"  mult={label}: dim={r_full}, closed={closed}")
    return r_full, closed

print("Testing angle multipliers in 4x4 v,w:")
test_mult(2, "2 (octave)")
test_mult(3, "3 (triple)")
test_mult(1.5, "3/2 (perfect fifth)")
test_mult(4, "4")
test_mult(5, "5")
