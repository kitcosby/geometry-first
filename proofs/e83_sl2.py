from sage.all import *
# E83: Negative n. 2x2 traceless from operators. Does it give sl(2)?
# M(x) = [[urx/2, sxp],[cxp, -urx/2]]. Traceless 2x2.
# Take 4 x's, commutators. Check for 3-dim closed (sl(2), rank 1).

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def M2(deg):
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    # Traceless 2x2.
    return np.array([[urx/2, sxp],[cxp, -urx/2]], dtype=float)

print("2x2 traceless M(x)=[[urx/2, sxp],[cxp, -urx/2]]")
print()
G=[]
for p in [15,30,45,60]:
    G.append(M2(p))
    # print(f"  x={p}: tr={np.trace(M2(p)):.2e}")

# Commutators.
G.append(G[0]@G[1]-G[1]@G[0])
G.append(G[0]@G[2]-G[2]@G[0])
G.append(G[0]@G[3]-G[3]@G[0])

M=np.array([g.flatten() for g in G]).T  # 4x7
r=np.linalg.matrix_rank(M, tol=1e-8)
print(f"Rank of 7 matrices: {r} (expect 3 for sl(2))")
print()

# If rank>3, find independent set and check closure.
# 2x2 traceless is at most 3-dim.
if r>=3:
    # Get 3 independent.
    U,sv,vh=np.linalg.svd(M)
    # Take first 3 with sv>tol.
    # Actually, let's just take G[0], G[1], and [G0,G1] if independent.
    G3=[]
    for g in G:
        if len(G3)==0:
            G3.append(g)
        else:
            Mall=np.array([h.flatten() for h in G3]).T
            r0=np.linalg.matrix_rank(Mall, tol=1e-8)
            Mall2=np.array([h.flatten() for h in G3+[g]]).T
            r1=np.linalg.matrix_rank(Mall2, tol=1e-8)
            if r1>r0 and r1<=3:
                G3.append(g)
        if len(G3)==3:
            break
    print(f"Found {len(G3)} independent 2x2 traceless.")
    if len(G3)==3:
        M3=np.array([g.flatten() for g in G3]).T
        # Closure: 3 commutators.
        closed=True
        for i in range(3):
            for j in range(i+1,3):
                comm=(G3[i]@G3[j]-G3[j]@G3[i]).flatten()
                coeffs=np.linalg.lstsq(M3, comm, rcond=None)[0]
                recon=M3@coeffs
                if np.linalg.norm(recon-comm)>1e-6:
                    closed=False
                    print(f"  [{i},{j}] not in span, residual {np.linalg.norm(recon-comm):.2e}")
        print(f"Closed: {closed}")
        if closed:
            print("  -> sl(2) (rank 1). Negative n works!")
        # Killing signature for sl(2,R) should be (2,1).
        # Compute structure constants.
        n=3
        f=np.zeros((n,n,n))
        for i in range(n):
            for j in range(n):
                comm=(G3[i]@G3[j]-G3[j]@G3[i]).flatten()
                coeffs=np.linalg.lstsq(M3, comm, rcond=None)[0]
                f[i,j,:]=coeffs
        Kmat=np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                s=0
                for k in range(n):
                    for l in range(n):
                        s+=f[i,k,l]*f[j,l,k]
                Kmat[i,j]=s
        ev=np.linalg.eigvalsh(Kmat)
        pos=np.sum(ev>1e-6); neg=np.sum(ev<-1e-6)
        print(f"Killing signature: ({pos},{neg}) (expect (2,1) for sl(2,R))")
