from sage.all import *
# E87: Small intervals. Test base b for b^{n-3}x.
# b = 1.01364 (Pythagorean comma), 1.05946 (semitone), 1.5 (fifth).

import math
import numpy as np
import random

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def max_n_for_base(b):
    # Returns highest n (3..11) that reaches full dim with base b.
    for n in range(3,12):
        expected=n*n-1
        # Build X_n.
        def Xn(deg, nn):
            if nn==3:
                srx,sxp,cxp,crx,urx,uxp=ops(deg)
                return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
            else:
                Xp=Xn(deg, nn-1); m=nn-1; pw=m-2
                ang=(b**pw)*math.radians(deg)
                s_k=math.sin(ang); c_k=math.cos(ang)
                v=np.zeros(m); v[0]=s_k; v[1]=c_k
                w=np.zeros(m); w[0]=s_k; w[1]=c_k
                X=np.zeros((nn,nn)); X[0:m,0:m]=Xp; X[0:m,m]=v; X[m,0:m]=w
                return X
        G=[Xn(p,n) for p in [15,30,45,60]]
        G.append(G[0]@G[1]-G[1]@G[0]); G.append(G[0]@G[2]-G[2]@G[0])
        G.append(G[0]@G[3]-G[3]@G[0]); A=G[0]@G[1]-G[1]@G[0]; G.append(A@G[2]-G[2]@A)
        Gfull=list(G)
        for _ in range(200):
            ch=False; pairs=[]
            if len(Gfull)<40:
                pairs=[(i,j) for i in range(len(Gfull)) for j in range(i+1,len(Gfull))]
            else:
                for __ in range(200):
                    i=random.randrange(len(Gfull)); j=random.randrange(len(Gfull))
                    if i<j: pairs.append((i,j))
            for (i,j) in pairs:
                comm=Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]
                Mall=np.array([g.flatten() for g in Gfull]).T
                r0=np.linalg.matrix_rank(Mall,tol=1e-7)
                if r0>=expected: break
                M2=np.array([g.flatten() for g in Gfull+[comm]]).T
                r1=np.linalg.matrix_rank(M2,tol=1e-7)
                if r1>r0:
                    Gfull.append(comm); ch=True
                    if r1>=expected: break
            if not ch: break
            Mall=np.array([g.flatten() for g in Gfull]).T
            if np.linalg.matrix_rank(Mall,tol=1e-7)>=expected: break
        Mall=np.array([g.flatten() for g in Gfull]).T
        r=np.linalg.matrix_rank(Mall,tol=1e-7)
        if r<expected:
            return n-1  # last good
    return 11

bases={
    "Pythagorean comma": (3/2)**12 / 2**7,
    "Chromatic semitone": 2**(1/12),
    "1.1": 1.1,
    "1.25": 1.25,
    "Perfect fifth (3/2)": 1.5,
}
print("Max n reached for each base:")
for name,b in bases.items():
    print(f"  {name} (b={b:.5f}): ", end="", flush=True)
    maxn=max_n_for_base(b)
    print(f"n={maxn} (sl({maxn}), rank {maxn-1})")
print("Done.")
