from sage.all import *
# E85: 3^n pattern. n x n uses angle 3^{n-3} x.
# Avoids 2^k aliasing. Test to 11.

import math
import numpy as np
import random

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Xn_mat_3pow(deg, n):
    # X_3 = T_a. X_{m+1}=[[X_m, v],[w^T,0]], v,w use 3^{m-2}x.
    # m=3->4: 3^1=3x. m=4->5: 3^2=9x. n: 3^{n-3}x.
    if n==3:
        srx,sxp,cxp,crx,urx,uxp=ops(deg)
        return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
    else:
        Xprev = Xn_mat_3pow(deg, n-1)
        m = n-1
        pw = m-2  # 3^pw x. m=3->pw=1 (3x). m=10->pw=8 (6561x).
        ang = (3**pw) * math.radians(deg)
        s_k = math.sin(ang); c_k = math.cos(ang)
        v=np.zeros(m); v[0]=s_k; v[1]=c_k
        w=np.zeros(m); w[0]=s_k; w[1]=c_k
        X=np.zeros((n,n))
        X[0:m,0:m]=Xprev
        X[0:m,m]=v
        X[m,0:m]=w
        return X

print("3^n pattern: n x n uses 3^{n-3} x.")
print()
for n in range(3, 12):
    expected=n*n-1
    G=[Xn_mat_3pow(p,n) for p in [15,30,45,60]]
    G.append(G[0]@G[1]-G[1]@G[0])
    G.append(G[0]@G[2]-G[2]@G[0])
    G.append(G[0]@G[3]-G[3]@G[0])
    A=G[0]@G[1]-G[1]@G[0]
    G.append(A@G[2]-G[2]@A)
    Gfull=list(G)
    # Grow.
    for _ in range(300):
        changed=False
        # Sample pairs.
        pairs=[]
        if len(Gfull)<40:
            pairs=[(i,j) for i in range(len(Gfull)) for j in range(i+1,len(Gfull))]
        else:
            for __ in range(300):
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
                Gfull.append(comm); changed=True
                if r1>=expected: break
        if not changed: break
        Mall=np.array([g.flatten() for g in Gfull]).T
        if np.linalg.matrix_rank(Mall,tol=1e-7)>=expected: break
    Mall=np.array([g.flatten() for g in Gfull]).T
    r=np.linalg.matrix_rank(Mall,tol=1e-7)
    # Closure sample.
    closed=True
    if r==expected:
        for (i,j) in [(a,b) for a in range(min(15,len(Gfull))) for b in range(a+1,min(15,len(Gfull)))][:80]:
            comm=(Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]).flatten()
            coeff=np.linalg.lstsq(Mall,comm,rcond=None)[0]
            if np.linalg.norm(Mall@coeff-comm)>1e-5:
                closed=False; break
    st="OK" if (r==expected and closed) else "FAIL"
    print(f"  n={n:2d}: dim={r:3d} (expect {expected:3d}), closed={closed}, {st}", flush=True)
print()
print("Done.")
