from sage.all import *
# E82: 2^n pattern to 11x11 (sl(11), rank 10, 120-dim).
# X_n uses 2^{n-3}x in v,w. Test closure for n=3..11.

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def sincos_pow2(deg, k):
    # Returns (sin(2^k x), cos(2^k x)) via iterated operator double-angle.
    # s_1 = sin(2x)/4, c_1 = cos(2x). Iterate: s_{j+1}=2*s_j*c_j, c_{j+1}=c_j^2-16*s_j^2.
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    s = 1/(urx+uxp)  # sin(2x)/4
    c = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    # s = sin(2^1 x)/4, c = cos(2^1 x). Need 2^k.
    for j in range(1, k):
        s_new = 2*s*c  # sin(2^{j+1}x)/4
        c_new = c*c - 16*s*s  # cos(2^{j+1}x)
        s, c = s_new, c_new
    # Return actual sin, cos (not /4).
    return 4*s, c

def Xn_mat(deg, n):
    # Recursive: X_3 = T_a. X_{m+1} = [[X_m, v],[w^T, 0]].
    # v,w use 2^{m-2}x? Let's check:
    # X_4 (m=3->4): uses 2^1 x. m-2 = 1. 2^1. ✓
    # X_5 (m=4->5): uses 2^2 x. m-2 = 2. 2^2. ✓
    # X_n: built from X_{n-1}, uses 2^{n-3}x.
    if n==3:
        srx,sxp,cxp,crx,urx,uxp=ops(deg)
        return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)
    else:
        Xprev = Xn_mat(deg, n-1)
        m = n-1  # X_m -> X_{m+1}
        k = m-2  # power: 2^k x. For m=3, k=1 (2x). For m=10, k=8 (256x).
        s_k, c_k = sincos_pow2(deg, k)
        v=np.zeros(m); v[0]=s_k; v[1]=c_k
        w=np.zeros(m); w[0]=s_k; w[1]=c_k
        X=np.zeros((n,n))
        X[0:m,0:m]=Xprev
        X[0:m,m]=v
        X[m,0:m]=w
        return X

print("Testing 2^n pattern for n=3..11.")
print("Expected dim: n^2-1. Checking closure.")
print()
for n in range(3, 12):
    expected = n*n-1
    # Build generators.
    G=[]
    for p in [15,30,45,60]:
        G.append(Xn_mat(p, n))
    # Add commutators until stable or max.
    # Start with 4, add [G0,Gi], then nested.
    G.append(G[0]@G[1]-G[1]@G[0])
    G.append(G[0]@G[2]-G[2]@G[0])
    G.append(G[0]@G[3]-G[3]@G[0])
    A=G[0]@G[1]-G[1]@G[0]
    G.append(A@G[2]-G[2]@A)
    # Grow.
    Gfull=list(G)
    # Limit iterations for large n.
    max_iter = 200 if n<=6 else 500
    it=0
    changed=True
    while changed and it<max_iter and len(Gfull)<expected+5:
        changed=False
        it+=1
        # Try commutators of existing.
        # For large n, sample rather than all pairs.
        import random
        pairs_to_try = []
        if len(Gfull) < 30:
            pairs_to_try = [(i,j) for i in range(len(Gfull)) for j in range(i+1,len(Gfull))]
        else:
            # Sample.
            for _ in range(200):
                i=random.randrange(len(Gfull)); j=random.randrange(len(Gfull))
                if i<j: pairs_to_try.append((i,j))
        for (i,j) in pairs_to_try:
            comm=Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]
            Mall=np.array([g.flatten() for g in Gfull]).T
            r0=np.linalg.matrix_rank(Mall, tol=1e-7)
            if r0>=expected: break
            Mall2=np.array([g.flatten() for g in Gfull+[comm]]).T
            r1=np.linalg.matrix_rank(Mall2, tol=1e-7)
            if r1>r0:
                Gfull.append(comm)
                changed=True
                if r1>=expected: break
        if len(Gfull)>=expected:
            # Check rank.
            Mall=np.array([g.flatten() for g in Gfull]).T
            if np.linalg.matrix_rank(Mall, tol=1e-7)>=expected:
                break
    Mall=np.array([g.flatten() for g in Gfull]).T
    r=np.linalg.matrix_rank(Mall, tol=1e-7)
    # Closure check (sample for large n).
    closed=True
    if r==expected:
        # Check a sample of commutators.
        test_pairs = [(i,j) for i in range(min(20,len(Gfull))) for j in range(i+1,min(20,len(Gfull)))]
        Mfull=Mall
        for (i,j) in test_pairs[:100]:
            comm=(Gfull[i]@Gfull[j]-Gfull[j]@Gfull[i]).flatten()
            coeffs=np.linalg.lstsq(Mfull, comm, rcond=None)[0]
            recon=Mfull@coeffs
            if np.linalg.norm(recon-comm)>1e-5:
                closed=False; break
    status = "OK" if (r==expected and closed) else ("RANK_OK" if r==expected else "FAIL")
    print(f"  n={n:2d}: dim={r:3d} (expect {expected:3d}), closed={closed}, {status}")
print()
print("Done.")
