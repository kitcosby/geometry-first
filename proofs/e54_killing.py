from sage.all import *
# E54: Killing form on our 8 sl(3,R) generators.
# Split into compact (so(3), 3-dim) vs non-compact (5-dim).
# Killing(X,Y) = Tr(ad_X ad_Y). Signature should be (5,3) for sl(3,R).

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta(deg):
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)

# 8 basis (from E47)
G=[]
for p in [15,30,45,60]: G.append(Ta(p))
G.append(Ta(15)@Ta(30)-Ta(30)@Ta(15))
G.append(Ta(15)@Ta(45)-Ta(45)@Ta(15))
G.append(Ta(15)@Ta(60)-Ta(60)@Ta(15))
A=Ta(15); B=Ta(30); C=Ta(45)
G.append((A@B-B@A)@C - C@(A@B-B@A))
print(f"8 generators defined.")
# Compute structure constants: [G_i,G_j] = f_ijk G_k.
# Solve for f via least squares (8 unknowns, 9 equations per (i,j)).
n=8
# Build matrix M where M[k,:] = G_k.flatten(). Then f_ijk = M^{-1} @ [G_i,G_j].flatten().
M=np.array([g.flatten() for g in G]).T  # 9x8
# Use pseudoinverse
Mpinv=np.linalg.pinv(M)
f=np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        comm=G[i]@G[j]-G[j]@G[i]
        f[i,j,:]=Mpinv@comm.flatten()
print(f"Structure constants computed.")
# Killing form: K_ij = sum_{k,l} f_ikl f_jlk  (Tr(ad_i ad_j))
K=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        s=0
        for k in range(n):
            for l in range(n):
                s+=f[i,k,l]*f[j,l,k]
        K[i,j]=s
print(f"\nKilling matrix K (8x8):")
print(np.array2string(K, precision=3, suppress_small=True))
# Eigenvalues of K
eig=np.linalg.eigvalsh(K)
print(f"\nKilling eigenvalues: {np.sort(eig)}")
n_pos=np.sum(eig>1e-6); n_neg=np.sum(eig<-1e-6); n_zero=np.sum(np.abs(eig)<=1e-6)
print(f"Signature: ({n_pos} positive, {n_neg} negative, {n_zero} zero)")
print(f"Expected for sl(3,R): (5,3,0).")
if n_pos==5 and n_neg==3:
    print("MATCH: Signature (5,3) confirms sl(3,R).")
else:
    print(f"MISMATCH: Got ({n_pos},{n_neg}), expected (5,3).")
# Find compact subalgebra: eigenvectors with negative Killing eigenvalue
# (compact generators have negative Killing form in standard convention)
print(f"\nWROTE /root/e54_killing.txt")
with open('/root/e54_killing.txt','w') as f:
    f.write(f"eig={sorted(eig)}\n")
    f.write(f"sig=({n_pos},{n_neg},{n_zero})\n")
