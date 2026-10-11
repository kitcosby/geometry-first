from sage.all import *
# E55: Extract explicit so(3) (compact) and non-compact generators.
# From E54: Killing has 3 negative (compact), 5 positive (non-compact) eigenvalues.
# New basis: H_a = sum_i v^{(a)}_i G_i, where v^{(a)} are Killing eigenvectors.

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

G=[]
for p in [15,30,45,60]: G.append(Ta(p))
G.append(Ta(15)@Ta(30)-Ta(30)@Ta(15))
G.append(Ta(15)@Ta(45)-Ta(45)@Ta(15))
G.append(Ta(15)@Ta(60)-Ta(60)@Ta(15))
A=Ta(15); B=Ta(30); C=Ta(45)
G.append((A@B-B@A)@C - C@(A@B-B@A))
n=8
# Structure constants and Killing (from E54)
M=np.array([g.flatten() for g in G]).T
Mpinv=np.linalg.pinv(M)
f=np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        comm=G[i]@G[j]-G[j]@G[i]
        f[i,j,:]=Mpinv@comm.flatten()
K=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        s=0
        for k in range(n):
            for l in range(n):
                s+=f[i,k,l]*f[j,l,k]
        K[i,j]=s
# Eigendecomposition of K
eigvals, eigvecs=np.linalg.eigh(K)  # eigvecs[:,a] is eigenvector for eigvals[a]
idx=np.argsort(eigvals)
eigvals=eigvals[idx]; eigvecs=eigvecs[:,idx]
print(f"Killing eigenvalues (sorted): {eigvals}")
# Negative: compact (so(3)). Positive: non-compact.
compact_idx=[i for i in range(n) if eigvals[i]<-1e-6]
noncompact_idx=[i for i in range(n) if eigvals[i]>1e-6]
print(f"Compact (neg): {len(compact_idx)} gens, indices {compact_idx}")
print(f"Non-compact (pos): {len(noncompact_idx)} gens, indices {noncompact_idx}")
print()
# Build explicit matrices: H_a = sum_i v_i^{(a)} G_i
def build_H(a):
    v=eigvecs[:,a]
    H=np.zeros((3,3))
    for i in range(n):
        H+=v[i]*G[i]
    return H
print("=== COMPACT so(3) TRIPLE (3 generators) ===")
for a in compact_idx:
    H=build_H(a)
    print(f"\nH_compact[{a}] (Killing eig={eigvals[a]:.2f}):")
    print(np.array2string(H, precision=4, suppress_small=True))
    # Verify: should be anti-symmetric? For so(3) in defining rep, generators are anti-symmetric?
    # Actually sl(3,R) defining rep: compact generators are anti-symmetric (so(3) subset).
    asym=np.max(np.abs(H+H.T))
    print(f"  Max |H+H^T| (0 if anti-symmetric): {asym:.2e}")
print()
print("=== NON-COMPACT (5 generators) ===")
for a in noncompact_idx:
    H=build_H(a)
    print(f"\nH_noncompact[{a}] (Killing eig={eigvals[a]:.2f}):")
    print(np.array2string(H, precision=4, suppress_small=True))
print()
print("WROTE /root/e55_split.txt")
with open('/root/e55_split.txt','w') as f:
    f.write(f"compact_idx={compact_idx}\n")
    f.write(f"noncompact_idx={noncompact_idx}\n")
    for a in compact_idx:
        f.write(f"\nCompact[{a}]:\n{build_H(a)}\n")
    for a in noncompact_idx:
        f.write(f"\nNoncompact[{a}]:\n{build_H(a)}\n")
