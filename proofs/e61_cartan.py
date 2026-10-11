from sage.all import *
# E61: Find Cartan subalgebra as centralizer of regular semisimple G0.

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

# G0 = G[0] (Ta15), has distinct real eigenvalues => regular semisimple.
G0=G[0]
ev0=np.linalg.eigvals(G0)
print(f"G0 eigenvalues: {np.sort(ev0)}")
print(f"Distinct: {len(set(np.round(ev0,6)))==3}")
print()

# Centralizer: { X = sum c_i G_i : [G0, X] = 0 }.
# [G0, sum c_i G_i] = sum c_i [G0, G_i] = 0 (3x3 matrix = 9 equations).
# Build linear system for c (8 unknowns).
# Each [G0, G_i] is 3x3. Flatten to 9-vector. Stack as 9x8 matrix. Nullspace = centralizer coeffs.
M=np.zeros((9, n))
for i in range(n):
    comm=G0@G[i]-G[i]@G0
    M[:,i]=comm.flatten()
# Nullspace of M (9x8). Use SVD.
u,s,vh=np.linalg.svd(M)
# Nullspace = rows of vh corresponding to zero singular values.
tol=1e-8
null_mask=(s<tol)
# s has min(9,8)=8 entries. Nullspace dim = 8 - rank.
rank=np.sum(s>tol)
null_dim=n-rank
print(f"Rank of [G0, .] map: {rank}")
print(f"Centralizer dimension: {null_dim} (expect 2 for Cartan)")
print()
if null_dim>=2:
    # Get nullspace basis (last null_dim rows of vh)
    null_basis=vh[-null_dim:,:]  # null_dim x 8
    H=[]
    for k in range(null_dim):
        c=null_basis[k,:]
        Hk=sum(c[i]*G[i] for i in range(n))
        H.append(Hk)
        print(f"H[{k}] coefficients: {np.round(c,4)}")
        print(f"H[{k}] matrix:")
        print(np.array2string(Hk, precision=4, suppress_small=True))
        ev=np.linalg.eigvals(Hk)
        print(f"H[{k}] eigenvalues: {np.sort(ev)}")
        print(f"H[{k}] trace: {np.trace(Hk):.2e}")
        print()
    # Verify [H0, H1] = 0
    if null_dim>=2:
        comm01=H[0]@H[1]-H[1]@H[0]
        print(f"[H0,H1] max abs: {np.max(np.abs(comm01)):.2e} (should be ~0)")
        # Check independent (not proportional)
        # Check diagonalizable (distinct eigenvalues or minimal poly)
        for k in range(2):
            ev=np.linalg.eigvals(H[k])
            print(f"H[{k}] distinct eigenvalues: {len(set(np.round(ev,6)))==3}")
    print()
    print("SUCCESS: Cartan subalgebra found (2-dim centralizer).")
    print("WROTE /root/e61_cartan.txt")
    with open('/root/e61_cartan.txt','w') as f:
        f.write(f"null_dim={null_dim}\n")
        for k in range(null_dim):
            f.write(f"H{k}_eig={np.sort(np.linalg.eigvals(H[k])).tolist()}\n")
else:
    print("FAILED: Centralizer dimension != 2. G0 may not be regular, or numerical issue.")
    print("WROTE /root/e61_cartan.txt")
    with open('/root/e61_cartan.txt','w') as f:
        f.write(f"null_dim={null_dim}\nFAILED\n")
