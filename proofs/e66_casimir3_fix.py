from sage.all import *
# E66: Fix E59 -- cubic Casimir by SOLVING (not guessing normalization).
# Find symmetric c_{ijk} such that C3 = sum c_{ijk} G_i G_j G_k commutes with all G_l.
# Then check if C3 is scalar.

import math
import numpy as np
from itertools import combinations_with_replacement

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
# Symmetric triples (i<=j<=k). Index map.
triples=list(combinations_with_replacement(range(n),3))
N=len(triples)
print(f"Symmetric triples: {N}")
# For each triple, the symmetrized product S_{ijk} = sum over perms of G_i G_j G_k.
# C3 = sum_t c_t * S_t. Want [C3, G_l]=0 for all l.
# [S_t, G_l] is 3x3. Flatten to 9. Stack: (8*9) x N matrix. Nullspace = solutions.
def sym_prod(i,j,k):
    # Sum over all distinct permutations.
    import itertools
    perms=set(itertools.permutations([i,j,k]))
    S=np.zeros((3,3))
    for (a,b,c_) in perms:
        S+=G[a]@G[b]@G[c_]
    return S

Eq=np.zeros((8*9, N))
for t,(i,j,k) in enumerate(triples):
    S=sym_prod(i,j,k)
    for l in range(n):
        comm=S@G[l]-G[l]@S
        Eq[l*9:(l+1)*9, t]=comm.flatten()
# Nullspace.
u,s,vh=np.linalg.svd(Eq)
rank=np.sum(s>1e-7)
null_dim=N-rank
print(f"Rank: {rank}, null dim: {null_dim} (expect >=1 for cubic Casimir)")
if null_dim>=1:
    # Take first null vector.
    c=vh[-1,:]
    # Build C3.
    C3=np.zeros((3,3))
    for t,(i,j,k) in enumerate(triples):
        C3+=c[t]*sym_prod(i,j,k)
    print("C3:")
    print(np.array2string(C3, precision=4, suppress_small=True))
    c3=np.trace(C3)/3
    dev=np.max(np.abs(C3-c3*np.eye(3)))
    print(f"Tr(C3)/3 = {c3:.6f}")
    print(f"Max |C3-c3*I|: {dev:.2e}")
    if dev<1e-4:
        print(f"SUCCESS: C3 = {c3:.6f} * I. Cubic Casimir eigenvalue = {c3:.6f}.")
    else:
        print("C3 not scalar. May need different symmetrization or it's not central.")
    print()
    print("WROTE /root/e66_casimir3_fix.txt")
    with open('/root/e66_casimir3_fix.txt','w') as f:
        f.write(f"c3={c3:.6f}\ndev={dev:.2e}\n")
else:
    print("No cubic Casimir found (null dim 0).")
    with open('/root/e66_casimir3_fix.txt','w') as f:
        f.write("FAILED\n")
