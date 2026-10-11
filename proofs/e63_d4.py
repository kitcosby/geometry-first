from sage.all import *
# E63: D4 quarter-turn action on the 5 non-compact (symmetric) generators.
# D4 swaps: srx<->crx, sxp<->cxp, hence urx->-urx, uxp->-uxp.
# For T_a = [[urx,srx,0],[cxp,-urx,uxp],[0,-uxp,0]],
# D4(T_a) = [[-urx,crx,0],[sxp,urx,-uxp],[0,uxp,0]].

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

def D4_of_Ta(deg):
    # D4 swaps srx<->crx, sxp<->cxp; urx->-urx, uxp->-uxp.
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    # New entries: srx'=crx, sxp'=cxp, cxp'=sxp, crx'=srx, urx'=-urx, uxp'=-uxp.
    # D4(T_a) = [[urx', srx', 0],[cxp', -urx', uxp'],[0, -uxp', 0]]
    #         = [[-urx, crx, 0],[sxp, urx, -uxp],[0, uxp, 0]].
    return np.array([[-urx, crx, 0],[sxp, urx, -uxp],[0, uxp, 0]], dtype=float)

# Build G_i and D4(G_i).
G=[]
D4G=[]
for p in [15,30,45,60]:
    G.append(Ta(p))
    D4G.append(D4_of_Ta(p))
# Commutators: D4([A,B]) = [D4(A), D4(B)].
def comm(A,B): return A@B-B@A
G.append(comm(G[0],G[1])); D4G.append(comm(D4G[0],D4G[1]))
G.append(comm(G[0],G[2])); D4G.append(comm(D4G[0],D4G[2]))
G.append(comm(G[0],G[3])); D4G.append(comm(D4G[0],D4G[3]))
G.append(comm(comm(G[0],G[1]),G[2]))
D4G.append(comm(comm(D4G[0],D4G[1]),D4G[2]))
n=8
# Symmetric parts (non-compact directions).
S=[(g+g.T)/2 for g in G]
D4S=[(g+g.T)/2 for g in D4G]
# Check: is D4S[i] in span(S)?
# Build matrix from flattened S, project D4S.
SM=np.array([s.flatten() for s in S]).T  # 9x8
SMpinv=np.linalg.pinv(SM)
print("Checking if D4 preserves the symmetric (non-compact) span:")
print("i | residual |D4S[i] - proj|")
print("--|-------------------")
preserved=True
for i in range(n):
    coeff=SMpinv@D4S[i].flatten()
    proj=(SM@coeff).reshape(3,3)
    res=np.max(np.abs(D4S[i]-proj))
    print(f"{i} | {res:.2e}")
    if res>1e-6:
        preserved=False
print()
if preserved:
    print("D4 PRESERVES the non-compact span. Computing action...")
    # Get 8x8 matrix of D4 in S-basis (but S may be rank-deficient; use 5-dim basis).
    # Find rank and basis.
    u,s,vh=np.linalg.svd(SM)
    rank=np.sum(s>1e-8)
    print(f"Rank of symmetric span: {rank} (expect 5)")
    # Use first 'rank' left singular vectors as orthonormal basis for the space.
    # Project D4 action.
    # For each basis vector, express D4.
    # Simpler: find 5x5 matrix via least squares on independent set.
    # Get indices of 5 independent S_i.
    # Use QR with pivoting? Simpler: use SVD basis.
    # Basis: U_rank (9 x rank). D4 acts as linear map on R^9, restrict.
    # Represent D4 as 9x9 (on flattened symmetric matrices).
    # D4 on entries: need 9x9 matrix.
    # Actually, easier: work in the 5-dim coordinate via pinv.
    # Find 5 independent S_i.
    # Find 5 independent S_i via greedy rank check.
    idx=[]
    cur=np.zeros((9,0))
    for i in range(n):
        trial=np.column_stack([cur, S[i].flatten()]) if cur.shape[1]>0 else S[i].flatten().reshape(9,1)
        if np.linalg.matrix_rank(trial, tol=1e-8)>len(idx):
            idx.append(i)
            cur=trial
        if len(idx)==5: break
    idx=np.array(idx)
    print(f"Independent indices: {idx}")
    B5=np.array([S[i].flatten() for i in idx]).T  # 9x5
    B5pinv=np.linalg.pinv(B5)
    # D4 action on these 5: D4S[i] = sum_j M_D4[j,i] S[idx[j]].
    M_D4=np.zeros((5,5))
    for i, ii in enumerate(idx):
        coeff=B5pinv@D4S[ii].flatten()
        M_D4[:,i]=coeff
    print("D4 action matrix (5x5) on non-compact basis:")
    print(np.array2string(M_D4, precision=3, suppress_small=True))
    # Eigenvalues of D4 action.
    ev=np.linalg.eigvals(M_D4)
    print(f"D4 eigenvalues on 5-dim: {np.sort(ev)}")
    print()
    print("Interpretation:")
    print("  Eigenvalue +1: D4-invariant (candidate: gravitational/time).")
    print("  Eigenvalue -1: D4-odd.")
    print("  Complex pair: D4 swaps a 2-dim subspace (candidate: electric/magnetic).")
    print()
    print("WROTE /root/e63_d4.txt")
    with open('/root/e63_d4.txt','w') as f:
        f.write(f"preserved=True\n")
        f.write(f"D4_eig={np.sort(ev).tolist()}\n")
else:
    print("D4 DOES NOT preserve the non-compact span.")
    print("The D4 symmetry does not act within our algebra.")
    print("WROTE /root/e63_d4.txt")
    with open('/root/e63_d4.txt','w') as f:
        f.write(f"preserved=False\n")
