from sage.all import *
# E65: Fix E62 -- 2D root vectors with EXACT centralizer (no rounding).
# Recompute everything in one script.

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
# Structure constants
M=np.array([g.flatten() for g in G]).T
Mpinv=np.linalg.pinv(M)
f=np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        comm=G[i]@G[j]-G[j]@G[i]
        f[i,j,:]=Mpinv@comm.flatten()

# EXACT centralizer of G0 (no hardcoding, no rounding).
G0=G[0]
CM=np.zeros((9,n))
for i in range(n):
    CM[:,i]=(G0@G[i]-G[i]@G0).flatten()
u,s,vh=np.linalg.svd(CM)
rank=np.sum(s>1e-8)
null_dim=n-rank
print(f"Centralizer dim: {null_dim}")
null_basis=vh[-null_dim:,:]
# H0 = first null vector, H1 = G0 (which is in centralizer).
# Get exact coefficients.
c0=null_basis[0,:]
# Ensure c0 is not proportional to [1,0,...] (i.e., not just G0).
# If it is, take the other.
if np.max(np.abs(c0-np.array([1,0,0,0,0,0,0,0])))<1e-6 or np.max(np.abs(c0+np.array([1,0,0,0,0,0,0,0])))<1e-6:
    c0=null_basis[1,:]
    print("Took second null vector for H0 (first was G0).")
H0=sum(c0[i]*G[i] for i in range(n))
H1=G0
# Verify commute (exact).
comm_H=np.max(np.abs(H0@H1-H1@H0))
print(f"[H0,H1] max abs: {comm_H:.2e} (exact, no rounding)")
print()

# Adjoint.
def ad_of(coeffs):
    ad=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                ad[k,j]+=coeffs[i]*f[i,j,k]
    return ad
ad_H0=ad_of(c0)
ad_H1=ad_of(np.array([1,0,0,0,0,0,0,0]))
comm_ad=np.max(np.abs(ad_H0@ad_H1-ad_H1@ad_H0))
print(f"[ad_H0,ad_H1] max abs: {comm_ad:.2e} (should be ~0)")
print()
# Joint eigenvalues.
ev0, V0=np.linalg.eig(ad_H0)
print(f"ad_H0 eigenvalues: {np.sort(ev0)}")
try:
    V0inv=np.linalg.inv(V0)
    Bm=V0inv@ad_H1@V0
    print()
    print("Joint (alpha(H0), alpha(H1)):")
    roots=[]
    for idx in range(n):
        if abs(ev0[idx])>1e-6:
            off=np.max(np.abs(np.delete(Bm[idx,:],idx)))
            a0=ev0[idx].real; a1=Bm[idx,idx].real
            print(f"  ({a0:+.6f}, {a1:+.6f})  off-diag {off:.1e}")
            roots.append((a0,a1))
    print(f"\nFound {len(roots)} roots.")
    # Check A2: largest = sum of two smaller (positive).
    pos=sorted([r for r in roots if r[0]>0])
    if len(pos)>=3:
        print(f"Positive alpha(H0): {[p[0] for p in pos]}")
        # Check additive
        s01=pos[0][0]+pos[1][0]
        print(f"Sum of two smaller: {s01:.6f}, largest: {pos[2][0]:.6f}, diff: {abs(s01-pos[2][0]):.2e}")
    print()
    print("WROTE /root/e65_roots_fix.txt")
    with open('/root/e65_roots_fix.txt','w') as f:
        f.write(f"comm_ad={comm_ad:.2e}\n")
        for r in roots:
            f.write(f"{r[0]:.6f},{r[1]:.6f}\n")
except np.linalg.LinAlgError as e:
    print(f"Failed: {e}")
    with open('/root/e65_roots_fix.txt','w') as f:
        f.write("FAILED\n")
