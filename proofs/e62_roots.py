from sage.all import *
# E62: Root system of sl(3,R) from our Cartan (H0, H1).
# Compute ad_{H0}, ad_{H1} (8x8), find joint eigenvalues => roots.

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
# Cartan basis from E61: H1=G0=Ta15, H0=linear combo.
# Recompute H0 coefficients from E61 output.
c0=np.array([0.0, 0.5407, -0.7909, 0.2788, 0.0306, -0.0553, 0.0217, 0.0001])
H0=sum(c0[i]*G[i] for i in range(n))
H1=G[0]  # Ta15
# Express H0, H1 in G-basis (for adjoint).
# H0 = sum c0[i] G_i. H1 = G_0.
# ad_H = sum coeff_i * ad_{G_i}.
# ad_{G_i} as 8x8: (ad_i)[k,j] = f_{i,j,k}.
def ad_of(coeffs):
    # coeffs: 8-vector, H = sum coeffs[i] G_i. Return 8x8 ad_H.
    ad=np.zeros((n,n))
    for i in range(n):
        if abs(coeffs[i])<1e-15: continue
        for j in range(n):
            for k in range(n):
                ad[k,j]+=coeffs[i]*f[i,j,k]
    return ad

ad_H0=ad_of(c0)
ad_H1=ad_of(np.array([1,0,0,0,0,0,0,0]))
# Verify they commute
comm_ad=ad_H0@ad_H1-ad_H1@ad_H0
print(f"[ad_H0, ad_H1] max abs: {np.max(np.abs(comm_ad)):.2e}")
print()
# Eigenvalues of ad_H0
ev0, V0=np.linalg.eig(ad_H0)
print("ad_H0 eigenvalues:")
print(np.sort(ev0))
print()
# Transform ad_H1 to ad_H0 eigenbasis
# V0 columns are eigenvectors. B = V0^{-1} @ ad_H1 @ V0.
try:
    V0inv=np.linalg.inv(V0)
    B=V0inv@ad_H1@V0
    # B should be block-diagonal. Extract diagonal (for 1-dim blocks).
    # Find indices where ev0 is (near) zero vs non-zero.
    print("Joint eigenvalues (alpha(H0), alpha(H1)) for non-zero roots:")
    print("Root | alpha(H0) | alpha(H1)")
    print("-----|-----------|----------")
    roots=[]
    for idx in range(n):
        if abs(ev0[idx])>1e-6:  # non-zero root
            # B[idx,idx] should be alpha(H1). Check off-diagonal small.
            off=np.max(np.abs(np.delete(B[idx,:], idx)))
            aH0=ev0[idx]
            aH1=B[idx,idx]
            # Only real parts (should be real for split Cartan)
            print(f"  {len(roots)}  | {aH0.real:+.4f}  | {aH1.real:+.4f}  (off-diag max {off:.1e})")
            roots.append((aH0.real, aH1.real))
    print()
    print(f"Found {len(roots)} non-zero roots (expect 6).")
    # Check ± pairs
    print()
    print("Root vectors (alpha(H0), alpha(H1)):")
    for r in roots:
        print(f"  ({r[0]:+.4f}, {r[1]:+.4f})")
    print()
    print("WROTE /root/e62_roots.txt")
    with open('/root/e62_roots.txt','w') as f:
        for r in roots:
            f.write(f"{r[0]:.6f},{r[1]:.6f}\n")
except np.linalg.LinAlgError as e:
    print(f"LinAlgError: {e}")
    print("V0 not invertible; ad_H0 may not be diagonalizable.")
    print("WROTE /root/e62_roots.txt")
    with open('/root/e62_roots.txt','w') as f:
        f.write("FAILED\n")
