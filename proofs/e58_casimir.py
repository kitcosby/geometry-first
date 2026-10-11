from sage.all import *
# E58: Casimir eigenvalues for our sl(3,R) 3-dim rep.
# Quadratic Casimir C2 = sum_{ij} (K^{-1})^{ij} G_i G_j.
# Should be c2 * I_3. Find c2.

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
# Structure constants, Killing (from E54)
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
# Inverse Killing
Kinv=np.linalg.inv(K)
# Quadratic Casimir: C2 = sum_{ij} Kinv^{ij} G_i G_j (3x3 matrix)
C2=np.zeros((3,3))
for i in range(n):
    for j in range(n):
        C2+=Kinv[i,j]*(G[i]@G[j])
print("Quadratic Casimir C2 (3x3):")
print(np.array2string(C2, precision=6, suppress_small=True))
# Check if proportional to identity
# C2 should be c2*I. Compute c2 = Tr(C2)/3.
c2=np.trace(C2)/3
print(f"\nTr(C2)/3 = {c2:.6f}")
# Check deviation from c2*I
dev=np.max(np.abs(C2 - c2*np.eye(3)))
print(f"Max |C2 - c2*I|: {dev:.2e}")
if dev<1e-6:
    print(f"SUCCESS: C2 = {c2:.6f} * I. Casimir eigenvalue is {c2:.6f}.")
else:
    print(f"NOT proportional to identity (dev={dev:.2e}).")
    print("May need better basis or normalization.")
print()
# Also: eigenvalues of individual G_i (3x3 matrices)
print("Eigenvalues of the 8 generators (3x3 each):")
for i, g in enumerate(G):
    ev=np.linalg.eigvals(g)
    print(f"  G[{i}]: {np.sort(ev)}")
print()
print("WROTE /root/e58_casimir.txt")
with open('/root/e58_casimir.txt','w') as f:
    f.write(f"c2={c2:.6f}\n")
    f.write(f"dev={dev:.2e}\n")
