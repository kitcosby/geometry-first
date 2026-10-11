from sage.all import *
# E74: Symmetric space SL(3,R)/SO(3).
# p = 5-dim symmetric traceless (tangent at identity).
# Metric from trace form. Sectional curvature via R(X,Y)Z = -[[X,Y],Z].

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

# p = symmetric parts.
S=[(g+g.T)/2 for g in G]
# Orthonormal basis for span(S) w.r.t. <X,Y>=Tr(XY).
SM=np.array([s.flatten() for s in S]).T  # 9x8
U,sv,vh=np.linalg.svd(SM)
rank=np.sum(sv>1e-8)
print(f"p dimension: {rank}")
U5=U[:,:rank]  # 9x5 orthonormal columns
# Basis matrices E_a (3x3).
E=[U5[:,a].reshape(3,3) for a in range(rank)]
# Verify orthonormal: Tr(E_a E_b) = delta.
print("Orthonormality check (Tr(E_a E_b)):")
Gram=np.array([[np.trace(E[a]@E[b]) for b in range(rank)] for a in range(rank)])
print(f"  Max |Gram - I|: {np.max(np.abs(Gram-np.eye(rank))):.2e}")
print()

# Sectional curvature: K(E_a, E_b) = <R(E_a,E_b)E_b, E_a> / (|Ea|^2|Eb|^2 - <Ea,Eb>^2)
# R(X,Y)Z = -[[X,Y],Z]. For orthonormal basis, denominator = 1.
# K(a,b) = -Tr( [[E_a,E_b],E_b] @ E_a ).
print("Sectional curvatures K(E_a, E_b):")
Ksec=np.zeros((rank,rank))
for a in range(rank):
    for b in range(rank):
        if a==b:
            Ksec[a,b]=0
            continue
        Xa, Xb = E[a], E[b]
        comm_ab = Xa@Xb - Xb@Xa  # in so(3)
        # [[Xa,Xb],Xb]
        R = (comm_ab@Xb - Xb@comm_ab)
        # <R, Xa> = Tr(R @ Xa). K = -Tr(R Xa) / 1.
        Ksec[a,b] = -np.trace(R @ Xa)
print(np.array2string(Ksec, precision=3, suppress_small=True))
print()
print(f"Max sectional curvature: {np.max(Ksec):.4f}")
print(f"Min sectional curvature: {np.min(Ksec):.4f}")
print("(Expect <= 0 for SL(3,R)/SO(3).)")
print()
# Find flat directions (K=0).
print("Pairs with |K| < 1e-8 (flat):")
flat=[]
for a in range(rank):
    for b in range(a+1, rank):
        if abs(Ksec[a,b])<1e-8:
            flat.append((a,b))
            print(f"  ({a},{b})")
if not flat:
    print("  None exactly flat (may be near-zero).")
    # Find smallest.
    min_val=np.min(np.abs(Ksec+np.eye(rank)*1e9))  # exclude diagonal
    print(f"  Smallest |K| off-diagonal: {min_val:.2e}")
print()
# Scalar curvature (sum over pairs).
# For orthonormal basis, scalar = 2*sum_{a<b} K(a,b).
scalar=2*np.sum([Ksec[a,b] for a in range(rank) for b in range(a+1,rank)])
print(f"Scalar curvature: {scalar:.4f}")
print()
print("WROTE /root/e74_symspace.txt")
with open('/root/e74_symspace.txt','w') as f:
    f.write(f"maxK={np.max(Ksec):.6f}\nminK={np.min(Ksec):.6f}\nscalar={scalar:.6f}\n")
