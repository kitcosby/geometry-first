from sage.all import *
# E79: Symmetric space SL(4,R)/SO(4).
# 9-dimensional. Metric from trace form. Sectional curvature.

import numpy as np

# Standard sl(4,R) basis: E_{ij} (i!=j) + Cartan.
# so(4) = antisymmetric (6-dim). p = symmetric traceless (9-dim).
# Basis for p (orthonormal w.r.t. Tr(XY)):
# Diagonal: D1=diag(1,-1,0,0)/sqrt2, D2=diag(1,1,-2,0)/sqrt6, D3=diag(1,1,1,-3)/sqrt12
# Off-diagonal symmetric: (E_{ij}+E_{ji})/sqrt2 for i<j (6 of them).
# Total: 3 + 6 = 9.

def Eij(i,j):
    M=np.zeros((4,4)); M[i,j]=1; return M

p_basis=[]
# Diagonal.
p_basis.append(np.diag([1,-1,0,0])/np.sqrt(2))
p_basis.append(np.diag([1,1,-2,0])/np.sqrt(6))
p_basis.append(np.diag([1,1,1,-3])/np.sqrt(12))
# Off-diagonal.
for i in range(4):
    for j in range(i+1,4):
        p_basis.append((Eij(i,j)+Eij(j,i))/np.sqrt(2))

assert len(p_basis)==9
print(f"p dimension: {len(p_basis)}")
# Verify orthonormal.
Gram=np.array([[np.trace(p_basis[a]@p_basis[b]) for b in range(9)] for a in range(9)])
print(f"Max |Gram - I|: {np.max(np.abs(Gram-np.eye(9))):.2e}")
print()

# Sectional curvature: K(a,b) = -Tr([[Ea,Eb],Eb] @ Ea).
print("Sectional curvature matrix (9x9, symmetric, zero diagonal):")
Ksec=np.zeros((9,9))
for a in range(9):
    for b in range(9):
        if a==b: continue
        Xa, Xb = p_basis[a], p_basis[b]
        comm_ab = Xa@Xb - Xb@Xa
        R = comm_ab@Xb - Xb@comm_ab
        Ksec[a,b] = -np.trace(R @ Xa)

# Print compactly.
np.set_printoptions(precision=2, suppress=True)
print(Ksec)
print()
print(f"Max K: {np.max(Ksec):.4f} (expect 0 on diagonal)")
print(f"Min K (off-diag): {np.min(Ksec+np.eye(9)*1e9):.4f} (expect <0)")
# Scalar curvature.
scalar=2*np.sum([Ksec[a,b] for a in range(9) for b in range(a+1,9)])
print(f"Scalar curvature: {scalar:.4f}")
print()
# Rank of symmetric space (max flat dimension).
# For SL(n,R)/SO(n), rank = n-1 = 3.
print("Rank of SL(4,R)/SO(4) = 3 (max dimension of flat totally geodesic submanifold).")
print("Our Cartan (diagonal matrices) gives a 3-dim flat.")
# Verify: K(Di, Dj) should be 0 for diagonal Di, Dj (they commute).
print("K(D1,D2), K(D1,D3), K(D2,D3) (should be 0, flat Cartan):")
print(f"  {Ksec[0,1]:.2e}, {Ksec[0,2]:.2e}, {Ksec[1,2]:.2e}")
