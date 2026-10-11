from sage.all import *
# E56: Find so(3) triple directly.
# so(3) = {X in sl(3,R) | X^T = -X} (anti-symmetric).
# Project our 8 generators onto anti-symmetric subspace.

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
# Anti-symmetrize each: X_anti = (X - X^T)/2
G_anti=[(g-g.T)/2 for g in G]
# What is the rank of the anti-symmetric parts?
vecs=np.array([g.flatten() for g in G_anti])
r=np.linalg.matrix_rank(vecs, tol=1e-8)
print(f"Rank of anti-symmetric parts: {r}")
print(f"Expected for so(3): 3.")
if r==3:
    print("SUCCESS: Anti-symmetric parts span so(3).")
    # Find 3 independent
    # Use SVD to get orthonormal basis for the 3-dim space
    u,s,vh=np.linalg.svd(vecs, full_matrices=False)
    print(f"Singular values: {s[:5]}")
    # The top 3 right-singular vectors (rows of vh) give orthonormal basis in R^9.
    # Reshape to 3x3.
    print(f"\nOrthonormal so(3) basis (from SVD):")
    for i in range(3):
        J=vh[i].reshape(3,3)
        print(f"\nJ_{i}:")
        print(np.array2string(J, precision=4, suppress_small=True))
        # Verify anti-symmetric
        print(f"  Max |J+J^T|: {np.max(np.abs(J+J.T)):.2e}")
        # Verify so(3) brackets: [J0,J1] should be proportional to J2, etc.
    # Check brackets
    J0=vh[0].reshape(3,3); J1=vh[1].reshape(3,3); J2=vh[2].reshape(3,3)
    c01=J0@J1-J1@J0; c12=J1@J2-J2@J1; c20=J2@J0-J0@J2
    print(f"\nBracket checks (should be proportional to basis):")
    # [J0,J1] = a*J0 + b*J1 + c*J2. Solve.
    M=np.array([J0.flatten(), J1.flatten(), J2.flatten()]).T
    for name, comm in [("[J0,J1]",c01), ("[J1,J2]",c12), ("[J2,J0]",c20)]:
        coeff, *_=np.linalg.lstsq(M, comm.flatten(), rcond=None)
        print(f"  {name} = {coeff[0]:.3f}*J0 + {coeff[1]:.3f}*J1 + {coeff[2]:.3f}*J2")
else:
    print(f"MISMATCH: rank {r}, not 3.")
print()
print("WROTE /root/e56_so3.txt")
