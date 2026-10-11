from sage.all import *
# E76: Corner parameter test.
# T_a(x; lam, mu) = [[urx, srx, lam],[cxp, -urx, uxp],[mu, -uxp, 0]]
# For which (lam, mu) do we still get rank 8 + closure?

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta(deg, lam, mu):
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    return np.array([[urx, srx, lam],[cxp, -urx, uxp],[mu, -uxp, 0]], dtype=float)

def test_closure(lam, mu):
    # Build 8 generators.
    G=[]
    for p in [15,30,45,60]:
        G.append(Ta(p, lam, mu))
    G.append(G[0]@G[1]-G[1]@G[0])
    G.append(G[0]@G[2]-G[2]@G[0])
    G.append(G[0]@G[3]-G[3]@G[0])
    A=G[0]@G[1]-G[1]@G[0]
    G.append(A@G[2]-G[2]@A)
    # Rank.
    M=np.array([g.flatten() for g in G]).T  # 9x8
    r=np.linalg.matrix_rank(M, tol=1e-8)
    if r != 8:
        return r, False
    # Closure: 28 commutators in span?
    # Use least squares.
    U,sv,vh=np.linalg.svd(M)
    # Project each commutator onto span, check residual.
    for i in range(8):
        for j in range(i+1, 8):
            comm=(G[i]@G[j]-G[j]@G[i]).flatten()
            # Coefficients via pseudoinverse.
            coeffs=np.linalg.lstsq(M, comm, rcond=None)[0]
            recon=M@coeffs
            if np.linalg.norm(recon-comm) > 1e-6:
                return r, False
    return r, True

print("Grid search over (lam, mu):")
print("  (rank, closed?)")
vals=[-2,-1,-0.5,0,0.5,1,2]
good=[]
for lam in vals:
    row=[]
    for mu in vals:
        r, closed = test_closure(lam, mu)
        tag="Y" if (r==8 and closed) else "."
        row.append(tag)
        if r==8 and closed:
            good.append((lam,mu))
    print(f"  lam={lam:5.1f}: {' '.join(row)}")
print()
print(f"mu vals:        {' '.join([f'{m:5.1f}' for m in vals])}")
print()
print(f"Good (rank 8 + closed): {good}")
print()
# Refine around (0,0).
print("Refined near (0,0):")
for lam in [0, 0.01, 0.1]:
    for mu in [0, 0.01, 0.1]:
        r, closed = test_closure(lam, mu)
        print(f"  ({lam},{mu}): rank={r}, closed={closed}")
print()
print("WROTE /root/e76_corner.txt")
with open('/root/e76_corner.txt','w') as f:
    f.write(f"good={good}\n")
