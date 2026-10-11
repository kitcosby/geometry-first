from sage.all import *
# E69: Axiom 0 -- opposite handedness as DUAL representation.
# W0: G_i (3-dim, C3<0). Z: H_i = -G_i^T (dual 3bar, C3>0).
# What real Lie algebra do {G_i, H_i} generate?

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
# Dual: H_i = -G_i^T.
H=[-g.T for g in G]
print("Verifying H_i are traceless:")
for i,h in enumerate(H):
    if abs(np.trace(h))>1e-10:
        print(f"  H[{i}] trace={np.trace(h):.2e} NOT traceless!")
print("  All traceless.")
print()
# Combined.
All=G+H  # 16 matrices
Mflat=np.array([m.flatten() for m in All]).T
rank=np.linalg.matrix_rank(Mflat, tol=1e-8)
print(f"Rank of {{G_i, H_i}} (16 matrices): {rank}")
print()
# Close under brackets.
extended=list(All)
for i in range(len(All)):
    for j in range(i+1, len(All)):
        comm=All[i]@All[j]-All[j]@All[i]
        M2=np.array([m.flatten() for m in extended]).T
        coeff=np.linalg.lstsq(M2, comm.flatten(), rcond=None)[0]
        proj=(M2@coeff).reshape(3,3)
        if np.max(np.abs(comm-proj))>1e-6:
            extended.append(comm)
            if len(extended)>40: break
    if len(extended)>40: break
Mext=np.array([m.flatten() for m in extended]).T
rank_ext=np.linalg.matrix_rank(Mext, tol=1e-8)
print(f"After brackets: {len(extended)} matrices, rank {rank_ext}")
print()
# Is it still sl(3,R) (rank 8) or larger?
if rank_ext==8:
    print("CLOSED: {G_i, H_i} still span sl(3,R) (8-dim).")
    print("The dual rep is already contained; no new dimensions.")
    print("Opposite handedness does NOT expand the algebra.")
elif rank_ext>8:
    print(f"EXPANDED: rank {rank_ext} > 8. New structure!")
else:
    print(f"Unexpected rank {rank_ext}.")
print()
print("WROTE /root/e69_dual.txt")
with open('/root/e69_dual.txt','w') as f:
    f.write(f"rank={rank}\nrank_ext={rank_ext}\n")
