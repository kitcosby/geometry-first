from sage.all import *
# E53: Coincidence audit - Lie closure test.
# E46 found 8 matrices spanning sl(3,R) (rank 8).
# But span != algebra. Test: is the span CLOSED under [.,.]?
# If yes: it's (isomorphic to) sl(3,R). If no: just a vector space.

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

# The 8 basis from E47
phases=[15,30,45,60]
G=[]
for p in phases: G.append(Ta(p))
# Add commutators [Ta15,Ta30], [Ta15,Ta45], [Ta15,Ta60]
G.append(Ta(15)@Ta(30)-Ta(30)@Ta(15))
G.append(Ta(15)@Ta(45)-Ta(45)@Ta(15))
G.append(Ta(15)@Ta(60)-Ta(60)@Ta(15))
# Add double commutator
A=Ta(15); B=Ta(30); C=Ta(45)
G.append((A@B-B@A)@C - C@(A@B-B@A))
print(f"8 basis matrices defined.")
# Verify rank 8
vecs=np.array([g.flatten() for g in G])
r0=np.linalg.matrix_rank(vecs, tol=1e-8)
print(f"Rank: {r0} (should be 8)")
# Test closure: for each i<j, is [G_i,G_j] in span?
print(f"\nTesting {28} commutators for closure...")
failures=[]
for i in range(8):
    for j in range(i+1,8):
        comm=G[i]@G[j]-G[j]@G[i]
        # Check if comm is in span: rank should not increase
        test_vecs=np.vstack([vecs, comm.flatten()])
        r=np.linalg.matrix_rank(test_vecs, tol=1e-8)
        if r>r0:
            failures.append((i,j))
            #print(f"  [G{i},G{j}] NOT in span (rank {r0}->{r})")
if not failures:
    print("SUCCESS: All 28 commutators in span. CLOSED under bracket.")
    print("The 8 matrices form a Lie algebra (8-dim => isomorphic to sl(3,R)).")
    print("Coincidence #1 is REAL, not numerology.")
else:
    print(f"FAILURE: {len(failures)} commutators outside span.")
    print(f"Examples: {failures[:5]}")
    print("The 8 matrices are just a vector space, NOT a Lie algebra.")
    print("Coincidence #1 is NUMEROLOGY (no algebraic structure).")
print()
print("WROTE /root/e53_closure.txt")
with open('/root/e53_closure.txt','w') as f:
    f.write(f"closed={len(failures)==0}\n")
    f.write(f"failures={len(failures)}\n")
