from sage.all import *
# E45: Generator search for sl(3,R).
# Start: T_a(x) = [[urx, srx, 0],[cxp,-urx,uxp],[0,-uxp,0]] (E30, in sl(3,R)).
# Generate: T_a(x_i) at 5 phases, dT_a/dx, commutators [T_a(x_i),T_a(x_j)].
# Question: does the span have dimension 8 (= dim sl(3,R))?

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

def dTa(deg, h=0.01):
    # Numerical derivative dT_a/dx (x in degrees)
    return (Ta(deg+h)-Ta(deg-h))/(2*h*math.pi/180)  # d/dx (x in radians)

# Build generator set
phases=[15,30,45,60,75]
mats=[]
for p in phases:
    mats.append(Ta(p))
for p in phases:
    mats.append(dTa(p))
# Commutators
for i in range(len(phases)):
    for j in range(i+1, len(phases)):
        A=Ta(phases[i]); B=Ta(phases[j])
        mats.append(A@B - B@A)

print(f"Generated {len(mats)} matrices (5 Ta + 5 dTa + 10 commutators).")
# Flatten to vectors in R^9, compute rank
vecs=np.array([m.flatten() for m in mats])
rank=np.linalg.matrix_rank(vecs, tol=1e-8)
print(f"Rank of span: {rank}")
print(f"dim sl(3,R) = 8.")
if rank==8:
    print("SUCCESS: Generators span sl(3,R)!")
elif rank<8:
    print(f"INCOMPLETE: span is {rank}-dim, need {8-rank} more independent generators.")
else:
    print(f"UNEXPECTED: rank {rank} > 8 (should be <=8 for 3x3 traceless).")
# Check: are all in sl(3,R)? (trace zero)
traces=[np.trace(m) for m in mats]
print(f"Max |trace|: {max(abs(t) for t in traces):.2e} (should be ~0)")
print("WROTE /root/e45_gen.txt")
with open('/root/e45_gen.txt','w') as f:
    f.write(f"rank={rank}\n")
