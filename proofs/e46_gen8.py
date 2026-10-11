from sage.all import *
# E46: Find the 8th generator for sl(3,R).
# E45 got rank 7 from 20 matrices. Try: more phases, 2nd derivatives,
# double commutators, and different operator combinations.

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
    return (Ta(deg+h)-Ta(deg-h))/(2*h*math.pi/180)

def d2Ta(deg, h=0.01):
    return (Ta(deg+h)-2*Ta(deg)+Ta(deg-h))/((h*math.pi/180)**2)

# Start with E45's set
phases=[15,30,45,60,75]
mats=[]
for p in phases: mats.append(Ta(p))
for p in phases: mats.append(dTa(p))
for i in range(len(phases)):
    for j in range(i+1,len(phases)):
        mats.append(Ta(phases[i])@Ta(phases[j])-Ta(phases[j])@Ta(phases[i]))

def rank_of(mats):
    vecs=np.array([m.flatten() for m in mats])
    return np.linalg.matrix_rank(vecs, tol=1e-8)

print(f"Base rank: {rank_of(mats)}")
# Try adding more phases
for p in [10, 20, 40, 50, 70, 80, 85]:
    new_mats=mats+[Ta(p)]
    r=rank_of(new_mats)
    if r>rank_of(mats):
        print(f"  Adding Ta({p}deg): rank {rank_of(mats)} -> {r} *IMPROVEMENT*")
        mats=new_mats
    #else: print(f"  Ta({p}): no improvement")
print(f"After more phases: rank={rank_of(mats)}")
# Try 2nd derivatives
for p in phases:
    new_mats=mats+[d2Ta(p)]
    r=rank_of(new_mats)
    if r>rank_of(mats):
        print(f"  Adding d2Ta({p}): rank -> {r} *IMPROVEMENT*")
        mats=new_mats
print(f"After 2nd derivs: rank={rank_of(mats)}")
# Try double commutators
base_Ta=[Ta(p) for p in phases]
for i in range(3):
    for j in range(3):
        for k in range(3):
            if len({i,j,k})<3: continue
            C1=base_Ta[i]@base_Ta[j]-base_Ta[j]@base_Ta[i]
            C2=C1@base_Ta[k]-base_Ta[k]@C1
            new_mats=mats+[C2]
            r=rank_of(new_mats)
            if r>rank_of(mats):
                print(f"  Adding [[Ta{i},Ta{j}],Ta{k}]: rank -> {r} *IMPROVEMENT*")
                mats=new_mats
                break
        if rank_of(mats)==8: break
    if rank_of(mats)==8: break
print(f"Final rank: {rank_of(mats)}")
if rank_of(mats)==8:
    print("SUCCESS: Full sl(3,R) generators found!")
else:
    print(f"STILL INCOMPLETE: rank {rank_of(mats)}, need {8-rank_of(mats)} more.")
print("WROTE /root/e46_gen8.txt")
with open('/root/e46_gen8.txt','w') as f:
    f.write(f"rank={rank_of(mats)}\n")
