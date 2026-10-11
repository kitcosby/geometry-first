from sage.all import *
# E47: Explicit 8 sl(3,R) basis + D4 quarter-turn intertwiner search.
# Part 1: Extract 8 independent matrices from E46's set.
# Part 2: Find X (4x4) s.t. X*[srx,sxp,cxp,crx](x) = [srx,sxp,cxp,crx](x+pi/2).

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

# --- Part 1: Explicit 8 basis ---
# Use FULL E46 set: 5 Ta + 5 dTa + 10 commutators + double commutator
phases=[15,30,45,60,75]
mats=[]
labels=[]
for p in phases:
    mats.append(Ta(p)); labels.append(f"Ta({p})")
def dTa_num(deg, h=0.01):
    r1=Ta(deg+h); r2=Ta(deg-h)
    return (r1-r2)/(2*h*math.pi/180)
for p in phases:
    mats.append(dTa_num(p)); labels.append(f"dTa({p})")
for i in range(len(phases)):
    for j in range(i+1, len(phases)):
        A=Ta(phases[i]); B=Ta(phases[j])
        mats.append(A@B-B@A); labels.append(f"[Ta{phases[i]},Ta{phases[j]}]")
# Add double commutator (the 8th)
A=Ta(15); B=Ta(30); C=Ta(45)
D=(A@B-B@A)@C - C@(A@B-B@A)
mats.append(D); labels.append("[[Ta15,Ta30],Ta45]")
# Find 8 independent via QR with pivoting
vecs=np.array([m.flatten() for m in mats])
# Use SVD to find rank and select independent
u,s,vh=np.linalg.svd(vecs, full_matrices=False)
rank=np.sum(s>1e-8)
print(f"Part 1: Rank={rank} from {len(mats)} matrices.")
# Greedy select 8 independent
selected=[]; sel_labels=[]
remaining=list(range(len(mats)))
for _ in range(8):
    best_i=None; best_r=0
    curr_rank=len(selected)
    for i in remaining:
        test=selected+[i]
        r=np.linalg.matrix_rank(vecs[test], tol=1e-8)
        if r>best_r:
            best_r=r; best_i=i
    if best_i is not None and best_r>len(selected):
        selected.append(best_i); sel_labels.append(labels[best_i])
        remaining.remove(best_i)
print(f"Selected {len(selected)} independent: {sel_labels}")
print("Explicit basis matrices (flattened):")
for i, lab in zip(selected, sel_labels):
    print(f"\n{lab}:")
    print(np.array2string(mats[i], precision=4, suppress_small=True))
print()
with open('/root/e47_basis.txt','w') as f:
    f.write(f"Selected: {sel_labels}\n")
    for i, lab in zip(selected, sel_labels):
        f.write(f"\n{lab}:\n{np.array2string(mats[i], precision=6)}\n")

# --- Part 2: D4 quarter-turn intertwiner ---
# Want X (4x4, constant) s.t. X * v(x) = v(x+90deg) for all x in QI,
# where v(x)=[srx,sxp,cxp,crx]^T.
# Sample at multiple x, solve for X via least squares.
print("Part 2: D4 quarter-turn intertwiner search.")
print("Want X*v(x) = v(x+90deg).")
xs=[10,20,30,40,50,60,70,80]
# Build: for each x, v(x) -> v(x+90). Stack as linear system for X.
# X is 4x4 (16 unknowns). Each x gives 4 equations (one per output component).
A_rows=[]; b_rows=[]
for x in xs:
    r1=math.radians(x); r2=math.radians(x+90)
    # v(x)
    s1=math.sin(r1); c1=math.cos(r1)
    srx1=(1+c1)/s1; sxp1=1/srx1; cxp1=(1+s1)/c1; crx1=1/cxp1
    v1=np.array([srx1,sxp1,cxp1,crx1])
    # v(x+90)
    s2=math.sin(r2); c2=math.cos(r2)
    # Handle signs: use absolute values as in definitions
    # srx=|csc|+cot, etc. At x+90, sin/cos swap with signs.
    # Compute directly with formulas, handling quadrants.
    # For x in (0,90), x+90 in (90,180) (QII).
    # In QII: sin>0, cos<0.
    # srx=|csc|+cot = 1/|sin| + cos/sin.
    # Let's just compute numerically with abs.
    srx2=abs(1/s2)+c2/s2 if abs(s2)>1e-9 else 0
    sxp2=1/srx2 if srx2!=0 else 0
    cxp2=abs(1/c2)+s2/c2 if abs(c2)>1e-9 else 0
    crx2=1/cxp2 if cxp2!=0 else 0
    v2=np.array([srx2,sxp2,cxp2,crx2])
    # X*v1 = v2. For each output i: sum_j X[i,j]*v1[j] = v2[i].
    for i in range(4):
        row=np.zeros(16)
        for j in range(4):
            row[i*4+j]=v1[j]
        A_rows.append(row); b_rows.append(v2[i])
A_mat=np.array(A_rows); b_vec=np.array(b_rows)
# Least squares
X_flat, residuals, rank_ls, sv = np.linalg.lstsq(A_mat, b_vec, rcond=None)
X_mat=X_flat.reshape(4,4)
print(f"Least-squares X (4x4):")
print(np.array2string(X_mat, precision=4, suppress_small=True))
# Check residual
pred=A_mat@X_flat
err=np.max(np.abs(pred-b_vec))
print(f"Max |X*v(x)-v(x+90)|: {err:.2e}")
if err<1e-6:
    print("SUCCESS: Exact intertwiner found!")
else:
    print(f"INCOMPLETE: No exact constant intertwiner (err={err:.2e}).")
    print("The D4 action may not be linear with constant coefficients.")
print("WROTE /root/e47_intertwiner.txt")
with open('/root/e47_intertwiner.txt','w') as f:
    f.write(f"X=\n{X_mat}\nerr={err:.2e}\n")
