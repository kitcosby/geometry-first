from sage.all import *
# E64: Do other octants give new generators beyond our 8?
# T_a(x) with full |.| handling. Check if T_a(105), T_a(195), etc. in span.

import math
import numpy as np

def ops_full(deg):
    # Full definitions with |.|: srx=|csc x|+cot x, cxp=|sec x|+tan x.
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-12 or abs(c)<1e-12:
        return None
    srx=abs(1/s)+c/s; sxp=1/srx
    cxp=abs(1/c)+s/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta_full(deg):
    v=ops_full(deg)
    if v is None: return None
    srx,sxp,cxp,crx,urx,uxp=v
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)

# Our 8 G_i (QI only).
def ops_qi(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp
def Ta_qi(deg):
    srx,sxp,cxp,crx,urx,uxp=ops_qi(deg)
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)

G=[]
for p in [15,30,45,60]: G.append(Ta_qi(p))
G.append(Ta_qi(15)@Ta_qi(30)-Ta_qi(30)@Ta_qi(15))
G.append(Ta_qi(15)@Ta_qi(45)-Ta_qi(45)@Ta_qi(15))
G.append(Ta_qi(15)@Ta_qi(60)-Ta_qi(60)@Ta_qi(15))
A=Ta_qi(15); B=Ta_qi(30); C=Ta_qi(45)
G.append((A@B-B@A)@C - C@(A@B-B@A))
M=np.array([g.flatten() for g in G]).T
Mpinv=np.linalg.pinv(M)

print("Testing T_a(x) for x in other octants:")
print("x   | in span? | residual")
print("----|----------|---------")
new_gens=[]
for deg in [105, 120, 135, 150, 195, 210, 225, 240, 285, 300, 315, 330]:
    T=Ta_full(deg)
    if T is None:
        print(f"{deg} | undefined| -")
        continue
    coeff=Mpinv@T.flatten()
    proj=(M@coeff).reshape(3,3)
    res=np.max(np.abs(T-proj))
    in_span=res<1e-6
    print(f"{deg} | {str(in_span):8s} | {res:.2e}")
    if not in_span:
        new_gens.append((deg, T))
print()
if new_gens:
    print(f"Found {len(new_gens)} NEW generators outside our 8-dim span.")
    print("These would expand the algebra beyond sl(3,R).")
    # Check: do they generate su(2)? Look at antisymmetric parts.
    print()
    print("Antisymmetric parts of new generators (candidate su(2)):")
    for deg, T in new_gens[:3]:
        Asym=(T-T.T)/2
        print(f"  T({deg}) antisym:")
        print(np.array2string(Asym, precision=3, suppress_small=True))
    print()
    print("WROTE /root/e64_octants.txt")
    with open('/root/e64_octants.txt','w') as f:
        f.write(f"new_count={len(new_gens)}\n")
        for deg, T in new_gens:
            f.write(f"deg={deg}\n")
else:
    print("All other-octant T_a(x) are IN our 8-dim span.")
    print("No new generators. The algebra is closed across octants.")
    print()
    print("WROTE /root/e64_octants.txt")
    with open('/root/e64_octants.txt','w') as f:
        f.write("new_count=0\n")
