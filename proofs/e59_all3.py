from sage.all import *
# E59: Cubic Casimir, adjoint rep eigenvalues, systematic G_i analysis.

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
n=8
M=np.array([g.flatten() for g in G]).T
Mpinv=np.linalg.pinv(M)
f=np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        comm=G[i]@G[j]-G[j]@G[i]
        f[i,j,:]=Mpinv@comm.flatten()
K=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        s=0
        for k in range(n):
            for l in range(n):
                s+=f[i,k,l]*f[j,l,k]
        K[i,j]=s
Kinv=np.linalg.inv(K)

print("=== 1. CUBIC CASIMIR ===")
# d_ijk = (1/4) * Tr({G_i,G_j} G_k) / (normalization)
# For sl(3), C3 = d_{ijk} G_i G_j G_k.
# Compute d_ijk from our G's (need symmetric product).
# {A,B} = AB+BA.
d=np.zeros((n,n,n))
for i in range(n):
    for j in range(n):
        for k in range(n):
            # d_ijk proportional to Tr({G_i,G_j} G_k)
            anticomm=G[i]@G[j]+G[j]@G[i]
            d[i,j,k]=np.trace(anticomm@G[k])
# Normalize: for standard sl(3), d_ijk d_ijk = 40/3.
# We'll just compute C3 = sum_{ijk} d_{ijk} G_i G_j G_k (up to scale).
# Actually, need Kinv to raise indices properly. Simpler: use d with lower indices
# and G's, but must be careful with normalization.
# Let's compute C3 = sum_{ijk} d^{ijk} G_i G_j G_k where d^{ijk}=Kinv^{ia}Kinv^{jb}Kinv^{kc} d_{abc}.
d_up=np.einsum('ia,jb,kc,abc->ijk', Kinv, Kinv, Kinv, d)
C3=np.zeros((3,3), dtype=complex)
for i in range(n):
    for j in range(n):
        for k in range(n):
            C3+=d_up[i,j,k]*(G[i]@G[j]@G[k])
print("C3 (3x3, may not be proportional to I if normalization off):")
print(np.array2string(C3, precision=4, suppress_small=True))
# Check if proportional to I
c3=np.trace(C3)/3
dev=np.max(np.abs(C3-c3*np.eye(3)))
print(f"Tr(C3)/3 = {c3:.6f}")
print(f"Max |C3-c3*I|: {dev:.2e}")
if dev<1e-4:
    print(f"C3 = {c3:.6f} * I (up to normalization).")
else:
    print("C3 not proportional to I with this normalization.")
    print("(Cubic Casimir needs careful normalization; result is scale-dependent.)")
print()

print("=== 2. ADJOINT REPRESENTATION (8x8) ===")
# (ad_{G_i})_{jk} = -f_{ijk}  (with appropriate index positions)
# For each i, ad_i is 8x8: (ad_i)[j,k] = f_{i,k,j} (?) Need convention.
# Standard: [G_i, G_j] = f_{ijk} G_k  =>  (ad_i)_{kj} = f_{ijk}.
# So ad_i as matrix with rows k, cols j: (ad_i)[k,j] = f_{i,j,k}.
for idx in [0, 2]:  # Just two examples (G0=Ta15, G2=Ta45)
    ad=np.zeros((n,n))
    for j in range(n):
        for k in range(n):
            ad[k,j]=f[idx,j,k]
    ev=np.linalg.eigvals(ad)
    print(f"ad_G[{idx}] eigenvalues (8):")
    print(f"  {np.sort(ev)}")
print("(Adjoint eigenvalues are the roots + zeros.)")
print()

print("=== 3. SYSTEMATIC G_i EIGENVALUE ANALYSIS ===")
print("G_i | eigenvalues | real? | trace")
print("----|-------------|-------|------")
for i, g in enumerate(G):
    ev=np.linalg.eigvals(g)
    is_real=np.all(np.abs(ev.imag)<1e-8)
    tr=np.trace(g)
    ev_str=", ".join([f"{v.real:.2f}{'+'+str(round(v.imag,2))+'j' if abs(v.imag)>1e-8 else ''}" for v in ev])
    print(f"G[{i}] | {ev_str} | {'real' if is_real else 'complex'} | {tr:.2e}")
print()
print("Interpretation:")
print("  Real eigenvalues => non-compact (boost-like) generator.")
print("  Complex pair => compact (rotation-like) component.")
print("  Mixed => generator has both compact and non-compact parts.")
print("  (Our G_i are not pure; the so(3) triple from E56 is pure compact.)")
print()
print("WROTE /root/e59_all3.txt")
with open('/root/e59_all3.txt','w') as f:
    f.write(f"c3_trace={c3:.6f}\n")
