from sage.all import *
# E68: Axiom 0 -- imaginary sector via x -> i*y.
# srx(iy) = |csc(iy)| + cot(iy). Compute T_a(iy), get real/imag parts,
# see what real Lie algebra they generate.

import math
import cmath
import numpy as np

def ops_complex(z):
    # z complex. srx=|csc z|+cot z, etc. |.| = modulus.
    s=cmath.sin(z); c=cmath.cos(z)
    if abs(s)<1e-12 or abs(c)<1e-12:
        return None
    srx=abs(1/s)+(c/s); sxp=1/srx
    cxp=abs(1/c)+(s/c); crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta_complex(z):
    v=ops_complex(z)
    if v is None: return None
    srx,sxp,cxp,crx,urx,uxp=v
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=complex)

# Real sector: our 8 G_i (3x3 real).
def ops_real(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp
def Ta_real(deg):
    srx,sxp,cxp,crx,urx,uxp=ops_real(deg)
    return np.array([[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]], dtype=float)

G_real=[]
for p in [15,30,45,60]: G_real.append(Ta_real(p))
G_real.append(Ta_real(15)@Ta_real(30)-Ta_real(30)@Ta_real(15))
G_real.append(Ta_real(15)@Ta_real(45)-Ta_real(45)@Ta_real(15))
G_real.append(Ta_real(15)@Ta_real(60)-Ta_real(60)@Ta_real(15))
A=Ta_real(15); B=Ta_real(30); C=Ta_real(45)
G_real.append((A@B-B@A)@C - C@(A@B-B@A))

# Imaginary sector: T_a(i*y) for y in degrees->radians.
print("T_a(i*y) for y (complex 3x3):")
for y_deg in [15, 30, 45]:
    y=math.radians(y_deg)
    T=Ta_complex(1j*y)
    print(f"  y={y_deg}deg: T_a(i*y) =")
    print(np.array2string(T, precision=3, suppress_small=True))
    print(f"  Trace: {np.trace(T):.3f}")
    print()

# Get real and imaginary parts as REAL 3x3 matrices.
# For a complex 3x3 M, Re(M) and Im(M) are real 3x3.
print("Real and imaginary parts (as real 3x3):")
Mats=[]
for y_deg in [15,30,45,60]:
    y=math.radians(y_deg)
    T=Ta_complex(1j*y)
    Re=T.real; Im=T.imag
    # Check traceless.
    print(f"  y={y_deg}: Re trace={np.trace(Re):.2e}, Im trace={np.trace(Im):.2e}")
    Mats.append(Re); Mats.append(Im)
print()
# What real Lie algebra do {Re, Im} generate?
# Combine with G_real. Compute span of all, plus commutators.
print("Combined real+imaginary span:")
All=G_real+Mats  # 8 + 8 = 16 matrices (3x3 real)
Mflat=np.array([m.flatten() for m in All]).T  # 9x16
rank=np.linalg.matrix_rank(Mflat, tol=1e-8)
print(f"  Rank of 16 matrices: {rank}")
# Add commutators to see if it closes or grows.
# Take commutators of Re/Im with each other and with G_real.
print("  Adding commutators...")
extended=list(All)
for i in range(len(All)):
    for j in range(i+1, len(All)):
        comm=All[i]@All[j]-All[j]@All[i]
        # Check if in span.
        M2=np.array([m.flatten() for m in extended]).T
        coeff=np.linalg.lstsq(M2, comm.flatten(), rcond=None)[0]
        proj=(M2@coeff).reshape(3,3)
        if np.max(np.abs(comm-proj))>1e-6:
            extended.append(comm)
            if len(extended)>30: break
    if len(extended)>30: break
Mext=np.array([m.flatten() for m in extended]).T
rank_ext=np.linalg.matrix_rank(Mext, tol=1e-8)
print(f"  After commutators: {len(extended)} matrices, rank {rank_ext}")
print(f"  (sl(3,C) as real is 16-dim; sl(3,R)+sl(3,R) is 16-dim)")
print()
print("WROTE /root/e68_axiom0.txt")
with open('/root/e68_axiom0.txt','w') as f:
    f.write(f"rank16={rank}\nrank_ext={rank_ext}\ncount_ext={len(extended)}\n")
