from sage.all import *
# E70 pilot: EXACT symbolic matrices and rank.
# Angles 15,30,45,60 have algebraic sin/cos in QQbar.

K = QQbar
sqrt2 = K(2).sqrt()
sqrt3 = K(3).sqrt()
sqrt6 = sqrt2*sqrt3

# Exact sin/cos.
# 15deg: sin=(sqrt6-sqrt2)/4, cos=(sqrt6+sqrt2)/4
# 30deg: sin=1/2, cos=sqrt3/2
# 45deg: sin=sqrt2/2, cos=sqrt2/2
# 60deg: sin=sqrt3/2, cos=1/2
trig = {
    15: ((sqrt6-sqrt2)/4, (sqrt6+sqrt2)/4),
    30: (K(1)/2, sqrt3/2),
    45: (sqrt2/2, sqrt2/2),
    60: (sqrt3/2, K(1)/2),
}

def ops_exact(deg):
    s,c = trig[deg]
    # QI: s>0, c>0, so |csc|=1/s, |sec|=1/c.
    srx = (1+c)/s
    sxp = 1/srx
    cxp = (1+s)/c
    crx = 1/cxp
    urx = srx-crx
    uxp = cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta_exact(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    M = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    return M

print("Exact T_a matrices:")
for p in [15,30,45,60]:
    T = Ta_exact(p)
    print(f"  T_a({p}): trace={T.trace()}, det={T.det()}")
    # Verify traceless exactly.
    assert T.trace() == 0, f"Not traceless at {p}!"
print("  All traceless (EXACT).")
print()

# Build 8 generators.
G = []
for p in [15,30,45,60]:
    G.append(Ta_exact(p))
G.append(G[0]*G[1]-G[1]*G[0])
G.append(G[0]*G[2]-G[2]*G[0])
G.append(G[0]*G[3]-G[3]*G[0])
G.append((G[0]*G[1]-G[1]*G[0])*G[2] - G[2]*(G[0]*G[1]-G[1]*G[0]))

# Exact rank: stack vectorized as rows, compute rank over QQbar.
Mrows = []
for g in G:
    Mrows.append([g[i,j] for i in range(3) for j in range(3)])
Mmat = matrix(K, Mrows)  # 8x9
r = Mmat.rank()
print(f"Exact rank of 8 generators: {r}")
if r == 8:
    print("PROVED: 8 linearly independent traceless 3x3 matrices => span = sl(3,R).")
else:
    print(f"Rank {r} != 8. NOT full.")
print()

# Exact quadratic Casimir (need structure constants; do via linear algebra).
# For now, just verify [G_i,G_j] in span (exact).
print("Exact Lie closure check (all 28 commutators in span):")
# Build span basis: use Mmat row space. For each commutator, check if in row space.
# Solve Mmat^T * c = comm_flat (least squares won't work exactly; use linear solve on independent rows).
# Get 8 independent rows (all of them, since rank 8).
# For comm in span iff rank doesn't increase.
ok = True
for i in range(8):
    for j in range(i+1, 8):
        comm = G[i]*G[j]-G[j]*G[i]
        cflat = [comm[a,b] for a in range(3) for b in range(3)]
        # Augment matrix and check rank.
        Maug = matrix(K, Mrows + [cflat])
        if Maug.rank() != 8:
            print(f"  [{i},{j}] NOT in span! rank {Maug.rank()}")
            ok = False
if ok:
    print("  All 28 commutators in span (EXACT). Lie algebra PROVED.")
print()
print("WROTE /root/e70_exact.txt")
with open('/root/e70_exact.txt','w') as f:
    f.write(f"rank={r}\nclosure={ok}\n")
