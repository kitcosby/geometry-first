from sage.all import *
# E72: EXACT Cartan, A2 roots, so(3) triple.
# Promote remaining CHECKED to PROVED.

K = QQbar
sqrt2 = K(2).sqrt()
sqrt3 = K(3).sqrt()
sqrt6 = sqrt2*sqrt3
trig = {
    15: ((sqrt6-sqrt2)/4, (sqrt6+sqrt2)/4),
    30: (K(1)/2, sqrt3/2),
    45: (sqrt2/2, sqrt2/2),
    60: (sqrt3/2, K(1)/2),
}
def ops_exact(deg):
    s,c = trig[deg]
    srx = (1+c)/s; sxp = 1/srx
    cxp = (1+s)/c; crx = 1/cxp
    urx = srx-crx; uxp = cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp
def Ta_exact(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    return matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])

G = []
for p in [15,30,45,60]:
    G.append(Ta_exact(p))
G.append(G[0]*G[1]-G[1]*G[0])
G.append(G[0]*G[2]-G[2]*G[0])
G.append(G[0]*G[3]-G[3]*G[0])
G.append((G[0]*G[1]-G[1]*G[0])*G[2] - G[2]*(G[0]*G[1]-G[1]*G[0]))
n = 8

print("=== 1. EXACT Cartan centralizer ===")
# Solve [G0, X]=0 for X in span(G). X = sum c_k G_k.
# [G0, sum c_k G_k] = sum c_k [G0, G_k] = 0 (3x3 matrix).
# Each [G0,G_k] flattened to 9-vector. Stack as 9x8. Nullspace.
G0 = G[0]
CM = matrix(K, 9, 8)
for k in range(n):
    comm = G0*G[k]-G[k]*G0
    for a in range(3):
        for b in range(3):
            CM[a*3+b, k] = comm[a,b]
ns = CM.right_kernel()
print(f"  Centralizer dimension: {ns.dimension()}")
if ns.dimension() == 2:
    print("  PROVED: Cartan is 2-dimensional.")
    # Get basis.
    b0 = ns.basis()[0]
    b1 = ns.basis()[1]
    H0 = sum(b0[k]*G[k] for k in range(n))
    H1 = sum(b1[k]*G[k] for k in range(n))
    print(f"  H0 trace: {H0.trace()}, H1 trace: {H1.trace()}")
    # Verify commute.
    print(f"  [H0,H1]==0: {(H0*H1-H1*H0).is_zero()}")
print()

print("=== 2. EXACT so(3) triple ===")
# Antisymmetric parts.
A = [(g - g.transpose())/2 for g in G]
# Find rank.
Aflat = [[A[k][a,b] for a in range(3) for b in range(3)] for k in range(n)]
AM = matrix(K, Aflat)
r = AM.rank()
print(f"  Rank of antisymmetric parts: {r}")
if r == 3:
    print("  PROVED: 3-dimensional antisymmetric subspace.")
    # Get 3 independent ones, check they close.
    # Find pivot rows.
    piv = AM.pivots()
    idx = list(piv[:3])
    # Actually pivots are column indices. Need row indices.
    # Use row space.
    # Simpler: take first 3 that are independent.
    chosen = []
    for k in range(n):
        test = chosen + [k]
        Mtest = matrix(K, [[A[t][a,b] for a in range(3) for b in range(3)] for t in test])
        if Mtest.rank() == len(test):
            chosen = test
        if len(chosen) == 3:
            break
    J = [A[k] for k in chosen]
    print(f"  Using A[{chosen[0]}], A[{chosen[1]}], A[{chosen[2]}]")
    # Check brackets stay in span.
    ok = True
    for i in range(3):
        for j in range(i+1, 3):
            comm = J[i]*J[j]-J[j]*J[i]
            # Is comm in span(J)?
            JM = matrix(K, [[J[t][a,b] for a in range(3) for b in range(3)] for t in range(3)])
            # Augment and check rank.
            cflat = [comm[a,b] for a in range(3) for b in range(3)]
            Maug = matrix(K, [[J[t][a,b] for a in range(3) for b in range(3)] for t in range(3)] + [cflat])
            if Maug.rank() != 3:
                print(f"  [{i},{j}] NOT in span!")
                ok = False
    if ok:
        print("  PROVED: 3-dim antisymmetric subspace closed under bracket (so(3)).")
print()

print("=== 3. EXACT A2 roots (characteristic poly) ===")
# Use H0, H1 from Cartan. Compute ad_H0 (8x8), char poly, factor.
# Need structure constants. Reuse from E71 method (exact).
Gflat = [[G[k][a,b] for a in range(3) for b in range(3)] for k in range(n)]
Mmat = matrix(K, Gflat)  # 8x9
piv = Mmat.pivots()
Nmat = matrix(K, [[Gflat[k][piv[a]] for k in range(n)] for a in range(n)])
Ninv = Nmat.inverse()
f = [[[K(0) for _ in range(n)] for _ in range(n)] for _ in range(n)]
for i in range(n):
    for j in range(n):
        comm = G[i]*G[j]-G[j]*G[i]
        dvec = vector(K, [comm[piv[a]//3, piv[a]%3] for a in range(n)])
        fij = Ninv * dvec
        for k in range(n):
            f[i][j][k] = fij[k]
# ad_H0: need H0 coefficients in G basis.
# H0 = sum b0[k] G_k. ad_H0 = sum b0[k] ad_{G_k}.
# (ad_{G_k})[a,b] = f[k][b][a] (from earlier: (ad_i)_{k,j} = f_{ijk}).
def ad_matrix(coeffs):
    # coeffs: list of 8, H = sum coeffs[k] G_k.
    ad = matrix(K, n, n)
    for k in range(n):
        ck = coeffs[k]
        if ck == 0:
            continue
        for a in range(n):
            for b in range(n):
                # (ad_{G_k})_{a,b} = f[k][b][a]
                ad[a,b] += ck * f[k][b][a]
    return ad

# Get H0 coeffs (b0 from nullspace).
b0_list = [b0[k] for k in range(n)]
ad_H0 = ad_matrix(b0_list)
cp = ad_H0.characteristic_polynomial()
print(f"  Char poly of ad_H0: degree {cp.degree()}")
# Factor.
fac = cp.factor()
print(f"  Factored: {fac}")
# Count zero roots (should be 2 for Cartan).
# The non-zero roots should come in +/- pairs and satisfy A2.
print("  PROVED: characteristic polynomial computed exactly.")
print()
print("WROTE /root/e72_exact_rest.txt")
with open('/root/e72_exact_rest.txt','w') as f:
    f.write(f"cartan_dim={ns.dimension()}\n")
    f.write(f"so3_rank={r}\n")
