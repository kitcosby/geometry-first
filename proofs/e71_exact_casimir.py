from sage.all import *
# E71: EXACT quadratic Casimir and Killing signature.
# Uses QQbar throughout. C2 should be exactly (4/9)*I.

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

print("Computing exact structure constants...")
# [G_i,G_j] = sum_k f_{ijk} G_k. For each (i,j), solve N * f = d
# where N[a,k] = (G_k) at flat position piv[a], d[a] = comm at piv[a].
Gflat = [[G[k][a,b] for a in range(3) for b in range(3)] for k in range(n)]
Mmat = matrix(K, Gflat)  # 8x9, rows=G_k flattened
piv = Mmat.pivots()  # pivot cols = independent flat positions
print(f"  Pivot positions: {piv}")
# N[a,k] = Gflat[k][piv[a]], a=0..7, k=0..7. 8x8.
Nmat = matrix(K, [[Gflat[k][piv[a]] for k in range(n)] for a in range(n)])
Ninv = Nmat.inverse()
f = [[[K(0) for _ in range(n)] for _ in range(n)] for _ in range(n)]
for i in range(n):
    for j in range(n):
        comm = G[i]*G[j]-G[j]*G[i]
        dvec = vector(K, [comm[piv[a]//3, piv[a]%3] for a in range(n)])
        f_ij = Ninv * dvec
        for k in range(n):
            f[i][j][k] = f_ij[k]
# Verify: reconstruct one commutator fully (all 9 entries).
print("  Verifying f_ijk on [G0,G1] (all 9 entries)...")
recon = sum(f[0][1][k]*G[k] for k in range(n))
actual = G[0]*G[1]-G[1]*G[0]
verr = max(abs((recon-actual)[a,b]) for a in range(3) for b in range(3))
print(f"  Max reconstruction error: {verr} (should be 0)")
print("  Structure constants done.")

print("Computing exact Killing form...")
Kmat = matrix(K, n, n)
for i in range(n):
    for j in range(n):
        s = K(0)
        for k in range(n):
            for l in range(n):
                s += f[i][k][l]*f[j][l][k]
        Kmat[i,j] = s
print("  Killing matrix done.")

print("Exact Killing signature...")
# Eigenvalues in QQbar, count signs.
ev = Kmat.eigenvalues()
pos = sum(1 for e in ev if e > 0)
neg = sum(1 for e in ev if e < 0)
zero = sum(1 for e in ev if e == 0)
print(f"  Positive: {pos}, Negative: {neg}, Zero: {zero}")
if pos==5 and neg==3 and zero==0:
    print("  PROVED: signature (5,3).")
print()

print("Exact inverse Killing and quadratic Casimir...")
Kinv = Kmat.inverse()
C2 = matrix(K, 3, 3, [[K(0)]*3 for _ in range(3)])
for i in range(n):
    for j in range(n):
        C2 += Kinv[i,j] * (G[i]*G[j])
print("  C2:")
print(C2)
# Check if C2 == (4/9)*I.
target = (K(4)/K(9)) * matrix.identity(K, 3)
diff = C2 - target
maxdiff = max(abs(diff[a,b]) for a in range(3) for b in range(3))
print(f"  Max |C2 - (4/9)I|: {maxdiff}")
if maxdiff == 0:
    print("  PROVED: C2 = (4/9)*I exactly.")
else:
    # Check if it's scalar (maybe different value).
    # C2 should be c*I. Compute c = trace/3.
    c = C2.trace()/3
    dev = max(abs((C2 - c*matrix.identity(K,3))[a,b]) for a in range(3) for b in range(3))
    print(f"  Trace/3 = {c}, deviation from scalar: {dev}")
    if dev == 0:
        print(f"  PROVED: C2 = ({c})*I exactly (scalar, value {c}).")
print()
print("WROTE /root/e71_exact_casimir.txt")
with open('/root/e71_exact_casimir.txt','w') as f:
    f.write(f"sig_pos={pos}\nsig_neg={neg}\n")
    f.write(f"C2_exact_4_9={maxdiff==0}\n")
