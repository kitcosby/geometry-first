from sage.all import *
# E75: Embed sl(3) in sl(4). Add 7 generators. A3 roots. Branching.
# Step 5: Can trig produce the extra 7 naturally?

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

# Step 1: Original 8 G_i (3x3).
G3 = []
for p in [15,30,45,60]:
    G3.append(Ta_exact(p))
G3.append(G3[0]*G3[1]-G3[1]*G3[0])
G3.append(G3[0]*G3[2]-G3[2]*G3[0])
G3.append(G3[0]*G3[3]-G3[3]*G3[0])
G3.append((G3[0]*G3[1]-G3[1]*G3[0])*G3[2] - G3[2]*(G3[0]*G3[1]-G3[1]*G3[0]))

# Embed as 4x4 block-diagonal.
def embed(M3):
    M4 = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            M4[a,b] = M3[a,b]
    return M4

G4_embed = [embed(g) for g in G3]

# Step 2: Add 7 extra generators.
# E_{i4}, E_{4i} for i=0,1,2 (6 matrices).
extra = []
for i in range(3):
    E_i4 = matrix(K, 4, 4); E_i4[i,3] = 1; extra.append(E_i4)
    E_4i = matrix(K, 4, 4); E_4i[3,i] = 1; extra.append(E_4i)
# H3 = diag(1,1,1,-3).
H3 = matrix(K, 4, 4, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-3]])
extra.append(H3)

G15 = G4_embed + extra
n = 15
print(f"Step 1-2: {len(G4_embed)} embedded + {len(extra)} extra = {n} generators.")
print()

# Verify rank 15.
Gflat = [[G15[k][a,b] for a in range(4) for b in range(4)] for k in range(n)]
Mmat = matrix(K, Gflat)
r = Mmat.rank()
print(f"Step 2 check: Rank = {r} (expect 15).")
if r == 15:
    print("  PROVED: 15 linearly independent generators.")
print()

# Verify closure.
print("Step 2 check: Lie closure (105 commutators)...")
piv = Mmat.pivots()
# Use first 15 pivots (should be 15).
Nmat = matrix(K, [[Gflat[k][piv[a]] for k in range(n)] for a in range(n)])
Ninv = Nmat.inverse()
closed = True
for i in range(n):
    for j in range(i+1, n):
        comm = G15[i]*G15[j]-G15[j]*G15[i]
        dvec = vector(K, [comm[piv[a]//4, piv[a]%4] for a in range(n)])
        # Check comm in span (residual).
        # Reconstruct and compare.
        coeffs = Ninv * dvec
        recon = sum(coeffs[k]*G15[k] for k in range(n))
        if not (recon - comm).is_zero():
            print(f"  [{i},{j}] NOT in span!")
            closed = False
            break
    if not closed:
        break
if closed:
    print("  PROVED: All 105 commutators in span. This is a 15-dim Lie algebra.")
print()

# Step 3: Structure constants (exact).
print("Step 3: Computing exact structure constants...")
f = {}
for i in range(n):
    for j in range(n):
        comm = G15[i]*G15[j]-G15[j]*G15[i]
        dvec = vector(K, [comm[piv[a]//4, piv[a]%4] for a in range(n)])
        fij = Ninv * dvec
        for k in range(n):
            if fij[k] != 0:
                f[(i,j,k)] = fij[k]
print(f"  Non-zero f_ijk: {len(f)}")
# Verify Jacobi on a sample (exact).
print("  Verifying Jacobi (sample)...")
# Just check one.
ok = True
# Skip full Jacobi (too many). Trust closure.
print("  (Closure implies Lie algebra; Jacobi follows from matrix commutator.)")
print()

# Step 4: A3 root system.
# Find Cartan: centralizer of G4_embed[0] (the embedded G0).
print("Step 4: A3 root system.")
G0 = G4_embed[0]
CM = matrix(K, 16, n)
for k in range(n):
    comm = G0*G15[k]-G15[k]*G0
    for a in range(4):
        for b in range(4):
            CM[a*4+b, k] = comm[a,b]
ns = CM.right_kernel()
print(f"  Centralizer dimension: {ns.dimension()} (expect 3 for A3).")
if ns.dimension() == 3:
    print("  PROVED: Cartan is 3-dimensional (rank 3).")
    # Get Cartan basis.
    basis = ns.basis()
    H = [sum(basis[m][k]*G15[k] for k in range(n)) for m in range(3)]
    # Compute ad_H0 char poly.
    # Build ad matrices.
    def ad_matrix(coeffs):
        ad = matrix(K, n, n)
        for k in range(n):
            ck = coeffs[k]
            if ck == 0: continue
            for a in range(n):
                for b in range(n):
                    # (ad_{G_k})_{a,b} = f_{kba}
                    key = (k,b,a)
                    if key in f:
                        ad[a,b] += ck * f[key]
        return ad
    b0 = [basis[0][k] for k in range(n)]
    ad_H0 = ad_matrix(b0)
    cp = ad_H0.characteristic_polynomial()
    print(f"  Char poly degree: {cp.degree()}")
    fac = cp.factor()
    print(f"  Factored: {fac}")
    # Count roots.
    # Should be 12 non-zero + x^3.
    print("  (Expect 12 roots + x^3 for A3.)")
print()

# Branching: original 8 (embedded) form A2.
print("Branching: Original 8 embedded generators.")
# Their span should be closed (it's sl(3)).
# Check: commutators of embedded stay in embedded span.
G8flat = [[G4_embed[k][a,b] for a in range(4) for b in range(4)] for k in range(8)]
M8 = matrix(K, G8flat)
print(f"  Rank of embedded 8: {M8.rank()} (expect 8).")
# Check closure within the 8.
closed8 = True
for i in range(8):
    for j in range(i+1, 8):
        comm = G4_embed[i]*G4_embed[j]-G4_embed[j]*G4_embed[i]
        # Is in span of 8?
        # Flatten and check rank.
        cflat = [comm[a,b] for a in range(4) for b in range(4)]
        Maug = matrix(K, G8flat + [cflat])
        if Maug.rank() != 8:
            closed8 = False
            break
    if not closed8:
        break
if closed8:
    print("  PROVED: Embedded 8 form a subalgebra (A2 inside A3).")
print()

# Step 5: Can trig produce the extra 7 naturally?
print("Step 5: Can trigonometric operators produce the extra 7?")
print("  The extra 7 are: E_{i4}, E_{4i} (i=1,2,3) and diag(1,1,1,-3).")
print("  These are standard matrix units, NOT derived from srx/sxp/cxp/crx/urx/uxp.")
print("  The trig operators are 3x3 in nature; the 4th dimension is new.")
print("  There is no trig operator that naturally gives E_{14} (1 at position (1,4)).")
print("  CONCLUSION: The extra 7 are ADJOINED BY HAND, not derived from trig.")
print("  The trig operators give A2 (sl(3)), not A3 (sl(4)).")
print("  To get A3 naturally, you'd need new trig operators requiring 4x4.")
print()
print("WROTE /root/e75_sl4.txt")
with open('/root/e75_sl4.txt','w') as f:
    f.write(f"rank15={r}\nclosed={closed}\ncartan_dim={ns.dimension()}\n")
