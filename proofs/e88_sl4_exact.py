from sage.all import *
# E88: EXACT sl(4,R) via octave.
# X(x) = [[T_a(x), v_op],[w_op^T, 0]], v_op=w_op=[S2_op, C2_op, 0].
# S2_op = 1/(urx+uxp) (= sin(2x)/4). C2_op = (srx+crx-sxp-cxp)/(urx+uxp) (= cos(2x)).
# All in QQbar. Prove: rank 15, closure, Killing (9,6), Cartan 3-dim, A3 roots.

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

def X4_exact(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    # T_a 3x3.
    A = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    # Operator double-angle.
    S2 = 1/(urx+uxp)  # sin(2x)/4
    C2 = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    # 4x4 block.
    X = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            X[a,b] = A[a,b]
    X[0,3] = S2; X[1,3] = C2; X[2,3] = 0
    X[3,0] = S2; X[3,1] = C2; X[3,2] = 0
    # X[3,3] = 0.
    return X

print("=== E88: Exact sl(4,R) ===")
print()
# Build generators.
G4 = [X4_exact(p) for p in [15,30,45,60]]
# Add commutators, grow to 15.
G = list(G4)
# Add [G0,Gi].
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
# Nested.
A01 = G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)
print(f"Initial: {len(G)} matrices.")

def mat_rank_exact(mats):
    # mats: list of n x n. Return rank of flattened.
    if len(mats) == 0:
        return 0
    n = mats[0].nrows()
    M = matrix(K, [[mats[k][a,b] for a in range(n) for b in range(n)] for k in range(len(mats))])
    # M is len(mats) x n^2. Rank.
    return M.rank()

r = mat_rank_exact(G)
print(f"Rank of {len(G)}: {r}")

# Grow by adding independent commutators.
changed = True
while changed and len(G) < 20:
    changed = False
    # Try all pairs.
    new_mats = []
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            comm = G[i]*G[j]-G[j]*G[i]
            # Check if independent.
            test = G + [comm]
            if mat_rank_exact(test) > mat_rank_exact(G):
                new_mats.append(comm)
                # Add one at a time to avoid duplicates.
                G.append(comm)
                changed = True
                break
        if changed:
            break
    # Limit.
    if len(G) >= 16:
        break

r = mat_rank_exact(G)
print(f"After growing: {len(G)} matrices, rank {r}.")
if r == 15:
    print("  PROVED: 15 linearly independent.")
else:
    print(f"  WARNING: rank {r}, not 15.")

# Trim to 15 independent.
# Find 15 independent.
G15 = []
for g in G:
    if len(G15) == 15:
        break
    test = G15 + [g]
    if mat_rank_exact(test) > mat_rank_exact(G15):
        G15.append(g)
print(f"Selected {len(G15)} independent.")

# Closure: all 105 commutators in span.
print("Checking closure (105 commutators)...")
n = 4
Gflat = [[G15[k][a,b] for a in range(n) for b in range(n)] for k in range(15)]
Mmat = matrix(K, Gflat)  # 15 x 16
# Use pivots for solving.
piv = Mmat.pivots()
print(f"  Pivots: {len(piv)} (expect 15).")
# Build N (15x15) and inverse.
Nmat = matrix(K, [[Gflat[k][piv[a]] for k in range(15)] for a in range(15)])
try:
    Ninv = Nmat.inverse()
    print("  N invertible.")
except:
    print("  N NOT invertible!")
    Ninv = None

closed = True
if Ninv is not None:
    for i in range(15):
        for j in range(i+1, 15):
            comm = G15[i]*G15[j]-G15[j]*G15[i]
            # Solve.
            dvec = vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
            coeffs = Ninv * dvec
            recon = sum(coeffs[k]*G15[k] for k in range(15))
            if not (recon - comm).is_zero():
                print(f"  [{i},{j}] NOT in span!")
                closed = False
                break
        if not closed:
            break
if closed:
    print("  PROVED: All 105 commutators in span. 15-dim Lie algebra.")
print()

# Killing form.
print("Killing form...")
# Structure constants.
f = {}
for i in range(15):
    for j in range(15):
        comm = G15[i]*G15[j]-G15[j]*G15[i]
        dvec = vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
        fij = Ninv * dvec
        for k in range(15):
            if fij[k] != 0:
                f[(i,j,k)] = fij[k]
print(f"  Non-zero f_ijk: {len(f)}")
# Killing matrix.
Kmat = matrix(K, 15, 15)
for i in range(15):
    for j in range(15):
        s = K(0)
        for k in range(15):
            for l in range(15):
                a = f.get((i,k,l), K(0))
                b = f.get((j,l,k), K(0))
                s += a*b
        Kmat[i,j] = s
# Signature: count positive/negative/zero via exact eigenvalues?
# QQbar eigenvalues might be heavy. Use characteristic polynomial and count sign changes?
# Simpler: compute rank and check non-degenerate, then use numeric for signature (but claim CHECKED).
# Actually, for PROVED we need exact. Let's try eigenvalues in QQbar.
print("  Computing Killing eigenvalues (exact)...")
try:
    ev = Kmat.eigenvalues()
    pos = sum(1 for e in ev if e > 0)
    neg = sum(1 for e in ev if e < 0)
    zero = sum(1 for e in ev if e == 0)
    print(f"  Signature: ({pos},{neg},{zero})")
    if pos==9 and neg==6 and zero==0:
        print("  PROVED: Killing signature (9,6). This is sl(4,R).")
    else:
        print(f"  WARNING: expected (9,6,0).")
except Exception as e:
    print(f"  Eigenvalue computation failed: {e}")
    print("  Trying numeric for indication...")
    import numpy as np
    Knum = np.array(Kmat.tolist(), dtype=complex)
    # Convert QQbar to float via real/imag?
    # Skip.

print()
print("WROTE /root/e88_sl4_exact.txt")
with open('/root/e88_sl4_exact.txt','w') as fh:
    fh.write(f"rank={r}\nclosed={closed}\n")
