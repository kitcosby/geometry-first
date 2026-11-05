# blending_odd_solve.py — Sage (toe), targeted imports only.
# Blending intertwiner hunt, steps (a)/(b) under the 2026-10-16 addendum:
# odd-candidate solve. Realization identical to proofs/d2_structure.py:
#   e = (sx x I, Y x I, sz x sx, sz x Y) in Cl(2,2); E_mu = e_mu x I4 (R sector),
#   E_{4+nu} = Om4 x e_nu (J sector), Om4 = e0 e1 e2 e3. Signature (+,-,+,-, +, -, +, -).
# Method: coefficient space over the 256 Clifford monomials E_S (S subset of 0..7).
# Left/right multiplication by a gamma is a signed permutation of coefficients.
#   V  = {X : X E_mu = E_{4+mu} X}      (intertwiner R->J, d2's V)
#   W  = {X : X E_{4+mu} = E_mu X}      (intertwiner J->R; with C^2=I, V+W = condition (2))
#   O  = {X : X Om + Om X = 0}          (odd part; Om = E_0...E_7)
# All arithmetic exact over QQ. Prints dims, parity splits, and the V+W basis.
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix, identity_matrix
from itertools import combinations

I2 = matrix(QQ, [[1, 0], [0, 1]])
Xm = matrix(QQ, [[0, 1], [1, 0]])
Z = matrix(QQ, [[1, 0], [0, -1]])
Y = matrix(QQ, [[0, -1], [1, 0]])
e = [Xm.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(Xm), Z.tensor_product(Y)]
I4 = identity_matrix(QQ, 4)
Om4 = e[0] * e[1] * e[2] * e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Om4.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ, 16)
Z16 = matrix(QQ, 16, 16)

print("=== gates ===")
sig = [1, -1, 1, -1, 1, -1, 1, -1]
print("gamma squares ok:", all(E[a] * E[a] == sig[a] * I16 for a in range(8)))
print("gamma anticommute ok:", all(E[a] * E[b] + E[b] * E[a] == Z16
                                  for a in range(8) for b in range(8) if a != b))
Om = I16
for a in range(8):
    Om = Om * E[a]
print("Om^2 = I:", Om * Om == I16)
print("Om anticommutes with each gamma:", all(Om * E[a] + E[a] * Om == Z16 for a in range(8)))
# E6 gate (record reproduction, cf. blending_step_d_screen.md)
S = [[E[a] * E[b] for b in range(8)] for a in range(8)]
def comm(A, B):
    return A * B - B * A
print("E6 gate [S04,S15]=[S05,S14]=0:", comm(S[0][4], S[1][5]) == Z16 and comm(S[0][5], S[1][4]) == Z16)
print("E6 gate other mixed pairs nonzero:", all(comm(S[a][b], S[c][d]) != Z16
      for (a, b, c, d) in [(0, 4, 0, 5), (0, 4, 1, 4), (1, 5, 0, 5), (1, 5, 1, 4)]))
T = S[0][4] * S[1][5] - S[0][5] * S[1][4]
print("E6 gate T^2=4I, tr T=0:", T * T == 4 * I16 and T.trace() == 0)
print("E6 gate T = -2 E0E1E4E5:", T == -2 * E[0] * E[1] * E[4] * E[5])

print("=== monomial coefficient machinery ===")
subs = []
for k in range(9):
    subs += [tuple(c) for c in combinations(range(8), k)]
MON = []
for Ss in subs:
    M = I16
    for a in Ss:
        M = M * E[a]
    MON.append(M)
key2 = {}
for i, M in enumerate(MON):
    key2[tuple(M.list())] = (i, 1)
    key2[tuple((-M).list())] = (i, -1)

def mult_maps(A):
    # returns (Lft, Rgt): 256x256 ZZ matrices, columns = coords of A*M_i, M_i*A
    Lft = matrix(ZZ, 256, 256)
    Rgt = matrix(ZZ, 256, 256)
    for i, M in enumerate(MON):
        j, s = key2[tuple((A * M).list())]
        Lft[j, i] = s
        j2, s2 = key2[tuple((M * A).list())]
        Rgt[j2, i] = s2
    return Lft, Rgt

L = []
R = []
for a in range(8):
    la, ra = mult_maps(E[a])
    L.append(la)
    R.append(ra)
# parity gate: monomial anticommutes with Om iff odd grade; commutes iff even grade
par_ok = True
for i, Ss in enumerate(subs):
    anti = (MON[i] * Om + Om * MON[i] == Z16)
    commutes = (MON[i] * Om - Om * MON[i] == Z16)
    if len(Ss) % 2 == 1 and not anti:
        par_ok = False
    if len(Ss) % 2 == 0 and not commutes:
        par_ok = False
print("parity gate (odd monomials anticommute with Om, even commute):", par_ok)
print("dim O (odd coordinate subspace):", sum(1 for Ss in subs if len(Ss) % 2 == 1))

print("=== linear solves over QQ ===")
def stack4(mats):
    M = mats[0]
    for m in mats[1:]:
        M = M.stack(m)
    return M

CV = stack4([R[mu] - L[4 + mu] for mu in range(4)]).change_ring(QQ)
CW = stack4([R[4 + mu] - L[mu] for mu in range(4)]).change_ring(QQ)
kerV = CV.right_kernel()
kerW = CW.right_kernel()
BV = kerV.basis_matrix().transpose()   # 256 x dimV, columns = coeff vectors
BW = kerW.basis_matrix().transpose()
print("dim V =", kerV.dimension(), " dim W =", kerW.dimension())

# V+W via restriction of W-conditions to V basis (small system), cross-checked below
kerVW = (CW * BV).right_kernel()
BVW = BV * kerVW.basis_matrix().transpose()
dVW = kerVW.dimension()
kerVW2 = (CV * BW).right_kernel()
print("dim V+W (via V) =", dVW, " dim V+W (via W) =", kerVW2.dimension())

even_idx = [i for i, Ss in enumerate(subs) if len(Ss) % 2 == 0]
odd_idx = [i for i, Ss in enumerate(subs) if len(Ss) % 2 == 1]

def parity_split(B, name):
    # B: 256 x d coefficient basis. Returns (dim odd part, dim even part).
    d = B.ncols()
    if d == 0:
        print(name, ": empty space")
        return (0, 0)
    SELo = matrix(QQ, len(even_idx), 256)
    for r, i in enumerate(even_idx):
        SELo[r, i] = 1
    SELe = matrix(QQ, len(odd_idx), 256)
    for r, i in enumerate(odd_idx):
        SELe[r, i] = 1
    do = (SELo * B).right_kernel().dimension()   # coeffs with zero even part = odd part
    de = (SELe * B).right_kernel().dimension()   # coeffs with zero odd part = even part
    print(name, ": dim odd part =", do, " dim even part =", de, " (sum vs dim:", do + de, "vs", d, ")")
    return (do, de)

parity_split(BV, "V")
parity_split(BW, "W")
parity_split(BVW, "V+W")

# direct verification: each V+W basis element (as a matrix) satisfies both intertwining systems
def coeffs_to_mat(x):
    M = Z16
    for i in range(256):
        if x[i] != 0:
            M = M + x[i] * MON[i]
    return M

print("=== V+W basis (sparse coefficients) ===")
for j in range(dVW):
    x = BVW.column(j)
    X = coeffs_to_mat(x)
    okV = all(X * E[mu] == E[4 + mu] * X for mu in range(4))
    okW = all(X * E[4 + mu] == E[mu] * X for mu in range(4))
    terms = [(subs[i], str(x[i])) for i in range(256) if x[i] != 0]
    print("basis", j, "intertwines V:", okV, " W:", okW, " nterms:", len(terms))
    print("  ", terms)
print("DONE-PART1")
