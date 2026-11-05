# blending_odd_solve2.py — Sage (toe), targeted imports only.
# Part 2 of the odd-candidate solve: the unique double intertwiner B0 of V+W,
# its involution C = B0/sqrt(c) (B0^2 = c I), full condition table for C,
# condition-(4) bracket/anticommutator supports, and spec step (c):
# the Lie closure of {C} u {E_0..E_7} in M(16, QQ), exact, with the
# saturation pass as certificate. Run on toe as a .sage copy (Sage init
# quirk: plain .py under `sage file.py` fails on rational_field import;
# interactive stdin mangles compound blocks).
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
Om = I16
for a in range(8):
    Om = Om * E[a]

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

def stack4(mats):
    M = mats[0]
    for m in mats[1:]:
        M = M.stack(m)
    return M

CV = stack4([R[mu] - L[4 + mu] for mu in range(4)]).change_ring(QQ)
CW = stack4([R[4 + mu] - L[mu] for mu in range(4)]).change_ring(QQ)
kerV = CV.right_kernel()
BV = kerV.basis_matrix().transpose()
kerVW = (CW * BV).right_kernel()
BVW = BV * kerVW.basis_matrix().transpose()
print("dim V+W recomputed:", kerVW.dimension())
x0 = BVW.column(0)
B0 = Z16
for i in range(256):
    if x0[i] != 0:
        B0 = B0 + x0[i] * MON[i]
print("B0 intertwines both ways:",
      all(B0 * E[mu] == E[4 + mu] * B0 for mu in range(4)),
      all(B0 * E[4 + mu] == E[mu] * B0 for mu in range(4)))
B0sq = B0 * B0
c = B0sq[0, 0]
print("B0^2 = c I, c =", c, " scalar check:", B0sq == c * I16)
cn = ZZ(c.numerator())
cd = ZZ(c.denominator())
print("c positive:", c > 0, " num/den perfect squares:", cn.is_square(), cd.is_square())
s = QQ(cn.sqrt()) / QQ(cd.sqrt())
C = B0 / s
print("normalizer s =", s)
print("=== condition table for C = B0/s ===")
print("(1) C^2 = I:", C * C == I16)
tab = []
for a in range(8):
    tgt = E[a + 4] if a < 4 else E[a - 4]
    tab.append(C * E[a] * C == tgt)
print("(2) Ad_C(E_a) = sector partner for all a:", tab)
print("(3) C Om C^-1 = +Om:", C * Om * C == Om, " ; = -Om:", C * Om * C == -Om)
print("    [C, Om] = 0:", C * Om - Om * C == Z16, " ; {C, Om} = 0:", C * Om + Om * C == Z16)
print("    tr C =", C.trace())

# coefficient extraction: PHI columns = flattened monomials; x = PHI^-1 flat(X)
PHI = matrix(QQ, 256, 256)
for i, M in enumerate(MON):
    fl = M.list()
    for r in range(256):
        if fl[r] != 0:
            PHI[r, i] = fl[r]
PHIinv = PHI.inverse()

def coeffs(X):
    fl = matrix(QQ, 256, 1, X.list())
    v = PHIinv * fl
    return [(subs[i], str(v[i, 0])) for i in range(256) if v[i, 0] != 0]

print("=== (4) brackets / anticommutators of C with sector generators (monomial supports) ===")
for a in range(8):
    br = C * E[a] - E[a] * C
    ac = C * E[a] + E[a] * C
    print("a =", a, " [C,E_a] terms:", coeffs(br))
    print("a =", a, " {C,E_a} terms:", coeffs(ac))
print("{C,C} = 2I:", C * C + C * C == 2 * I16)

print("=== Lie closure in M(16, QQ) ===")
from sage.modules.free_module_element import vector

def flat(M):
    return vector(QQ, M.list())

def unflat(v):
    return matrix(QQ, 16, 16, list(v))

class Span:
    def __init__(self):
        self.piv = []   # list of pivot columns
        self.rows = []  # RREF vectors
        self.mats = []  # parallel matrix forms
    def reduce(self, v):
        v = vector(QQ, list(v))
        for pc, row in zip(self.piv, self.rows):
            if v[pc] != 0:
                v = v - v[pc] * row
        return v
    def try_add(self, M):
        v = self.reduce(flat(M))
        if v.is_zero():
            return False
        pc = min(i for i in range(256) if v[i] != 0)
        v = v / v[pc]
        for k in range(len(self.rows)):
            if self.rows[k][pc] != 0:
                self.rows[k] = self.rows[k] - self.rows[k][pc] * v
        ins = 0
        while ins < len(self.piv) and self.piv[ins] < pc:
            ins += 1
        self.piv.insert(ins, pc)
        self.rows.insert(ins, v)
        self.mats.insert(ins, unflat(v))
        return True
    def dim(self):
        return len(self.piv)

def lie_closure(gens, label):
    sp = Span()
    for g in gens:
        sp.try_add(g)
    changed = True
    passes = 0
    while changed:
        changed = False
        passes += 1
        cur = list(sp.mats)
        for b in cur:
            for g in gens:
                h = b * g - g * b
                if h != Z16 and sp.try_add(h):
                    changed = True
    # certificate: one explicit full pass, all brackets reduce to zero
    cert = 0
    ok = True
    for b in list(sp.mats):
        for g in gens:
            h = b * g - g * b
            cert += 1
            if not sp.reduce(flat(h)).is_zero():
                ok = False
    print(label, ": dim =", sp.dim(), " passes =", passes,
          " certificate brackets checked =", cert, " all in span:", ok)
    return sp

sp0 = lie_closure(E, "closure of gammas alone")
sp1 = lie_closure(E + [C], "closure of gammas + C")

# parity split of the big closure (coefficient space)
def parity_of_span(sp, label):
    k = sp.dim()
    BC = matrix(QQ, 256, k)
    for j, M in enumerate(sp.mats):
        v = PHIinv * matrix(QQ, 256, 1, M.list())
        for i in range(256):
            if v[i, 0] != 0:
                BC[i, j] = v[i, 0]
    even_idx = [i for i, Ss in enumerate(subs) if len(Ss) % 2 == 0]
    odd_idx = [i for i, Ss in enumerate(subs) if len(Ss) % 2 == 1]
    SELo = matrix(QQ, len(even_idx), 256)
    for r, i in enumerate(even_idx):
        SELo[r, i] = 1
    SELe = matrix(QQ, len(odd_idx), 256)
    for r, i in enumerate(odd_idx):
        SELe[r, i] = 1
    do = (SELo * BC).right_kernel().dimension()
    de = (SELe * BC).right_kernel().dimension()
    print(label, ": odd part dim =", do, " even part dim =", de, " total =", sp.dim())

parity_of_span(sp1, "closure(gammas + C)")
print("DONE-PART2")
