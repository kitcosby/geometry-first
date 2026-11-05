# blending_odd_check.py — desktop cross-check (exact integer arithmetic + mod-p ranks).
# Independent rebuild of the d2_structure.py realization; verifies the toe results:
#   dim V = 16, dim V+W = 1, V+W purely even (odd part 0), the double intertwiner B0,
#   B0^2 = c*I, and (if c is a rational square) the involution C = B0/sqrt(c):
#   C^2 = I, exact sector-swap conjugation table, C Om C^-1 = +Om (condition (3) fails).
import numpy as np
from fractions import Fraction
from itertools import combinations

I2 = np.array([[1, 0], [0, 1]], dtype=object)
Xm = np.array([[0, 1], [1, 0]], dtype=object)
Z = np.array([[1, 0], [0, -1]], dtype=object)
Y = np.array([[0, -1], [1, 0]], dtype=object)
kron = np.kron
e = [kron(Xm, I2), kron(Y, I2), kron(Z, Xm), kron(Z, Y)]
I4 = np.eye(4, dtype=object)
Om4 = e[0] @ e[1] @ e[2] @ e[3]
E = [kron(e[mu], I4) for mu in range(4)] + [kron(Om4, e[nu]) for nu in range(4)]
I16 = np.eye(16, dtype=object)
Z16 = np.zeros((16, 16), dtype=object)
Om = I16.copy()
for a in range(8):
    Om = Om @ E[a]

subs = []
for k in range(9):
    subs += [tuple(c) for c in combinations(range(8), k)]
MON = []
for Ss in subs:
    M = I16.copy()
    for a in Ss:
        M = M @ E[a]
    MON.append(M)

# B0 coefficients transcribed from the toe Sage run (blending_odd_solve.py output).
B0_TERMS = [((0, 1, 2, 3), 1), ((0, 1, 2, 7), -1), ((0, 1, 3, 6), 1), ((0, 1, 6, 7), 1),
            ((0, 2, 3, 5), -1), ((0, 2, 5, 7), -1), ((0, 3, 5, 6), 1), ((0, 5, 6, 7), -1),
            ((1, 2, 3, 4), 1), ((1, 2, 4, 7), 1), ((1, 3, 4, 6), -1), ((1, 4, 6, 7), 1),
            ((2, 3, 4, 5), 1), ((2, 4, 5, 7), -1), ((3, 4, 5, 6), 1), ((4, 5, 6, 7), 1)]
idx = {s: i for i, s in enumerate(subs)}
B0 = Z16.copy()
for Ss, cf in B0_TERMS:
    B0 = B0 + cf * MON[idx[Ss]]

def eq(A, B):
    return np.array_equal(A, B)

print("B0 intertwines R->J (V):", all(eq(B0 @ E[mu], E[4 + mu] @ B0) for mu in range(4)))
print("B0 intertwines J->R (W):", all(eq(B0 @ E[4 + mu], E[mu] @ B0) for mu in range(4)))
B0sq = B0 @ B0
c = B0sq[0, 0]
print("B0^2 = c*I with c =", c, ":", eq(B0sq, c * I16))
print("tr B0 =", sum(B0[i, i] for i in range(16)))
print("[B0, Om] = 0 (B0 even):", eq(B0 @ Om, Om @ B0))
print("{B0, Om} = 0? (would mean odd):", eq(B0 @ Om + Om @ B0, Z16))

import math
ci = int(c)
s = math.isqrt(ci)
print("c is a perfect square:", s * s == ci, " sqrt =", s)
if s * s == ci and ci > 0:
    C = [[Fraction(int(B0[i, j]), s) for j in range(16)] for i in range(16)]
    def mm(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(16)) for j in range(16)] for i in range(16)]
    Cf = [[Fraction(int(B0[i, j]), s) for j in range(16)] for i in range(16)]
    Ef = [[[Fraction(int(E[a][i, j])) for j in range(16)] for i in range(16)] for a in range(8)]
    Omf = [[Fraction(int(Om[i, j])) for j in range(16)] for i in range(16)]
    C2 = mm(Cf, Cf)
    print("C^2 = I:", all(C2[i][j] == (1 if i == j else 0) for i in range(16) for j in range(16)))
    table = []
    for a in range(8):
        conj = mm(mm(Cf, Ef[a]), Cf)  # C E_a C^{-1} = C E_a C since C^2 = I
        tgt = Ef[a + 4] if a < 4 else Ef[a - 4]
        table.append(all(conj[i][j] == tgt[i][j] for i in range(16) for j in range(16)))
    print("Ad_C swaps every gamma with its sector partner:", table)
    COmC = mm(mm(Cf, Omf), Cf)
    print("C Om C^-1 = +Om:", all(COmC[i][j] == Omf[i][j] for i in range(16) for j in range(16)),
          " = -Om:", all(COmC[i][j] == -Omf[i][j] for i in range(16) for j in range(16)))

# ---- mod-p rank cross-checks of the dimension claims (independent of Sage) ----
MONi = [np.array(M, dtype=np.int64) for M in MON]
Ei = [np.array(M, dtype=np.int64) for M in E]
key2 = {}
for i, M in enumerate(MONi):
    key2[M.tobytes()] = (i, 1)
    key2[(-M).tobytes()] = (i, -1)
L = []
R = []
for a in range(8):
    LA = np.zeros((256, 256), dtype=np.int64)
    RA = np.zeros((256, 256), dtype=np.int64)
    for i, M in enumerate(MONi):
        j, sg = key2[(Ei[a] @ M).tobytes()]
        LA[j, i] = sg
        j2, sg2 = key2[(M @ Ei[a]).tobytes()]
        RA[j2, i] = sg2
    L.append(LA)
    R.append(RA)
CV = np.vstack([R[mu] - L[4 + mu] for mu in range(4)])
CW = np.vstack([R[4 + mu] - L[mu] for mu in range(4)])
SEL_even = np.zeros((128, 256), dtype=np.int64)
for r, i in enumerate([i for i, Ss in enumerate(subs) if len(Ss) % 2 == 0]):
    SEL_even[r, i] = 1

def rank_mod(A, p):
    A = (A % p).astype(np.int64).copy()
    rows, cols = A.shape
    r = 0
    for col in range(cols):
        piv = -1
        for rr in range(r, rows):
            if A[rr, col] % p != 0:
                piv = rr
                break
        if piv < 0:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, col]), -1, p)
        A[r] = (A[r] * inv) % p
        for rr in range(rows):
            if rr != r and A[rr, col] % p != 0:
                A[rr] = (A[rr] - A[rr, col] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r

for p in (10007, 10009):
    rv = rank_mod(CV, p)
    rvw = rank_mod(np.vstack([CV, CW]), p)
    rvwo = rank_mod(np.vstack([CV, CW, SEL_even]), p)
    print(f"p={p}: rank CV={rv} (dim V={256-rv}), rank CV+CW={rvw} (dim V+W={256-rvw}), "
          f"rank CV+CW+odd-forcing={rvwo} (odd part of V+W dim={256-rvwo})")
print("DONE-CHECK")
