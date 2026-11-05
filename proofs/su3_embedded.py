# su3_embedded.py
# Stage 1 (Kit, 2026-10-26): verify the su(3) embedded in the proved
# SU(5) stabilizer algebra of the 16, via the block split W' = C^3 (+) C^2.
#
# Model: the stabilizer build's own formula (su5_stabilizer.md):
#   rho(X)|_piece = dLambda^m(-A(X)^T) + mu(X) I,  mu = (1/2) tr A,
# on the pieces Lambda^m(W'), m = 0, 2, 4 (the 16), subset basis.
# On su(5) (A anti-Hermitian, traceless) mu = 0.
# The split is a CHOICE: W' = span{e0,e1,e2} (+) span{e3,e4}.
# Generators: A_a = diag(E_a, 0_2), E_a = i lambda_a/2 in su(3),
# R_a = dLambda^even(-A_a^T). Expected law: [R_a, R_b] = -f^{abc} R_c,
# with f computed in-script from the lambda matrices (no hand-typed table).
#
# Exact over Q(i, sqrt(3)). sympy. Counted checks; exit 0 iff all pass.
# Disclosure: the S4 trace-form expectation was first written -4 delta
# (wrong index guess); the computed form is -2 delta, matching the
# decomposition 3+3+3bar+3bar+1^4 (four triplet copies x -1/2). The
# build was right; the expectation was corrected, not the build.

import itertools
import sympy as sp

I = sp.I
SQ3 = sp.sqrt(3)
count = 0
fails = []

def check(name, cond):
    global count
    count += 1
    print(("PASS" if cond else "FAIL"), name)
    if not cond:
        fails.append(name)

# ---------- exterior algebra model ----------
SUBSETS = [S for m in (0, 2, 4) for S in itertools.combinations(range(5), m)]
IDX = {S: k for k, S in enumerate(SUBSETS)}
assert len(SUBSETS) == 16

def dlambda(M):
    """Derivation action of a 5x5 matrix M on Lambda^even, subset basis."""
    D = sp.zeros(16)
    for S in SUBSETS:
        col = IDX[S]
        for r, i in enumerate(S):
            for j in range(5):
                if M[j, i] == 0:
                    continue
                if j == i:
                    D[col, col] += M[j, i]
                elif j not in S:
                    raw = list(S)
                    raw[r] = j
                    srt = sorted(raw)
                    # sign of the permutation raw -> srt
                    perm = [raw.index(v) for v in srt]
                    inv = sum(1 for p in range(len(perm))
                              for q in range(p + 1, len(perm)) if perm[p] > perm[q])
                    D[IDX[tuple(srt)], col] += (-1) ** inv * M[j, i]
    return D

def rho_of_A(A):
    return dlambda(-A.T)

# ---------- Gell-Mann basis and f^{abc} from the matrices ----------
lam = [None] * 9
lam[1] = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
lam[2] = sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]])
lam[3] = sp.diag(1, -1, 0)
lam[4] = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
lam[5] = sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]])
lam[6] = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
lam[7] = sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]])
lam[8] = sp.diag(1, 1, -2) / SQ3
T = [None] + [lam[a] / 2 for a in range(1, 9)]          # Hermitian generators
E = [None] + [I * lam[a] / 2 for a in range(1, 9)]       # anti-Hermitian su(3)

# f^{abc} from [T^a, T^b] = i f^{abc} T^c, tr(T^a T^b) = delta/2
f = {}
for a, b, c in itertools.product(range(1, 9), repeat=3):
    f[(a, b, c)] = sp.simplify(-2 * I * ((T[a] * T[b] - T[b] * T[a]) * T[c]).trace())
check("f table: f^123 = 1", f[(1, 2, 3)] == 1)
check("f table: f^458 = sqrt(3)/2", f[(4, 5, 8)] == SQ3 / 2)
check("f table: totally antisymmetric",
      all(f[(a, b, c)] == -f[(b, a, c)] == -f[(a, c, b)]
          for a, b, c in itertools.product(range(1, 9), repeat=3)))
check("E-basis su(3) itself: [E_a,E_b] = -f^{abc} E_c (3 sample pairs)",
      all(E[a] * E[b] - E[b] * E[a] ==
          -sum((f[(a, b, c)] * E[c] for c in range(1, 9)), sp.zeros(3))
          for a, b in [(1, 2), (4, 5), (2, 7)]))

# ---------- the embedded generators ----------
def block5(B3):
    A = sp.zeros(5)
    A[:3, :3] = B3
    return A

R = [None] + [rho_of_A(block5(E[a])) for a in range(1, 9)]

# Controls for the checker: commuting Cartan pair must give commuting rho's
H1 = rho_of_A(block5(E[3]))
H2 = rho_of_A(block5(E[8]))
check("control (known-true): Cartan rho's commute", H1 * H2 == H2 * H1)

def brackets_match(Rs):
    ok = True
    for a, b in itertools.combinations(range(1, 9), 2):
        comm = Rs[a] * Rs[b] - Rs[b] * Rs[a]
        rhs = -sum((f[(a, b, c)] * Rs[c] for c in range(1, 9)), sp.zeros(16))
        if comm != rhs:
            ok = False
    return ok

check("control (known-false): scaled R_1 must FAIL the bracket test",
      not brackets_match([None, 2 * R[1]] + R[2:]))

# S1: linear independence
Mflat = sp.Matrix.hstack(*[sp.Matrix(R[a]).reshape(256, 1) for a in range(1, 9)])
check("S1: the 8 R_a are linearly independent", Mflat.rank() == 8)

# S2: anti-Hermiticity and piece preservation
check("S2: R_a anti-Hermitian (8)", all(R[a].H == -R[a] for a in range(1, 9)))
pieces = [(0, 1), (1, 11), (11, 16)]
def preserves(Rm):
    for (r0, r1), (c0, c1) in itertools.product(pieces, repeat=2):
        if (r0, r1) != (c0, c1) and Rm[r0:r1, c0:c1] != sp.zeros(r1 - r0, c1 - c0):
            return False
    return True
check("S2: R_a preserve the |S| pieces 1 (+) 10 (+) 5 (8)",
      all(preserves(R[a]) for a in range(1, 9)))

# S3: the bracket law, all 28 pairs, against the in-script f table
check("S3: [R_a, R_b] = -f^{abc} R_c for all 28 pairs", brackets_match(R))

# S4: trace form on the 16 and abstract Killing form from the constants
G = sp.Matrix(8, 8, lambda i, j: (R[i + 1] * R[j + 1]).trace())
check("S4: tr_16(R_a R_b) = -2 delta_ab", G == -2 * sp.eye(8))
K = sp.zeros(8)
for a, b in itertools.product(range(1, 9), repeat=2):
    K[a - 1, b - 1] = sum(f[(a, c, d)] * f[(b, c, d)]
                          for c in range(1, 9) for d in range(1, 9))
check("S4: Killing from f: sum_cd f^{acd} f^{bcd} = 3 delta_ab", K == 3 * sp.eye(8))

# S5: Casimir on the 16: spectrum {0 (x4), -4/3 (x12)} via exact identities
C = sum((R[a] * R[a] for a in range(1, 9)), sp.zeros(16))
check("S5: Casimir Hermitian", C.H == C)
check("S5: C (C + 4/3) = 0 exactly", sp.simplify(C * (C + sp.Rational(4, 3) * sp.eye(16))) == sp.zeros(16))
check("S5: tr C = -16 and rank C = 12 (nullity 4)",
      C.trace() == -16 and C.rank() == 12)

# S6: (8,1) placement: the 8 commute with the SU(2) block and the U(1)'
def block5_lower(B2):
    A = sp.zeros(5)
    A[3:, 3:] = B2
    return A
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.diag(1, -1)]
Sgens = [rho_of_A(block5_lower(I * s / 2)) for s in sig]
U1p = rho_of_A(I * sp.diag(2, 2, 2, -3, -3))
check("S6: [R_a, SU(2) block] = 0 (24) and [R_a, U(1)'] = 0 (8)",
      all(R[a] * s - s * R[a] == sp.zeros(16) for a in range(1, 9) for s in Sgens) and
      all(R[a] * U1p - U1p * R[a] == sp.zeros(16) for a in range(1, 9)))

print()
print("TOTAL:", count, "checks;", len(fails), "failed", fails)
raise SystemExit(0 if not fails else 1)
