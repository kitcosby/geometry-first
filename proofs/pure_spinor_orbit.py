#!/usr/bin/env python3
"""
Pure spinor orbit for the Spin(10) Weyl spinor 16.

Claim under test: the Spin(10,C)-orbit of a highest-weight (Fock vacuum)
spinor in P(16) is the pure spinor variety OG+(5,10): dim 10, stabilizer
the maximal parabolic with Levi GL(5,C); the other Weyl spinor 16' gives
the opposite component OG-(5,10).

Method: exact arithmetic over Q(i) with the corpus Cl(10) gamma matrices
(Pauli recursion, as in spin10_spinor.md). Everything is a finite exact
matrix computation: Clifford relations, chirality split, annihilator
(purity) of the vacuum and of the top vector, stabilizer of the spinor
line (nullspace dimension + block form in the W (+) W' decomposition),
Casimir on both Weyl blocks.

Checker controls (run first, per symbolic_checker_validation):
  TRUE  control: {gamma_1, gamma_2} = 0 must be accepted.
  FALSE control: gamma_1 gamma_2 - gamma_2 gamma_1 must be REJECTED as 0
                 (it equals 2 gamma_1 gamma_2 != 0).
  FALSE control (purity): psi_c = e_0 + e_{11110} (vacuum plus a 4-fold
                 occupation, no 2-fold part) must NOT have annihilator
                 dimension 5 -- its degree-4 part is not the square of
                 any degree-2 part, violating the pure-spinor quadrics.
"""
import sys
from sympy import Matrix, eye, zeros, I, Rational, kronecker_product

PASS = 0
FAIL = 0
FAILURES = []

def check(name, cond):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS  {name}")
    else:
        FAIL += 1
        FAILURES.append(name)
        print(f"  FAIL  {name}")

def is_zero(M):
    return all(x == 0 for x in M)

def eq(A, B):
    return is_zero(A - B)

# ---------------------------------------------------------------- gammas
s1 = Matrix([[0, 1], [1, 0]])
s2 = Matrix([[0, -I], [I, 0]])
s3 = Matrix([[1, 0], [0, -1]])
gam = [s1, s2]
while gam[0].rows < 32:
    n = gam[0].rows
    gam = [kronecker_product(g, s3) for g in gam] + \
          [kronecker_product(eye(n), s1), kronecker_product(eye(n), s2)]
G = gam  # G[a], a = 0..9, 32x32, entries in {0, +-1, +-i}
I32 = eye(32)

print("== checker controls ==")
check("TRUE control: {g1,g2} = 0 accepted", eq(G[0]*G[1] + G[1]*G[0], zeros(32)))
check("FALSE control: g1 g2 - g2 g1 rejected as 0",
      not is_zero(G[0]*G[1] - G[1]*G[0]))

print("== Clifford relations and chirality ==")
ok = True
for a in range(10):
    for b in range(10):
        want = 2 * I32 if a == b else zeros(32)
        if not eq(G[a]*G[b] + G[b]*G[a], want):
            ok = False
check("{gamma_a, gamma_b} = 2 delta_ab I, all 100 pairs", ok)

prod = I32
for a in range(10):
    prod = prod * G[a]
Gamma = I * prod
check("Gamma^2 = I", eq(Gamma*Gamma, I32))
check("Gamma diagonal", all(Gamma[i, j] == 0 for i in range(32)
                            for j in range(32) if i != j))
check("tr Gamma = 0", Gamma.trace() == 0)
# Gamma = -Z^{x5}: eigenvalue on basis vector e_n is -(-1)^{|n|}
check("Gamma e_0 = -e_0 (vacuum sector is Gamma = -1)",
      Gamma[0, 0] == -1)
check("Gamma e_31 = +e_31 (top vector in Gamma = +1)",
      Gamma[31, 31] == 1)

# spin generators and their Weyl blocks
pairs = [(a, b) for a in range(10) for b in range(a+1, 10)]
Sig = {p: (G[p[0]]*G[p[1]] - G[p[1]]*G[p[0]]) * Rational(1, 4) for p in pairs}
minus_idx = [n for n in range(32) if Gamma[n, n] == -1]
plus_idx = [n for n in range(32) if Gamma[n, n] == 1]
check("both chirality sectors have dim 16",
      len(minus_idx) == 16 and len(plus_idx) == 16)

def block(M, idx):
    return M.extract(idx, idx)

SigM = {p: block(Sig[p], minus_idx) for p in pairs}   # the 16 (vacuum sector)
SigP = {p: block(Sig[p], plus_idx) for p in pairs}    # the 16'

print("== Casimir on both Weyl blocks (ties to layer 1) ==")
for name, S in (("16 (Gamma=-1)", SigM), ("16' (Gamma=+1)", SigP)):
    C = zeros(16)
    for p in pairs:
        C += S[p] * S[p]
    check(f"Casimir on {name} = -45/4 I", eq(C, Rational(-45, 4) * eye(16)))

# ---------------------------------------------------------------- purity
# L_psi : C^10 -> C^32, v |-> (sum_a v_a gamma_a) psi.  W_psi = ker L_psi.
def annihilator(psi):
    M = Matrix.hstack(*[G[a] * psi for a in range(10)])  # 32x10
    return M.nullspace()  # basis of ker, as column vectors in C^10

def iso_gram(W):
    B = zeros(len(W))
    for i, w in enumerate(W):
        for j, z in enumerate(W):
            B[i, j] = sum(w[k] * z[k] for k in range(10))
    return B

e0 = zeros(32, 1); e0[0] = 1
e31 = zeros(32, 1); e31[31] = 1

print("== purity of the vacuum (the 16) ==")
W = annihilator(e0)
check("dim W_{e0} = 5 (e0 is a pure spinor)", len(W) == 5)
check("W_{e0} is isotropic: Gram = 0", is_zero(iso_gram(W)))
# independent second route: the a_k coefficient vectors
def avec(k, sign):
    v = zeros(10, 1)
    v[2*k] = Rational(1, 2)
    v[2*k+1] = sign * I * Rational(1, 2)
    return v
a_vecs = [avec(k, +1) for k in range(5)]    # a_k  = (g_{2k} + i g_{2k+1})/2
ad_vecs = [avec(k, -1) for k in range(5)]   # a_k+ = (g_{2k} - i g_{2k+1})/2
MW = Matrix.hstack(*W)
check("each a_k coefficient vector lies in W_{e0}",
      all(MW.rank() == Matrix.hstack(MW, v).rank() for v in a_vecs))
check("each a_k annihilates e0 directly",
      all(is_zero(sum((v[a] * G[a] for a in range(10)), zeros(32)) * e0)
          for v in a_vecs))

print("== purity of the top vector (the 16') and the family split ==")
Wp = annihilator(e31)
check("dim W_{e31} = 5 (e31 is a pure spinor)", len(Wp) == 5)
check("W_{e31} is isotropic: Gram = 0", is_zero(iso_gram(Wp)))
MWp = Matrix.hstack(*Wp)
check("each a_k+ coefficient vector lies in W_{e31}",
      all(MWp.rank() == Matrix.hstack(MWp, v).rank() for v in ad_vecs))
gram_cross = zeros(5)
for i, w in enumerate(W):
    for j, z in enumerate(Wp):
        gram_cross[i, j] = sum(w[k] * z[k] for k in range(10))
check("pairing W x W' nondegenerate (det Gram != 0, nullspace bases)",
      gram_cross.det() != 0)
gram_av = zeros(5)
for i, v in enumerate(a_vecs):
    for j, z in enumerate(ad_vecs):
        gram_av[i, j] = sum(v[k] * z[k] for k in range(10))
check("explicit a_k x a_k+ pairing Gram = (1/2) I",
      eq(gram_av, Rational(1, 2) * eye(5)))
both = Matrix.hstack(MW, MWp)
check("W (+) W' = C^10 (rank 10, intersection 0)",
      both.rank() == 10)
# parity: dim(W n W') = 0 is EVEN, n = 5 is odd -> opposite families.
print("== same-family parity sanity (graph of a skew map) ==")
# W_S = { w + S w }, S: W -> W', matrix E_01 - E_10 in the two bases
Smat = zeros(5); Smat[0, 1] = 1; Smat[1, 0] = -1
WS = []
for k in range(5):
    v = MW[:, k] + MWp * Smat[:, k]
    WS.append(v)
MWS = Matrix.hstack(*WS)
inter = Matrix.hstack(MW, -MWS).nullspace()
# dim(W n W_S) = 5 + 5 - rank[W | W_S]
dim_inter = 10 - Matrix.hstack(MW, MWS).rank()
check("dim(W n W_graph) = 3 (odd, same family as W)", dim_inter == 3)

print("== FALSE control: a non-pure spinor ==")
psi_c = zeros(32, 1); psi_c[0] = 1; psi_c[30] = 1   # e_0 + e_{11110}
Wc = annihilator(psi_c)
print(f"       annihilator dim of psi_c = {len(Wc)}")
check("psi_c is NOT pure (annihilator dim != 5)", len(Wc) != 5)

# ------------------------------------------------------------- stabilizer
print("== stabilizer of the vacuum line (spin side) ==")
# X = sum c_ab SigM_ab ; condition: X e0 in C e0, e0 = first minus-index (0)
pos0 = minus_idx.index(0)
rows = [r for r in range(16) if r != pos0]
Aeq = zeros(15, 45)
for j, p in enumerate(pairs):
    col = SigM[p][:, pos0]
    for i, r in enumerate(rows):
        Aeq[i, j] = col[r]
stab_spin = Aeq.nullspace()
check("spin-side stabilizer dim = 35", len(stab_spin) == 35)
check("45 - 35 = 10 = dim OG(5,10)", 45 - len(stab_spin) == 10)

print("== stabilizer of W (vector side) and Levi form ==")
# E_ab basis of so(10,C); adapted basis P = [W | W']; condition: lower-left
# 5x5 block of P^{-1} Y P vanishes.
P = both
Pinv = P.inv()
E = {}
for (a, b) in pairs:
    M = zeros(10); M[a, b] = 1; M[b, a] = -1
    E[(a, b)] = M
AeqV = zeros(25, 45)
for j, p in enumerate(pairs):
    Yh = Pinv * E[p] * P
    for i in range(25):
        AeqV[i, j] = Yh[5 + i % 5, i // 5]
stab_vec = AeqV.nullspace()
check("vector-side stabilizer dim = 35", len(stab_vec) == 35)

def blocks_of(cvec):
    Y = sum((cvec[j] * E[p] for j, p in enumerate(pairs)), zeros(10))
    Yh = Pinv * Y * P
    return Yh[:5, :5], Yh[:5, 5:], Yh[5:, 5:]

okA = okB = okD = True
Arows, Brows = [], []
for cvec in stab_vec:
    A, B, D = blocks_of(cvec)
    if not is_zero(D + A.T):
        okD = False
    if not is_zero(B + B.T):
        okB = False
    Arows.append(list(A))
    Brows.append([B[i, j] for i in range(5) for j in range(i+1, 5)])
check("Levi block form: D = -A^T on all 35 basis elements", okD)
check("unipotent block: B = -B^T (skew) on all 35 basis elements", okB)
check("A-projection surjective onto gl(5): rank 25",
      Matrix(Arows).rank() == 25)
check("B-projection surjective onto Lambda^2: rank 10",
      Matrix(Brows).rank() == 10)
# Levi part = kernel of B-projection: dim 35 - 10 = 25 = dim GL(5) algebra

print("== rho identifies the two stabilizers ==")
ok = True
for cvec in stab_vec:
    Y = sum((cvec[j] * Sig[p] for j, p in enumerate(pairs)), zeros(32))
    img = Y * e0
    if any(img[m] != 0 for m in range(32) if m != 0):
        ok = False
check("rho(stab_vec) preserves the line C e0 (all 35)", ok)
ok = True
for cvec in stab_spin:
    X = sum((cvec[j] * SigM[p] for j, p in enumerate(pairs)), zeros(16))
    img = X[:, pos0]
    if any(img[r] != 0 for r in range(16) if r != pos0):
        ok = False
check("spin stabilizer basis indeed fixes the line (self-consistency)", ok)

print("== line vs vector stabilizer: the character on the vacuum line ==")
# For X in the line stabilizer, rho(X) e0 = mu(X) e0. On the vector side,
# with blocks (A, B, -A^T) in the W (+) W' decomposition, the computation
# gives mu = +(1/2) tr(A) -- the Fock vacuum charge, with A acting on the
# ANNIHILATOR space W (the more common -1/2 sign is the same statement
# parametrized by the action on the creation space W'). The VECTOR
# stabilizer is ker mu inside the 35-dim line stabilizer: dim 34, Levi
# sl(5) (traceless A) + the same 10-dim unipotent radical.
mus = []
okmu = True
for cvec in stab_vec:
    A, B, D = blocks_of(cvec)
    Y = sum((cvec[j] * Sig[p] for j, p in enumerate(pairs)), zeros(32))
    mu = (Y * e0)[0]
    mus.append(mu)
    if mu != Rational(1, 2) * A.trace():
        okmu = False
check("character identity: mu = +(1/2) tr(A) on all 35 basis elements", okmu)
# the functional on the spin-side basis, for the kernel dimension:
mus_spin = []
for cvec in stab_spin:
    X = sum((cvec[j] * SigM[p] for j, p in enumerate(pairs)), zeros(16))
    mus_spin.append(X[pos0, pos0])
null_mu = Matrix([mus_spin]).nullspace()
check("vector stabilizer = ker mu has dim 34 (Levi sl(5): 24, + unipotent 10)",
      len(null_mu) == 34)

print()
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
