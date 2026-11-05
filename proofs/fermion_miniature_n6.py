#!/usr/bin/env python3
"""
n=6 fermion miniature check (exact, sympy over Q(i)).

Construction: Cl(6) gamma matrices on C^8 by the Jordan-Wigner / Pauli
monomial construction (the same Pauli-string engine as the so-ladder
spinor work). Sections:

  A. Clifford relations, chirality Gamma, Sigma_ab blocks
  B. Per-block analysis: su(4) rank / closure / Killing / Casimir,
     4 vs 4bar conjugacy via the symmetric triple trace (d-symbol)
  C. Lambda^2 of the even block: so(6,C) w.r.t. the wedge form Q,
     anti-Hermiticity, invariant real structure (Hodge conjugation)
  D. Fock / CAR: a_k, a_k^dagger from gamma pairs; occupation numbers;
     Gamma = (-1)^N; vacuum; Fock basis = chirality sectors
  E. Double cover: exp(2*pi*Sigma_12) = -I on the spinor and on the
     4-block, +I on Lambda^2(4) and on the vector 6

Every check is an exact matrix identity over Q(i); no sampling.
"""
from fractions import Fraction
from itertools import combinations
from sympy import (Matrix, I, eye, zeros, kronecker_product, Rational,
                   conjugate, re as sy_re, im as sy_im, sin, cos, pi)

FAIL = []
COUNT = [0]

def check(cond, label):
    COUNT[0] += 1
    if not cond:
        FAIL.append(label)
        print("FAIL:", label)

def is_zero(M):
    return all(x == 0 for x in (M))

def eq(A, B):
    return is_zero(A - B)

# ---------------------------------------------------------------- A
sx = Matrix([[0, 1], [1, 0]])
sy = Matrix([[0, -I], [I, 0]])
sz = Matrix([[1, 0], [0, -1]])
id2 = eye(2)

gam = []
for k in range(6):
    p = k // 2                      # pair index 0,1,2
    fac = []
    for f in range(3):
        if f < p:
            fac.append(sz)
        elif f == p:
            fac.append(sx if k % 2 == 0 else sy)
        else:
            fac.append(id2)
    gam.append(kronecker_product(fac[0], fac[1], fac[2]))

I8 = eye(8)
for a in range(6):
    check(eq(gam[a].H, gam[a]), f"gamma{a+1} Hermitian")
    check(eq(gam[a] * gam[a], I8), f"gamma{a+1}^2 = I")
    for b in range(a + 1, 6):
        check(eq(gam[a] * gam[b] + gam[b] * gam[a], zeros(8)),
              f"{{gamma{a+1},gamma{b+1}}} = 0")

prod = gam[0]
for a in range(1, 6):
    prod = prod * gam[a]
GAM = kronecker_product(sz, sz, sz)
check(eq(I * prod, GAM), "Gamma = i*gamma1*gamma2*gamma3*gamma4*gamma5*gamma6 = sz(x)sz(x)sz")
check(eq(GAM.H, GAM), "Gamma Hermitian")
check(eq(GAM * GAM, I8), "Gamma^2 = I")
check(GAM.trace() == 0, "tr Gamma = 0")
for a in range(6):
    check(eq(GAM * gam[a] + gam[a] * GAM, zeros(8)), f"{{Gamma,gamma{a+1}}} = 0")

even_idx = [i for i in range(8) if GAM[i, i] == 1]
odd_idx = [i for i in range(8) if GAM[i, i] == -1]
check(len(even_idx) == 4 and len(odd_idx) == 4, "Gamma sectors 4+4")
print("even (Gamma=+1) basis indices:", even_idx, " odd:", odd_idx)

pairs6 = list(combinations(range(6), 2))          # (a,b), a<b, 15 pairs
Sig = [ (gam[a] * gam[b]) / 2 for (a, b) in pairs6 ]   # Sigma_ab on C^8
for s, (a, b) in zip(Sig, pairs6):
    check(eq(s.H, -s), f"Sigma{a+1}{b+1} anti-Hermitian")
    check(s.trace() == 0, f"Sigma{a+1}{b+1} traceless")
    check(eq(s * GAM - GAM * s, zeros(8)), f"[Sigma{a+1}{b+1},Gamma] = 0")
cas_dirac = sum((s * s for s in Sig), zeros(8))
check(eq(cas_dirac, Rational(-15, 4) * I8), "Dirac Casimir sum Sigma^2 = -15/4 I8")

def block(M, idx):
    return M.extract(idx, idx)

Sig_even = [block(s, even_idx) for s in Sig]
Sig_odd = [block(s, odd_idx) for s in Sig]

# ---------------------------------------------------------------- B
def vec_real(M):
    es = list(M)
    return Matrix([sy_re(x) for x in es] + [sy_im(x) for x in es])

def analyze(mats, name):
    n = mats[0].rows
    F = Matrix.hstack(*[vec_real(M) for M in mats])       # (2n^2) x 15
    r = F.rank()
    check(r == 15, f"{name}: real rank of the 15 block generators = 15")
    piv = list(F.T.rref()[1])                              # independent rows
    S = F.extract(piv, list(range(15)))
    Sinv = S.inv()
    check(S.det() != 0, f"{name}: feature minor invertible")
    def coeff(M):
        return Sinv * vec_real(M).extract(piv, [0])
    f = {}
    for a in range(15):
        for b in range(15):
            br = mats[a] * mats[b] - mats[b] * mats[a]
            c = coeff(br)
            check(eq(F * c, vec_real(br)), f"{name}: closure [{a},{b}]")
            f[(a, b)] = [Fraction(int(x.p), int(x.q)) for x in c]
    K = zeros(15)
    for a in range(15):
        for b in range(15):
            s = Fraction(0)
            for c in range(15):
                for d in range(15):
                    s += f[(a, c)][d] * f[(b, d)][c]
            K[a, b] = Rational(s.numerator, s.denominator)
    lam = K[0, 0]
    check(eq(K, lam * eye(15)), f"{name}: Killing = lambda * I, lambda = {lam}")
    cas = sum((M * M for M in mats), zeros(n))
    cval = cas[0, 0]
    check(eq(cas, cval * eye(n)), f"{name}: Casimir scalar, value {cval}")
    return {"killing_lambda": lam, "casimir": cval, "rank": r}

res_even = analyze(Sig_even, "even block")
res_odd = analyze(Sig_odd, "odd block")
print("even block:", res_even)
print("odd block:", res_odd)

def dsym(mats, a, b, c):
    return (mats[a] * (mats[b] * mats[c] + mats[c] * mats[b])).trace()

found = None
for a in range(15):
    for b in range(15):
        for c in range(15):
            v = dsym(Sig_even, a, b, c)
            if v != 0:
                found = (a, b, c, v, dsym(Sig_odd, a, b, c))
                break
        if found:
            break
    if found:
        break
check(found is not None, "d-symbol: a nonzero symmetric triple trace exists")
a, b, c, v_even, v_odd = found
print(f"d-symbol triple (Sigma indices {a},{b},{c}): even = {v_even}, odd = {v_odd}")
check(v_odd == -v_even and v_even != 0,
      "d-symbol: odd block = -(even block), so the blocks are conjugate, non-equivalent")

# vector Casimir control (record: -5)
cas_vec = zeros(6)
for (a, b) in pairs6:
    J = zeros(6)
    J[a, b] = 1
    J[b, a] = -1
    cas_vec += J * J
check(eq(cas_vec, -5 * eye(6)), "vector Casimir sum J^2 = -5 I6 (record value)")

# ---------------------------------------------------------------- C
pairs4 = list(combinations(range(4), 2))
pidx = {p: i for i, p in enumerate(pairs4)}

def perm_sign(t):
    inv = sum(1 for i in range(4) for j in range(i + 1, 4) if t[i] > t[j])
    return -1 if inv % 2 else 1

def wedge_vec(u, v):
    w = zeros(6, 1)
    for (p, q) in pairs4:
        w[pidx[(p, q)], 0] = u[p] * v[q] - u[q] * v[p]
    return w

def lam2(X):
    L = zeros(6)
    for (i, j) in pairs4:
        ei = zeros(4, 1); ei[i] = 1
        ej = zeros(4, 1); ej[j] = 1
        col = wedge_vec(X * ei, ej) + wedge_vec(ei, X * ej)
        for r in range(6):
            L[r, pidx[(i, j)]] = col[r]
    return L

Q = zeros(6)
Jm = zeros(6)
for (i, j) in pairs4:
    comp = tuple(x for x in range(4) if x not in (i, j))
    s = perm_sign((i, j, comp[0], comp[1]))
    Q[pidx[(i, j)], pidx[comp]] = s
    Jm[pidx[comp], pidx[(i, j)]] = s
check(eq(Jm * Jm, eye(6)), "Hodge conjugation J: J^2 = +I on Lambda^2")

Lam = [lam2(X) for X in Sig_even]
for L, (a, b) in zip(Lam, pairs6):
    check(eq(L.H, -L), f"Lambda2 Sigma{a+1}{b+1} anti-Hermitian")
    check(eq(L.T * Q + Q * L, zeros(6)), f"Lambda2 Sigma{a+1}{b+1} in so(6,C) wrt Q")
    check(eq(L * Jm, Jm * L.conjugate()), f"Lambda2 Sigma{a+1}{b+1} commutes with real structure J")
FL = Matrix.hstack(*[vec_real(L) for L in Lam])
check(FL.rank() == 15, "Lambda2 image has rank 15")

# ---------------------------------------------------------------- D
aa = [(gam[2 * k] + I * gam[2 * k + 1]) / 2 for k in range(3)]
ad = [(gam[2 * k] - I * gam[2 * k + 1]) / 2 for k in range(3)]
for k in range(3):
    check(eq(ad[k], aa[k].H), f"a{k+1}^dagger is the adjoint of a{k+1}")
    for l in range(3):
        d = I8 if k == l else zeros(8)
        check(eq(aa[k] * ad[l] + ad[l] * aa[k], d), f"{{a{k+1},a{l+1}^dagger}} = delta")
        check(eq(aa[k] * aa[l] + aa[l] * aa[k], zeros(8)), f"{{a{k+1},a{l+1}}} = 0")
        check(eq(ad[k] * ad[l] + ad[l] * ad[k], zeros(8)), f"{{a{k+1}^dagger,a{l+1}^dagger}} = 0")
N = [ad[k] * aa[k] for k in range(3)]
for k in range(3):
    check(eq(N[k] * N[k], N[k]), f"N{k+1}^2 = N{k+1} (occupation projector)")
    for l in range(3):
        check(eq(N[k] * N[l] - N[l] * N[k], zeros(8)), f"[N{k+1},N{l+1}] = 0")
P = (I8 - 2 * N[0]) * (I8 - 2 * N[1]) * (I8 - 2 * N[2])
check(eq(P, GAM), "prod_k (1 - 2 N_k) = Gamma, i.e. Gamma = (-1)^(N_1+N_2+N_3)")

Astack = Matrix.vstack(aa[0], aa[1], aa[2])
ns = Astack.nullspace()
check(len(ns) == 1, "common kernel of a1,a2,a3 is 1-dimensional (unique vacuum)")
vac = ns[0]
check(eq(vac.H * vac, Matrix([[1]])), "vacuum normalized")
states = []
for mask in range(8):
    w = vac
    for k in range(3):
        if mask & (1 << k):
            w = ad[k] * w
    states.append(w)
    pc = bin(mask).count("1")
    check(eq(GAM * w, ((-1) ** pc) * w),
          f"Fock state mask {mask}: Gamma eigenvalue (-1)^{pc}")
U = Matrix.hstack(*states)
check(eq(U.H * U, I8), "the 8 Fock states are an orthonormal basis of C^8")

# ---------------------------------------------------------------- E
S12 = Sig[0]                                    # Sigma_12, pairs6[0] = (0,1)
check(eq(S12 * S12, Rational(-1, 4) * I8), "Sigma12^2 = -1/4 I8")
check(eq((2 * S12) * (2 * S12), -I8),
      "exp(2 pi Sigma12) = (exp(pi Sigma12))^2 = (2 Sigma12)^2 = -I8 on the spinor")
B12 = Sig_even[0]
check(eq(B12 * B12, Rational(-1, 4) * eye(4)), "block Sigma12^2 = -1/4 I4")
check(eq((2 * B12) * (2 * B12), -eye(4)), "exp(2 pi Sigma12) = -I4 on the 4-block")
G12 = Lam[0]
check(eq(G12 * G12 * G12 + G12, zeros(6)), "Lambda2(Sigma12): G^3 + G = 0 (eigenvalues 0, +/- i)")
expL = eye(6) + sin(2 * pi) * G12 + (1 - cos(2 * pi)) * (G12 * G12)
check(eq(expL, eye(6)), "exp(2 pi * Lambda2(Sigma12)) = +I6 on Lambda^2(4)")
J12 = zeros(6); J12[0, 1] = 1; J12[1, 0] = -1
check(eq(J12 * J12 * J12 + J12, zeros(6)), "vector J12: J^3 + J = 0")
expJ = eye(6) + sin(2 * pi) * J12 + (1 - cos(2 * pi)) * (J12 * J12)
check(eq(expJ, eye(6)), "exp(2 pi J12) = +I6 on the vector 6")

# ---------------------------------------------------------------- summary
print()
print(f"CHECKS RUN: {COUNT[0]}, FAILURES: {len(FAIL)}")
if FAIL:
    print("FAILED LABELS:", FAIL)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
