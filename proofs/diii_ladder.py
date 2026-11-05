#!/usr/bin/env python3
"""
The DIII ladder: the Spin(2n)/U(n) layer stack for general n,
derived once and verified rung by rung (n = 2..6).

Per rung, exact over Q(i):
  A. gamma construction, Weyl sector (dim 2^{n-1}), spinor
     Casimir = -n(2n-1)/4 on the sector.
  B. Parabolic (vector-side solve): line stabilizer of the pure
     spinor, dim n^2 + n(n-1)/2, block form [[A, B], [0, -A^T]],
     character mu = (1/2) tr A on the spinor line; orbit dim
     n(n-1)/2 = dim OG+(n, 2n).
  C. Compact stabilizer (real solve): dim n^2; A-blocks
     anti-Hermitian spanning u(n); centre dim 1; derived =
     ker mu, dim n^2 - 1, A-blocks spanning su(n).
  D. Branching: pieces |S| = m (m even) have dim C(n, m),
     charge m - n/2 under Z (A(Z) = -i I_n); operator identity
     rho(X) = dLambda^m(-A^T) + mu I on the sparse u(n) basis;
     weights test with a distinct-weight Cartan element:
     piece-m eigenvalues are exactly {-i * sum_{k in S} w_k}.
  E. Bundle/Kostant algebra: m (A-block = 0) real dim n(n-1);
     [h, m] subset m, [m, m] subset h; kappa = -tr nondegenerate
     on m; kappa(h, m) = 0; cubic phi = 0 on all triples.
  F. Monomial engine (n <= 5 only): gamma relations on
     Lambda*(n) (dim 2^nu, nu = n(n-1)/2), spin lift sigma
     intertwining + trace 0, gamma-contraction = nu * I.
     At n = 6 (2^15 = 32768 monomials) the evaluation is not
     run; those statements are the n-uniform CAR derivation
     of diii_ladder.md, whose per-rung matrix inputs (ad on n
     = dLambda^2(A), checked in E via the radical commutators)
     ARE verified at n = 6.

Controls per rung: checker TRUE/FALSE pair on the gamma
anticommutators; the centre candidate Z = -sum_k E_{2k,2k+1}
must be central with A(Z) = -i I (a derived identification,
not an assumption).
"""
import sys
from fractions import Fraction
from itertools import combinations
from sympy import (Matrix, eye, zeros, I, Rational, kronecker_product,
                   Symbol, re as Re, im as Im)

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

def real_nullspace_from_complex_rows(rows):
    rr = []
    for row in rows:
        rr.append([Re(x) for x in row])
        rr.append([Im(x) for x in row])
    return Matrix(rr).nullspace()

def complex_nullity(rows, ncols):
    """nullity over Q(i) of a complex matrix given as row lists,
    via the doubled real system."""
    rr = []
    for row in rows:
        rr.append([Re(x) for x in row] + [-Im(x) for x in row])
        rr.append([Im(x) for x in row] + [Re(x) for x in row])
    M = Matrix(rr)
    return ncols - (M.rank() // 2)

# ------------------------------------------------------- sparse engine
Z2 = (Fraction(0), Fraction(0))
def fp(z):
    r, i = Re(z), Im(z)
    return (Fraction(r.p, r.q), Fraction(i.p, i.q))
def cadd(u, v):
    return (u[0] + v[0], u[1] + v[1])
def cmul(u, v):
    return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def dacc(acc, d, scale):
    for m, c in d.items():
        v = cmul(c, scale)
        acc[m] = cadd(acc[m], v) if m in acc else v
    return acc
def clean(d):
    return {m: c for m, c in d.items() if c != Z2}
def popc(x):
    return bin(x).count("1")

def run_rung(n, engine):
    print(f"===== rung n = {n}: so({2*n}), spinor {2**(n-1)} =====")
    N = 2 * n
    D = 2 ** n
    s1 = Matrix([[0, 1], [1, 0]])
    s2 = Matrix([[0, -I], [I, 0]])
    s3 = Matrix([[1, 0], [0, -1]])
    gam = [s1, s2]
    while gam[0].rows < D:
        d0 = gam[0].rows
        gam = [kronecker_product(g, s3) for g in gam] + \
              [kronecker_product(eye(d0), s1),
               kronecker_product(eye(d0), s2)]
    G = gam
    check(f"n={n} TRUE control: {{g1,g2}} = 0",
          eq(G[0]*G[1] + G[1]*G[0], zeros(D)))
    check(f"n={n} FALSE control: [g1,g2] != 0",
          not is_zero(G[0]*G[1] - G[1]*G[0]))
    prod = eye(D)
    for a in range(N):
        prod = prod * G[a]
    # (g_1...g_{2n})^2 = (-1)^n: Gamma = iP for n odd, P for n even
    Gamma = I * prod if n % 2 == 1 else prod
    # the rung's spinor is the Weyl sector containing the vacuum e_0;
    # its Gamma sign alternates with n (period 4), so take it from
    # the matrix, not from a fixed convention
    g00 = Gamma[0, 0]
    minus_idx = [k for k in range(D) if Gamma[k, k] == g00]
    pos0 = minus_idx.index(0)
    sec = len(minus_idx)
    check(f"n={n} Weyl sector dim {2**(n-1)}", sec == 2**(n-1))
    pairs = [(a, b) for a in range(N) for b in range(a+1, N)]
    NP = len(pairs)
    Sig = {p: (G[p[0]]*G[p[1]] - G[p[1]]*G[p[0]]) * Rational(1, 4)
           for p in pairs}
    SigM = {p: Sig[p].extract(minus_idx, minus_idx) for p in pairs}
    # A. Casimir
    Cas = sum((SigM[p] * SigM[p] for p in pairs), zeros(sec))
    target = -Rational(n * (2*n - 1), 4)
    check(f"n={n} spinor Casimir = {target} on sector",
          eq(Cas, target * eye(sec)))
    # adapted basis
    def avec(k, sign):
        v = zeros(N, 1)
        v[2*k] = Rational(1, 2)
        v[2*k+1] = sign * I * Rational(1, 2)
        return v
    Pa = Matrix.hstack(*([avec(k, +1) for k in range(n)] +
                         [avec(k, -1) for k in range(n)]))
    PaInv = Pa.inv()
    E = {}
    for (a, b) in pairs:
        M = zeros(N); M[a, b] = 1; M[b, a] = -1
        E[(a, b)] = M
    Yhat = {p: PaInv * E[p] * Pa for p in pairs}
    SK = list(combinations(range(n), 2))
    nu = len(SK)
    U = {}
    for t, (a, b) in enumerate(SK):
        M = zeros(n); M[a, b] = 1; M[b, a] = -1
        U[t] = M
    # sparse u(n) basis from A-units
    Aunits = []
    for p in range(n):
        A = zeros(n); A[p, p] = I; Aunits.append(A)
    for p in range(n):
        for q in range(p+1, n):
            A = zeros(n); A[p, q] = 1; A[q, p] = -1; Aunits.append(A)
            A = zeros(n); A[p, q] = I; A[q, p] = I; Aunits.append(A)
    Hs, Ahat, rhoS, muS = [], [], [], []
    ok_sparse = True
    for A in Aunits:
        Yh = zeros(N); Yh[:n, :n] = A; Yh[n:, n:] = -A.T
        Y = Pa * Yh * PaInv
        if any(Im(Y[r, c]) != 0 for r in range(N) for c in range(N)):
            ok_sparse = False
        Yr = Matrix(N, N, lambda r, c: Re(Y[r, c]))
        # E-coords for the spin action
        cf = None
        Hs.append(Yr); Ahat.append(A)
    # E-coordinate solver (shared)
    VE = Matrix.hstack(*[E[p].reshape(N*N, 1) for p in pairs])
    pivE = VE.T.rref()[1]
    VpE = VE.extract(list(pivE), list(range(NP)))
    VpEInv = VpE.inv()
    def ecoords(Y):
        return VpEInv * Y.reshape(N*N, 1).extract(list(pivE), [0])
    for Yr, A in zip(Hs, Ahat):
        cf = ecoords(Yr)
        if not eq(sum((cf[j] * E[p] for j, p in enumerate(pairs)),
                      zeros(N)), Yr):
            ok_sparse = False
        R = sum((cf[j] * SigM[p] for j, p in enumerate(pairs)),
                zeros(sec))
        col = R[:, pos0]
        if any(col[r] != 0 for r in range(sec) if r != pos0):
            ok_sparse = False
        rhoS.append(R); muS.append(col[pos0])
    check(f"n={n} sparse u(n) basis: real, line-preserving", ok_sparse)
    check(f"n={n} mu = (1/2) tr A on sparse basis",
          all(muS[j] == Rational(1, 2) * Ahat[j].trace()
              for j in range(n * n)))
    # B. parabolic (vector side): Y W subset W  <=>  C-block = 0
    rowsC = [[Yhat[p][r, s] for p in pairs]
             for r in range(n, N) for s in range(n)]
    dim_par = complex_nullity(rowsC, NP)
    check(f"n={n} parabolic dim_C = n^2 + nu = {n*n + nu}",
          dim_par == n * n + nu)
    check(f"n={n} orbit dim_C = nu = {nu}", NP - dim_par == nu)
    # C. compact stabilizer: also B-block = 0
    rowsBC = rowsC + [[Yhat[p][r, s] for p in pairs]
                      for r in range(n) for s in range(n, N)]
    stab = real_nullspace_from_complex_rows(rowsBC)
    check(f"n={n} compact stabilizer dim_R = n^2 = {n*n}",
          len(stab) == n * n)
    Ys, As, mus = [], [], []
    okAH = True
    for cv in stab:
        Y = sum((cv[j] * E[p] for j, p in enumerate(pairs)), zeros(N))
        A = (PaInv * Y * Pa)[:n, :n]
        if not is_zero(A.conjugate().T + A):
            okAH = False
        Ys.append(Y); As.append(A)
        cf = ecoords(Y)
        R = sum((cf[j] * SigM[p] for j, p in enumerate(pairs)),
                zeros(sec))
        mus.append(R[pos0, pos0])
    check(f"n={n} stabilizer A-blocks anti-Hermitian", okAH)
    def real_rank(mats_entries):
        return Matrix(mats_entries).rank()
    Avec = [[Re(A[r, c]) for r in range(n) for c in range(n)] +
            [Im(A[r, c]) for r in range(n) for c in range(n)]
            for A in As]
    check(f"n={n} A-span = u(n) (real rank {n*n})",
          real_rank(Avec) == n * n)
    # centre + derived via brackets in E-coords of Ys
    Yc = [ecoords(Y) for Y in Ys]
    VY = Matrix.hstack(*[Y.reshape(N*N, 1) for Y in Ys])
    pivY = VY.T.rref()[1]
    VpY = VY.extract(list(pivY), list(range(n * n)))
    VpYInv = VpY.inv()
    def ycoords(Y):
        return VpYInv * Y.reshape(N*N, 1).extract(list(pivY), [0])
    brackets = {}
    for i in range(n * n):
        for j in range(i+1, n * n):
            brackets[(i, j)] = Ys[i]*Ys[j] - Ys[j]*Ys[i]
    Bvec = Matrix.hstack(*[b.reshape(N*N, 1)
                           for b in brackets.values()])
    check(f"n={n} derived dim = n^2 - 1 = {n*n - 1}",
          Bvec.rank() == n * n - 1)
    murow = Matrix([[m for m in mus]])
    check(f"n={n} ker mu dim = n^2 - 1, derived subset ker mu",
          len(murow.nullspace()) == n * n - 1 and
          all(murow * ycoords(b) == zeros(1, 1)
              for b in brackets.values()))
    zrows = []
    for j in range(n * n):
        Mj = zeros(N*N, n * n)
        for i in range(n * n):
            Mj[:, i] = (Ys[i]*Ys[j] - Ys[j]*Ys[i]).reshape(N*N, 1)
        zrows.append(Mj)
    zspace = Matrix.vstack(*zrows).nullspace()
    check(f"n={n} centre dim 1", len(zspace) == 1)
    # Z candidate
    zcand = zeros(NP, 1)
    for k in range(n):
        zcand[pairs.index((2*k, 2*k+1))] = -1
    YZ = sum((zcand[j] * E[p] for j, p in enumerate(pairs)), zeros(N))
    AZ = (PaInv * YZ * Pa)[:n, :n]
    check(f"n={n} Z = -sum E_(2k,2k+1): A(Z) = -i I, central",
          eq(AZ, -I * eye(n)) and
          all(is_zero(YZ * Ys[j] - Ys[j] * YZ) for j in range(n * n)))
    rhoZ = sum((zcand[j] * SigM[p] for j, p in enumerate(pairs)),
               zeros(sec))
    def subset_of(idx):
        return frozenset(k for k in range(n)
                         if (idx >> (n - 1 - k)) & 1)
    pieces = {}
    for i2, idx in enumerate(minus_idx):
        pieces.setdefault(len(subset_of(idx)), []).append(i2)
    okch = True
    for m, idxs in pieces.items():
        if m % 2 == 1:
            okch = False
        for i2 in idxs:
            if rhoZ[i2, i2] != (m - Rational(n, 2)) * I:
                okch = False
    check(f"n={n} charges: piece m has rho(Z)/i = m - n/2, even m only",
          okch)
    check(f"n={n} piece dims = C(n, m)",
          all(len(pieces[m]) == len(list(combinations(range(n), m)))
              for m in pieces))
    # D. operator identity + weights on sparse basis
    eps = Symbol('eps')
    unit_cache = {}
    def dLam(m, A):
        if m not in unit_cache:
            subs = list(combinations(range(n), m))
            units = {}
            for p in range(n):
                for q in range(n):
                    Mx = eye(n) + eps * Matrix(
                        n, n, lambda a, b: 1 if (a, b) == (p, q) else 0)
                    Mm = zeros(len(subs))
                    for si, S in enumerate(subs):
                        for ti, T in enumerate(subs):
                            Mm[si, ti] = Mx.extract(
                                list(S), list(T)).det().coeff(eps, 1)
                    units[(p, q)] = Mm
            unit_cache[m] = (subs, units)
        subs, units = unit_cache[m]
        Mm = zeros(len(subs))
        for p in range(n):
            for q in range(n):
                if A[p, q] != 0:
                    Mm += A[p, q] * units[(p, q)]
        return subs, Mm
    def perm_for(m, subs):
        num = [subset_of(minus_idx[i2]) for i2 in pieces[m]]
        P = zeros(len(subs))
        for r, S in enumerate(num):
            P[r, subs.index(tuple(sorted(S)))] = 1
        return P
    okid = True
    for m in sorted(pieces):
        subs, _ = dLam(m, Ahat[0])
        P = perm_for(m, subs)
        for j in range(n * n):
            _, Mm = dLam(m, -Ahat[j].T)
            model = P * Mm * P.T + muS[j] * eye(len(subs))
            R = rhoS[j].extract(pieces[m], pieces[m])
            if not eq(R, model):
                okid = False
    check(f"n={n} rho(X) = dLambda^m(-A^T) + mu I, sparse basis, all pieces",
          okid)
    # weights test
    w = [2*k - (n - 1) for k in range(n)]
    Aw = zeros(n)
    for k in range(n):
        Aw[k, k] = I * w[k]
    Yhw = zeros(N); Yhw[:n, :n] = Aw; Yhw[n:, n:] = -Aw.T
    Yw = Pa * Yhw * PaInv
    okw = all(Im(Yw[r, c]) == 0 for r in range(N) for c in range(N))
    if okw:
        Ywr = Matrix(N, N, lambda r, c: Re(Yw[r, c]))
        cfw = ecoords(Ywr)
        Rw = sum((cfw[j] * SigM[p] for j, p in enumerate(pairs)),
                 zeros(sec))
        def skey(z):
            return (Re(z), Im(z))
        for m in sorted(pieces):
            R = Rw.extract(pieces[m], pieces[m])
            if any(R[a, b] != 0 for a in range(R.rows)
                   for b in range(R.cols) if a != b):
                okw = False
                break
            got = sorted([R[a, a] for a in range(R.rows)], key=skey)
            want = sorted([-I * sum(w[k] for k in S)
                           for S in combinations(range(n), m)], key=skey)
            if got != want:
                okw = False
    check(f"n={n} weights: piece m eigenvalues {{-i sum_S w}}", okw)
    # E. m, brackets, kappa, cubic
    rowsA = [[Yhat[p][r, s] for p in pairs]
             for r in range(n) for s in range(n)]
    mspace = real_nullspace_from_complex_rows(rowsA)
    Ms = [sum((cv[j] * E[p] for j, p in enumerate(pairs)), zeros(N))
          for cv in mspace]
    check(f"n={n} dim_R m = n(n-1) = {n*(n-1)}",
          len(Ms) == n * (n - 1))
    Vall = Matrix.hstack(*[Y.reshape(N*N, 1) for Y in Ys + Ms])
    pivA = Vall.T.rref()[1]
    VpA = Vall.extract(list(pivA), list(range(NP)))
    VpAInv = VpA.inv()
    def coords(Y):
        return VpAInv * Y.reshape(N*N, 1).extract(list(pivA), [0])
    nth = n * n
    ok_hm = ok_mm = True
    for i in range(nth):
        for j in range(len(Ms)):
            cf = coords(Ys[i]*Ms[j] - Ms[j]*Ys[i])
            if any(cf[t] != 0 for t in range(nth)):
                ok_hm = False
    for i in range(len(Ms)):
        for j in range(i+1, len(Ms)):
            cf = coords(Ms[i]*Ms[j] - Ms[j]*Ms[i])
            if any(cf[t] != 0 for t in range(nth, NP)):
                ok_mm = False
    check(f"n={n} [h, m] subset m", ok_hm)
    check(f"n={n} [m, m] subset h (symmetric)", ok_mm)
    def tprod(P, Q):
        return sum(x * y for x, y in zip(P, Q.T))
    Gk = Matrix(len(Ms), len(Ms),
                lambda j, k: -tprod(Ms[j], Ms[k]))
    check(f"n={n} kappa on m nondegenerate", Gk.det() != 0)
    check(f"n={n} kappa(h, m) = 0",
          all(-tprod(Ys[i], Ms[j]) == 0 for i in range(nth)
              for j in range(len(Ms))))
    okphi = True
    MM = [[Ms[j] * Ms[k] for k in range(len(Ms))] for j in range(len(Ms))]
    for j in range(len(Ms)):
        for k in range(j+1, len(Ms)):
            for l in range(k+1, len(Ms)):
                if -tprod(MM[j][k], Ms[l]) + tprod(MM[k][j], Ms[l]) != 0:
                    okphi = False
    check(f"n={n} cubic phi = 0 on all triples", okphi)
    # radical commutator: ad on n is dLambda^2(A) (engine input, all n)
    okrad = True
    for j in range(n * n):
        Dj = Hs[j]
        Aj = Ahat[j]
        for t in range(nu):
            Yh = zeros(N); Yh[:n, n:] = U[t]
            NB = Pa * Yh * PaInv
            comm = Dj * NB - NB * Dj
            Yh2 = PaInv * comm * Pa
            if not is_zero(Yh2[:n, :n]) or not is_zero(Yh2[n:, :]) or \
               not eq(Yh2[:n, n:], Aj * U[t] + U[t] * Aj.T):
                okrad = False
    check(f"n={n} ad on radical: B -> A B + B A^T (sparse basis)", okrad)
    if not engine:
        print(f"  SKIP  n={n} monomial engine (2^{nu} monomials) -- "
              f"covered by the n-uniform CAR derivation; inputs above")
        return
    # F. monomial engine
    print(f"  -- monomial engine: nu = {nu}, S dim {2**nu} --")
    NMASK = 2 ** nu
    def eps_op(t, mask):
        if (mask >> t) & 1:
            return None
        sgn = -1 if popc(mask & ((1 << t) - 1)) % 2 else 1
        return (mask | (1 << t), sgn)
    def iota_op(t, mask):
        if not ((mask >> t) & 1):
            return None
        sgn = -1 if popc(mask & ((1 << t) - 1)) % 2 else 1
        return (mask & ~(1 << t), sgn)
    def apply_lin(terms, mask):
        out = {}
        for (kind, t, c) in terms:
            r = eps_op(t, mask) if kind == 0 else iota_op(t, mask)
            if r is not None:
                m2, sgn = r
                v = c if sgn == 1 else (-c[0], -c[1])
                out[m2] = cadd(out[m2], v) if m2 in out else v
        return out
    rad = []
    for t in range(nu):
        Yh = zeros(N); Yh[:n, n:] = U[t]
        rad.append(Pa * Yh * PaInv)
    for t in range(nu):
        Yh = zeros(N); Yh[n:, :n] = U[t] * Rational(1, 2)
        rad.append(Pa * Yh * PaInv)
    GramC = Matrix(2*nu, 2*nu,
                   lambda a, b: -(rad[a] * rad[b]).trace())
    want = zeros(2*nu); want[:nu, nu:] = eye(nu); want[nu:, :nu] = eye(nu)
    check(f"n={n} kappa(m_C basis) = [[0, I], [I, 0]]",
          eq(GramC, want))
    gterms = [[(0, t, (Fraction(1), Fraction(0)))] for t in range(nu)] + \
             [[(1, t, (Fraction(1), Fraction(0)))] for t in range(nu)]
    def compose(F, terms, mask):
        acc = {}
        for m1, c1 in apply_lin(terms, mask).items():
            dacc(acc, F(m1), c1)
        return acc
    okcliff = True
    for a in range(2 * nu):
        for b in range(a, 2 * nu):
            kval = fp(GramC[a, b])
            for mask in range(NMASK):
                acc = compose(lambda m: apply_lin(gterms[a], m),
                              gterms[b], mask)
                dacc(acc, compose(lambda m: apply_lin(gterms[b], m),
                                  gterms[a], mask),
                     (Fraction(1), Fraction(0)))
                if clean(acc) != ({mask: kval} if kval != Z2 else {}):
                    okcliff = False
    check(f"n={n} Clifford relations on Lambda*(n), all monomials",
          okcliff)
    # sigma
    def dLam2(A):
        Mx = zeros(nu)
        for s, (a, b) in enumerate(SK):
            for t in range(nu):
                Mx[s, t] = (A * U[t] + U[t] * A.T)[a, b]
        return Mx
    def sigma_eval(Mx, trM, mask):
        acc = {}
        ht = (trM[0] / 2, trM[1] / 2)
        acc[mask] = (-ht[0], -ht[1])
        for t in range(nu):
            r = iota_op(t, mask)
            if r is None:
                continue
            m1, s1 = r
            for s in range(nu):
                c = Mx[s][t]
                if c == Z2:
                    continue
                r2 = eps_op(s, m1)
                if r2 is None:
                    continue
                m2, s2 = r2
                v = c if s1 * s2 == 1 else (-c[0], -c[1])
                acc[m2] = cadd(acc[m2], v) if m2 in acc else v
        return acc
    oklift = oktr = True
    for X, A in zip(Hs, Ahat):
        MxM = dLam2(A)
        Mx = [[fp(MxM[s, t]) for t in range(nu)] for s in range(nu)]
        trM = fp(MxM.trace())
        NM = zeros(nu)
        for s, (a, b) in enumerate(SK):
            for t in range(nu):
                NM[s, t] = (-A.T * U[t] - U[t] * A)[a, b]
        if not eq(NM, -MxM.T):
            oklift = False
        Nf = [[fp(NM[s, t]) for t in range(nu)] for s in range(nu)]
        sig = [sigma_eval(Mx, trM, mask) for mask in range(NMASK)]
        trsum = Z2
        for mask in range(NMASK):
            if mask in sig[mask]:
                trsum = cadd(trsum, sig[mask][mask])
        if trsum != Z2:
            oktr = False
        for t in range(nu):
            for kind in (0, 1):
                MxUse = Mx if kind == 0 else Nf
                gt = [(kind, t, (Fraction(1), Fraction(0)))]
                for mask in range(NMASK):
                    lhs = {}
                    for m1, c1 in apply_lin(gt, mask).items():
                        dacc(lhs, sig[m1], c1)
                    for m1, c1 in sig[mask].items():
                        d = apply_lin(gt, m1)
                        d = {m: (-v[0], -v[1]) for m, v in d.items()}
                        dacc(lhs, d, c1)
                    rhs = {}
                    for s in range(nu):
                        c = MxUse[s][t]
                        if c != Z2:
                            dacc(rhs, apply_lin(
                                [(kind, s, (Fraction(1), Fraction(0)))],
                                mask), c)
                    if clean(lhs) != clean(rhs):
                        oklift = False
    check(f"n={n} sigma intertwining, all sparse X, all monomials",
          oklift)
    check(f"n={n} tr sigma(X) = 0", oktr)
    # real-basis gamma + contraction
    Bterms = []
    for j in range(len(Ms)):
        Yh = PaInv * Ms[j] * Pa
        terms = []
        for t, (a, b) in enumerate(SK):
            cb = fp(Yh[a, n + b])
            if cb != Z2:
                terms.append((0, t, cb))
            cc = fp(Yh[n + a, b])
            if cc != Z2:
                terms.append((1, t, (2 * cc[0], 2 * cc[1])))
        Bterms.append(terms)
    GkInv = Gk.inv()
    Gf = [[fp(GkInv[k, j]) for j in range(len(Ms))]
          for k in range(len(Ms))]
    okcontr = True
    nuval = (Fraction(nu), Fraction(0))
    for mask in range(NMASK):
        acc = {}
        for j in range(len(Ms)):
            dual = {}
            for k in range(len(Ms)):
                if Gf[k][j] != Z2:
                    dacc(dual, apply_lin(Bterms[k], mask), Gf[k][j])
            for m1, c1 in dual.items():
                dacc(acc, apply_lin(Bterms[j], m1), c1)
        if clean(acc) != {mask: nuval}:
            okcontr = False
    check(f"n={n} contraction sum gamma gamma* = nu I = {nu} I",
          okcontr)

for n in range(2, 7):
    run_rung(n, engine=(n <= 5))

print()
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
