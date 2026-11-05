#!/usr/bin/env python3
"""
Kostant's cubic Dirac operator for Spin(10)/U(5), built exactly.

Conventions fixed once: kappa(X, Y) = -tr(XY) on so(10); gamma is
the Clifford action for kappa|_m on the spinor module
S = Lambda*(n), n the B-block radical (dim_C n = 10, dim S = 1024);
D = sum_j X_j (x) gamma(X^{j*}), X_j a basis of m, X^{j*} its
kappa-dual basis. Kostant's D^2 formula (Casimir difference) is the
cited theorem; what is COMPUTED here, exactly, is every finite
input of the construction:

  1. kappa on m: Gram matrix nondegenerate; dual basis computed;
     kappa(h, m) = 0; kappa invariant under h on m.
  2. Clifford model: gamma on the m_C basis {e_t, e-bar_t}
     (wedge / contraction on Lambda*(n)) satisfies
     {gamma_a, gamma_b} = kappa_{ab} I on all 1024 monomials,
     with the Gram of that basis verified = [[0, I], [I, 0]] from
     the radical matrices. Diagonal Clifford relations for the
     real basis gamma(M_j) checked the same way (off-diagonal
     real relations follow by C-linearity from the C-basis ones).
  3. Spin lift sigma: h -> End(S), sigma(X) = sum M_pq a+_p a_q
     - (1/2) tr(M), M = ad(X)|_n = dLambda^2(A). Verified: the
     intertwining [sigma(X), gamma(v)] = gamma([X, v]) for all
     sparse-basis X and all m_C basis v, by evaluation on all
     1024 monomials; and tr sigma(X) = 0. (A lift satisfying the
     intertwining is a homomorphism: the defect commutes with all
     of Cl(m), hence is scalar -- Cl(m_C) = End(S) in even
     dimension, standard -- and trace 0 kills the scalar.)
  4. gamma-contraction: sum_j gamma(M_j) gamma(M^{j*}) = 10 I
     on all monomials (the canonical metric contraction).
  5. The cubic term VANISHES: phi(X, Y, Z) = kappa([X, Y], Z) = 0
     for all m basis triples j<k<l (the rest by antisymmetry).
     Reason visible in the computation: [m, m] subset h (bundle
     build) and kappa(h, m) = 0 (item 1). So for this Hermitian
     symmetric pair Kostant's operator is the plain Dirac sum.

Controls (per symbolic_checker_validation):
  TRUE  control: {gamma_1, gamma_2} = 0 accepted; kappa(e_t,
                 e-bar_t) = 1 for the radical pairing.
  FALSE control: gamma_1 gamma_2 - gamma_2 gamma_1 rejected as 0;
                 kappa(e_1, e_2) rejected as 1 (it is 0).
"""
import sys
from fractions import Fraction
from itertools import combinations
from sympy import (Matrix, eye, zeros, I, Rational, kronecker_product,
                   re as Re, im as Im)

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
G = gam

print("== checker controls ==")
check("TRUE control: {g1,g2} = 0 accepted", eq(G[0]*G[1] + G[1]*G[0], zeros(32)))
check("FALSE control: g1 g2 - g2 g1 rejected as 0",
      not is_zero(G[0]*G[1] - G[1]*G[0]))

pairs = [(a, b) for a in range(10) for b in range(a+1, 10)]
def avec(k, sign):
    v = zeros(10, 1)
    v[2*k] = Rational(1, 2)
    v[2*k+1] = sign * I * Rational(1, 2)
    return v
Pa = Matrix.hstack(*([avec(k, +1) for k in range(5)] +
                     [avec(k, -1) for k in range(5)]))
PaInv = Pa.inv()
E = {}
for (a, b) in pairs:
    M = zeros(10); M[a, b] = 1; M[b, a] = -1
    E[(a, b)] = M
Yhat = {p: PaInv * E[p] * Pa for p in pairs}

def real_nullspace(eqs):
    rows = []
    for row in eqs:
        rows.append([Re(x) for x in row])
        rows.append([Im(x) for x in row])
    return Matrix(rows).nullspace()

print("== h (solved) and the sparse u(5) basis ==")
eqs = []
for (r, s) in [(r, s) for r in range(5) for s in range(5, 10)] + \
              [(r, s) for r in range(5, 10) for s in range(5)]:
    eqs.append([Yhat[p][r, s] for p in pairs])
stab = real_nullspace(eqs)
Ys = [sum((cv[j] * E[p] for j, p in enumerate(pairs)), zeros(10))
      for cv in stab]
check("dim h = 25 (solved)", len(stab) == 25)
Vh = Matrix.hstack(*[Y.reshape(100, 1) for Y in Ys])
def in_span_h(Y):
    return Vh.rank() == Matrix.hstack(Vh, Y.reshape(100, 1)).rank()
# sparse basis: A-units of u(5)
Aunits = []
for p in range(5):
    A = zeros(5); A[p, p] = I; Aunits.append(A)
for p in range(5):
    for q in range(p+1, 5):
        A = zeros(5); A[p, q] = 1; A[q, p] = -1; Aunits.append(A)
        A = zeros(5); A[p, q] = I; A[q, p] = I; Aunits.append(A)
Hs, Ahat = [], []
ok_sparse = True
for A in Aunits:
    Yh = zeros(10); Yh[:5, :5] = A; Yh[5:, 5:] = -A.T
    Y = Pa * Yh * PaInv
    if any(Im(Y[r, c]) != 0 for r in range(10) for c in range(10)):
        ok_sparse = False
    Yr = Matrix(10, 10, lambda r, c: Re(Y[r, c]))
    if not in_span_h(Yr):
        ok_sparse = False
    if not eq((PaInv * Yr * Pa)[:5, :5], A):
        ok_sparse = False
    Hs.append(Yr); Ahat.append(A)
check("25 sparse A-unit elements: real, in h, A recovered", ok_sparse)

print("== m, kappa, duals ==")
eqsA = [[Yhat[p][r, s] for p in pairs] for r in range(5) for s in range(5)]
mspace = real_nullspace(eqsA)
Ms = [sum((cv[j] * E[p] for j, p in enumerate(pairs)), zeros(10))
      for cv in mspace]
check("dim_R m = 20", len(Ms) == 20)
MM = [[Ms[j] * Ms[k] for k in range(20)] for j in range(20)]
def tprod(P, Q):
    """tr(P Q) entrywise, without forming the product."""
    return sum(x * y for x, y in zip(P, Q.T))
Gk = Matrix(20, 20, lambda j, k: -tprod(Ms[j], Ms[k]))
check("kappa Gram on m nondegenerate (det != 0)", Gk.det() != 0)
GkInv = Gk.inv()
check("dual basis: kappa(M_j, M^{k*}) = delta_jk",
      eq(Gk * GkInv, eye(20)))
check("kappa(h, m) = 0 on solved bases (25 x 20)",
      all(-tprod(Ys[i], Ms[j]) == 0 for i in range(25)
          for j in range(20)))
okinv = True
for X in Hs:
    for j in range(20):
        for k in range(j+1, 20):
            # tr([X,M_j] M_k) + tr(M_j [X,M_k])
            #   = -tr(X M_j M_k) + tr(M_j M_k X), the negative of
            #   what is tested; vanishing is tr(X P) = tr(P X),
            #   P = M_j M_k (cyclicity)
            if tprod(X, MM[j][k]) - tprod(MM[j][k], X) != 0:
                okinv = False
check("kappa invariant under h on m (sparse basis, all pairs)", okinv)

print("== cubic form vanishes ==")
okphi = True
for j in range(20):
    for k in range(j+1, 20):
        for l in range(k+1, 20):
            v = -tprod(MM[j][k], Ms[l]) + tprod(MM[k][j], Ms[l])
            if v != 0:
                okphi = False
check("phi(M_j, M_k, M_l) = kappa([M_j,M_k], M_l) = 0, all j<k<l",
      okphi)

print("== the m_C basis and its kappa Gram ==")
SK = list(combinations(range(5), 2))
U = {}
for t, (a, b) in enumerate(SK):
    M = zeros(5); M[a, b] = 1; M[b, a] = -1
    U[t] = M
rad = []   # 20 radical matrices: e_t (B-type), e-bar_t (C = U/2,
           # the kappa-dual of e_t: tr(U * U/2) = -1, kappa = +1)
for t in range(10):
    Yh = zeros(10); Yh[:5, 5:] = U[t]
    rad.append(Pa * Yh * PaInv)
for t in range(10):
    Yh = zeros(10); Yh[5:, :5] = U[t] * Rational(1, 2)
    rad.append(Pa * Yh * PaInv)
GramC = Matrix(20, 20, lambda a, b: -(rad[a] * rad[b]).trace())
want = zeros(20); want[:10, 10:] = eye(10); want[10:, :10] = eye(10)
check("TRUE control + Gram: kappa(m_C basis) = [[0, I], [I, 0]]",
      eq(GramC, want))
check("FALSE control: kappa(e_1, e_2) = 0, not 1", GramC[0, 1] == 0)

# ------------------------------------------------- sparse evaluation engine
Z2 = (Fraction(0), Fraction(0))
def fp(z):
    """sympy gaussian rational -> (Fraction, Fraction)."""
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
    """terms: list of (kind, t, coeff-pair); kind 0 = eps, 1 = iota."""
    out = {}
    for (kind, t, c) in terms:
        r = eps_op(t, mask) if kind == 0 else iota_op(t, mask)
        if r is not None:
            m2, sgn = r
            v = c if sgn == 1 else (-c[0], -c[1])
            out[m2] = cadd(out[m2], v) if m2 in out else v
    return out
def compose(F, Gterms, mask):
    """F: dict-valued evaluator; apply G then F."""
    acc = {}
    for m1, c1 in apply_lin(Gterms, mask).items():
        dacc(acc, F(m1), c1)
    return acc
NMASK = 1024
check("monomial count 1024, even/odd split 512/512",
      NMASK == 1024 and sum(1 for m in range(1024)
                             if popc(m) % 2 == 0) == 512)

print("== Clifford relations on the m_C basis (evaluation) ==")
gterms = []
for t in range(10):
    gterms.append([(0, t, (Fraction(1), Fraction(0)))])       # gamma(e_t)
for t in range(10):
    gterms.append([(1, t, (Fraction(1), Fraction(0)))])       # gamma(ebar_t)
def anticomm_ok(a, b, kval):
    for mask in range(NMASK):
        acc = compose(lambda m: apply_lin(gterms[a], m), gterms[b], mask)
        dacc(acc, compose(lambda m: apply_lin(gterms[b], m),
                          gterms[a], mask), (Fraction(1), Fraction(0)))
        want = {mask: kval} if kval != Z2 else {}
        if clean(acc) != want:
            return False
    return True
okcliff = True
for a in range(20):
    for b in range(a, 20):
        kval = fp(GramC[a, b])
        if not anticomm_ok(a, b, kval):
            okcliff = False
check("{gamma_a, gamma_b} = kappa_ab I, all a <= b, all monomials",
      okcliff)

print("== real-basis gamma: diagonal Clifford ==")
Bterms = []
for j in range(20):
    Yh = PaInv * Ms[j] * Pa
    terms = []
    for t, (a, b) in enumerate(SK):
        cb = fp(Yh[a, 5 + b])
        if cb != Z2:
            terms.append((0, t, cb))
        cc = fp(Yh[5 + a, b])
        if cc != Z2:
            terms.append((1, t, (2 * cc[0], 2 * cc[1])))
    Bterms.append(terms)
okdiag = True
for j in range(20):
    kval = fp(Gk[j, j])
    for mask in range(NMASK):
        acc = compose(lambda m: apply_lin(Bterms[j], m), Bterms[j], mask)
        acc = {m: (2 * c[0], 2 * c[1]) for m, c in acc.items()}
        want = {mask: kval} if kval != Z2 else {}
        if clean(acc) != want:
            okdiag = False
check("{gamma(M_j), gamma(M_j)} = 2 gamma^2 = kappa_jj I (diagonal)",
      okdiag)

print("== spin lift sigma: intertwining + trace ==")
def dLam2(A):
    Mx = zeros(10)
    for s, (a, b) in enumerate(SK):
        for t, (c, d) in enumerate(SK):
            Mx[s, t] = (A * U[t] + U[t] * A.T)[a, b]
    return Mx
def sigma_eval(Mx, trM, mask):
    acc = {}
    half_tr = (trM[0] / 2, trM[1] / 2)
    acc[mask] = (-half_tr[0], -half_tr[1])
    for t in range(10):
        r = iota_op(t, mask)
        if r is None:
            continue
        m1, s1 = r
        for s in range(10):
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
oklift = True
oktr = True
for X, A in zip(Hs, Ahat):
    MxM = dLam2(A)
    Mx = [[fp(MxM[s, t]) for t in range(10)] for s in range(10)]
    trM = fp(MxM.trace())
    # N = action on e-bar basis; consistency N = -M^T
    NM = zeros(10)
    for s, (a, b) in enumerate(SK):
        for t, (c, d) in enumerate(SK):
            NM[s, t] = (-A.T * U[t] - U[t] * A)[a, b]
    if not eq(NM, -MxM.T):
        oklift = False
    Nf = [[fp(NM[s, t]) for t in range(10)] for s in range(10)]
    sig = [sigma_eval(Mx, trM, mask) for mask in range(NMASK)]
    trsum = Z2
    for mask in range(NMASK):
        if mask in sig[mask]:
            trsum = cadd(trsum, sig[mask][mask])
    if trsum != Z2:
        oktr = False
    for t in range(10):
        for kind in (0, 1):
            MxUse = Mx if kind == 0 else Nf
            gt = [(kind, t, (Fraction(1), Fraction(0)))]
            for mask in range(NMASK):
                # [sigma, gamma_t] (mask) vs gamma([X, v_t]) (mask)
                lhs = {}
                for m1, c1 in apply_lin(gt, mask).items():
                    dacc(lhs, sig[m1], c1)
                for m1, c1 in sig[mask].items():
                    d = apply_lin(gt, m1)
                    d = {m: (-v[0], -v[1]) for m, v in d.items()}
                    dacc(lhs, d, c1)
                rhs = {}
                for s in range(10):
                    c = MxUse[s][t]
                    if c != Z2:
                        dacc(rhs, apply_lin(
                            [(kind, s, (Fraction(1), Fraction(0)))], mask), c)
                if clean(lhs) != clean(rhs):
                    oklift = False
check("[sigma(X), gamma(v)] = gamma([X, v]), all sparse X, all v",
      oklift)
check("tr sigma(X) = 0 on S for all sparse X", oktr)

print("== gamma-contraction ==")
def gamma_real_eval(j, mask):
    return apply_lin(Bterms[j], mask)
okcontr = True
ten = (Fraction(10), Fraction(0))
Gf = [[fp(GkInv[k, j]) for j in range(20)] for k in range(20)]
for mask in range(NMASK):
    acc = {}
    for j in range(20):
        dual = {}
        for k in range(20):
            if Gf[k][j] != Z2:
                dacc(dual, gamma_real_eval(k, mask), Gf[k][j])
        for m1, c1 in dual.items():
            dacc(acc, gamma_real_eval(j, m1), c1)
    if clean(acc) != {mask: ten}:
        okcontr = False
check("sum_j gamma(M_j) gamma(M^{j*}) = 10 I on all monomials",
      okcontr)

print()
print("D assembled as: D = sum_j M_j (x) gamma(M^{j*}) in "
      "U(g) (x) Cl(m); cubic term 0 (verified); D^2 Casimir "
      "formula is Kostant's theorem (cited).")
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
