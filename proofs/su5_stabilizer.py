#!/usr/bin/env python3
"""
The SU(5) x U(1) stabilizer of the Spin(10) 16, built exactly.

Claim under test: the stabilizer in the COMPACT real form Spin(10) of
the pure spinor line C e_0 (the Fock vacuum of pure_spinor_orbit.py)
has Lie algebra u(5) = su(5) (+) u(1), dim 25 -- the compact form of
the parabolic's Levi GL(5, C) of the orbit build. Under it the 16
branches as Lambda^even of the CREATION space W' (which carries the
contragredient action -A^T, i.e. a 5bar): Lambda^0 (+) Lambda^2 (+)
Lambda^4 = 1 (+) 10bar (+) 5 as SU(5) representations (the last via
Lambda^4(5bar) = 5 (x) det, det trivial on SU(5)), with U(1) charges
(for the centre generator normalised by A = -i I on W) of -5/2, -1/2,
+3/2 on the three pieces (doubled: -5, -1, +3). The stabilizer of
the vector e_0 itself (not the line) is the su(5) part, dim 24.
Note the chirality sensitivity: the conjugate sector 16' branches as
the conjugate 1 (+) 10 (+) 5bar with charges +5/2, +1/2, -3/2 read
in the same normalisation (verified below on the odd pieces).

Method: exact arithmetic over Q and Q(i). The real stabilizer is
SOLVED (real nullspace of the W- and W'-preservation equations in the
adapted basis [W | W']), not assumed. Everything else -- block form,
anti-Hermiticity, brackets, centre, derived algebra, spin action,
branching operator identity, epsilon intertwining Lambda^4 = 5bar,
commutants -- is a finite exact matrix computation.

Controls (per symbolic_checker_validation):
  TRUE  control: {gamma_1, gamma_2} = 0 accepted; the plane-rotation
                 sum Y0 = sum_k E_{2k,2k+1} IS in the stabilizer.
  FALSE control: gamma_1 gamma_2 - gamma_2 gamma_1 rejected as 0;
                 E_{0,2} is NOT in the stabilizer (it moves W).
"""
import sys
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

def conjM(M):
    return Matrix(M.rows, M.cols,
                  lambda r, c: Re(M[r, c]) - I * Im(M[r, c]))

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
I32 = eye(32)

print("== checker controls ==")
check("TRUE control: {g1,g2} = 0 accepted", eq(G[0]*G[1] + G[1]*G[0], zeros(32)))
check("FALSE control: g1 g2 - g2 g1 rejected as 0",
      not is_zero(G[0]*G[1] - G[1]*G[0]))

prod = I32
for a in range(10):
    prod = prod * G[a]
Gamma = I * prod
minus_idx = [n for n in range(32) if Gamma[n, n] == -1]
pos0 = minus_idx.index(0)
pairs = [(a, b) for a in range(10) for b in range(a+1, 10)]
Sig = {p: (G[p[0]]*G[p[1]] - G[p[1]]*G[p[0]]) * Rational(1, 4) for p in pairs}
SigM = {p: Sig[p].extract(minus_idx, minus_idx) for p in pairs}
e0 = zeros(32, 1); e0[0] = 1

# W and W' bases (annihilator / creation directions), adapted basis
def avec(k, sign):
    v = zeros(10, 1)
    v[2*k] = Rational(1, 2)
    v[2*k+1] = sign * I * Rational(1, 2)
    return v
Wa = [avec(k, +1) for k in range(5)]
Wp = [avec(k, -1) for k in range(5)]
Pa = Matrix.hstack(*(Wa + Wp))
PaInv = Pa.inv()

# vector-rep generators
E = {}
for (a, b) in pairs:
    M = zeros(10); M[a, b] = 1; M[b, a] = -1
    E[(a, b)] = M
Yhat = {p: PaInv * E[p] * Pa for p in pairs}   # in adapted basis

print("== solve the real stabilizer ==")
# Unknowns c in Q^45 (REAL). Conditions: B-block and C-block of
# Yhat(c) vanish (preserve W' and W). 100 real equations (Re/Im).
rows = []
for (r, s) in [(r, s) for r in range(5) for s in range(5, 10)] + \
              [(r, s) for r in range(5, 10) for s in range(5)]:
    rows.append([Re(Yhat[p][r, s]) for p in pairs])
    rows.append([Im(Yhat[p][r, s]) for p in pairs])
Meq = Matrix(rows)
stab = Meq.nullspace()
check("real stabilizer (line, compact form) has dim 25", len(stab) == 25)

def in_span(basis_vecs, v):
    M = Matrix.hstack(*basis_vecs)
    return M.rank() == Matrix.hstack(M, v).rank()

y0 = zeros(45, 1)
for k in range(5):
    y0[pairs.index((2*k, 2*k+1))] = 1
check("TRUE control: Y0 = sum_k E_{2k,2k+1} is in the stabilizer",
      in_span(stab, y0))
e02 = zeros(45, 1); e02[pairs.index((0, 2))] = 1
check("FALSE control: E_{0,2} is NOT in the stabilizer",
      not in_span(stab, e02))

# per-element data
Ys, Ahat, rho16, mu = [], [], [], []
okB = okC = okD = okAH = okmu = okline = True
for cvec in stab:
    Y = sum((cvec[j] * E[p] for j, p in enumerate(pairs)), zeros(10))
    Yh = PaInv * Y * Pa
    A = Yh[:5, :5]
    if not is_zero(Yh[:5, 5:]):
        okB = False
    if not is_zero(Yh[5:, :5]):
        okC = False
    if not eq(Yh[5:, 5:], -A.T):
        okD = False
    if not is_zero(conjM(A).T + A):
        okAH = False
    R = sum((cvec[j] * SigM[p] for j, p in enumerate(pairs)), zeros(16))
    col = R[:, pos0]
    if any(col[r] != 0 for r in range(16) if r != pos0):
        okline = False
    m = col[pos0]
    if m != Rational(1, 2) * A.trace():
        okmu = False
    Ys.append(Y); Ahat.append(A); rho16.append(R); mu.append(m)
check("all 25 preserve W and W' (B = C = 0 blocks)", okB and okC)
check("D-block = -A^T on all 25 (so(10) form)", okD)
check("A-blocks anti-Hermitian: A* + A = 0 (the u(5) condition)", okAH)
check("spin action preserves the line C e_0 on all 25", okline)
check("mu = (1/2) tr(A) on all 25 (orbit-build character)", okmu)
# A-projection: real rank 25 inside the 25-real-dim space u(5)
Avec = []
for A in Ahat:
    Avec.append([Re(A[r, c]) for r in range(5) for c in range(5)] +
                [Im(A[r, c]) for r in range(5) for c in range(5)])
check("A-projection has real rank 25 => span_R{A} = u(5)",
      Matrix(Avec).rank() == 25)

print("== centre and derived algebra ==")
V = Matrix.hstack(*[Y.reshape(100, 1) for Y in Ys])  # 100x25, rank 25
piv = V.T.rref()[1]   # pivot cols of V^T = independent ROWS of V
Vp = V.extract(list(piv), list(range(25)))
VpInv = Vp.inv()
def express(Y):
    return VpInv * Y.reshape(100, 1).extract(list(piv), [0])
brackets = {}
for i in range(25):
    for j in range(i+1, 25):
        brackets[(i, j)] = Ys[i]*Ys[j] - Ys[j]*Ys[i]
Bvec = Matrix.hstack(*[b.reshape(100, 1) for b in brackets.values()])
check("derived algebra (bracket span) has dim 24", Bvec.rank() == 24)
# centre: z with [Y(z), Y_j] = 0 for all j
zrows = []
for j in range(25):
    Mj = zeros(100, 25)
    for i in range(25):
        comm = Ys[i]*Ys[j] - Ys[j]*Ys[i]
        Mj[:, i] = comm.reshape(100, 1)
    zrows.append(Mj)
Zmat = Matrix.vstack(*zrows)
zspace = Zmat.nullspace()
check("centre has dim 1", len(zspace) == 1)
# ker mu and equality with derived
murow = Matrix([[m for m in mu]])
kern = murow.nullspace()
check("ker mu has dim 24", len(kern) == 24)
der_in_ker = True
for b in brackets.values():
    if murow * express(b) != zeros(1, 1):
        der_in_ker = False
check("every bracket lies in ker mu", der_in_ker)
ker_vecs = Matrix.hstack(*[sum((k[t] * Ys[t] for t in range(25)),
                               zeros(10)).reshape(100, 1) for k in kern])
check("derived algebra = ker mu (rank of union = 24)",
      Matrix.hstack(Bvec, ker_vecs).rank() == 24)
# derived A-blocks: traceless + span sl(5) inside su(5)
derA = []
for b in brackets.values():
    cf = express(b)
    Ab = sum((cf[t] * Ahat[t] for t in range(25)), zeros(5))
    derA.append(Ab)
check("derived A-blocks all traceless (inside sl(5))",
      all(A.trace() == 0 for A in derA))
derAvec = []
for A in derA:
    derAvec.append([Re(A[r, c]) for r in range(5) for c in range(5)] +
                   [Im(A[r, c]) for r in range(5) for c in range(5)])
check("derived A-blocks have real rank 24 => they are su(5)",
      Matrix(derAvec).rank() == 24)

print("== the centre generator and its charges ==")
z0 = zspace[0]
Yz = sum((z0[t] * Ys[t] for t in range(25)), zeros(10))
Az = sum((z0[t] * Ahat[t] for t in range(25)), zeros(5))
lam = Az[0, 0]
check("centre A-block is a scalar matrix lam * I, lam != 0",
      eq(Az, lam * eye(5)) and lam != 0)
r = Im(lam)  # lam = i r, r rational != 0
scale = -1 / r
Zc = [scale * z0[t] for t in range(25)]       # now A(Z) = -i I
YZ = scale * Yz
AZ = -I * eye(5)
rhoZ16 = sum((Zc[t] * rho16[t] for t in range(25)), zeros(16))
rhoZ32 = sum((scale * z0[t] * sum((stab[t][j] * Sig[p]
             for j, p in enumerate(pairs)), zeros(32))
             for t in range(25)), zeros(32))
check("[rho(Y_j), rho(Z)] = 0 on the 16 for all j (spin-side centrality)",
      all(eq(rho16[j]*rhoZ16 - rhoZ16*rho16[j], zeros(16)) for j in range(25)))
check("rho(Z) is diagonal on the 32",
      all(rhoZ32[a, b] == 0 for a in range(32) for b in range(32) if a != b))
# charges by popcount (a-index k <-> bit (4-k) of the 32-index)
def subset_of(n):
    return frozenset(k for k in range(5) if (n >> (4 - k)) & 1)
pieces = {m: [i for i, n in enumerate(minus_idx)
              if len(subset_of(n)) == m] for m in (0, 2, 4)}
want = {0: Rational(-5, 2), 2: Rational(-1, 2), 4: Rational(3, 2)}
okcharge = True
for m in (0, 2, 4):
    for i in pieces[m]:
        n = minus_idx[i]
        if rhoZ32[n, n] != want[m] * I:
            okcharge = False
check("charges of rho(Z)/i on pieces m=0,2,4: -5/2, -1/2, +3/2", okcharge)
check("piece dims 1, 10, 5",
      (len(pieces[0]), len(pieces[2]), len(pieces[4])) == (1, 10, 5))

print("== branching operator identity ==")
eps = Symbol('eps')
# dLambda^m of gl(5) matrix units via epsilon-minors, LEX subset order
def dLam_units(m):
    subs = list(combinations(range(5), m))
    out = {}
    for p in range(5):
        for q in range(5):
            U = eye(5) + eps * Matrix(5, 5, lambda a, b: 1 if (a, b) == (p, q) else 0)
            M = zeros(len(subs))
            for si, S in enumerate(subs):
                for ti, T in enumerate(subs):
                    d = U.extract(list(S), list(T)).det()
                    M[si, ti] = d.coeff(eps, 1)
            out[(p, q)] = M
    return subs, out
def dLam(m, A, cache={}):
    if m not in cache:
        cache[m] = dLam_units(m)
    subs, units = cache[m]
    M = zeros(len(subs))
    for p in range(5):
        for q in range(5):
            if A[p, q] != 0:
                M += A[p, q] * units[(p, q)]
    return subs, M
# permutation: lex order -> numeric (sector) order within each piece
def perm_for(m, subs):
    num = [subset_of(minus_idx[i]) for i in pieces[m]]
    P = zeros(len(subs))
    for r, S in enumerate(num):
        P[r, subs.index(tuple(sorted(S)))] = 1
    return P
okid = True
for m in (0, 2, 4):
    subs, _ = dLam(m, Ahat[0])
    P = perm_for(m, subs)
    for j in range(25):
        _, Mm = dLam(m, -Ahat[j].T)
        model = P * Mm * P.T + mu[j] * eye(len(subs))
        R = rho16[j].extract(pieces[m], pieces[m])
        if not eq(R, model):
            okid = False
check("rho(X)|_piece = dLambda^m(-A(X)^T) + mu(X) I, all 25, m = 0,2,4",
      okid)
# control: dLambda^2 minors vs the explicit Lie formula on a sample A
A0 = Ahat[0]
subs2, M2 = dLam(2, A0)
okf = True
for si, (i1, i2) in enumerate(subs2):
    for ti, (k1, k2) in enumerate(subs2):
        f = (A0[i1, k1] * (1 if i2 == k2 else 0)
             - A0[i1, k2] * (1 if i2 == k1 else 0)
             - A0[i2, k1] * (1 if i1 == k2 else 0)
             + A0[i2, k2] * (1 if i1 == k1 else 0))
        if M2[si, ti] != f:
            okf = False
check("control: epsilon-minor dLambda^2 = explicit wedge formula", okf)

print("== Lambda^4 = 5 (not 5bar): epsilon intertwining ==")
subs4, _ = dLam(4, Ahat[0])
P4 = perm_for(4, subs4)
Peps = zeros(5)
for k in range(5):
    S = tuple(t for t in range(5) if t != k)
    perm = [k] + list(S)
    inv = sum(1 for a in range(5) for b in range(a+1, 5) if perm[a] > perm[b])
    Peps[k, subs4.index(S)] = (-1) ** inv
Peps_num = Peps * P4.T   # rows: dual index k; cols: numeric piece order
okE = True
for j in range(25):
    R = rho16[j].extract(pieces[4], pieces[4])
    if not eq(Peps_num * R,
              (Ahat[j] - Ahat[j].trace() * eye(5) + mu[j] * eye(5)) * Peps_num):
        okE = False
check("eps intertwines Lambda^4 with A + (mu - tr A) I: the defining "
      "action, scalar-twisted (on su(5): exactly the 5)", okE)

print("== weights: which of 5/5bar, 10/10bar the pieces are ==")
# Cartan element with A = diag(i, 2i, 3i, 4i, -10i): distinct weights,
# traceless, anti-Hermitian. Its spin eigenvalues on a piece decide
# the isomorphism type with no convention freedom.
aval = [I, 2*I, 3*I, 4*I, -10*I]
Ac = zeros(5)
for k in range(5):
    Ac[k, k] = aval[k]
Yhc = zeros(10); Yhc[:5, :5] = Ac; Yhc[5:, 5:] = -Ac.T
Yc = Pa * Yhc * PaInv
check("Cartan element is real (lies in the compact algebra)",
      all(Im(Yc[r, c]) == 0 for r in range(10) for c in range(10)))
VE = Matrix.hstack(*[E[p].reshape(100, 1) for p in pairs])
pivE = VE.T.rref()[1]
cfE = VE.extract(list(pivE), list(range(45))).inv() * \
    Yc.reshape(100, 1).extract(list(pivE), [0])
Rc = sum((cfE[j] * SigM[p] for j, p in enumerate(pairs)), zeros(16))
def eigs_on(m):
    R = Rc.extract(pieces[m], pieces[m])
    if any(R[a, b] != 0 for a in range(R.rows) for b in range(R.cols)
           if a != b):
        return None
    return sorted([R[a, a] for a in range(R.rows)],
                  key=lambda z: (Re(z), Im(z)))
from itertools import combinations as _comb
check("piece 4 weights = {a_k}: Lambda^4 is the defining 5",
      eigs_on(4) == sorted(aval, key=lambda z: (Re(z), Im(z))))
check("piece 2 weights = {-(a_i + a_j)}: Lambda^2 is the 10bar",
      eigs_on(2) == sorted([-(aval[i] + aval[j])
                            for i, j in _comb(range(5), 2)],
                           key=lambda z: (Re(z), Im(z))))
# conjugate sector 16': charges on the odd pieces under the same Z
plus_idx = [n for n in range(32) if Gamma[n, n] == 1]
def subset_of_n(n):
    return frozenset(k for k in range(5) if (n >> (4 - k)) & 1)
okodd = True
wantodd = {1: Rational(-3, 2), 3: Rational(1, 2), 5: Rational(5, 2)}
for n in plus_idx:
    m = len(subset_of_n(n))
    if rhoZ32[n, n] != wantodd[m] * I:
        okodd = False
check("16' charges under same Z: m=1,3,5 give -3/2, +1/2, +5/2 "
      "(conjugate branching)", okodd)

print("== irreducibility: Burnside span ==")
der_coefs = []
seen = []
for b in brackets.values():
    cf = express(b)
    der_coefs.append(cf)
# basis of derived in Y-coords: take independent rows of the bracket vecs
Bmat = Matrix.hstack(*[b.reshape(100, 1) for b in brackets.values()])
pivcols = list(Bmat.rref()[1])   # independent bracket columns
der_basis_coefs = []
keys = list(brackets.keys())
for pc in pivcols:
    i, j = keys[pc]
    der_basis_coefs.append(express(brackets[(i, j)]))
def rho_of(cf, m=None):
    R = sum((cf[t] * rho16[t] for t in range(25)), zeros(16))
    return R if m is None else R.extract(pieces[m], pieces[m])
def A_of(cf):
    return sum((cf[t] * Ahat[t] for t in range(25)), zeros(5))
def burnside_dim(mats):
    # Span of {I} u {M_t} u {M_t M_s}: if this is all of End(piece)
    # (dim d^2), the generated algebra is End(piece), which forces
    # irreducibility (any invariant subspace is End-stable, hence
    # 0 or everything). Length-2 products suffice when they span.
    d = mats[0].rows
    vecs = [eye(d).reshape(d * d, 1)]
    vecs += [M.reshape(d * d, 1) for M in mats]
    vecs += [(mats[t] * mats[s]).reshape(d * d, 1)
             for t in range(len(mats)) for s in range(len(mats))]
    return Matrix.hstack(*vecs).rank()
check("derived algebra generates End(W): defining 5 irreducible",
      burnside_dim([A_of(cf) for cf in der_basis_coefs]) == 25)
check("derived algebra generates End(10 piece): irreducible",
      burnside_dim([rho_of(cf, 2) for cf in der_basis_coefs]) == 100)
check("derived algebra generates End(5 piece): irreducible",
      burnside_dim([rho_of(cf, 4) for cf in der_basis_coefs]) == 25)
check("no invariant vector in the 10 piece (no singlet inside)",
      Matrix.vstack(*[rho_of(cf, 2) for cf in der_basis_coefs]).nullspace()
      .__len__() == 0)

print()
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
