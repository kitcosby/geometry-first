#!/usr/bin/env python3
"""
Phi: the intertwiner between the two models of the 27 of E6 --
the representation-theoretic model of e6_rep27.py (V27 = S + V + 1,
Pauli/Fock spinor basis) and the Jordan model of e6_jordan.py /
e6_jordan_e6.py (J = H_3(O_C), Peirce pieces J1, J1/2, J0).
All arithmetic exact over Q(i) (Fraction pairs), except the rep27
front-end rebuild, which follows e6_rep27.py in sympy and is then
converted.

Construction outline:
  P1  both e6 actions assembled as 78 dict-operators.
  P2  a regular torus element h: on the rep side
      h_R = sum a_k Y_(2k,2k+1), a = (1, 2, 4, 8, -15) (distinct
      weights on V and on S, S-block diagonal); on the Jordan
      side h_J built from a Peirce-frame torus (boost on (d1,d2),
      rotations on four C-planes, solved inside the derived Levi
      by J0-restriction).  Spinor spectra compared as multisets;
      the parity flip (a4 -> -a4) applied if needed.
  P3  Phi0: eigenspace matching for h (Phi_V by explicit
      eigenvectors, Phi_S by exact nullspaces, Phi_1: e26 -> E).
  P4  beta0 on the Levi: Z' -> Zhat_n (charges), torus by weight
      matching, root vectors by [torus, X] = delta X solves in
      the derived algebra; their scalars t_g fixed by the bracket
      relations (multiplicative propagation); Phi's weight
      scalars tau fixed by the root intertwining ratios.
  P5  p-stage: alpha(T_x0), alpha(T_y0) and the remaining piece
      scales from the intertwining systems (bilinear products
      substituted to linear ones); gauges tau^1 = 1, mu = 1;
      alpha on the rest of p+- by Levi-bracket generation.
  P6  FULL VERIFICATION: Phi rho_R(g) = rho_J(alpha(g)) Phi for
      all 78 generators as operators; alpha a Lie isomorphism on
      all 3081 bracket pairs; det o Phi = c * C as cubic tensors;
      Phi invertible.
"""
from fractions import Fraction
from itertools import combinations
import sympy as sp

FAIL = []
NCHK = [0]

def ok(name, cond, detail=""):
    NCHK[0] += 1
    if not cond:
        FAIL.append((name, detail))
        print(f"FAIL {name} {detail}", flush=True)
    else:
        print(f"  ok {name}", flush=True)

# ================================================================ Q(i) pairs
Z2 = (Fraction(0), Fraction(0))
ONE2 = (Fraction(1), Fraction(0))
def padd(a, b): return (a[0] + b[0], a[1] + b[1])
def pneg(a): return (-a[0], -a[1])
def psub(a, b): return (a[0] - b[0], a[1] - b[1])
def pmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def pscale(a, f): return (a[0] * f, a[1] * f)
def pinv(a):
    d = a[0] * a[0] + a[1] * a[1]
    return (a[0] / d, -a[1] / d)
def pdiv(a, b): return pmul(a, pinv(b))
def phas(a): return a != Z2
def pkey(a): return (a[0], a[1])

# ================================================================ octonions
LINES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6),
         (2, 5, 7), (3, 4, 7), (3, 6, 5)]
OT = {}
for i in range(8):
    OT[(0, i)] = (i, 1)
    OT[(i, 0)] = (i, 1)
for i in range(1, 8):
    OT[(i, i)] = (0, -1)
for (a, b, c) in LINES:
    OT[(a, b)] = (c, 1); OT[(b, c)] = (a, 1); OT[(c, a)] = (b, 1)
    OT[(b, a)] = (c, -1); OT[(c, b)] = (a, -1); OT[(a, c)] = (b, -1)

def omul(A, B):
    out = {}
    for i, ci in A.items():
        for j, cj in B.items():
            k, s = OT[(i, j)]
            v = pmul(ci, cj)
            if s < 0:
                v = pneg(v)
            out[k] = padd(out.get(k, Z2), v)
            if not phas(out[k]):
                del out[k]
    return out

def oconj(A):
    return {i: (c if i == 0 else pneg(c)) for i, c in A.items()}

# ================================================================ J coords
def entry(X, i, j):
    if i == j:
        v = X.get(i)
        return {0: v} if v is not None else {}
    if (i, j) == (0, 1):
        return {u: X[3 + u] for u in range(8) if (3 + u) in X}
    if (i, j) == (1, 0):
        return oconj(entry(X, 0, 1))
    if (i, j) == (0, 2):
        return {u: X[11 + u] for u in range(8) if (11 + u) in X}
    if (i, j) == (2, 0):
        return oconj(entry(X, 0, 2))
    if (i, j) == (1, 2):
        return {u: X[19 + u] for u in range(8) if (19 + u) in X}
    if (i, j) == (2, 1):
        return oconj(entry(X, 1, 2))
    raise ValueError

def coords_of(M):
    out = {}
    for i in range(3):
        d = M[i][i]
        if any(u != 0 for u in d):
            return None
        if 0 in d:
            out[i] = d[0]
    for base, (i, j) in ((3, (0, 1)), (11, (0, 2)), (19, (1, 2))):
        for u, c in M[i][j].items():
            out[base + u] = c
    return out

def jprod(X, Y):
    E = [[entry(X, i, j) for j in range(3)] for i in range(3)]
    F = [[entry(Y, i, j) for j in range(3)] for i in range(3)]
    half = Fraction(1, 2)
    M = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            acc = {}
            for k in range(3):
                for P, Q in ((E, F), (F, E)):
                    if P[i][k] and Q[k][j]:
                        for u, c in omul(P[i][k], Q[k][j]).items():
                            acc[u] = padd(acc.get(u, Z2), c)
            M[i][j] = {u: pscale(c, half) for u, c in acc.items()
                       if phas(c)}
    return coords_of(M)

def coord_add(X, Y):
    out = dict(X)
    for k, v in Y.items():
        out[k] = padd(out.get(k, Z2), v)
        if not phas(out[k]):
            del out[k]
    return out

BASIS = [{i: ONE2} for i in range(27)]

def slot_oct(X, base):
    return {u: X[base + u] for u in range(8) if (base + u) in X}

def onorm0(A):
    return omul(A, oconj(A)).get(0, Z2)

def detJ(X):
    d = [X.get(i, Z2) for i in range(3)]
    A = slot_oct(X, 3)
    B = slot_oct(X, 11)
    C = slot_oct(X, 19)
    v = pmul(pmul(d[0], d[1]), d[2])
    v = psub(v, pmul(d[0], onorm0(C)))
    v = psub(v, pmul(d[1], onorm0(B)))
    v = psub(v, pmul(d[2], onorm0(A)))
    t = omul(omul(A, C), oconj(B)).get(0, Z2)
    return padd(v, pscale(t, Fraction(2)))

def det3J(I, J, K):
    eI, eJ, eK = BASIS[I], BASIS[J], BASIS[K]
    s = detJ(coord_add(coord_add(eI, eJ), eK))
    s = psub(s, detJ(coord_add(eI, eJ)))
    s = psub(s, detJ(coord_add(eI, eK)))
    s = psub(s, detJ(coord_add(eJ, eK)))
    s = padd(s, detJ(eI))
    s = padd(s, detJ(eJ))
    s = padd(s, detJ(eK))
    return pscale(s, Fraction(1, 6))

print("building det tensor...", flush=True)
TRIPLES = []
for i in range(27):
    TRIPLES.append((i, i, i))
    for j in range(i + 1, 27):
        TRIPLES.append((i, i, j))
        TRIPLES.append((i, j, j))
for i, j, k in combinations(range(27), 3):
    TRIPLES.append((i, j, k))
TDET = {}
for t in TRIPLES:
    v = det3J(*t)
    if phas(v):
        TDET[t] = v
assert len(TDET) == 89, len(TDET)

def tdet3(a, b, c):
    return TDET.get(tuple(sorted((a, b, c))), Z2)

PAIRL = {}
for pat, v in TDET.items():
    for L in set(pat):
        rem = list(pat)
        rem.remove(L)
        PAIRL.setdefault(tuple(sorted(rem)), []).append((L, v))

CHARGE = {}
CHARGE[0] = 4
for idx in range(3, 19):
    CHARGE[idx] = 1
for idx in [1, 2] + list(range(19, 27)):
    CHARGE[idx] = -2

# ================================================================ sparse GE
def nullspace(eqs, nvars):
    occ = [0] * nvars
    for e in eqs:
        for v in e:
            occ[v] += 1
    pivots = {}
    for e0 in eqs:
        e = dict(e0)
        while True:
            pv = [v for v in e if v in pivots]
            if not pv:
                break
            p = min(pv, key=lambda v: occ[v])
            f = e[p]
            pe = pivots[p]
            for v, c in pe.items():
                nv = psub(e.get(v, Z2), pmul(f, c))
                if phas(nv):
                    e[v] = nv
                elif v in e:
                    del e[v]
        if not e:
            continue
        p = min(e, key=lambda v: occ[v])
        inv = pinv(e[p])
        e = {v: pmul(c, inv) for v, c in e.items()}
        for q, qe in list(pivots.items()):
            if p in qe:
                f = qe[p]
                ne = dict(qe)
                for v, c in e.items():
                    nv = psub(ne.get(v, Z2), pmul(f, c))
                    if phas(nv):
                        ne[v] = nv
                    elif v in ne:
                        del ne[v]
                pivots[q] = ne
        pivots[p] = e
    free = [v for v in range(nvars) if v not in pivots]
    basis = []
    for f in free:
        x = {f: ONE2}
        for p, pe in pivots.items():
            if f in pe:
                x[p] = pneg(pe[f])
        basis.append(x)
    return basis

def rank_of_vectors(vecs, nvars):
    pivots = {}
    for v0 in vecs:
        e = dict(v0)
        while True:
            pv = [v for v in e if v in pivots]
            if not pv:
                break
            p = pv[0]
            f = e[p]
            pe = pivots[p]
            for v, c in pe.items():
                nv = psub(e.get(v, Z2), pmul(f, c))
                if phas(nv):
                    e[v] = nv
                elif v in e:
                    del e[v]
        if e:
            p = min(e)
            inv = pinv(e[p])
            pivots[p] = {v: pmul(c, inv) for v, c in e.items()}
    return len(pivots), sorted(pivots)

# ================================================================ sectors
SECTORS = {}
for qT in (0, 3, -3, 6, -6):
    varlist = [(r, c) for r in range(27) for c in range(27)
               if CHARGE[r] - CHARGE[c] == qT]
    vidx = {v: n for n, v in enumerate(varlist)}
    eqs = []
    for (I, J, K) in TRIPLES:
        if CHARGE[I] + CHARGE[J] + CHARGE[K] != -qT:
            continue
        e = {}
        for (col, x, y) in ((I, J, K), (J, I, K), (K, I, J)):
            for (L, v) in PAIRL.get(tuple(sorted((x, y))), []):
                if (L, col) in vidx:
                    key = vidx[(L, col)]
                    e[key] = padd(e.get(key, Z2), v)
        e = {k: v for k, v in e.items() if phas(v)}
        if e:
            eqs.append(e)
    sol = nullspace(eqs, len(varlist))
    SECTORS[qT] = (varlist, sol)
    print(f"  sector {qT:+d}: nullity {len(sol)}", flush=True)
assert [len(SECTORS[q][1]) for q in (0, 3, -3, 6, -6)] == \
    [46, 16, 16, 0, 0]

OPS = []
for qT in (0, 3, -3):
    varlist, sol = SECTORS[qT]
    for x in sol:
        M = {}
        for vi, val in x.items():
            M[varlist[vi]] = val
        OPS.append((M, qT))
NB = len(OPS)
assert NB == 78
LEVI = [M for (M, q) in OPS if q == 0]
NL = len(LEVI)
PmJ = [M for (M, q) in OPS if q == -3]   # module charge -3: J1->J1/2->J0
PpJ = [M for (M, q) in OPS if q == 3]
SEC0_IDX = [a for a in range(NB) if OPS[a][1] == 0]

# ================================================================ matrices
def matmul(A, B):
    rows = {}
    for (r, c), v in B.items():
        rows.setdefault(r, []).append((c, v))
    out = {}
    for (r, m), v in A.items():
        if m in rows:
            for c, w in rows[m]:
                key = (r, c)
                out[key] = padd(out.get(key, Z2), pmul(v, w))
                if not phas(out[key]):
                    del out[key]
    return out

def commutator(A, B):
    AB = matmul(A, B)
    BA = matmul(B, A)
    keys = set(AB) | set(BA)
    return {k: psub(AB.get(k, Z2), BA.get(k, Z2)) for k in keys
            if phas(psub(AB.get(k, Z2), BA.get(k, Z2)))}

def mat_scale(A, c):
    if not phas(c):
        return {}
    return {k: pmul(v, c) for k, v in A.items()}

def mat_add(*Ms):
    out = {}
    for M in Ms:
        for k, v in M.items():
            out[k] = padd(out.get(k, Z2), v)
            if not phas(out[k]):
                del out[k]
    return out

def lincomb(ops, cf):
    out = {}
    for a, c in cf.items():
        if phas(c):
            out = mat_add(out, mat_scale(ops[a], c))
    return out

def apply_op(T, v):
    """T dict op, v coord dict (27) -> coord dict."""
    out = {}
    for (r, c), val in T.items():
        if c in v:
            out[r] = padd(out.get(r, Z2), pmul(val, v[c]))
    return {k: x for k, x in out.items() if phas(x)}

def to_sympy(x):
    return sp.Rational(x[0].numerator, x[0].denominator) + \
        sp.I * sp.Rational(x[1].numerator, x[1].denominator)

def from_sympy(z):
    z = sp.expand(z)
    re = sp.re(z)
    im = sp.im(z)
    return (Fraction(re.p, re.q), Fraction(im.p, im.q))

def op_vec(M):
    return {r * 27 + c: v for (r, c), v in M.items()}

def make_express(ops):
    n = len(ops)
    vecs = [op_vec(M) for M in ops]
    _, pos = rank_of_vectors(vecs, 729)
    assert len(pos) == n
    MINOR = sp.Matrix(n, n, lambda k, a: to_sympy(
        ops[a].get((pos[k] // 27, pos[k] % 27), Z2)))
    MINV = MINOR.inv()
    def express(M):
        rhs = sp.Matrix(n, 1, lambda k, _: to_sympy(
            M.get((pos[k] // 27, pos[k] % 27), Z2)))
        sol = MINV * rhs
        coeffs = [from_sympy(sol[a]) for a in range(n)]
        rec = {}
        for a, c in enumerate(coeffs):
            if phas(c):
                rec = mat_add(rec, mat_scale(ops[a], c))
        return coeffs if rec == M else None
    return express

print("building minors (78 and Levi)...", flush=True)
express78 = make_express([M for (M, q) in OPS])

# Levi centre and derived-45 (as in e6_jordan_e6)
print("Levi commutators...", flush=True)
SCOMM = {}
for a in range(NL):
    for b in range(a + 1, NL):
        SCOMM[(a, b)] = commutator(LEVI[a], LEVI[b])

def scomm(a, b):
    if a == b:
        return {}
    if a < b:
        return SCOMM[(a, b)]
    return {k: pneg(v) for k, v in SCOMM[(b, a)].items()}

zeqs = []
for b in range(NL):
    entries = {}
    for a in range(NL):
        for pos, v in scomm(a, b).items():
            entries.setdefault(pos, {})[a] = v
    zeqs.extend(entries.values())
zcen = nullspace(zeqs, NL)
assert len(zcen) == 1, len(zcen)
ZHAT = {}
for ai, val in zcen[0].items():
    ZHAT = mat_add(ZHAT, mat_scale(LEVI[ai], val))
s_prime = pscale(ZHAT.get((0, 0), Z2), Fraction(1, 4))
assert phas(s_prime)
ZHATN = mat_scale(ZHAT, pinv(s_prime))     # block scalars (4, 1, -2)
ok("P1a Jordan side rebuilt: 78 ops; Zhat normalized to "
   "(4, 1, -2)", ZHATN.get((0, 0)) == pscale(ONE2, Fraction(4)))

piv_rows = []
pivots45 = {}
for a in range(NL):
    for b in range(a + 1, NL):
        vec = op_vec(SCOMM[(a, b)])
        e = dict(vec)
        while True:
            pv = [v for v in e if v in pivots45]
            if not pv:
                break
            p = pv[0]
            f = e[p]
            pe = pivots45[p]
            for v, c in pe.items():
                nv = psub(e.get(v, Z2), pmul(f, c))
                if phas(nv):
                    e[v] = nv
                elif v in e:
                    del e[v]
        if e:
            p = min(e)
            inv = pinv(e[p])
            pivots45[p] = {v: pmul(c, inv) for v, c in e.items()}
            piv_rows.append(SCOMM[(a, b)])
        if len(piv_rows) == 45:
            break
    if len(piv_rows) == 45:
        break
DERJ = piv_rows
assert len(DERJ) == 45
expressDER = make_express(DERJ)
print("  derived-45 basis and minor ready", flush=True)

# ================================================================ rep27 side
print("rebuilding rep27 front-end (sympy)...", flush=True)
from sympy import Matrix, eye, zeros, I as sI, Rational, \
    kronecker_product
s1 = Matrix([[0, 1], [1, 0]])
s2 = Matrix([[0, -sI], [sI, 0]])
s3 = Matrix([[1, 0], [0, -1]])
gam = [s1, s2]
while gam[0].rows < 32:
    d0 = gam[0].rows
    gam = [kronecker_product(g, s3) for g in gam] + \
          [kronecker_product(eye(d0), s1),
           kronecker_product(eye(d0), s2)]
G = gam
prod = eye(32)
for a in range(10):
    prod = prod * G[a]
Gamma = sI * prod
minus_idx = [k for k in range(32) if Gamma[k, k] == -1]
plus_idx = [k for k in range(32) if Gamma[k, k] == 1]
mpos = {v: i for i, v in enumerate(minus_idx)}
ppos = {v: i for i, v in enumerate(plus_idx)}
pairs = [(a, b) for a in range(10) for b in range(a + 1, 10)]
SigL = [(G[p[0]] * G[p[1]] - G[p[1]] * G[p[0]]) * Rational(1, 4)
        for p in pairs]
rho = [S.extract(minus_idx, minus_idx) for S in SigL]
rhop = [S.extract(plus_idx, plus_idx) for S in SigL]
Mchg = G[0] * G[2] * G[4] * G[6] * G[8]
Bmat = Mchg.extract(minus_idx, plus_idx)

def colmap_sec(Mlist, n=16):
    out = []
    for Mx in Mlist:
        cols = []
        for c in range(n):
            nz = [(r, Mx[r, c]) for r in range(n) if Mx[r, c] != 0]
            assert len(nz) == 1
            cols.append(nz[0])
        out.append(cols)
    return out
cm = colmap_sec(rho)
cmp_ = colmap_sec(rhop)
gcm = []
for a in range(10):
    cols = []
    for c in range(32):
        nz = [(r, G[a][r, c]) for r in range(32) if G[a][r, c] != 0]
        assert len(nz) == 1
        cols.append(nz[0])
    gcm.append(cols)
Wbas = [[[0] * 16 for _ in range(16)] for _ in range(10)]
for a in range(10):
    for t in range(16):
        tgt, ph = gcm[a][minus_idx[t]]
        u = ppos[tgt]
        for r in range(16):
            if Bmat[r, u] != 0:
                Wbas[a][r][t] = ph * Bmat[r, u]
E10 = {}
for (a, b) in pairs:
    Mx = zeros(10)
    Mx[a, b] = 1
    Mx[b, a] = -1
    E10[(a, b)] = Mx

WbasP = [[[from_sympy(Wbas[a][r][t]) if Wbas[a][r][t] != 0 else Z2
           for t in range(16)] for r in range(16)] for a in range(10)]
BP = [[from_sympy(Bmat[r, c]) if Bmat[r, c] != 0 else Z2
       for c in range(16)] for r in range(16)]

# assembled rep operators as dict ops on V27 (S 0..15, V 16..25,
# singlet 26)
RY = []
for j, p in enumerate(pairs):
    M = {}
    for c in range(16):
        t, ph = cm[j][c]
        M[(t, c)] = from_sympy(ph)
    Ej = E10[p]
    for bb in range(10):
        for b in range(10):
            if Ej[bb, b] != 0:
                M[(16 + bb, 16 + b)] = from_sympy(Ej[bb, b])
    RY.append(M)
RZp = {}
for t in range(16):
    RZp[(t, t)] = ONE2
for b in range(10):
    RZp[(16 + b, 16 + b)] = pscale(ONE2, Fraction(-2))
RZp[(26, 26)] = pscale(ONE2, Fraction(4))
RTx = []
for x in range(16):
    M = {(x, 26): ONE2}
    for t in range(16):
        for a in range(10):
            w = WbasP[a][x][t]
            if phas(w):
                M[(16 + a, t)] = w
    RTx.append(M)
RTy = []
a3 = (Fraction(-1, 16), Fraction(0))
a4 = (Fraction(-1, 8), Fraction(0))
for y in range(16):
    M = {}
    for b in range(10):
        tgt, ph = gcm[b][plus_idx[y]]
        tp = mpos[tgt]
        M[(tp, 16 + b)] = pmul(a3, from_sympy(ph))
    for r in range(16):
        if phas(BP[r][y]):
            M[(26, r)] = pmul(a4, BP[r][y])
    RTy.append(M)
REPB = RY + [RZp] + RTx + RTy
assert len(REPB) == 78
expressREP = make_express(REPB)
ok("P1b rep side rebuilt: 78 assembled operators on V27", True)

# rep cubic tensor from C(U) = W(psi,psi).v - s (v.v).
# ERRATUM (found in this build): e6_rep27.py's tensor-delta
# invariance check contracts a stored triple's three slots into
# one sorted bucket, double-counting the derivation action on
# repeated indices; its (2, -1) ratio is an artifact of that
# defect.  Under the derivation action the invariant ratio is
# (1, -1), re-verified here from scratch (delta = 0 for all 78
# generators, and det o Phi = 16i * C below).
def Cform(U):
    psi = {r: U[r] for r in range(16) if r in U}
    v = {a: U[16 + a] for a in range(10) if (16 + a) in U}
    s = U.get(26, Z2)
    tot = Z2
    for a in range(10):
        if a in v:
            w = Z2
            for r, pr in psi.items():
                for t, pt in psi.items():
                    if phas(WbasP[a][r][t]):
                        w = padd(w, pmul(WbasP[a][r][t], pmul(pr, pt)))
            tot = padd(tot, pmul(w, v[a]))
    if phas(s):
        vv = Z2
        for a, va in v.items():
            vv = padd(vv, pmul(va, va))
        tot = psub(tot, pmul(s, vv))
    return tot

def C3(I, J, K):
    eI = {I: ONE2}
    eJ = {J: ONE2}
    eK = {K: ONE2}
    s = Cform(coord_add(coord_add(eI, eJ), eK))
    s = psub(s, Cform(coord_add(eI, eJ)))
    s = psub(s, Cform(coord_add(eI, eK)))
    s = psub(s, Cform(coord_add(eJ, eK)))
    s = padd(s, Cform(eI))
    s = padd(s, Cform(eJ))
    s = padd(s, Cform(eK))
    return pscale(s, Fraction(1, 6))

TC = {}
for t in TRIPLES:
    v = C3(*t)
    if phas(v):
        TC[t] = v
pat = set()
for (i, j, k) in TC:
    def piece(n):
        return 0 if n < 16 else (1 if n < 26 else 2)
    pat.add(tuple(sorted((piece(i), piece(j), piece(k)))))
ok("P1c rep cubic tensor built; patterns are (S,S,V) and (V,V,1) "
   "only", pat == {(0, 0, 1), (1, 1, 2)}, str(pat))

# ================================================================ h elements
AVALS = [1, 2, 4, 8, -15]
plane_pos = [pairs.index((2 * k, 2 * k + 1)) for k in range(5)]
HR = {}
for k in range(5):
    HR = mat_add(HR, mat_scale(RY[plane_pos[k]],
                               (Fraction(AVALS[k]), Fraction(0))))
s_block_offdiag = [(r, c) for (r, c) in HR
                   if r < 16 and c < 16 and r != c]
ok("P2a h_R is diagonal on S (torus in the Fock basis)",
   not s_block_offdiag, str(s_block_offdiag[:3]))
lamS = [HR.get((t, t), Z2) for t in range(16)]

def uR(lam_sign, k):
    """rep V-eigenvector of h_R, eigenvalue lam_sign * i * a_k."""
    return {16 + 2 * k: ONE2,
            16 + 2 * k + 1: (Fraction(0), Fraction(lam_sign))}

def vec_scale(v, c):
    return {k: pmul(x, c) for k, x in v.items() if phas(pmul(x, c))}

def vec_add(*vs):
    out = {}
    for v in vs:
        for k, x in v.items():
            out[k] = padd(out.get(k, Z2), x)
    return {k: x for k, x in out.items() if phas(x)}

hr_ok = True
for k in range(5):
    for sgn in (1, -1):
        u = uR(sgn, k)
        lam = (Fraction(0), Fraction(sgn * AVALS[k]))
        if apply_op(HR, u) != vec_scale(u, lam):
            hr_ok = False
ok("P2b h_R V-eigenvectors/eigenvalues as constructed", hr_ok)

# Jordan torus by J0-restriction solves
J0ID = [1, 2] + list(range(19, 27))

def solve_rect(rows, b, n):
    """Solve M x = b exactly (rows: list of lists of pairs).
    Rref of the augmented matrix; asserts consistency and that all
    n columns are pivot columns (unique solution)."""
    Mm = sp.Matrix(len(rows), n,
                   lambda i, j: to_sympy(rows[i][j]))
    bv = sp.Matrix(len(rows), 1, lambda i, _: to_sympy(b[i]))
    A = Mm.row_join(bv)
    R, pivs = A.rref()
    assert list(pivs) == list(range(n)), \
        f"solve_rect: inconsistent or non-unique, pivots {pivs}"
    sol = [None] * n
    for ri, pc in enumerate(pivs):
        sol[pc] = from_sympy(R[ri, n])
    assert all(s is not None for s in sol)
    return sol

def derj_with_J0(target):
    """target: dict {(r, c): pair} on J0ID x J0ID.  Return the DERJ
    operator with that J0 block."""
    rows = []
    b = []
    for r in J0ID:
        for c in J0ID:
            rows.append([DERJ[a].get((r, c), Z2) for a in range(45)])
            b.append(target.get((r, c), Z2))
    cf = solve_rect(rows, b, 45)
    op = lincomb(DERJ, {a: cf[a] for a in range(45) if phas(cf[a])})
    for r in J0ID:
        for c in J0ID:
            if op.get((r, c), Z2) != target.get((r, c), Z2):
                return None, None
    return op, cf

KJ = []
KJ_COEF = []
_op, _cf = derj_with_J0({(1, 1): ONE2, (2, 2): pneg(ONE2)})
KJ.append(_op)
KJ_COEF.append(_cf)
for m in range(1, 5):
    c1 = 19 + 2 * (m - 1)
    c2 = c1 + 1
    _op, _cf = derj_with_J0({(c2, c1): ONE2, (c1, c2): pneg(ONE2)})
    KJ.append(_op)
    KJ_COEF.append(_cf)
ok("P2c Peirce-frame torus K_0..K_4 solved inside the derived "
   "Levi by J0 restriction", all(K is not None for K in KJ))

def build_hJ(a4sign):
    a = [1, 2, 4, 8, a4sign * 15]
    H = mat_scale(KJ[0], (Fraction(0), Fraction(a[0])))
    for m in range(1, 5):
        H = mat_add(H, mat_scale(KJ[m], (Fraction(a[m]),
                                         Fraction(0))))
    return H

_x = sp.Symbol('x')
_target_poly = 1
for _lam in lamS:
    _target_poly = _target_poly * (_x - to_sympy(_lam))
_target_poly = sp.expand(_target_poly)

def spinor_charpoly(H):
    M = sp.Matrix(16, 16, lambda i, j: to_sympy(
        H.get((3 + i, 3 + j), Z2)))
    return sp.expand(M.charpoly(_x).as_expr())

def spectrum_match(H):
    return sp.expand(spinor_charpoly(H) - _target_poly) == 0

HJ = build_hJ(-1)
FLIP = False
if not spectrum_match(HJ):
    HJ = build_hJ(1)
    FLIP = True
ok("P2d h_J spinor spectrum equals h_R spinor spectrum "
   "(characteristic polynomials)", spectrum_match(HJ))
print(f"  h_J parity flip applied: {FLIP}", flush=True)
# effective J-side torus weights: a^J_k = AVALS[k] for k < 4; the flip
# changes the sign of the a_4 weight only.  JSIGN[k] = a^J_k / a_k.
AJ = list(AVALS)
if FLIP:
    AJ[4] = -AJ[4]
JSIGN = [1 if AJ[k] == AVALS[k] else -1 for k in range(5)]

# ================================================================ Phi0
def uJ(lam_sign, k):
    """Jordan J0-eigenvector of h_J, eigenvalue lam_sign*i*a^J_k."""
    if k == 0:
        return {1: ONE2} if lam_sign == 1 else {2: ONE2}
    c1 = 19 + 2 * (k - 1)
    c2 = c1 + 1
    return {c1: ONE2, c2: (Fraction(0), Fraction(-lam_sign))}

hj_ok = True
for k in range(5):
    for sgn in (1, -1):
        u = uJ(sgn, k)
        lam = (Fraction(0), Fraction(sgn * AJ[k]))
        if apply_op(HJ, u) != vec_scale(u, lam):
            hj_ok = False
ok("P2e h_J J0-eigenvectors/eigenvalues as constructed", hj_ok)

PHI0 = {}
inv2 = (Fraction(1, 2), Fraction(0))
inv2i = (Fraction(0), Fraction(-1, 2))
for k in range(5):
    # Phi0(uR(sgn, k)) = uJ(sgn * JSIGN[k], k): same h-eigenvalue.
    up = uJ(JSIGN[k], k)
    um = uJ(-JSIGN[k], k)
    col_even = vec_scale(vec_add(up, um), inv2)
    col_odd = vec_scale(vec_add(up, vec_scale(um, pneg(ONE2))),
                        inv2i)
    for r, v in col_even.items():
        PHI0[(r, 16 + 2 * k)] = v
    for r, v in col_odd.items():
        PHI0[(r, 16 + 2 * k + 1)] = v
for t in range(16):
    lam = lamS[t]
    eqs = []
    for ri in range(16):
        e = {}
        for ci in range(16):
            v = HJ.get((3 + ri, 3 + ci), Z2)
            if ri == ci:
                v = psub(v, lam)
            if phas(v):
                e[ci] = v
        eqs.append(e)
    ns = nullspace(eqs, 16)
    assert len(ns) == 1, (t, len(ns))
    for ci, v in ns[0].items():
        PHI0[(3 + ci, t)] = v
PHI0[(0, 26)] = ONE2
PHI0M = PHI0
ok("P3 Phi0 intertwines h: Phi0 h_R = h_J Phi0",
   matmul(PHI0M, HR) == matmul(HJ, PHI0M))

def prop_ratio(T, w):
    """If T w = lam * w, return lam, else None."""
    Tw = apply_op(T, w)
    if not Tw:
        return Z2
    lam = None
    for k in set(w) | set(Tw):
        tv = Tw.get(k, Z2)
        wv = w.get(k, Z2)
        if not phas(wv):
            return None
        cand = pdiv(tv, wv) if phas(tv) else Z2
        if lam is None:
            lam = cand
        elif cand != lam:
            return None
    return lam

def vec_ratio(u, v):
    """If u = lam * v (coordinate dicts), return lam, else None."""
    if not u and not v:
        return Z2
    if not v:
        return None
    lam = None
    for k in set(u) | set(v):
        uv = u.get(k, Z2)
        vv = v.get(k, Z2)
        if not phas(vv):
            return None
        cand = pdiv(uv, vv) if phas(uv) else Z2
        if lam is None:
            lam = cand
        elif cand != lam:
            return None
    return lam

# ================================================================ P4: beta0
print("structure constants of the rep basis...", flush=True)
SC_REP = {}
for g in range(78):
    for h in range(g + 1, 78):
        cc = commutator(REPB[g], REPB[h])
        if cc:
            cf = expressREP(cc)
            assert cf is not None, (g, h)
            SC_REP[(g, h)] = {m: c for m, c in enumerate(cf)
                              if phas(c)}
print("structure constants of the derived Levi basis...", flush=True)
SC_DERJ = {}
for g in range(45):
    for h in range(g + 1, 45):
        cc = commutator(DERJ[g], DERJ[h])
        if cc:
            cf = expressDER(cc)
            assert cf is not None, (g, h)
            SC_DERJ[(g, h)] = {m: c for m, c in enumerate(cf)
                               if phas(c)}

def bracket_coords(SC, v, w):
    out = {}
    for g, vg in v.items():
        for h, wh in w.items():
            if g == h:
                continue
            if g < h:
                sc, sgn = SC.get((g, h), {}), 1
            else:
                sc, sgn = SC.get((h, g), {}), -1
            for m, c in sc.items():
                out[m] = padd(out.get(m, Z2),
                              pscale(pmul(pmul(vg, wh), c),
                                     Fraction(sgn)))
    return {m: c for m, c in out.items() if phas(c)}

def bcoords_rep(v, w):
    return bracket_coords(SC_REP, v, w)

def bcoords_derj(v, w):
    return bracket_coords(SC_DERJ, v, w)

def op_ratio(A, B):
    """A = lam * B as operators?  return lam else None."""
    lam = None
    for e, v in B.items():
        if e not in A:
            return None
        c = pdiv(A[e], v)
        if lam is None:
            lam = c
        elif c != lam:
            return None
    for e in A:
        if e not in B:
            return None
    return lam if lam is not None else Z2

HKR = [RY[plane_pos[k]] for k in range(5)]

# weight vectors of the 27 and their Phi0 images
WTS = []
for t in range(16):
    WTS.append(({t: ONE2},
                {r: PHI0M[(r, t)] for r in range(27)
                 if (r, t) in PHI0M}))
IMG_UR = {}
for k in range(5):
    for sgn in (1, -1):
        IMG_UR[(sgn, k)] = uJ(sgn * JSIGN[k], k)
wv_ok = True
for k in range(5):
    for sgn in (1, -1):
        w = uR(sgn, k)
        WTS.append((w, IMG_UR[(sgn, k)]))
        if apply_op(PHI0M, w) != IMG_UR[(sgn, k)]:
            wv_ok = False
WTS.append(({26: ONE2}, {0: ONE2}))
ok("P4a Phi0 maps the rep V-eigenvectors to the J0 eigenvectors "
   "of the same h-eigenvalue (flip-adjusted)", wv_ok)

centre_ok = True
for (w, img) in WTS:
    if prop_ratio(ZHATN, img) != prop_ratio(RZp, w):
        centre_ok = False
ok("P4b beta0(Z') = ZHATN: centre weights agree on all 27 weight "
   "vectors", centre_ok)

TORUSB = [ZHATN] + KJ
BK = []
BK_Z = []
BK_D = []
for k in range(5):
    rows = []
    b = []
    bas_imgs = {}
    for (w, img) in WTS:
        lam = prop_ratio(HKR[k], w)
        assert lam is not None
        ims = [apply_op(TB, img) for TB in TORUSB]
        coords = set(img)
        for im in ims:
            coords |= set(im)
        for r in coords:
            rows.append([im.get(r, Z2) for im in ims])
            b.append(pmul(lam, img.get(r, Z2)))
    x = solve_rect(rows, b, 6)
    BK.append(lincomb(TORUSB, {j: x[j] for j in range(6)
                               if phas(x[j])}))
    BK_Z.append(x[0])
    d = {}
    for m in range(5):
        if phas(x[m + 1]):
            for a in range(45):
                if phas(KJ_COEF[m][a]):
                    d[a] = padd(d.get(a, Z2),
                                pmul(x[m + 1], KJ_COEF[m][a]))
    BK_D.append({a: c for a, c in d.items() if phas(c)})
ok("P4c beta0 on the plane torus solved from the 27 weights; "
   "centre components (expect all zero: so(10) is perfect)",
   all(not phas(z) for z in BK_Z), f"z = {BK_Z}")

# root vectors of so(10) wrt the plane torus (rep side, by ansatz
# search over the four Sigma's mixing two planes)
roots = []
for k in range(5):
    for l in range(k + 1, 5):
        a, ap, bb, bp = 2 * k, 2 * k + 1, 2 * l, 2 * l + 1
        idxs = [pairs.index((a, bb)), pairs.index((ap, bb)),
                pairs.index((a, bp)), pairs.index((ap, bp))]
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    coef = {idxs[0]: ONE2,
                            idxs[1]: (Fraction(0), Fraction(s1)),
                            idxs[2]: (Fraction(0), Fraction(s2)),
                            idxs[3]: (Fraction(s3), Fraction(0))}
                    X = lincomb(RY, coef)
                    dk = op_ratio(commutator(HKR[k], X), X)
                    if dk is None:
                        continue
                    dl = op_ratio(commutator(HKR[l], X), X)
                    if dl is None:
                        continue
                    delta = [Z2] * 5
                    delta[k] = dk
                    delta[l] = dl
                    roots.append({"coef": coef,
                                  "delta": tuple(delta), "mat": X})
ok("P4d 40 root vectors of so(10) wrt the plane torus found on "
   "the rep side", len(roots) == 40, f"found {len(roots)}")

# J-side root representatives: solve [BK_k, X] = delta_k X in DERJ
CCk = [[bcoords_derj(BK_D[k], {a: ONE2}) for a in range(45)]
       for k in range(5)]
RJ = []
rj_ok = True
for gi, rt in enumerate(roots):
    eqs = []
    for k in range(5):
        dk = rt["delta"][k]
        for bb in range(45):
            e = {}
            for a in range(45):
                c = CCk[k][a].get(bb, Z2)
                if a == bb:
                    c = psub(c, dk)
                if phas(c):
                    e[a] = c
            if e:
                eqs.append(e)
    ns = nullspace(eqs, 45)
    if len(ns) != 1:
        rj_ok = False
        RJ.append(None)
    else:
        RJ.append(ns[0])
ok("P4e every root space on the Jordan side is exactly "
   "1-dimensional (45-coefficient solve)", rj_ok)
RJmat = [lincomb(DERJ, rj) if rj is not None else None for rj in RJ]

# change of basis: NEW = (5 plane torus, 40 roots) in RY coordinates
NB = sp.Matrix(45, 45, lambda i, j: 0)
for k in range(5):
    NB[plane_pos[k], k] = 1
for gi, rt in enumerate(roots):
    for ri, c in rt["coef"].items():
        NB[ri, 5 + gi] = to_sympy(c)
NBinv = NB.inv()
NBINV = [[from_sympy(NBinv[i, j]) if NBinv[i, j] != 0 else Z2
          for j in range(45)] for i in range(45)]

def new_coords(v):
    """RY-coefficient dict -> NEW-basis coefficient dict."""
    out = {}
    for m in range(45):
        s = Z2
        for c, vc in v.items():
            if phas(NBINV[m][c]):
                s = padd(s, pmul(NBINV[m][c], vc))
        if phas(s):
            out[m] = s
    return out

# t-relations from root-root brackets
rels = []
zchecks = []
shape_ok = True
for g in range(40):
    for h in range(g + 1, 40):
        bc = bcoords_rep(roots[g]["coef"], roots[h]["coef"])
        nbc = new_coords(bc)
        jc = bcoords_derj(RJ[g], RJ[h])
        root_terms = {m: c for m, c in nbc.items() if m >= 5}
        tor_terms = {m: c for m, c in nbc.items() if m < 5}
        zc = Z2
        for m, c in tor_terms.items():
            zc = padd(zc, pmul(c, BK_Z[m]))
        zchecks.append(zc)
        if not jc and not nbc:
            continue
        if len(root_terms) == 1 and not tor_terms:
            m = next(iter(root_terms))
            rm = RJ[m - 5]
            bst = next(iter(rm))
            if bst not in jc:
                shape_ok = False
                continue
            K = pdiv(pmul(root_terms[m], rm[bst]), jc[bst])
            rels.append(("lin", g, h, m - 5, K))
        elif not root_terms and tor_terms:
            bst = next(iter(jc)) if jc else None
            if bst is None:
                shape_ok = False
                continue
            Cnum = Z2
            for m, c in tor_terms.items():
                Cnum = padd(Cnum,
                            pmul(c, BK_D[m].get(bst, Z2)))
            rels.append(("prod", g, h, pdiv(Cnum, jc[bst])))
        elif not nbc and jc:
            shape_ok = False
        else:
            shape_ok = False
ok("P4f root-bracket shapes: every root-root bracket is a single "
   "root or a torus element, on both sides; centre components of "
   "the images vanish", shape_ok and all(not phas(z)
                                         for z in zchecks))

tval = {}
progress = True
while progress:
    progress = False
    for rel in rels:
        if rel[0] == "prod":
            _, g, h, C = rel
            if g in tval and h in tval:
                assert pmul(tval[g], tval[h]) == C, ("prod", rel)
            elif g in tval:
                tval[h] = pdiv(C, tval[g])
                progress = True
            elif h in tval:
                tval[g] = pdiv(C, tval[h])
                progress = True
        else:
            _, g, h, m, K = rel
            nk = sum(x in tval for x in (g, h, m))
            if nk == 3:
                assert pmul(tval[g], tval[h]) == pmul(K, tval[m]), \
                    ("lin", rel)
            elif nk == 2:
                if m not in tval:
                    tval[m] = pdiv(pmul(tval[g], tval[h]), K)
                elif g not in tval:
                    tval[g] = pdiv(pmul(K, tval[m]), tval[h])
                else:
                    tval[h] = pdiv(pmul(K, tval[m]), tval[g])
                progress = True
    if not progress and len(tval) < 40:
        g = next(g for g in range(40) if g not in tval)
        tval[g] = ONE2
        progress = True
ok("P4g root scale factors t_g solved from the bracket relations "
   "(seeded at 1 where the character freedom remains)",
   len(tval) == 40, f"{len(tval)} of 40")

# full vector verification of every root-root relation
rel_ok = True
for g in range(40):
    for h in range(g + 1, 40):
        bc = bcoords_rep(roots[g]["coef"], roots[h]["coef"])
        nbc = new_coords(bc)
        lhs = {a: pscale(c, 1) for a, c in
               bcoords_derj(RJ[g], RJ[h]).items()}
        lhs = {a: pmul(pmul(tval[g], tval[h]), c)
               for a, c in lhs.items()}
        rhs = {}
        for m, c in nbc.items():
            if m < 5:
                src = BK_D[m]
                cc = c
            else:
                src = RJ[m - 5]
                cc = pmul(c, tval[m - 5])
            for a, v in src.items():
                rhs[a] = padd(rhs.get(a, Z2), pmul(cc, v))
        rhs = {a: c for a, c in rhs.items() if phas(c)}
        if lhs != rhs:
            rel_ok = False
ok("P4h all 780 root-root bracket relations hold as full "
   "DERJ-coordinate vectors with the solved t_g", rel_ok)

# beta0 on the whole Levi, as matrices
def lincomb_mats(mats, cf):
    out = {}
    for i, c in cf.items():
        if phas(c):
            for e, v in mats[i].items():
                out[e] = padd(out.get(e, Z2), pmul(c, v))
    return {e: v for e, v in out.items() if phas(v)}

ALPHA = [None] * 78
for j in range(45):
    cj = {m: NBINV[m][j] for m in range(45) if phas(NBINV[m][j])}
    parts = []
    for m, c in cj.items():
        if m < 5:
            parts.append(mat_scale(BK[m], c))
        else:
            parts.append(mat_scale(RJmat[m - 5],
                                   pmul(c, tval[m - 5])))
    acc = {}
    for p in parts:
        acc = mat_add(acc, p)
    ALPHA[j] = acc
ALPHA[45] = ZHATN
ok("P4i beta0 assembled on all 46 Levi basis operators", True)

# ================================================================ P5: tau
def PHI0col(t):
    return {r: PHI0M[(r, t)] for r in range(27) if (r, t) in PHI0M}

SNODES = [("S", t) for t in range(16)]
VNODES = [("V", sgn, k) for k in range(5) for sgn in (1, -1)]

def node_rep(nd):
    if nd[0] == "S":
        return {nd[1]: ONE2}
    return uR(nd[1], nd[2])

def node_img(nd):
    if nd[0] == "S":
        return PHI0col(nd[1])
    return IMG_UR[(nd[1], nd[2])]

def find_node(y, piece):
    """y: rep coord dict.  Return (node, r) with y = r * node_rep."""
    if not y:
        return None
    if piece == "S":
        if len(y) != 1:
            return None
        t = next(iter(y))
        return ("S", t), y[t]
    supp = set(y)
    k0 = (min(supp) - 16) // 2
    if not supp <= {16 + 2 * k0, 16 + 2 * k0 + 1}:
        return None
    for sgn in (1, -1):
        r = vec_ratio(y, uR(sgn, k0))
        if r is not None:
            return ("V", sgn, k0), r
    return None

def tau_bfs(nodes, piece):
    tau = {nodes[0]: ONE2}
    queue = [nodes[0]]
    while queue:
        w = queue.pop()
        for gi in range(40):
            XG = roots[gi]["mat"]
            XJG = mat_scale(RJmat[gi], tval[gi])
            y = apply_op(XG, node_rep(w))
            fn = find_node(y, piece)
            if fn is None:
                continue
            w2, r = fn
            u = apply_op(XJG, node_img(w))
            R = vec_ratio(u, node_img(w2))
            assert R is not None, (w, gi)
            want = pdiv(pmul(tau[w], R), r)
            if w2 in tau:
                assert tau[w2] == want, (w, w2, gi)
            else:
                tau[w2] = want
                queue.append(w2)
    return tau

TAU_Sn = tau_bfs(SNODES, "S")
TAU_Vn = tau_bfs(VNODES, "V")
ok("P5a tau ratios: root-action graphs connect all 16 spinor "
   "weights and all 10 vector weights, cycle-consistent",
   len(TAU_Sn) == 16 and len(TAU_Vn) == 10,
   f"S: {len(TAU_Sn)}/16, V: {len(TAU_Vn)}/10")
TAU_S = {t: TAU_Sn[("S", t)] for t in range(16)}
TAU_V = {(sgn, k): TAU_Vn[("V", sgn, k)] for k in range(5)
         for sgn in (1, -1)}

def PHIUR(sgn, k):
    return vec_scale(IMG_UR[(sgn, k)], TAU_V[(sgn, k)])

def PhiRatV(b):
    k = b // 2
    up = PHIUR(1, k)
    um = PHIUR(-1, k)
    if b % 2 == 0:
        return vec_scale(vec_add(up, um), inv2)
    return vec_scale(vec_add(up, vec_scale(um, pneg(ONE2))), inv2i)

# ================================================================ P6: p
TJm = PmJ
TJp = PpJ
x0 = 0
y0s = [y for y in range(16) if phas(BP[0][y])]
assert len(y0s) == 1
y0 = y0s[0]

eqs = []
for r in range(3, 19):
    e = {}
    for a in range(16):
        v = TJm[a].get((r, 0), Z2)
        if phas(v):
            e[a] = v
    c0 = pmul(TAU_S[0], PHI0col(0).get(r, Z2))
    if phas(c0):
        e[16] = pneg(c0)
    if e:
        eqs.append(e)
for t in range(16):
    img_t = PHI0col(t)
    Aimg = [apply_op(TJm[a], img_t) for a in range(16)]
    Vhat = {}
    for av in range(10):
        w = WbasP[av][x0][t]
        if phas(w):
            for rr, vv in PhiRatV(av).items():
                Vhat[rr] = padd(Vhat.get(rr, Z2), pmul(w, vv))
    for r in [1, 2] + list(range(19, 27)):
        e = {}
        for a in range(16):
            v = Aimg[a].get(r, Z2)
            if phas(v):
                e[a] = pmul(v, TAU_S[t])
        if phas(Vhat.get(r, Z2)):
            e[17] = pneg(Vhat[r])
        if e:
            eqs.append(e)
ns = nullspace(eqs, 18)
assert len(ns) == 1, f"x-stage nullity {len(ns)}"
sol = ns[0]
Acoef = {a: sol[a] for a in range(16) if a in sol}
# unknown 16 is mu = sigma_S/sigma_1 (the singlet equation is
# A(E) = mu * TAU_S[0] * Phi0(e_0), from sigma_1 * A(E) =
# sigma_S * TAU_S[0] * Phi0(e_0)); unknown 17 is kap = sigma_V/sigma_S
mu = sol.get(16, Z2)
kap = sol.get(17, Z2)
assert phas(mu) and phas(kap), (mu, kap)
lam = pinv(mu)   # sigma_1/sigma_S
Aop = lincomb(TJm, Acoef)
ALPHA[46 + x0] = Aop
print(f"  x-stage solved: sigma_1/sigma_S = {lam}, "
      f"sigma_V/sigma_S = {kap}", flush=True)

eqs = []
for sr in range(16):
    img = PHI0col(sr)
    e = {}
    for a in range(16):
        v = apply_op(TJp[a], img).get(0, Z2)
        if phas(v):
            e[a] = pmul(v, TAU_S[sr])
    c0 = pmul(a4, BP[sr][y0])
    if phas(c0):
        e[16] = pneg(c0)
    if e:
        eqs.append(e)
ns = nullspace(eqs, 17)
assert len(ns) == 1, f"y-stage nullity {len(ns)}"
sol = ns[0]
lam2 = sol.get(16, Z2)
assert phas(lam2)
fac = pdiv(lam, lam2)
Ccoef = {a: pmul(sol[a], fac) for a in range(16) if a in sol}
Cop = lincomb(TJp, Ccoef)
ALPHA[62 + y0] = Cop
ycheck = True
for b in range(10):
    tgt, phs = gcm[b][plus_idx[y0]]
    tp = mpos[tgt]
    ph = from_sympy(phs)
    lhs = vec_scale(apply_op(Cop, PhiRatV(b)), kap)
    rhs = vec_scale(PHI0col(tp), pmul(a3, pmul(ph, TAU_S[tp])))
    if lhs != rhs:
        ycheck = False
ok("P6a y-stage: T_y0 solved from the S-input equations (scale "
   "locked to the x-stage through sigma_1); V-input intertwining "
   "verified directly on all 10 vector inputs", ycheck)

bc = SC_REP.get((46, 62 + y0), {})
assert all(m < 46 for m in bc), "p-p bracket leaves Levi+centre"
rhs = lincomb_mats(ALPHA, bc)
ok("P6b [alpha(T_x0), alpha(T_y0)] = beta0([T_x0, T_y0])",
   commutator(Aop, Cop) == rhs)

vis = {x0}
frontier = True
while len(vis) < 16 and frontier:
    frontier = False
    for xp in list(vis):
        for j in range(45):
            tgt, phs = cm[j][xp]
            if tgt in vis:
                continue
            s = op_ratio(commutator(REPB[j], REPB[46 + xp]),
                         REPB[46 + tgt])
            assert s is not None
            ALPHA[46 + tgt] = mat_scale(
                commutator(ALPHA[j], ALPHA[46 + xp]), pinv(s))
            vis.add(tgt)
            frontier = True
ok("P6c alpha on all 16 T_x by Levi bracket generation", len(vis) == 16)
vis = {y0}
frontier = True
while len(vis) < 16 and frontier:
    frontier = False
    for yp in list(vis):
        for j in range(45):
            tgt, phs = cmp_[j][yp]
            if tgt in vis:
                continue
            s = op_ratio(commutator(REPB[j], REPB[62 + yp]),
                         REPB[62 + tgt])
            assert s is not None
            ALPHA[62 + tgt] = mat_scale(
                commutator(ALPHA[j], ALPHA[62 + yp]), pinv(s))
            vis.add(tgt)
            frontier = True
ok("P6d alpha on all 16 T_y by Levi bracket generation", len(vis) == 16)
assert all(A is not None for A in ALPHA)

# final Phi
PHI = {}
for t in range(16):
    for r, v in PHI0col(t).items():
        PHI[(r, t)] = pmul(TAU_S[t], v)
for b in range(10):
    for r, v in PhiRatV(b).items():
        PHI[(r, 16 + b)] = pmul(kap, v)
PHI[(0, 26)] = lam

def PHIcol(i):
    return {r: PHI[(r, i)] for r in range(27) if (r, i) in PHI}

# ================================================================ verify
v1 = all(matmul(PHI, REPB[g]) == matmul(ALPHA[g], PHI)
         for g in range(78))
ok("V1 Phi intertwines all 78 generators: Phi rho_R(g) = "
   "rho_J(alpha(g)) Phi", v1)
v2 = True
for g in range(78):
    for h in range(g + 1, 78):
        lhs = commutator(ALPHA[g], ALPHA[h])
        rhs = lincomb_mats(ALPHA, SC_REP.get((g, h), {}))
        if lhs != rhs:
            v2 = False
ok("V2 alpha is a Lie algebra isomorphism: all 3003 bracket "
   "pairs", v2)
def det3v(X, Y, Z):
    """Polarized determinant on three coordinate dicts."""
    s = detJ(coord_add(coord_add(X, Y), Z))
    s = psub(s, detJ(coord_add(X, Y)))
    s = psub(s, detJ(coord_add(X, Z)))
    s = psub(s, detJ(coord_add(Y, Z)))
    s = padd(s, detJ(X))
    s = padd(s, detJ(Y))
    s = padd(s, detJ(Z))
    return pscale(s, Fraction(1, 6))

c_ratio = None
for trip, val in TC.items():
    d = det3v(PHIcol(trip[0]), PHIcol(trip[1]), PHIcol(trip[2]))
    if phas(d):
        c_ratio = pdiv(d, val)
        break
v3 = c_ratio is not None
if v3:
    for trip in TRIPLES:
        d = det3v(PHIcol(trip[0]), PHIcol(trip[1]),
                  PHIcol(trip[2]))
        if d != pmul(c_ratio, TC.get(trip, Z2)):
            v3 = False
            break
ok("V3 det o Phi = c * C as cubic tensors (all 3654 triples)",
   v3, f"c = {c_ratio}")
v4 = rank_of_vectors([PHIcol(i) for i in range(27)], 27)[0] == 27
ok("V4 Phi has rank 27 (isomorphism of vector spaces)", v4)

print(f"\nCHECKS: {NCHK[0] - len(FAIL)} passed, {len(FAIL)} "
      f"failed, {NCHK[0]} total", flush=True)
for name, d in FAIL:
    print(f"  FAILED: {name} {d}", flush=True)
