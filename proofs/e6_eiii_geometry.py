#!/usr/bin/env python3
"""
The EIII homogeneous geometry on the Jordan-side model --
the last layer of the E6 stack.

Front end (rebuilt here): J = H_3(O_C), its det tensor, the
determinant-preserving algebra e6 as 78 explicit operators with
all 3081 commutators in the basis (checks J0-J4c as in
e6_jordan_e6.py), plus the conjugation machinery of
e6_real_forms.py (trace-form adjoint and theta, the Peirce
grading involution psi, the Killing scalar Killing = 4 tr_27
verified on all 3081 pairs, real fixed bases).

New content (checks G): the Hermitian symmetric pair of type
EIII, at the basepoint given by the Peirce decomposition of J
with respect to E = diag(1, 0, 0).

  * Charges: every basis operator is ad(Z')-homogeneous for the
    grading operator Z' = 6L_E - 2id, with charges in {0, +3, -3}
    and counts (46, 16, 16) (G1).  Charges add under brackets;
    the +/-6 brackets vanish (G2).  This is the 3-grading
    e6 = p- (+) k (+) p+ of the Hermitian symmetric pair.
  * k (charge 0) is a subalgebra with 1-dimensional centre
    C Z' and 45-dimensional derived algebra (G3), and k is
    exactly the algebra of operators preserving the three
    Peirce pieces (J1, J1/2, J0) of the basepoint (G4) --
    concretely so(10) + C: the Levi of the 62-parabolic, here
    appearing as the (reductive) isotropy algebra.
  * Killing on p: p+ and p- are isotropic and dually paired
    (G5).  p+ as a k-module: Z' acts as the scalar +3, the
    action is irreducible (commutant = scalars), and its
    Casimir under the J9 prescription is -15/8, the spinor
    value (G6) -- the holomorphic tangent space of EIII is the
    16-dimensional spinor module of so(10).
  * The real form: sigma_EIII (built in e6_real_forms.py, F4)
    preserves k and p; its 78 real fixed operators split as
    46 in k (k0) and 32 in p (p0) (G7).  Ad_psi is the pair
    involution: +1 on k, -1 on p at operator level.  Killing
    is negative definite on k0 (0, 46) and positive definite
    on p0 (32, 0): the Cartan decomposition of e6(-14), whose
    noncompact dual is the bounded symmetric domain of type
    EIII, real dimension 32 = complex dimension 16.
  * Hermitian data: H = iZ'/3 lies in k0, is central there,
    and ad(H) preserves p0, is real in the fixed basis, squares
    to -id, and is Killing-orthogonal: the invariant complex
    structure of EIII (G8).
  * Compact dual: sigma_compact fixes k0 pointwise and fixes
    i.p0 elementwise, so the compact form's Cartan
    decomposition is k0 (+) i p0 with the same isotropy
    algebra k0 (G9) -- EIII as a compact Hermitian symmetric
    space is the compact dual of the e6(-14) space.
  * The real bracket relations [k0, k0] in k0, [k0, p0] in p0,
    [p0, p0] in k0 are verified on all fixed-basis pairs (G10).

Cited, not computed: the classification name "EIII" and the
group-level isotropy (Spin(10) x U(1))/Z_4; everything about
the pair, the real forms, the signatures, the tangent module
and the complex structure is computed exactly over Q(i).
"""
from fractions import Fraction
from itertools import combinations, permutations
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
def phas(a): return a != Z2

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

def coord_scale(X, f):
    return {} if f == 0 else {k: pscale(v, f) for k, v in X.items()}

BASIS = [{i: ONE2} for i in range(27)]

def slot_oct(X, base):
    return {u: X[base + u] for u in range(8) if (base + u) in X}

def onorm0(A):
    return omul(A, oconj(A)).get(0, Z2)

def det(X):
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

def det3(I, J, K):
    eI, eJ, eK = BASIS[I], BASIS[J], BASIS[K]
    s = det(coord_add(coord_add(eI, eJ), eK))
    s = psub(s, det(coord_add(eI, eJ)))
    s = psub(s, det(coord_add(eI, eK)))
    s = psub(s, det(coord_add(eJ, eK)))
    s = padd(s, det(eI))
    s = padd(s, det(eJ))
    s = padd(s, det(eK))
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
    v = det3(*t)
    if phas(v):
        TDET[t] = v
print(f"  det tensor: {len(TDET)} nonzero entries", flush=True)

def tdet3(a, b, c):
    return TDET.get(tuple(sorted((a, b, c))), Z2)

# pair index for delta contractions: PAIRL[(x, y)] = [(L, val)] with
# tdet3(L, x, y) = val != 0
PAIRL = {}
for pat, v in TDET.items():
    for L in set(pat):
        rem = list(pat)
        rem.remove(L)
        PAIRL.setdefault(tuple(sorted(rem)), []).append((L, v))

# ================================================================ charges
CHARGE = {}
CHARGE[0] = 4
for idx in range(3, 19):
    CHARGE[idx] = 1
for idx in [1, 2] + list(range(19, 27)):
    CHARGE[idx] = -2

# ================================================================ sparse GE
def nullspace(eqs, nvars, pivot_hint=None):
    """eqs: list of dict {var: pair}.  Returns list of basis vectors
    (dicts {var: pair}) spanning {x : A x = 0}."""
    occ = [0] * nvars
    for e in eqs:
        for v in e:
            occ[v] += 1
    pivots = {}   # pvar -> eq dict, normalized, no other pivot vars
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
        # eliminate p from existing pivot equations
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
    """rank of a set of vectors (dicts) over Q(i), via the same
    incremental pivot reduction; returns (rank, pivot_vars)."""
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
    print(f"sector q={qT:+d}: {len(varlist)} vars, {len(eqs)} eqs; "
          f"solving...", flush=True)
    sol = nullspace(eqs, len(varlist))
    print(f"  nullity {len(sol)}", flush=True)
    SECTORS[qT] = (varlist, sol)

dims = {q: len(SECTORS[q][1]) for q in SECTORS}
ok("J1 sector dims (0, +3, -3, +6, -6) = (46, 16, 16, 0, 0), total 78",
   dims == {0: 46, 3: 16, -3: 16, 6: 0, -6: 0} and
   sum(dims.values()) == 78, str(dims))

# assemble basis operators as dict matrices {(r, c): pair}
OPS = []        # (matrix, sector q)
for qT in (0, 3, -3):
    varlist, sol = SECTORS[qT]
    for x in sol:
        M = {}
        for vi, val in x.items():
            M[varlist[vi]] = val
        OPS.append((M, qT))
NB = len(OPS)
print(f"  assembled {NB} basis operators", flush=True)

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
    """Scale a dict matrix by a Q(i) scalar (a Fraction pair)."""
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

def trace_mat(A):
    t = Z2
    for i in range(27):
        if (i, i) in A:
            t = padd(t, A[(i, i)])
    return t

def delta_zero(T):
    """True iff delta_T Tdet = 0 (direct contraction, all triples)."""
    for (I, J, K) in TRIPLES:
        tot = Z2
        for (col, x, y) in ((I, J, K), (J, I, K), (K, I, J)):
            for (L, v) in PAIRL.get(tuple(sorted((x, y))), []):
                c = T.get((L, col))
                if c is not None:
                    tot = padd(tot, pmul(c, v))
        if phas(tot):
            return False
    return True

# J0 regression guard: the known solution Z_J = 6 L_E - 2 id
# (Jordan build, S8) must satisfy delta = 0 under THIS contraction.
LE_cols = [jprod({0: ONE2}, BASIS[j]) for j in range(27)]
ZJ = {}
for c in range(27):
    for r, v in LE_cols[c].items():
        ZJ[(r, c)] = pscale(v, Fraction(6))
    ZJ[(c, c)] = psub(ZJ.get((c, c), Z2), pscale(ONE2, Fraction(2)))
ok("J0 known solution: delta_{6 L_E - 2 id} det = 0 (direct)",
   delta_zero(ZJ))

# J2: every basis op preserves det (solver cross-check)
ok(f"J2 all {NB} basis operators satisfy delta_T det = 0 (direct)",
   all(delta_zero(M) for (M, q) in OPS))

# J3: closure spot-verification on commutators
fixed = [0, 1, 2, 46, 47, 48, 62, 63, 64, 77]
nspot = 0
spot_ok = True
for a in fixed:
    for b in range(NB):
        if not delta_zero(commutator(OPS[a][0], OPS[b][0])):
            spot_ok = False
        nspot += 1
ok(f"J3 delta_[Ta,Tb] det = 0 for {nspot} commutators "
   f"({len(fixed)} fixed ops x all 78)", spot_ok)

# ================================================================ minor
def op_vec(M):
    return {r * 27 + c: v for (r, c), v in M.items()}

print("building pivot minor for the 78...", flush=True)
vecs = [op_vec(M) for (M, q) in OPS]
_, minor_pos = rank_of_vectors(vecs, 729)
assert len(minor_pos) == NB

def to_sympy(x):
    return sp.Rational(x[0].numerator, x[0].denominator) + \
        sp.I * sp.Rational(x[1].numerator, x[1].denominator)

def from_sympy(z):
    z = sp.expand(z)
    re = sp.re(z)
    im = sp.im(z)
    return (Fraction(re.p, re.q), Fraction(im.p, im.q))

MINOR = sp.Matrix(NB, NB, lambda k, a: to_sympy(
    OPS[a][0].get((minor_pos[k] // 27, minor_pos[k] % 27), Z2)))
MINOR_INV = MINOR.inv()
print("  minor inverted", flush=True)

def express(M):
    """Coefficients of M in the OPS basis via the minor, with full
    matrix verification.  Returns coeff list or None."""
    rhs = sp.Matrix(NB, 1, lambda k, _: to_sympy(
        M.get((minor_pos[k] // 27, minor_pos[k] % 27), Z2)))
    sol = MINOR_INV * rhs
    coeffs = [from_sympy(sol[a]) for a in range(NB)]
    rec = {}
    for a, c in enumerate(coeffs):
        if phas(c):
            rec = mat_add(rec, mat_scale(OPS[a][0], c))
    return coeffs if rec == M else None

# all pairwise commutators, expressed
print("computing all pairwise commutators...", flush=True)
COMM = {}       # (a, b) -> coefficient list in OPS basis
comm_fail = None
for a in range(NB):
    for b in range(a + 1, NB):
        C = commutator(OPS[a][0], OPS[b][0])
        cf = express(C) if C else [Z2] * NB
        if cf is None:
            comm_fail = (a, b)
        COMM[(a, b)] = cf
ok("J4a all 3081 pairwise commutators lie in the span (closure, "
   "fully verified)", comm_fail is None, str(comm_fail))

def comm_coeffs(a, b):
    if a == b:
        return [Z2] * NB
    if a < b:
        return COMM[(a, b)]
    return [pneg(c) for c in COMM[(b, a)]]

# J4b: centre of the 78 = 0 ; [g, g] = g
ceqs = []
for b in range(NB):
    for r in range(NB):
        # coefficient equation: (sum_a c_a [Ta, Tb]) projected at
        # minor position r must vanish -- use full: build per-b,
        # per-minor-row equations (sufficient: minor separates)
        e = {}
        for a in range(NB):
            v = comm_coeffs(a, b)[r]
            if phas(v):
                e[a] = v
        if e:
            ceqs.append(e)
cen = nullspace(ceqs, NB)
ok("J4b centre of the 78-dim algebra is 0", len(cen) == 0,
   f"nullity {len(cen)}")
allc = [tuple(comm_coeffs(a, b)) for a in range(NB)
        for b in range(a + 1, NB)]
allc = [{a: v for a, v in enumerate(cf) if phas(v)} for cf in allc]
rk, _ = rank_of_vectors(allc, NB)
ok("J4c [g, g] spans the full 78 (perfect)", rk == 78, f"rank {rk}")

def pconj(a):
    return (a[0], -a[1])

def pinv_(a):
    return pinv(a)

# ---- Jordan product tensor, L operators, trace form
PROD = {}
for i in range(27):
    for j in range(i, 27):
        PROD[(i, j)] = jprod(BASIS[i], BASIS[j])
        PROD[(j, i)] = PROD[(i, j)]

def jprod_uv(u, v):
    """Jordan product of two coord dicts via the PROD tensor."""
    out = {}
    for i, ci in u.items():
        for j, cj in v.items():
            cij = pmul(ci, cj)
            for r, w in PROD[(i, j)].items():
                out[r] = padd(out.get(r, Z2), pmul(cij, w))
                if not phas(out[r]):
                    del out[r]
    return out

# trace bilinear form B(X, Y) = tr(X o Y), tr = sum of diagonal coords
GRAM = [[Z2] * 27 for _ in range(27)]
for i in range(27):
    for j in range(27):
        t = Z2
        for d in range(3):
            t = padd(t, PROD[(i, j)].get(d, Z2))
        GRAM[i][j] = t
ok("R1a trace form Gram is diagonal with entries (1,1,1, 2 x24)",
   all(GRAM[i][j] == (ONE2 if i < 3 else (Fraction(2), Fraction(0)))
       for i in range(27) for j in range(27) if i == j) and
   all(not phas(GRAM[i][j]) for i in range(27) for j in range(27)
       if i != j))

def adjoint(M):
    """B-adjoint: (T*)_{c,r} = (G_r / G_c) * T_{r,c}, G diagonal."""
    out = {}
    for (r, c), v in M.items():
        g = pmul(GRAM[r][r], pinv(GRAM[c][c]))
        out[(c, r)] = pmul(g, v)
    return out

def theta_op(M):
    """theta(T) = -T^{*B} at operator level."""
    return {k: pneg(v) for k, v in adjoint(M).items()}

ok("R1b adjoint is an involution on all 78 basis operators",
   all(adjoint(adjoint(M)) == M for (M, q) in OPS))

# Z' = 6 L_E - 2 id (Jordan build S8): Peirce block scalars (4, 1, -2)
LE_cols = [jprod({0: ONE2}, BASIS[j]) for j in range(27)]
ZP = {}
for c in range(27):
    for r, v in LE_cols[c].items():
        ZP[(r, c)] = pscale(v, Fraction(6))
    ZP[(c, c)] = psub(ZP.get((c, c), Z2), pscale(ONE2, Fraction(2)))
cfZP = express(ZP)
ZPCOORD = {a: v for a, v in enumerate(cfZP) if phas(v)} if cfZP else None
ok("R1c Z' = 6L_E - 2id lies in the algebra, block scalars (4, 1, -2)",
   cfZP is not None and ZP.get((0, 0)) == pscale(ONE2, Fraction(4))
   and ZP.get((3, 3)) == ONE2 and ZP.get((1, 1)) == pscale(ONE2, Fraction(-2)))

# ---- conjugations on J coords
def tau_vec(X, D):
    return {k: pscale(pconj(v), Fraction(D[k])) for k, v in X.items()
            if phas(pscale(pconj(v), Fraction(D[k])))}

def tau_op(M, D):
    return {(r, c): pscale(pconj(v), Fraction(D[r] * D[c]))
            for (r, c), v in M.items()}

D0 = [1] * 27
DS = []
for idx in range(27):
    if idx < 3:
        DS.append(1)
    else:
        u = (idx - 3) % 8
        DS.append(-1 if u % 2 == 1 else 1)

def tau_fixed_basis(D):
    """Real basis of the tau-fixed subspace of J (coord dicts)."""
    out = []
    for k in range(27):
        if D[k] == 1:
            out.append({k: ONE2})
        else:
            out.append({k: (Fraction(0), Fraction(1))})
    return out

# ---- exact signature of a symmetric rational matrix
def signature(Ain):
    n = len(Ain)
    A = [[Fraction(Ain[i][j]) for j in range(n)] for i in range(n)]
    npos = nneg = nzero = 0
    for k in range(n):
        if A[k][k] == 0:
            piv = None
            for i in range(k + 1, n):
                if A[i][i] != 0:
                    piv = i
                    break
            if piv is not None:
                A[k], A[piv] = A[piv], A[k]
                for r in range(n):
                    A[r][k], A[r][piv] = A[r][piv], A[r][k]
            else:
                src = None
                for i in range(k + 1, n):
                    if A[i][k] != 0:
                        src = i
                        break
                if src is None:
                    nzero += 1
                    continue
                oldcol = [A[r][src] for r in range(n)]
                for j in range(n):
                    A[k][j] += A[src][j]
                for r in range(n):
                    A[r][k] += oldcol[r]
        d = A[k][k]
        if d > 0:
            npos += 1
        elif d < 0:
            nneg += 1
        else:
            nzero += 1
            continue
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = A[i][k] / d
                for j in range(n):
                    A[i][j] -= f * A[k][j]
                for r in range(n):
                    A[r][i] -= f * A[r][k]
    return npos, nneg, nzero

ok("R1d signature routine controls: diag(1,-1) -> (1,1,0); "
   "[[0,1],[1,0]] -> (1,1,0); [[2,1],[1,2]] -> (2,0,0)",
   signature([[1, 0], [0, -1]]) == (1, 1, 0) and
   signature([[0, 1], [1, 0]]) == (1, 1, 0) and
   signature([[2, 1], [1, 2]]) == (2, 0, 0))

def combo_ops(cf):
    out = {}
    for a, c in cf.items():
        if phas(c):
            out = mat_add(out, mat_scale(OPS[a][0], c))
    return out

def apply_op(T, v):
    out = {}
    for (r, c), val in T.items():
        if c in v:
            out[r] = padd(out.get(r, Z2), pmul(val, v[c]))
    return {k: x for k, x in out.items() if phas(x)}

def trace_prod(A, B):
    t = Z2
    for (r, c), v in A.items():
        w = B.get((c, r))
        if w is not None:
            t = padd(t, pmul(v, w))
    return t


# ---- phi automorphisms of J (signed permutations on coords)
def phi_diag(signs):
    return {(i, i): pscale(ONE2, Fraction(s)) for i, s in enumerate(signs)}

def phi_perm_pi():
    PH = {}
    def put(r, c, s):
        PH[(r, c)] = pscale(ONE2, Fraction(s))
    put(0, 0, 1); put(1, 2, 1); put(2, 1, 1)
    for u in range(8):
        put(3 + u, 11 + u, 1)
        put(11 + u, 3 + u, 1)
    put(19, 19, 1)
    for u in range(1, 8):
        put(19 + u, 19 + u, -1)
    return PH

SD_psi = []
for idx in range(27):
    SD_psi.append(-1 if 3 <= idx <= 18 else 1)
SD_D = [1, 1, 1] + [-1] * 16 + [1] * 8
SD_D2 = [1, 1, 1] + [1] * 8 + [-1] * 8 + [-1] * 8

def matmul_full(A, B):
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

def ad_phi(M, PH):
    """Conjugation of an operator by phi: phi M phi^-1, phi^2 = I."""
    return matmul_full(matmul_full(PH, M), PH)

def phi_vec(X, PH):
    out = {}
    for (r, c), v in PH.items():
        if c in X:
            out[r] = padd(out.get(r, Z2), pmul(v, X[c]))
            if not phas(out[r]):
                del out[r]
    return out

PHIS = {
    "psi": phi_diag(SD_psi),
    "D": phi_diag(SD_D),
    "D2": phi_diag(SD_D2),
    "s": phi_diag(DS),
    "pi": phi_perm_pi(),
}
IDENT27 = {(i, i): ONE2 for i in range(27)}
ok("R2a each phi is an involution on J coords (phi^2 = I)",
   all(matmul_full(PH, PH) == IDENT27 for PH in PHIS.values()))

# phi is a Jordan automorphism + det-preserving, on basis pairs / tests
test_els = list(BASIS)
pair_sums = {}
for i in range(27):
    for j in range(i + 1, 27):
        pair_sums[(i, j)] = coord_add(BASIS[i], BASIS[j])
test_els += list(pair_sums.values())
triple_sums = []
for i, j, k in combinations(range(27), 3):
    triple_sums.append(coord_add(coord_add(BASIS[i], BASIS[j]), BASIS[k]))
test_els += triple_sums

def phi_checks(nm):
    PH = PHIS[nm]
    aut = True
    for i in range(27):
        for j in range(i, 27):
            lhs = phi_vec(PROD[(i, j)], PH)
            rhs = jprod_uv(phi_vec(BASIS[i], PH), phi_vec(BASIS[j], PH))
            if lhs != rhs:
                aut = False
    detp = all(det(phi_vec(X, PH)) == det(X) for X in test_els)
    return aut and detp

# ---- coordinate machinery for conjugations of the algebra
CD = {}
for a in range(NB):
    for b in range(a + 1, NB):
        CD[(a, b)] = {c: v for c, v in enumerate(COMM[(a, b)]) if phas(v)}


def bdict(x, a):
    if x == a:
        return {}
    if x < a:
        return CD[(x, a)]
    return {c: pneg(v) for c, v in CD[(a, x)].items()}

def bracket_cd(u, v):
    """Bracket of two OPS-coordinate dicts, via COMM. Returns dict."""
    out = {}
    for a, ca in u.items():
        for b, cb in v.items():
            if a == b:
                continue
            if a < b:
                cd, sgn = CD[(a, b)], 1
            else:
                cd, sgn = CD[(b, a)], -1
            cab = pmul(ca, cb)
            if sgn < 0:
                cab = pneg(cab)
            for c, w in cd.items():
                nv = padd(out.get(c, Z2), pmul(cab, w))
                if phas(nv):
                    out[c] = nv
                elif c in out:
                    del out[c]
    return out

def sigma_cols(transform):
    """Columns of sigma in the OPS basis, or None if span not preserved."""
    cols = []
    for a in range(NB):
        cf = express(transform(OPS[a][0]))
        if cf is None:
            return None
        cols.append({i: v for i, v in enumerate(cf) if phas(v)})
    return cols

def compose_cols(A, B):
    out = []
    for b in range(78):
        acc = {}
        for a, v in B[b].items():
            for r, w in A[a].items():
                nv = padd(acc.get(r, Z2), pmul(w, v))
                if phas(nv):
                    acc[r] = nv
                elif r in acc:
                    del acc[r]
        out.append(acc)
    return out

def conj_cols(A):
    return [{r: pconj(v) for r, v in col.items()} for col in A]

def cols_is_id(A):
    return all(A[i] == {i: ONE2} for i in range(78))

def sigma_vec(S, v):
    """sigma applied to a coordinate dict (antilinear)."""
    out = {}
    for a, ca in v.items():
        cca = pconj(ca)
        for r, w in S[a].items():
            nv = padd(out.get(r, Z2), pmul(cca, w))
            if phas(nv):
                out[r] = nv
            elif r in out:
                del out[r]
    return out


def sigma_checks(S):
    """Involution + bracket compatibility in OPS coordinates.
    Bracket compatibility is checked for all 3003 basis pairs using
    the one-sided expansion H[z][a] = [T_z, sigma(T_a)]:
      [sigma Ta, sigma Tb] = -sum_z S_b[z] * H[z][a],
    and sigma([Ta, Tb]) directly from the COMM tensor."""
    inv = cols_is_id(compose_cols(S, conj_cols(S)))
    if not inv:
        return False, False
    H = {}
    for z in range(78):
        for a in range(78):
            acc = {}
            for b, cb in S[a].items():
                cd = bdict(z, b)
                for c, w in cd.items():
                    nv = padd(acc.get(c, Z2), pmul(cb, w))
                    if phas(nv):
                        acc[c] = nv
                    elif c in acc:
                        del acc[c]
            H[(z, a)] = acc
    for a in range(78):
        for b in range(a + 1, 78):
            lhs = sigma_vec(S, CD[(a, b)])
            rhs = {}
            for z, cz in S[b].items():
                for c, w in H[(z, a)].items():
                    nv = psub(rhs.get(c, Z2), pmul(cz, w))
                    if phas(nv):
                        rhs[c] = nv
                    elif c in rhs:
                        del rhs[c]
            if lhs != rhs:
                return inv, False
    return inv, True

ok("R2f the grading involution psi and the split sign map s are "
   "Jordan automorphisms of J preserving det (378 basis pairs, "
   "3303 test elements)", phi_checks("psi") and phi_checks("s"))
# ---- real fixed basis from a conjugation's coordinate matrix S
def real_fixed_basis(S):
    """Candidates u_a = e_a + sigma(e_a), w_a = i(e_a - sigma(e_a));
    extract a real-independent subset (returned as OPS-coord dicts)."""
    cands = []
    for a in range(78):
        u = {a: ONE2}
        for r, v in S[a].items():
            u[r] = padd(u.get(r, Z2), v)
            if not phas(u[r]):
                del u[r]
        cands.append(u)
        w = {a: (Fraction(0), Fraction(1))}
        for r, v in S[a].items():
            nv = pmul((Fraction(0), Fraction(1)), pneg(v))
            w[r] = padd(w.get(r, Z2), nv)
            if not phas(w[r]):
                del w[r]
        if w.get(a) == (Fraction(0), Fraction(1)) and len(w) == 1:
            pass
        cands.append(w)
    # real pivot extraction
    pivots = {}
    chosen = []
    for cand in cands:
        e = {}
        for r, v in cand.items():
            e[2 * r] = v[0]
            e[2 * r + 1] = v[1]
        e = {k: v for k, v in e.items() if v != 0}
        while True:
            pv = [v for v in e if v in pivots]
            if not pv:
                break
            p = pv[0]
            f = e[p]
            pe = pivots[p]
            for v, c in pe.items():
                nv = e.get(v, Fraction(0)) - f * c
                if nv != 0:
                    e[v] = nv
                elif v in e:
                    del e[v]
        if e:
            p = min(e)
            inv = Fraction(1) / e[p]
            pivots[p] = {v: c * inv for v, c in e.items()}
            chosen.append(cand)
    return chosen

def combo_ops(cf):
    out = {}
    for a, c in cf.items():
        if phas(c):
            out = mat_add(out, mat_scale(OPS[a][0], c))
    return out

def apply_op(T, v):
    out = {}
    for (r, c), val in T.items():
        if c in v:
            out[r] = padd(out.get(r, Z2), pmul(val, v[c]))
    return {k: x for k, x in out.items() if phas(x)}

def trace_prod(A, B):
    t = Z2
    for (r, c), v in A.items():
        w = B.get((c, r))
        if w is not None:
            t = padd(t, pmul(v, w))
    return t

def bdict(x, a):
    if x == a:
        return {}
    if x < a:
        return CD[(x, a)]
    return {c: pneg(v) for c, v in CD[(a, x)].items()}

def killing_ad(a, b):
    t = Z2
    for x in range(78):
        ba = bdict(x, a)
        if not ba:
            continue
        for y, v1 in ba.items():
            v2 = bdict(y, b).get(x)
            if v2 is not None:
                t = padd(t, pmul(v1, v2))
    return t

# scalar c with Killing = c * tr_27, from Z' (support-2 coords)
_zp_pairs = [(a, b) for a in ZPCOORD for b in ZPCOORD]
_num = Z2
for a, ca in ZPCOORD.items():
    for b, cb in ZPCOORD.items():
        _num = padd(_num, pmul(pmul(ca, cb), killing_ad(a, b)))
_den = trace_prod(ZP, ZP)
CKILL = pmul(_num, pinv(_den))
print(f"  Killing scalar c = {CKILL} (Killing = c * tr_27)", flush=True)
_kok = True
for a in range(78):
    for b in range(a, 78):
        if killing_ad(a, b) != pmul(CKILL, trace_prod(OPS[a][0], OPS[b][0])):
            _kok = False
ok("R3a Killing(X, Y) = c * tr_27(XY) on all 3081 basis pairs "
   "(adjoint trace from COMM vs 27-trace)", _kok)

# ================================================================ EIII geometry
print("EIII geometry: charges, symmetric pair, tangent module...",
      flush=True)
I2 = (Fraction(0), Fraction(1))
NEG1 = (Fraction(-1), Fraction(0))
THREE = (Fraction(3), Fraction(0))
MTHREE = (Fraction(-3), Fraction(0))

def unitv(a):
    return {a: ONE2}

# ---- G1: Peirce charges of the basis operators under ad(Z')
CHARGE = {}
for a in range(78):
    bc = bracket_cd(ZPCOORD, unitv(a))
    if not bc:
        CHARGE[a] = Z2
    elif len(bc) == 1 and a in bc:
        CHARGE[a] = bc[a]
    else:
        CHARGE[a] = None
K0 = [a for a in range(78) if CHARGE[a] == Z2]
PP = [a for a in range(78) if CHARGE[a] == THREE]
PM = [a for a in range(78) if CHARGE[a] == MTHREE]
ok("G1 every basis operator is ad(Z')-homogeneous, charges in "
   "{0, +3, -3} with counts (46, 16, 16): k = charge 0, p+ / p- the "
   "two 16s", all(CHARGE[a] is not None for a in range(78))
   and (len(K0), len(PP), len(PM)) == (46, 16, 16),
   f"{(len(K0), len(PP), len(PM))}")
KSET = set(K0)
PSET = set(PP) | set(PM)
PALL = PP + PM

# ---- G2: bracket charges add (symmetric-pair relations, complex level)
_g2 = all(CHARGE[a] is not None for a in range(78))
for (a, b), cd in CD.items():
    if CHARGE[a] is None or CHARGE[b] is None:
        continue
    s = padd(CHARGE[a], CHARGE[b])
    if s == Z2 or s == THREE or s == MTHREE:
        if any(CHARGE[c] != s for c in cd):
            _g2 = False
    elif cd:
        _g2 = False
ok("G2 charges add under brackets: [T_a, T_b] is supported on "
   "charge q_a + q_b, and is zero when q_a + q_b = +/-6; in "
   "particular [k, k] in k, [k, p+/-] in p+/-, [p+, p-] in k, "
   "[p+, p+] = [p-, p-] = 0", _g2)

# ---- G3: centre and derived algebra of k
eqs = []
for a in K0:
    for c in K0:
        row = {}
        for j, b in enumerate(K0):
            w = bdict(b, a).get(c)
            if w is not None:
                row[j] = w
        if row:
            eqs.append(row)
_centre = nullspace(eqs, 46)
ok("G3a centre of k is 1-dimensional (spanned by Z' via G1: "
   "[Z', k] = 0)", len(_centre) == 1, f"dim {len(_centre)}")
_der = []
for i, a in enumerate(K0):
    for b in K0[i + 1:]:
        if CD[(min(a, b), max(a, b))]:
            _der.append(CD[(min(a, b), max(a, b))])
_derrk = rank_of_vectors(_der, 78)[0]
ok("G3b derived algebra [k, k] has dim 45", _derrk == 45, f"{_derrk}")

# ---- G4: k is exactly the Peirce-piece-preserving operators
BLOCK = {i: (0 if i == 0 else (1 if 3 <= i <= 18 else 2))
         for i in range(27)}
eqs = []
for r in range(27):
    for cc in range(27):
        if BLOCK[r] != BLOCK[cc]:
            row = {}
            for a in range(78):
                v = OPS[a][0].get((r, cc))
                if v is not None:
                    row[a] = v
            if row:
                eqs.append(row)
_ns = nullspace(eqs, 78)
_units = [{a: ONE2} for a in K0]
ok("G4 the operators preserving the three Peirce pieces "
   "(J1, J1/2, J0) form a 46-dim space, equal to k",
   len(_ns) == 46 and rank_of_vectors(_ns + _units, 78)[0] == 46,
   f"nullity {len(_ns)}")

# ---- G5: Killing form on p: isotropy and nondegenerate pairing
def killing2(A, B):
    return pmul(CKILL, trace_prod(A, B))
_g5 = True
for a in PALL:
    for b in PALL:
        if padd(CHARGE[a], CHARGE[b]) != Z2:
            if phas(killing2(OPS[a][0], OPS[b][0])):
                _g5 = False
ok("G5a Killing vanishes on p+ x p+ and on p- x p- "
   "(Killing = 4 tr_27)", _g5)
_rows = []
for a in PP:
    _rows.append({j: v for j, v in enumerate(
        [killing2(OPS[a][0], OPS[b][0]) for b in PM]) if phas(v)})
ok("G5b the Killing pairing p+ x p- is nondegenerate (rank 16)",
   rank_of_vectors(_rows, 16)[0] == 16)

# ---- G6: p+ as a k-module
def restr(a, src, dst):
    pos = {x: i for i, x in enumerate(dst)}
    M = [[Z2] * len(src) for _ in dst]
    for j, x in enumerate(src):
        for c, w in bdict(a, x).items():
            if c in pos:
                M[pos[c]][j] = w
    return M
RACT = {a: restr(a, PP, PP) for a in K0}

def mat16_add(*Ms):
    out = [[Z2] * 16 for _ in range(16)]
    for M in Ms:
        for i in range(16):
            for j in range(16):
                out[i][j] = padd(out[i][j], M[i][j])
    return out

def mat16_scale(M, cf):
    return [[pmul(cf, M[i][j]) for j in range(16)] for i in range(16)]

def mat16_mul(A, B):
    out = [[Z2] * 16 for _ in range(16)]
    for i in range(16):
        for k in range(16):
            if phas(A[i][k]):
                for j in range(16):
                    if phas(B[k][j]):
                        out[i][j] = padd(out[i][j], pmul(A[i][k], B[k][j]))
    return out

RZ = [[Z2] * 16 for _ in range(16)]
for a, cf in ZPCOORD.items():
    RZ = mat16_add(RZ, mat16_scale(RACT[a], cf))
ok("G6a Z' acts on p+ as the scalar +3 (the U(1) charge of the "
   "tangent module)", all(RZ[i][j] == (THREE if i == j else Z2)
                           for i in range(16) for j in range(16)))

def comm_dim(idxs):
    eqs = []
    for a in idxs:
        R = RACT[a]
        for i in range(16):
            for j in range(16):
                row = {}
                for k in range(16):
                    if phas(R[k][j]):
                        row[i * 16 + k] = padd(row.get(i * 16 + k, Z2),
                                               R[k][j])
                    if phas(R[i][k]):
                        row[k * 16 + j] = psub(row.get(k * 16 + j, Z2),
                                               R[i][k])
                if row:
                    eqs.append(row)
    return len(nullspace(eqs, 256))
_cd = None
for _take in (4, 8, 16, 46):
    _cd = comm_dim(K0[:_take])
    if _cd == 1:
        break
ok("G6b p+ is irreducible as a k-module (commutant of the "
   "restricted action is the scalars)", _cd == 1, f"commutant dim {_cd}")

# Casimir on p+ with the J9 prescription: Gram = -tr_27 on the
# derived 45, duals, C = sum rho(T_a) rho(T^a)
_indep = []
for v in _der:
    if rank_of_vectors(_indep + [v], 78)[0] == len(_indep) + 1:
        _indep.append(v)
    if len(_indep) == 45:
        break
_ops45 = [combo_ops(v) for v in _indep]
import sympy as _sp
_G = _sp.Matrix(45, 45, lambda i, j: to_sympy(
    pneg(trace_prod(_ops45[i], _ops45[j]))))
_Ginv = _G.inv()
_dual = []
for i in range(45):
    dv = {}
    for j in range(45):
        cf = from_sympy(_Ginv[i, j])
        if phas(cf):
            for c, w in _indep[j].items():
                nv = padd(dv.get(c, Z2), pmul(cf, w))
                if phas(nv):
                    dv[c] = nv
                elif c in dv:
                    del dv[c]
    _dual.append(dv)

def ract_of(cf):
    out = [[Z2] * 16 for _ in range(16)]
    for a, w in cf.items():
        out = mat16_add(out, mat16_scale(RACT[a], w))
    return out
_CP = [[Z2] * 16 for _ in range(16)]
for i in range(45):
    _CP = mat16_add(_CP, mat16_mul(ract_of(_indep[i]),
                                   ract_of(_dual[i])))
_sc = _CP[0][0]
ok("G6c Casimir of k' on p+ (J9 prescription) is the scalar "
   "-15/8 -- the same value as the spinor piece J1/2 of J: p+ is "
   "the 16-dimensional spinor module of so(10)",
   all(_CP[i][j] == (_sc if i == j else Z2)
       for i in range(16) for j in range(16))
   and _sc == (Fraction(-15, 8), Fraction(0)), f"{_sc}")

# ---- G7: the real form sigma_EIII and its Cartan decomposition
print("EIII geometry: real form and Cartan decomposition...",
      flush=True)
SE = sigma_cols(lambda M: tau_op(theta_op(ad_phi(M, PHIS["psi"])), D0))
ok("G7a sigma_EIII preserves the algebra (span)", SE is not None)
_inv, _brk = sigma_checks(SE)
ok("G7b sigma_EIII is an involution respecting all 3003 brackets "
   "in coordinates", _inv and _brk)
ok("G7c sigma preserves k and p (coordinate supports)",
   all(set(SE[a]) <= KSET for a in K0)
   and all(set(SE[a]) <= PSET for a in PALL))
FB = real_fixed_basis(SE)
FBK = [v for v in FB if set(v) <= KSET]
FBP = [v for v in FB if set(v) <= PSET]
ok("G7d the 78 real fixed operators split as 46 in k and 32 in p "
   "(every fixed vector is homogeneous)", len(FB) == 78
   and len(FBK) + len(FBP) == 78
   and (len(FBK), len(FBP)) == (46, 32),
   f"{(len(FBK), len(FBP))}")
ok("G7e Ad_psi is +1 on k and -1 on p at operator level "
   "(the symmetric-pair involution)",
   all(ad_phi(OPS[a][0], PHIS["psi"]) ==
       {k: pmul(v, ONE2 if a in KSET else NEG1)
        for k, v in OPS[a][0].items()} for a in range(78)))

def km_sig(vecs):
    ops = [combo_ops(v) for v in vecs]
    n = len(vecs)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = killing2(ops[i], ops[j])
            assert v[1] == 0
            M[i][j] = M[j][i] = v[0]
    return signature(M)
ok("G7f Killing on k0 (the real isotropy algebra) is negative "
   "definite: signature (0, 46)", km_sig(FBK) == (0, 46, 0),
   f"{km_sig(FBK)}")
ok("G7g Killing on p0 is positive definite: signature (32, 0) "
   "-- Cartan decomposition of e6(-14), character -14 = 32 - 46",
   km_sig(FBP) == (32, 0, 0), f"{km_sig(FBP)}")

# ---- G8: the complex structure on p0
HCOORD = {a: pscale(pmul(v, I2), Fraction(1, 3))
          for a, v in ZPCOORD.items()}

def solve_in_basis(basis, target):
    """basis: lin-indep coord dicts; target: coord dict.
    Returns coefficient list (pairs) or None."""
    n = len(basis)
    cols = []
    for v in basis:
        cols.append([v.get(r, Z2) for r in range(78)])
    tgt = [target.get(r, Z2) for r in range(78)]
    rows = []
    for r in range(78):
        row = [cols[j][r] for j in range(n)] + [tgt[r]]
        if any(phas(x) for x in row):
            rows.append(row)
    piv = []
    pr = 0
    for c in range(n):
        pivrow = None
        for i in range(pr, len(rows)):
            if phas(rows[i][c]):
                pivrow = i
                break
        if pivrow is None:
            continue
        rows[pr], rows[pivrow] = rows[pivrow], rows[pr]
        iv = pinv(rows[pr][c])
        rows[pr] = [pmul(x, iv) for x in rows[pr]]
        for i in range(len(rows)):
            if i != pr and phas(rows[i][c]):
                f = rows[i][c]
                rows[i] = [psub(rows[i][k], pmul(f, rows[pr][k]))
                           for k in range(n + 1)]
        piv.append(c)
        pr += 1
        if pr == len(rows):
            break
    for i in range(pr, len(rows)):
        if not any(phas(x) for x in rows[i][:n]) and phas(rows[i][n]):
            return None
    if len(piv) < n:
        return None
    out = [Z2] * n
    for i, c in enumerate(piv):
        out[c] = rows[i][n]
    return out

_Jcols = []
_jok = True
for v in FBP:
    cf = solve_in_basis(FBP, bracket_cd(HCOORD, v))
    if cf is None or any(x[1] != 0 for x in cf):
        _jok = False
    _Jcols.append(cf)
ok("G8a ad(iZ'/3) preserves p0 and is real in the fixed basis "
   "(a complex structure on the real tangent space)", _jok)
Jmat = [[_Jcols[j][i][0] for j in range(32)] for i in range(32)]
_J2 = [[sum(Jmat[i][k] * Jmat[k][j] for k in range(32))
        for j in range(32)] for i in range(32)]
ok("G8b the complex structure squares to -id on p0",
   all(_J2[i][j] == (Fraction(-1) if i == j else Fraction(0))
       for i in range(32) for j in range(32)))
_Pops = [combo_ops(v) for v in FBP]
_KMP = [[killing2(_Pops[i], _Pops[j])[0] for j in range(32)]
        for i in range(32)]
_JKJ = [[sum(Jmat[k][i] * _KMP[k][l] * Jmat[l][j]
             for k in range(32) for l in range(32))
         for j in range(32)] for i in range(32)]
ok("G8c the complex structure is Killing-orthogonal: "
   "K(JX, JY) = K(X, Y) on p0 (Hermitian data)",
   _JKJ == _KMP)
ok("G8d iZ'/3 is central in k0: [iZ'/3, k0] = 0",
   all(bracket_cd(HCOORD, v) == {} for v in FBK))

# ---- G9: the compact dual shares k0
SC = sigma_cols(lambda M: tau_op(theta_op(M), D0))
ok("G9a sigma_compact preserves the algebra and is an involution "
   "in coordinates (its bracket verification is F3 of "
   "e6_real_forms.py)", SC is not None
   and cols_is_id(compose_cols(SC, conj_cols(SC))))
ok("G9b sigma_compact fixes k0 pointwise in coordinates: the "
   "compact form has the same real isotropy algebra k0",
   all(sigma_vec(SC, v) == v for v in FBK))
ok("G9c sigma_compact fixes i.p0 elementwise: the compact dual "
   "Cartan decomposition is k0 + i p0",
   all(sigma_vec(SC, {a: pmul(c, I2) for a, c in v.items()})
       == {a: pmul(c, I2) for a, c in v.items()} for v in FBP))

# ---- G10: brackets of the real decomposition
def supp_in(v, S):
    return all(a in S for a in v)
ok("G10a [k0, k0] in k0 (all 1035 fixed-basis pairs)",
   all(supp_in(bracket_cd(u, v), KSET)
       for i, u in enumerate(FBK) for v in FBK[i + 1:]))
ok("G10b [k0, p0] in p0 (all 46 x 32 fixed-basis pairs)",
   all(supp_in(bracket_cd(u, v), PSET) for u in FBK for v in FBP))
ok("G10c [p0, p0] in k0 (all 496 fixed-basis pairs)",
   all(supp_in(bracket_cd(u, v), KSET)
       for i, u in enumerate(FBP) for v in FBP[i + 1:]))

print(f"\n{len(FAIL)} failures / {NCHK[0]} checks")
if FAIL:
    raise SystemExit(1)
print("ALL PASS")
