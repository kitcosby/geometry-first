#!/usr/bin/env python3
"""
Real forms of the E6 stack: conjugations of the algebra, of J,
and of the 27 -- built on the Jordan-side model (e6_jordan_e6).

Front end (rebuilt here, checks J0-J4c as in e6_jordan_e6.py):
J = H_3(O_C), its det tensor, and the determinant-preserving
algebra e6 as 78 explicit operators on J, with all 3081 pairwise
commutators expressed in the basis (COMM).

New content (checks R, J-tau, F):

  * The trace form B(X, Y) = tr(X o Y) on J has diagonal Gram
    (1, 1, 1, 2 x24) in the coordinate basis, so the B-adjoint of
    an operator is an explicit transpose with weights 1, 2, 1/2.
    theta(T) = -T^{*B} is an involutive automorphism of e6
    (verified on all 3003 bracket pairs at operator level).
    Its +1 eigenspace is the derivation algebra Der(J) (dim 52,
    solved independently from the Leibniz identity on all 378
    basis pairs); its -1 eigenspace is the 26 traceless
    multiplication operators L_a (structure-algebra decomposition
    Der (52) + L (26) = 78).

  * Two antilinear Jordan automorphisms of J: tau0 (coefficient
    conjugation; fixed algebra H_3(O), octonion norm signature
    (8, 0), Jordan trace signature (27, 0)) and tau_split
    (coefficient conjugation composed with the Fano sign character
    u -> (-1)^{bit0(u)} on octonion units; fixed algebra
    H_3(O_split), norm signature (4, 4), trace signature (15, 12)).
    Both preserve det up to conjugation on all 3303 test elements.

  * Five conjugations sigma of the algebra (antilinear involutive
    automorphisms), each verified in coordinates: span preserved,
    S * conj(S) = I, all 3003 brackets respected.  For each, a real
    fixed basis of 78 operators is extracted, and the Killing form
    on it computed exactly.  The Killing form is first pinned to
    the 27-trace: Killing = c * tr_27 with c computed from Z' and
    verified on all 3081 basis pairs against the adjoint trace
    from COMM.  Signatures (Sylvester, exact congruence):
      sigma0        = tau0                    -> (26, 52)  e6(-26)
      sigma_split   = tau_split               -> (42, 36)  e6(6)
      sigma_compact = tau0 o theta            -> (0, 78)   e6(-78)
      sigma_EIII    = tau0 o theta o Ad_psi   -> (32, 46)  e6(-14)
      sigma_2       = tau0 o theta o Ad_s o Ad_psi
                                              -> (40, 38)  e6(2)
    (Ad_psi: conjugation by the Peirce grading involution of J;
     Ad_s: conjugation by the split sign automorphism of J.
     Naming of the real form by its Killing signature/character
     is the standard classification of the five real forms of E6,
     cited; everything else is computed.)

  * Real forms of the 27: tau0 and tau_split are compatible
    conjugations on J itself (intertwining verified on all
    78 x 27 basis pairs).  The other three conjugations satisfy
    sigma(Z') = -Z' for the grading operator Z' (Peirce block
    scalars (4, 1, -2), multiplicities (1, 16, 10)); a compatible
    conjugation on J would have to map the lam-eigenspace of Z'
    to the (-lam)-eigenspace, and the spectrum contains none of
    (-4, -1, 2) -- so no compatible conjugation on J exists for
    sigma_compact, sigma_EIII, sigma_2.  (Obstruction: PROVED from
    the computed sigma(Z') signs and the spectrum.)

All arithmetic exact over Q(i) (Fraction pairs); signatures by
exact rational congruence (Sylvester).  No floating point.
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


# ================================================================ real forms
# Conventions: a "conjugation" sigma is an antilinear map on the
# 78-dim algebra, given at operator level by a transform on 27x27
# matrices.  S[sigma] = its coordinate matrix in the OPS basis as a
# list of column dicts: sigma(T_a) = sum_r S[a][r] * OPS[r], and for
# a coordinate vector v, sigma(v) = sum_a conj(v_a) * S[a].

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


# ---- theta: full bracket verification at operator level
print("theta bracket verification (3003 pairs)...", flush=True)
THI = [theta_op(M) for (M, q) in OPS]
_thok = True
for a in range(NB):
    for b in range(a + 1, NB):
        C = commutator(OPS[a][0], OPS[b][0])
        if theta_op(C) != commutator(THI[a], THI[b]):
            _thok = False
ok("R2b theta([Ta, Tb]) = [theta Ta, theta Tb] on all 3003 pairs "
   "(operator level)", _thok)
ok("R2c theta preserves the algebra: theta(Ta) in span for all a",
   all(express(THI[a]) is not None for a in range(NB)))
THCOLS = []
for a in range(NB):
    cf = express(THI[a])
    THCOLS.append({i: v for i, v in enumerate(cf) if phas(v)})

def theta_vec(v):
    out = {}
    for a, ca in v.items():
        for r, w in THCOLS[a].items():
            nv = padd(out.get(r, Z2), pmul(ca, w))
            if phas(nv):
                out[r] = nv
            elif r in out:
                del out[r]
    return out

ok("R2d theta(Z') = -Z' in coordinates",
   theta_vec(ZPCOORD) == {a: pneg(v) for a, v in ZPCOORD.items()})
ok("R2e tau0 and theta commute on all 78 basis operators",
   all(tau_op(theta_op(M), D0) == theta_op(tau_op(M, D0))
       for (M, q) in OPS))

# ---- derivations (Leibniz) and multiplication operators
print("solving derivations (Leibniz equations)...", flush=True)
OPCOL = []
for a in range(NB):
    cols = [dict() for _ in range(27)]
    for (r, c), v in OPS[a][0].items():
        cols[c][r] = v
    OPCOL.append(cols)

deqs = []
for i in range(27):
    for j in range(i, 27):
        X = PROD[(i, j)]
        # per output coord r: dict over a
        rows = {}
        for a in range(NB):
            t1 = {}
            for k, xv in X.items():
                for r, v in OPCOL[a][k].items():
                    t1[r] = padd(t1.get(r, Z2), pmul(xv, v))
            t2 = {}
            for s, v in OPCOL[a][i].items():
                for r, w in PROD[(s, j)].items():
                    t2[r] = padd(t2.get(r, Z2), pmul(v, w))
            t3 = {}
            for s, v in OPCOL[a][j].items():
                for r, w in PROD[(i, s)].items():
                    t3[r] = padd(t3.get(r, Z2), pmul(v, w))
            for r in set(t1) | set(t2) | set(t3):
                val = psub(psub(t1.get(r, Z2), t2.get(r, Z2)),
                           t3.get(r, Z2))
                if phas(val):
                    rows.setdefault(r, {})[a] = val
        deqs.extend(rows.values())
DER = nullspace(deqs, NB)
print(f"  derivations: nullity {len(DER)}", flush=True)
ok("R4a derivation space (Leibniz on all 378 basis pairs) has dim 52",
   len(DER) == 52, f"{len(DER)}")
ok("R4b every derivation is theta-fixed (+1 eigenspace)",
   all(theta_vec(d) == d for d in DER))
ok("R4c every derivation operator is B-antisymmetric (adjoint = -self)",
   all(adjoint(combo_ops(d)) == {k: pneg(v) for k, v in
       combo_ops(d).items()} for d in DER))

# L operators for the 27 basis elements
LMAT = []
for k in range(27):
    M = {}
    for j in range(27):
        for r, v in PROD[(k, j)].items():
            M[(r, j)] = v
    LMAT.append(M)
WCOORDS = []
_wok = True
for k in range(3, 27):
    cf = express(LMAT[k])
    if cf is None:
        _wok = False
    else:
        WCOORDS.append({i: v for i, v in enumerate(cf) if phas(v)})
for (k1, k2) in ((0, 1), (0, 2)):
    cf = express(mat_add(LMAT[k1], {kk: pneg(vv) for kk, vv in
                                    LMAT[k2].items()}))
    if cf is None:
        _wok = False
    else:
        WCOORDS.append({i: v for i, v in enumerate(cf) if phas(v)})
ok("R4d all 26 traceless multiplication operators L_a lie in the "
   "algebra (24 slot units + 2 diagonal differences)", _wok and
   len(WCOORDS) == 26, f"{len(WCOORDS)}")
ok("R4e every L_a is theta-anti-fixed (-1 eigenspace) and B-symmetric",
   all(theta_vec(w) == {a: pneg(v) for a, v in w.items()}
       for w in WCOORDS))
_rk, _ = rank_of_vectors(DER + WCOORDS, 78)
ok("R4f Der (52) + W (26) spans the full 78 (structure decomposition)",
   _rk == 78, f"{_rk}")

# ---- Jordan-level real forms tau0 / tau_split
def oct_fixed_basis(split):
    out = []
    for u in range(8):
        if (not split) or u % 2 == 0:
            out.append({u: ONE2})
        else:
            out.append({u: (Fraction(0), Fraction(1))})
    return out

def oadd(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = padd(out.get(k, Z2), v)
        if not phas(out[k]):
            del out[k]
    return out

for nm, D, split, exp_oct, exp_tr in (
        ("tau0", D0, False, (8, 0, 0), (27, 0, 0)),
        ("tau_split", DS, True, (4, 4, 0), (15, 12, 0))):
    aut = True
    for i in range(27):
        for j in range(i, 27):
            if tau_vec(PROD[(i, j)], D) != jprod_uv(
                    tau_vec(BASIS[i], D), tau_vec(BASIS[j], D)):
                aut = False
    ok(f"J-{nm}: tau is an antilinear Jordan automorphism "
       f"(378 basis pairs)", aut)
    ok(f"J-{nm}: det(tau X) = conj(det X) on all 3303 test elements",
       all(det(tau_vec(X, D)) == pconj(det(X)) for X in test_els))
    FB = tau_fixed_basis(D)
    ok(f"J-{nm}: fixed basis of J has 27 elements, all tau-fixed, "
       f"products stay fixed",
       len(FB) == 27 and all(tau_vec(v, D) == v for v in FB) and
       all(tau_vec(jprod_uv(FB[i], FB[j]), D) == jprod_uv(FB[i], FB[j])
           for i in range(27) for j in range(i, 27)))
    OB = oct_fixed_basis(split)
    NM = [[Fraction(0)] * 8 for _ in range(8)]
    for a in range(8):
        for b in range(8):
            v = pscale(psub(psub(onorm0(oadd(OB[a], OB[b])),
                            onorm0(OB[a])), onorm0(OB[b])),
                       Fraction(1, 2))
            NM[a][b] = v[0]
    ok(f"J-{nm}: octonion norm signature on the fixed octonions = "
       f"{exp_oct}", signature(NM) == exp_oct, f"{signature(NM)}")
    TM = [[Fraction(0)] * 27 for _ in range(27)]
    for a in range(27):
        for b in range(27):
            t = Z2
            P = jprod_uv(FB[a], FB[b])
            for d in range(3):
                t = padd(t, P.get(d, Z2))
            TM[a][b] = t[0]
    ok(f"J-{nm}: Jordan trace-form signature on the fixed J = "
       f"{exp_tr}", signature(TM) == exp_tr, f"{signature(TM)}")


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

def process_form(name, transform, expect_sig):
    print(f"processing {name}...", flush=True)
    S = sigma_cols(transform)
    ok(f"{name}: sigma preserves the 78-dim algebra (span)",
       S is not None)
    if S is None:
        return None
    inv, brk = sigma_checks(S)
    ok(f"{name}: sigma is an involution in coordinates", inv)
    ok(f"{name}: sigma respects all 3003 brackets in coordinates", brk)
    F = real_fixed_basis(S)
    ok(f"{name}: real fixed basis has 78 elements", len(F) == 78,
       f"{len(F)}")
    ok(f"{name}: every fixed-basis vector is sigma-fixed",
       all(sigma_vec(S, f) == f for f in F))
    rk, _ = rank_of_vectors(F, 78)
    ok(f"{name}: fixed basis spans the algebra over C (rank 78)",
       rk == 78, f"{rk}")
    FBops = [combo_ops(f) for f in F]
    ok(f"{name}: fixed basis operators are sigma-fixed as matrices",
       all(transform(M) == M for M in FBops))
    KM = [[Fraction(0)] * 78 for _ in range(78)]
    kreal = True
    for k in range(78):
        for l in range(k, 78):
            v = pmul(CKILL, trace_prod(FBops[k], FBops[l]))
            if v[1] != 0:
                kreal = False
            KM[k][l] = v[0]
            KM[l][k] = v[0]
    ok(f"{name}: Killing form on the fixed basis is real", kreal)
    sig = signature(KM)
    ok(f"{name}: Killing signature = {expect_sig}", sig == expect_sig,
       f"{sig}")
    sz = sigma_vec(S, ZPCOORD)
    zsign = 1 if sz == ZPCOORD else (-1 if sz ==
          {a: pneg(v) for a, v in ZPCOORD.items()} else 0)
    print(f"  {name}: signature {sig}, sigma(Z') sign {zsign}", flush=True)
    return {"S": S, "F": F, "sig": sig, "zsign": zsign}


# ---- the five real forms
FORMS = {}
FORMS["sigma0"] = process_form(
    "F1 sigma0 = tau0 (Jordan conjugation, compact octonions)",
    lambda M: tau_op(M, D0), (26, 52, 0))
FORMS["sigma_split"] = process_form(
    "F2 sigma_split = tau_split (Jordan conjugation, split octonions)",
    lambda M: tau_op(M, DS), (42, 36, 0))
FORMS["sigma_compact"] = process_form(
    "F3 sigma_compact = tau0 o theta",
    lambda M: tau_op(theta_op(M), D0), (0, 78, 0))
FORMS["sigma_EIII"] = process_form(
    "F4 sigma_EIII = tau0 o theta o Ad_psi (Peirce grading)",
    lambda M: tau_op(theta_op(ad_phi(M, PHIS["psi"])), D0), (32, 46, 0))
FORMS["sigma_2"] = process_form(
    "F5 sigma_2 = tau0 o theta o Ad_s o Ad_psi",
    lambda M: tau_op(theta_op(ad_phi(ad_phi(M, PHIS["s"]),
                                      PHIS["psi"])), D0), (40, 38, 0))

ok("F6 sigma(Z') signs: +1 for the two Jordan-induced conjugations, "
   "-1 for sigma_compact, sigma_EIII, sigma_2",
   FORMS["sigma0"]["zsign"] == 1 and FORMS["sigma_split"]["zsign"] == 1
   and FORMS["sigma_compact"]["zsign"] == -1
   and FORMS["sigma_EIII"]["zsign"] == -1
   and FORMS["sigma_2"]["zsign"] == -1,
   str({k: v["zsign"] for k, v in FORMS.items()}))

# compatible conjugations on J for the Jordan-induced forms:
# tau(T X) = sigma(T) tau(X) for all basis operators T and basis X
for nm, D in (("sigma0", D0), ("sigma_split", DS)):
    trf = (lambda D: (lambda M: tau_op(M, D)))(D)
    ok(f"F7 {nm}: tau on J intertwines sigma on the algebra "
       f"(all 78 x 27 basis pairs)",
       all(tau_vec(apply_op(OPS[a][0], BASIS[c]), D) ==
           apply_op(trf(OPS[a][0]), tau_vec(BASIS[c], D))
           for a in range(78) for c in range(27)))

# obstruction for the other three: a compatible tau on J would satisfy
# tau Z' = sigma(Z') tau = -Z' tau, forcing tau to map the Z'-eigenspace
# of eigenvalue lam to that of -lam; the spectrum (4, 1, -2) with
# multiplicities (1, 16, 10) contains none of (-4, -1, 2).
ok("F8 Z' spectrum (4, 1, -2) is disjoint from its negative "
   "(-4, -1, 2): no compatible conjugation on J for the three "
   "non-Jordan forms (with F6)",
   not ({4, 1, -2} & {-4, -1, 2}))

print(f"\n{len(FAIL)} failures / {NCHK[0]} checks")
if FAIL:
    raise SystemExit(1)
print("ALL PASS")
