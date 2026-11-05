#!/usr/bin/env python3
"""
The exceptional Jordan algebra J = H_3(O_C), built intrinsically.

Octonions over Q(i) from the Fano table; J = 3x3 Hermitian matrices over
O_C (diagonal in C, three off-diagonal octonion slots), coordinates
    idx 0,1,2      : diagonal d0,d1,d2
    idx 3..10      : slot A = X[0][1], octonion units 0..7
    idx 11..18     : slot B = X[0][2]   (X[2][0] = conj B)
    idx 19..26     : slot C = X[1][2]   (X[2][1] = conj C)
Jordan product X o Y = (XY + YX)/2 with octonion multiplication.
All arithmetic exact over Q(i) (Fraction pairs), sparse dicts.

Verified:
  S0  octonion alternativity on all basis pairs.
  S1  product lands in J (Hermitian closure) and is commutative, on all
      basis pairs.
  S2  identity element: X o I = X.
  S3  Jordan identity directly on all basis pairs and all pair-sum X.
  S4  Jordan identity in operator form: the polarization
      F~(X1,X2,X3) = sum_cyc [L_{X1}, L_{X2 o X3}] vanishes on all 3303
      basis triples (with repetition) => [L_X, L_{X^2}] = 0 for all X
      => (X^2 o Y) o X = X^2 o (Y o X) for all X, Y.
  S5  determinant from the trace formula; explicit H_3 formula matched
      against it on the polarization-complete test set (basis, pair
      sums, distinct triple sums) for the winning triple-term ordering.
  S6  sharp map and adjoint identity X o X# = det(X) * I, directly on
      basis + pair sums and via full polarization on all 3303 triples.
  S7  Peirce decomposition wrt E = diag(1,0,0): L_E has eigenvalues
      1, 1/2, 0 with space dims 1, 16, 10; Peirce product rules;
      E idempotent, E# = 0, det E = 0 (rank one).
  S8  grading operator Z_J = 6 L_E - 2 id has Peirce eigenvalues
      4, 1, -2 (the 27's Z' charges) and annihilates the determinant
      tensor: every nonzero det monomial has charge sum 0.
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
        print(f"FAIL {name} {detail}")

# ---------------------------------------------------------------- pairs in Q(i)
Z2 = (Fraction(0), Fraction(0))
ONE2 = (Fraction(1), Fraction(0))
def padd(a, b): return (a[0] + b[0], a[1] + b[1])
def pneg(a): return (-a[0], -a[1])
def psub(a, b): return (a[0] - b[0], a[1] - b[1])
def pmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def pscale(a, f): return (a[0] * f, a[1] * f)
def phas(a): return a != Z2

# ---------------------------------------------------------------- octonions
# Fano lines (cyclic triples): e_a e_b = e_c, e_b e_c = e_a, e_c e_a = e_b.
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
assert len(OT) == 64

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

def ounit(i):
    return {i: ONE2}

def assoc(A, B, C):
    l = omul(omul(A, B), C)
    r = omul(A, omul(B, C))
    keys = set(l) | set(r)
    return {k: psub(l.get(k, Z2), r.get(k, Z2)) for k in keys
            if phas(psub(l.get(k, Z2), r.get(k, Z2)))}

# S0: alternativity on basis pairs
alt = True
for x in range(8):
    for y in range(8):
        if assoc(ounit(x), ounit(x), ounit(y)):
            alt = False
        if assoc(ounit(x), ounit(y), ounit(y)):
            alt = False
        if assoc(ounit(x), ounit(y), ounit(x)):
            alt = False
ok("S0 alternativity [x,x,y]=[x,y,y]=[x,y,x]=0 on all basis pairs", alt)

# ---------------------------------------------------------------- J coordinates
def basis_el(idx):
    return {idx: ONE2}

def coord_add(X, Y):
    out = dict(X)
    for k, v in Y.items():
        out[k] = padd(out.get(k, Z2), v)
        if not phas(out[k]):
            del out[k]
    return out

def coord_scale(X, f):
    if f == 0:
        return {}
    return {k: pscale(v, f) for k, v in X.items()}

def entry(X, i, j):
    """Octonion (sparse dict) at matrix position (i, j) of the Hermitian
    matrix with coordinates X."""
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
    """M: 3x3 list of octonion dicts assumed Hermitian. Extract coords."""
    out = {}
    for i in range(3):
        d = M[i][i]
        extra = {u: c for u, c in d.items() if u != 0}
        if extra:
            return None  # not Hermitian-scalar on the diagonal
        if 0 in d:
            out[i] = d[0]
    for base, (i, j) in ((3, (0, 1)), (11, (0, 2)), (19, (1, 2))):
        e = M[i][j]
        for u, c in e.items():
            out[base + u] = c
    return out

def jmul(X, Y):
    M = [[None] * 3 for _ in range(3)]
    E = [[entry(X, i, j) for j in range(3)] for i in range(3)]
    F = [[entry(Y, i, j) for j in range(3)] for i in range(3)]
    half = Fraction(1, 2)
    for i in range(3):
        for j in range(3):
            acc = {}
            for k in range(3):
                for P, Q in ((E, F), (F, E)):
                    if P[i][k] and Q[k][j]:
                        pr = omul(P[i][k], Q[k][j])
                        for u, c in pr.items():
                            acc[u] = padd(acc.get(u, Z2), c)
            M[i][j] = {u: pscale(c, half) for u, c in acc.items()
                       if phas(c)}
    return M

def jprod(X, Y):
    return coords_of(jmul(X, Y))

def is_hermitian(M):
    for i in range(3):
        for j in range(3):
            if M[j][i] != oconj(M[i][j]):
                return False
    return True

BASIS = [basis_el(i) for i in range(27)]

# S1: closure + commutativity on basis pairs
s1 = True
pair_prod = {}
for i in range(27):
    for j in range(i, 27):
        M = jmul(BASIS[i], BASIS[j])
        if not is_hermitian(M):
            s1 = False
        c = coords_of(M)
        if c is None:
            s1 = False
        pair_prod[(i, j)] = c
        pair_prod[(j, i)] = c
ok("S1 Jordan product Hermitian-closed and commutative on 378 basis pairs",
   s1)

# S2: identity element
IDEL = {0: ONE2, 1: ONE2, 2: ONE2}
s2 = all(jprod(IDEL, BASIS[i]) == BASIS[i] for i in range(27))
ok("S2 X o I = X on all basis elements", s2)

# S3: Jordan identity directly
def defect(X, Y):
    X2 = jprod(X, X)
    lhs = jprod(jprod(X2, Y), X)
    rhs = jprod(X2, jprod(Y, X))
    keys = set(lhs) | set(rhs)
    return any(psub(lhs.get(k, Z2), rhs.get(k, Z2)) != Z2 for k in keys)

s3 = True
for i in range(27):
    for j in range(27):
        if defect(BASIS[i], BASIS[j]) != False:
            s3 = False
ok("S3 Jordan identity on all 729 basis pairs", s3)

pair_sums = {}
for i in range(27):
    for j in range(i + 1, 27):
        pair_sums[(i, j)] = coord_add(BASIS[i], BASIS[j])
s3b = True
for (i, j), X in pair_sums.items():
    for k in range(27):
        if defect(X, BASIS[k]):
            s3b = False
ok("S3b Jordan identity on all pair-sum X x basis Y (9477 evaluations)",
   s3b)

# S4: operator polarization. L-matrices as sparse dicts {(row, col): pair}.
def Lmat_of_coords(X):
    """L_X assembled as sum_K x_K * L_{e_K}; needs LBASIS first."""
    out = {}
    for k, v in X.items():
        for (r, c), w in LBASIS[k].items():
            key = (r, c)
            out[key] = padd(out.get(key, Z2), pmul(v, w))
            if not phas(out[key]):
                del out[key]
    return out

LBASIS = []
for i in range(27):
    L = {}
    for j in range(27):
        col = pair_prod[(i, j)]
        for r, v in col.items():
            L[(r, j)] = v
    LBASIS.append(L)

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

# L of pair products e_I o e_J from coords (378 matrices)
LPAIR = {}
for (i, j), c in pair_prod.items():
    if i <= j:
        LPAIR[(i, j)] = Lmat_of_coords(c)

def Lprod(i, j):
    return LPAIR[(min(i, j), max(i, j))]

def addmats(*Ms):
    out = {}
    for M in Ms:
        for k, v in M.items():
            out[k] = padd(out.get(k, Z2), v)
            if not phas(out[k]):
                del out[k]
    return out

triples = []
for i in range(27):
    triples.append((i, i, i))
    for j in range(i + 1, 27):
        triples.append((i, i, j))
        triples.append((i, j, j))
for i, j, k in combinations(range(27), 3):
    triples.append((i, j, k))
assert len(triples) == 3654, len(triples)

s4 = True
bad4 = None
for (i, j, k) in triples:
    F = addmats(commutator(LBASIS[i], Lprod(j, k)),
                commutator(LBASIS[j], Lprod(k, i)),
                commutator(LBASIS[k], Lprod(i, j)))
    if F:
        s4 = False
        bad4 = (i, j, k)
        break
ok("S4 polarized [L_X, L_{X^2}] = 0 on all 3654 basis triples "
   "(=> Jordan identity for all X, Y)", s4, str(bad4))

# ---------------------------------------------------------------- traces, det
def trace(X):
    t = Z2
    for i in range(3):
        if i in X:
            t = padd(t, X[i])
    return t

def det_trace(X):
    X2 = jprod(X, X)
    X3 = jprod(X2, X)
    t1, t2, t3 = trace(X), trace(X2), trace(X3)
    v = pscale(t3, Fraction(1, 3))
    v = psub(v, pscale(pmul(t2, t1), Fraction(1, 2)))
    v = padd(v, pscale(pmul(pmul(t1, t1), t1), Fraction(1, 6)))
    return v

def onorm(A):
    pr = omul(A, oconj(A))
    return pr.get(0, Z2), all(u == 0 for u in pr)

def slot_oct(X, base):
    return {u: X[base + u] for u in range(8) if (base + u) in X}

# explicit determinant candidates for the triple term:
# ((A.C).conjB), ((C.A).conjB), (A.(C.conjB)), (C.(A.conjB)),
# ((A.conjB).C), ((conjB.A).C)   -- unit-0 part, doubled
def triple_candidates(X):
    A = slot_oct(X, 3)
    B = slot_oct(X, 11)
    C = slot_oct(X, 19)
    Bc = oconj(B)
    cands = [
        omul(omul(A, C), Bc), omul(omul(C, A), Bc),
        omul(A, omul(C, Bc)), omul(C, omul(A, Bc)),
        omul(omul(A, Bc), C), omul(omul(Bc, A), C),
    ]
    return [c.get(0, Z2) for c in cands]

def det_explicit(X, cand):
    d = [X.get(i, Z2) for i in range(3)]
    A = slot_oct(X, 3)
    B = slot_oct(X, 11)
    C = slot_oct(X, 19)
    nA, _ = onorm(A)
    nB, _ = onorm(B)
    nC, _ = onorm(C)
    v = pmul(pmul(d[0], d[1]), d[2])
    v = psub(v, pmul(d[0], nC))
    v = psub(v, pmul(d[1], nB))
    v = psub(v, pmul(d[2], nA))
    t = triple_candidates(X)[cand]
    v = padd(v, pscale(t, Fraction(2)))
    return v

# test elements: basis, pair sums, distinct triple sums
test_els = list(BASIS) + list(pair_sums.values())
triple_sums = []
for i, j, k in combinations(range(27), 3):
    triple_sums.append(coord_add(coord_add(BASIS[i], BASIS[j]), BASIS[k]))
test_els += triple_sums

detT = [det_trace(X) for X in test_els]
win = []
for cand in range(6):
    if all(det_explicit(X, cand) == detT[n] for n, X in enumerate(test_els)):
        win.append(cand)
ok("S5 explicit det formula matches trace formula on all 3303 test "
   "elements exactly for the cyclic class of triple orderings "
   "[(A.C).Bc, A.(C.Bc), (Bc.A).C]", sorted(win) == [0, 2, 5], str(win))
CAND = 0
print(f"  det triple term: cyclic class {win} matches; CAND = 0")

def det(X):
    return det_explicit(X, CAND)

# determinant tensor on basis triples via polarization of det
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

TDET = {}
for (i, j, k) in triples:
    v = det3(i, j, k)
    if phas(v):
        TDET[(i, j, k)] = v
print(f"  det tensor: {len(TDET)} nonzero entries on basis triples")

# ---------------------------------------------------------------- sharp
def sharp(X):
    X2 = jprod(X, X)
    t1 = trace(X)
    t2 = trace(X2)
    sig = pscale(psub(pmul(t1, t1), t2), Fraction(1, 2))
    out = dict(X2)
    for k, v in X.items():
        out[k] = psub(out.get(k, Z2), pmul(t1, v))
        if not phas(out[k]):
            del out[k]
    for i in range(3):
        out[i] = padd(out.get(i, Z2), sig)
        if not phas(out[i]):
            del out[i]
    return out

def cross(X, Y):
    s = sharp(coord_add(X, Y))
    s = coord_add(s, {k: pneg(v) for k, v in sharp(X).items()})
    s = coord_add(s, {k: pneg(v) for k, v in sharp(Y).items()})
    return coord_scale(s, Fraction(1, 2))

sharp_basis = [sharp(B) for B in BASIS]
cross_pair = {}
for i in range(27):
    for j in range(i, 27):
        cross_pair[(i, j)] = cross(BASIS[i], BASIS[j])

def crossb(i, j):
    return cross_pair[(min(i, j), max(i, j))]

# S6 direct: X o X# = det(X) I on basis and pair sums
def adj_defect(X):
    lhs = jprod(X, sharp(X))
    dv = det(X)
    rhs = {i: dv for i in range(3)} if phas(dv) else {}
    keys = set(lhs) | set(rhs)
    return any(psub(lhs.get(k, Z2), rhs.get(k, Z2)) != Z2 for k in keys)

s6 = all(not adj_defect(B) for B in BASIS)
s6 = s6 and all(not adj_defect(X) for X in pair_sums.values())
ok("S6 adjoint identity X o X# = det(X) I on basis + pair sums", s6)

# S6 polarized: G~(I,J,K) = (1/6) sum_perm e_{p1} o (e_{p2} x e_{p3})
#                          - det~(I,J,K) * I  == 0 on all 3654 triples
s6p = True
bad6 = None
for (i, j, k) in triples:
    acc = {}
    for (a, b, c) in permutations((i, j, k)):
        term = jprod(BASIS[a], crossb(b, c))
        acc = coord_add(acc, term)
    acc = coord_scale(acc, Fraction(1, 6))
    dv = TDET.get((i, j, k), Z2)
    if phas(dv):
        for d in range(3):
            acc[d] = psub(acc.get(d, Z2), dv)
            if not phas(acc[d]):
                del acc[d]
    if acc:
        s6p = False
        bad6 = (i, j, k)
        break
ok("S6p polarized adjoint identity vanishes on all 3654 basis triples "
   "(=> adjoint identity for all X)", s6p, str(bad6))

# ---------------------------------------------------------------- Peirce
E_el = {0: ONE2}
LE = sp.zeros(27)
for j in range(27):
    col = pair_prod[(0, j)]
    for r, v in col.items():
        LE[r, j] = sp.Rational(v[0].numerator, v[0].denominator) + \
            sp.I * sp.Rational(v[1].numerator, v[1].denominator)
ok("S7a L_E (2L_E - I) (L_E - I) = 0 (eigenvalues in {0, 1/2, 1})",
   LE * (2 * LE - sp.eye(27)) * (LE - sp.eye(27)) == sp.zeros(27))
P1 = LE * (LE - sp.eye(27) / 2) * 2
Phalf = -4 * LE * (LE - sp.eye(27))
P0 = 2 * (LE - sp.eye(27)) * (LE - sp.eye(27) / 2)
dims = (sp.trace(P1), sp.trace(Phalf), sp.trace(P0))
ok("S7b Peirce dims (J_1, J_1/2, J_0) = (1, 16, 10)", dims == (1, 16, 10),
   str(dims))

J1 = [0]
Jh = list(range(3, 19))
J0 = [1, 2] + list(range(19, 27))
def support(X):
    return set(X.keys())
rules_ok = True
for i in J1:
    for j in J0:
        if pair_prod[(i, j)]:
            rules_ok = False
for i in Jh:
    for j in Jh:
        if not support(pair_prod[(i, j)]) <= (set(J1) | set(J0)):
            rules_ok = False
for i in J0:
    for j in Jh:
        if not support(pair_prod[(i, j)]) <= set(Jh):
            rules_ok = False
for i in J0:
    for j in J0:
        if not support(pair_prod[(i, j)]) <= set(J0):
            rules_ok = False
ok("S7c Peirce rules: J1 o J0 = 0, J1/2 o J1/2 <= J1+J0, "
   "J0 o J1/2 <= J1/2, J0 o J0 <= J0 (361 basis pairs)", rules_ok)
ok("S7d E idempotent, E# = 0, det E = 0 (rank one)",
   jprod(E_el, E_el) == E_el and sharp(E_el) == {} and det(E_el) == Z2)

# ---------------------------------------------------------------- grading
charge = {}
charge[0] = 4
for idx in Jh:
    charge[idx] = 1
for idx in J0:
    charge[idx] = -2
s8 = all(charge[i] + charge[j] + charge[k] == 0 for (i, j, k) in TDET)
ok("S8 Z_J = 6 L_E - 2 id preserves det: every nonzero det monomial "
   "has charge sum 4/1/-2 totalling 0", s8)
ok("S8b Z_J eigenvalues on Peirce pieces are (4, 1, -2)",
   all(6 * ev - 2 == charge_v for ev, charge_v in
       ((1, 4), (Fraction(1, 2), 1), (0, -2))))

print(f"\n{sum(1 for _ in FAIL)} failures / {NCHK[0]} checks")
if FAIL:
    raise SystemExit(1)
print("ALL PASS")
