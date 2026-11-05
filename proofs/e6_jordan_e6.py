#!/usr/bin/env python3
"""
e6 on the Jordan side: the determinant-preserving algebra of
J = H_3(O_C) (built in e6_jordan.py), solved exactly over Q(i).

The det tensor Tdet (89 nonzero entries on basis triples) is
recomputed from the Jordan front-end.  An operator T in End(J)
preserves det infinitesimally iff delta_T Tdet = 0.  Because Z_J =
6 L_E - 2 id grades End(J) by Peirce-charge differences
{0, +-3, +-6}, the system delta_T Tdet = 0 splits into five
independent charge sectors, solved separately by sparse exact
Gaussian elimination over Q(i) (Fraction pairs).

Verified:
  J1  sector dimensions: charge 0 -> 46, +3 -> 16, -3 -> 16,
      +6 -> 0, -6 -> 0; total 78.
  J2  every assembled basis operator satisfies delta_T det = 0 by
      direct tensor contraction (solver cross-check).
  J3  closure spot-verification: delta_[Ta,Tb] det = 0 for 780
      commutators (10 fixed basis ops against all 78).  (General
      closure is the delta-homomorphism identity, proved in the
      note; the structurally distinguished brackets are verified
      in full in J6-J8.)
  J4  centre of the full 78 is 0; [g, g] spans the full 78
      (perfectness), from all 3081 commutators via a pivot minor.
  J5  the line stabilizer of E is the maximal parabolic of dim
      62 = Levi (charge-0 sector, 46) (+) nilradical P+ (16); the
      Levi's centre is 1-dimensional, spanned by Zhat whose Peirce
      block scalars are proportional to (4, 1, -2) -- the rep27
      Z' charges; the Levi's derived algebra has dim 45; P+ is an
      ideal of the parabolic (J6b).
  J6  Levi brackets: all 1035 commutators expressed in the
      78-basis via the pivot minor and verified as matrices.
  J7  the derived 45 act on J_0 antisymmetrically wrt the det_2
      bilinear form B_Q (so(10) on the 10, concretely).
  J8  p-brackets: [P+, P+] = 0 and [P-, P-] = 0 as matrices (all
      pairs); [P+, P-] lies in the Levi for all 256 pairs
      (pivot-minor membership, fully verified); [Zhat, P+-] has
      the single eigenvalue +-3s on each p-sector; [Levi, P+-]
      stays in the P+- sectors (J8d).
  J9  Casimir cross-check: with Gram = -tr_27(T_a T_b) on the
      derived 45 and dual bases, the Casimir operator is scalar
      on each Peirce piece; the scalars agree with the same
      prescription applied to the rep27 so(10) action (Sig on S,
      E on V, 0 on the singlet) -- the two h-modules coincide at
      the level of this complete invariant.
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

# ================================================================ parabolic
# The line stabilizer of E: T in span(OPS) with (T e_0)_k = 0 for
# all k != 0 (T e_0 = column 0 of T).  This is the maximal parabolic
# of dim 62 = 46 (Levi) + 16 (nilradical); the Levi is the charge-0
# sector, the nilradical the charge-+3 sector.
LEVI = [M for (M, q) in OPS if q == 0]
NL = len(LEVI)
assert NL == 46
Pp_all = [M for (M, q) in OPS if q == 3]

seqs = []
for k in range(1, 27):
    e = {}
    for a in range(NB):
        v = OPS[a][0].get((k, 0))
        if v is not None:
            e[a] = v
    if e:
        seqs.append(e)
stab_basis = nullspace(seqs, NB)     # coefficient vectors
ok("J5a line stabilizer of E (maximal parabolic) has dim 62",
   len(stab_basis) == 62, f"dim {len(stab_basis)}")

# sector-0 ops all satisfy the equations (column 0 has support in
# {0} by charge), and every P+ op has empty column 0 (charge +3 from
# a charge-4 column lands in charge 7, which does not exist).
sec0_in = all(all(r == 0 for (r, c) in M if c == 0) for M in LEVI)
pp_kill = all(not any(c == 0 for (r, c) in M) for M in Pp_all)
ok("J5b parabolic = Levi (+) P+ as vector spaces "
   "(sector-0 ops stabilize the line; P+ kills e_0; "
   "46 + 16 = 62, independent by charge)", sec0_in and pp_kill)

def combo(cf):
    """cf: dict {OPS index: Q(i) pair} -> the combined operator."""
    out = {}
    for a, c in cf.items():
        if phas(c):
            out = mat_add(out, mat_scale(OPS[a][0], c))
    return out

STAB = [combo(cf) for cf in stab_basis]
NS = len(STAB)

# centre of the Levi: coefficients d with [sum d_a LEVI_a, LEVI_b]
# = 0 for all b.
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
ok("J5c centre of the Levi is 1-dimensional",
   len(zcen) == 1, f"dim {len(zcen)}")

# ZHAT as a matrix: zcen vectors are over LEVI indices
ZHAT = {}
for ai, val in zcen[0].items():
    ZHAT = mat_add(ZHAT, mat_scale(LEVI[ai], val))
# block scalars: J1 = {0}, Jhalf = 3..18, J0 = {1,2} + 19..26
def block_scalar(M, idxs):
    vals = set()
    for r in idxs:
        for c in idxs:
            v = M.get((r, c), Z2)
            if r == c:
                vals.add(v)
            elif phas(v):
                return None
    return vals.pop() if len(vals) == 1 else None

lam1 = block_scalar(ZHAT, [0])
lamh = block_scalar(ZHAT, list(range(3, 19)))
lam0 = block_scalar(ZHAT, [1, 2] + list(range(19, 27)))
cross_ok = all(not phas(v) for (r, c), v in ZHAT.items()
               if (r in [0] and c not in [0]) or
                  (r in range(3, 19) and c not in range(3, 19)) or
                  (r in [1, 2] + list(range(19, 27)) and
                   c not in [1, 2] + list(range(19, 27))))
ratio_ok = False
if lam1 is not None and lamh is not None and lam0 is not None:
    s = pscale(lam1, Fraction(1, 4))
    ratio_ok = (phas(s) and lamh == s and
                lam0 == pscale(s, Fraction(-2)))
ok("J5d centre generator Zhat is block-scalar with Peirce "
   "eigenvalues proportional to (4, 1, -2)",
   cross_ok and ratio_ok,
   f"lam=({lam1}, {lamh}, {lam0})")
S_SCALE = pscale(lam1, Fraction(1, 4)) if lam1 is not None else Z2
print(f"  Zhat block scalars: (4, 1, -2) * s with s = {S_SCALE}",
      flush=True)

# derived algebra of the Levi: rank of its commutator matrices
# flattened as vectors in End(J).
sflat = [{r * 27 + c: v for (r, c), v in SCOMM[(a, b)].items()}
         for a in range(NL) for b in range(a + 1, NL)]
rk_s, _ = rank_of_vectors(sflat, 729)
ok("J5e derived algebra of the Levi has dim 45",
   rk_s == 45, f"rank {rk_s}")

# J6: Levi brackets in the 78-span, fully verified by express.
j6 = all(express(SCOMM[(a, b)]) is not None
         for a in range(NL) for b in range(a + 1, NL))
ok("J6 all 1035 Levi commutators verified in the 78-span", j6)

# J6b: the nilradical P+ is an ideal of the parabolic:
# [parabolic, P+] stays in the P+ sector for all 62 x 16 pairs.
j6b = True
for T in STAB:
    for P in Pp_all:
        C = commutator(T, P)
        cf = express(C) if C else [Z2] * NB
        if cf is None or any(phas(cf[t]) for t in range(NB)
                             if OPS[t][1] != 3):
            j6b = False
ok("J6b nilradical ideal: [parabolic, P+] in P+ sector "
   "(all 62 x 16 pairs, minor membership)", j6b)

# ================================================================ J7: so(10) on J0
# B_Q on J0 coords {1, 2} + 19..26 (C slot): for
# X = (d1, d2, C), Y = (d1', d2', C'):
# B_Q(X, Y) = (d1 d2' + d2 d1')/2 - bN(C, C'),
# bN(u, v) = ((u . conj v)[0] + (v . conj u)[0]) / 2
J0ID = [1, 2] + list(range(19, 27))
def bN(u, v):
    a = omul(u, oconj(v)).get(0, Z2)
    b = omul(v, oconj(u)).get(0, Z2)
    return pscale(padd(a, b), Fraction(1, 2))

def BQ(X, Y):
    v = pscale(padd(pmul(X.get(1, Z2), Y.get(2, Z2)),
                    pmul(X.get(2, Z2), Y.get(1, Z2))),
               Fraction(1, 2))
    return psub(v, bN(slot_oct(X, 19), slot_oct(Y, 19)))

def apply_op(T, X):
    out = {}
    for (r, c), v in T.items():
        if c in X:
            out[r] = padd(out.get(r, Z2), pmul(v, X[c]))
    return {k: v for k, v in out.items() if phas(v)}

# derived-45 basis: pivot commutators of the Levi
piv_rows = []
pivots45 = {}
for a in range(NL):
    for b in range(a + 1, NL):
        vec = {r * 27 + c: v for (r, c), v in SCOMM[(a, b)].items()}
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
DER45 = piv_rows
print(f"  derived-45 basis: {len(DER45)} commutators", flush=True)

j7 = True
for T in DER45:
    for ix in J0ID:
        for iy in J0ID:
            X = {ix: ONE2}
            Y = {iy: ONE2}
            TX = apply_op(T, X)
            TY = apply_op(T, Y)
            if padd(BQ(TX, Y), BQ(X, TY)) != Z2:
                j7 = False
ok("J7 derived 45 preserve B_Q on J0 (antisymmetric, 45 x 55 "
   "basis pairs)", j7)

# ================================================================ J8: p-brackets
Pp = [M for (M, q) in OPS if q == 3]
Pm = [M for (M, q) in OPS if q == -3]
j8a = all(commutator(Pp[a], Pp[b]) == {}
          for a in range(16) for b in range(a + 1, 16)) and \
      all(commutator(Pm[a], Pm[b]) == {}
          for a in range(16) for b in range(a + 1, 16))
ok("J8a [P+, P+] = 0 and [P-, P-] = 0 as matrices (all pairs)", j8a)

# coefficient indices belonging to sector 0:
sec0_idx = {a for a in range(NB) if OPS[a][1] == 0}
j8b = True
for a in range(16):
    for b in range(16):
        C = commutator(Pp[a], Pm[b])
        cf = express(C) if C else [Z2] * NB
        if cf is None or any(phas(cf[t]) for t in range(NB)
                             if t not in sec0_idx):
            j8b = False
ok("J8b [P+, P-] lies in the Levi (charge-0 sector) for all "
   "256 pairs (minor membership, fully verified)", j8b)

def sector_eig(T, Plist, mu):
    return all(commutator(T, P) == mat_scale(P, mu) for P in Plist)

muP = pscale(S_SCALE, Fraction(3))
muM = pscale(S_SCALE, Fraction(-3))
ok("J8c [Zhat, P] = (+3s) P on P+ and (-3s) P on P- "
   "(grading eigenvalues)", sector_eig(ZHAT, Pp, muP) and
   sector_eig(ZHAT, Pm, muM))

# [Levi, p] membership in the p sectors (the parabolic as a whole
# does NOT preserve P-: [P+, P-] lands in the Levi, J8b)
j8d = True
for T in LEVI:
    for Plist, qsec in ((Pp, 3), (Pm, -3)):
        for P in Plist:
            C = commutator(T, P)
            cf = express(C) if C else [Z2] * NB
            if cf is None or any(phas(cf[t]) for t in range(NB)
                                 if OPS[t][1] != qsec):
                j8d = False
ok("J8d [Levi, P+-] stays in the P+- sector (all 46 x 32 "
   "pairs, minor membership)", j8d)

# ================================================================ J9 Casimir
def trace_prod(A, B):
    t = Z2
    for (r, m), v in A.items():
        w = B.get((m, r))
        if w is not None:
            t = padd(t, pmul(v, w))
    return t

def casimir(ops45):
    n = len(ops45)
    G = sp.Matrix(n, n, lambda a, b: to_sympy(
        pneg(trace_prod(ops45[a], ops45[b]))))
    Ginv = G.inv()
    duals = []
    for a in range(n):
        D = {}
        for b in range(n):
            c = from_sympy(Ginv[a, b])
            if phas(c):
                D = mat_add(D, mat_scale(ops45[b], c))
        duals.append(D)
    C = {}
    for a in range(n):
        C = mat_add(C, matmul(ops45[a], duals[a]))
    return C

def piece_scalars(M, pieces):
    owner = {}
    for pi, idxs in enumerate(pieces):
        for i in idxs:
            owner[i] = pi
    vals = [None] * len(pieces)
    for (r, c), v in M.items():
        if owner.get(r) != owner.get(c):
            return None
        pi = owner[r]
        if r == c:
            if vals[pi] is None:
                vals[pi] = v
            elif vals[pi] != v:
                return None
        else:
            return None
    return [v if v is not None else Z2 for v in vals]

print("Casimir on the Jordan side...", flush=True)
CJ = casimir(DER45)
scJ = piece_scalars(CJ, [[0], list(range(3, 19)),
                         [1, 2] + list(range(19, 27))])
print(f"  J-side Casimir scalars (J1, J1/2, J0): {scJ}", flush=True)

# rep27 side: rebuild the so(10) action on the 27 as in e6_rep27.py
print("rebuilding rep27 so(10) action...", flush=True)
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
prod = eye(32)
for a in range(10):
    prod = prod * gam[a]
Gamma = sI * prod
minus_idx = [k for k in range(32) if Gamma[k, k] == -1]
pairs10 = [(a, b) for a in range(10) for b in range(a + 1, 10)]
SigL = [(gam[p[0]] * gam[p[1]] - gam[p[1]] * gam[p[0]]) *
        Rational(1, 4) for p in pairs10]
rho27 = [S.extract(minus_idx, minus_idx) for S in SigL]

def sym_to_pair(z):
    z = sp.expand(z)
    re = sp.re(z)
    im = sp.im(z)
    return (Fraction(re.p, re.q), Fraction(im.p, im.q))

REP45 = []
for j, p in enumerate(pairs10):
    M = {}
    R = rho27[j]
    for r in range(16):
        for c in range(16):
            if R[r, c] != 0:
                M[(r, c)] = sym_to_pair(R[r, c])
    a, b = p
    M[(16 + a, 16 + b)] = ONE2
    M[(16 + b, 16 + a)] = pneg(ONE2)
    REP45.append(M)

print("Casimir on the rep27 side...", flush=True)
CR = casimir(REP45)
scR = piece_scalars(CR, [list(range(16)), list(range(16, 26)), [26]])
print(f"  rep-side Casimir scalars (S, V, 1): {scR}", flush=True)

ok("J9 Casimir scalars agree across the two models "
   "(J1<->1, J1/2<->S, J0<->V)",
   scJ is not None and scR is not None and
   scJ[0] == scR[2] and scJ[1] == scR[0] and scJ[2] == scR[1],
   f"J={scJ} rep={scR}")

print(f"\n{len(FAIL)} failures / {NCHK[0]} checks")
if FAIL:
    raise SystemExit(1)
print("ALL PASS")
