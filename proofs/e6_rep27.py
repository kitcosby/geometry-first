#!/usr/bin/env python3
"""
The 27 of E6, built on the seed (e6_seed.py).

V27 = S + V + 1: S = vacuum-sector 16, V = C^10, 1 = singlet.
Under k = so(10) + Z6 the charges are FORCED by the bracket
[Z6, T_x] = +T_x (x in p+): with charges (al, be, ga) on
(S, V, 1) proportional to (1, -2, 4), the nonzero components
of T_x must raise charge by +1, giving -3c = 1, c = -1/3:
(al, be, ga) = (-1/3, 2/3, -4/3). With Z' = -3 Z6 these are
the textbook charges (1, -2, 4).

The p-actions, all so(10)-equivariant maps (Hom spaces
1-dimensional or 0 by the tensor decompositions; the zero
components are forced both by Hom-vanishing and by charge):
  T_x (x in p+):  1 -> S : s |-> a1 * s * x
                  S -> V : psi |-> a2 * W(x, psi),
                    W(x, psi)_a = B(x, gamma_a psi)
  T_y (y in p-):  V -> S : v |-> a3 * gamma(v) y
                  S -> 1 : psi |-> a4 * B(psi, y)
Only the products a1 a4 and a2 a3 are determined (piece-basis
rescaling gauge); they are SOLVED from the homomorphism
conditions, not chosen:
  [T_x, T_y]|_1 = lam * ga * B(x, y)   =>  a1 a4 = -1/8
      (lam = -3/32 is the seed's solved value, seed normalisations)
  [T_x, T_y]|_S : spinor Fierz-type identity  =>  a2 a3 = k2
  [T_x, T_y]|_V : vector identity, must give the SAME k2.
[T_x, T_x'] = 0 is W-symmetry; [T_y, T_y'] = 0 is the
gamma-B symmetry on 16'. Gauge fixed: a1 = a2 = 1,
a3 = k2, a4 = -1/8. The assembled action of all 78 generators
is then verified bracket by bracket on the 27 basis vectors.

Cubic invariant: the k-invariant symmetric cubics on V27 are
exactly span{ W(psi,psi).v , s (v.v) } (piece-content +
Z6-charge: every other content has nonzero charge or
vanishing Hom). p-invariance imposes one linear condition:
the ratio is solved and invariance verified for all 78
generators on all 27 basis vectors (via polarisation).
Irreducibility: one basis vector per piece is cyclic under
the 78 operators (complete reducibility cited).
"""
import sys
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

print("== infrastructure: spinors, B, A-table (as in the seed) ==")
s1 = Matrix([[0, 1], [1, 0]])
s2 = Matrix([[0, -I], [I, 0]])
s3 = Matrix([[1, 0], [0, -1]])
gam = [s1, s2]
while gam[0].rows < 32:
    d0 = gam[0].rows
    gam = [kronecker_product(g, s3) for g in gam] + \
          [kronecker_product(eye(d0), s1), kronecker_product(eye(d0), s2)]
G = gam
check("TRUE control: {g1,g2} = 0", eq(G[0]*G[1] + G[1]*G[0], zeros(32)))
check("FALSE control: [g1,g2] != 0", not is_zero(G[0]*G[1] - G[1]*G[0]))
prod = eye(32)
for a in range(10):
    prod = prod * G[a]
Gamma = I * prod
minus_idx = [k for k in range(32) if Gamma[k, k] == -1]
plus_idx = [k for k in range(32) if Gamma[k, k] == +1]
mpos = {v: i for i, v in enumerate(minus_idx)}
ppos = {v: i for i, v in enumerate(plus_idx)}
pairs = [(a, b) for a in range(10) for b in range(a+1, 10)]
SigL = [(G[p[0]]*G[p[1]] - G[p[1]]*G[p[0]]) * Rational(1, 4)
        for p in pairs]
rho = [S.extract(minus_idx, minus_idx) for S in SigL]
rhop = [S.extract(plus_idx, plus_idx) for S in SigL]
Mchg = G[0]*G[2]*G[4]*G[6]*G[8]
B = Mchg.extract(minus_idx, plus_idx)
check("B invariant (45 generators, as in the seed)",
      all(eq(S.T * Mchg + Mchg * S, zeros(32)) for S in SigL))
Gram = Matrix(45, 45, lambda j, l: -(SigL[j] * SigL[l]).trace())
Ginv = Gram.inv()
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
# gamma_a on the 32: monomial, sector-swapping
gcm = []
for a in range(10):
    cols = []
    for c in range(32):
        nz = [(r, G[a][r, c]) for r in range(32) if G[a][r, c] != 0]
        assert len(nz) == 1
        cols.append(nz[0])
    gcm.append(cols)
Atable = {}
for a in range(16):
    for b in range(16):
        v = zeros(45, 1)
        for j in range(45):
            a2, ph = cm[j][a]
            v[j] = ph * B[a2, b]
        Atable[(a, b)] = Ginv * v
lam = Rational(-3, 32)   # the seed's solved value

print("== the intertwiner W: (16 x 16) -> 10 ==")
# W(x, psi)_a = B(x, gamma_a psi); basis tensor Wbas[a][r][t]
Wbas = [[[0]*16 for _ in range(16)] for _ in range(10)]
for a in range(10):
    for t in range(16):
        tgt, ph = gcm[a][minus_idx[t]]       # gamma_a psi_t in plus
        u = ppos[tgt]
        for r in range(16):
            if B[r, u] != 0:
                Wbas[a][r][t] = ph * B[r, u]
check("W nonzero (some basis entry; W(e_vac, e_vac) = 0 by the "
      "subset-complement pairing, as it must be)",
      any(Wbas[a][r][t] != 0 for a in range(10)
          for r in range(16) for t in range(16)))
# symmetry as matrix identities on the 32:
#   W symmetric  <=>  block_--(M g_a) symmetric for all a
#   B(g_a ., .) symmetric on 16' <=> block_++(g_a^T M) symmetric
oksy1 = True
for a in range(10):
    blk = (Mchg * G[a]).extract(minus_idx, minus_idx)
    if not eq(blk, blk.T):
        oksy1 = False
check("W symmetric in its spinor arguments (all a)", oksy1)
oksy2 = True
for a in range(10):
    blk = (G[a].T * Mchg).extract(plus_idx, plus_idx)
    if not eq(blk, blk.T):
        oksy2 = False
check("B(gamma_a ., .) symmetric on 16' (all a)", oksy2)
# W equivariance: W(rho_j x, psi) + W(x, rho_j psi) = (E_j W)_vec
E = {}
for (a, b) in pairs:
    Mx = zeros(10); Mx[a, b] = 1; Mx[b, a] = -1
    E[(a, b)] = Mx
okW = True
for j, p in enumerate(pairs):
    Ej = E[p]
    for r in range(16):
        r2, phr = cm[j][r]
        for t in range(16):
            t2, pht = cm[j][t]
            for a in range(10):
                lhs = phr * Wbas[a][r2][t] + pht * Wbas[a][r][t2]
                rhs = sum(Ej[a, bb] * Wbas[bb][r][t] for bb in range(10))
                if lhs != rhs:
                    okW = False
check("W equivariant under so(10) (45 x 256 x 10)", okW)
# gamma intertwining: Sig_j gamma_b - gamma_b Sig_j acts as
# gamma_{E_j e_b} across sectors (checked on all 32 basis vectors)
okgi = True
for j, p in enumerate(pairs):
    Ej = E[p]
    for b in range(10):
        for c in range(32):
            # LHS on e_c: Sig_j (gamma_b e_c) - gamma_b (Sig_j e_c)
            t1, ph1 = gcm[b][c]
            # Sig_j on the 32 as sector colmaps:
            def sig32(j, idx):
                if idx in mpos:
                    t, ph = cm[j][mpos[idx]]
                    return minus_idx[t], ph
                t, ph = cmp_[j][ppos[idx]]
                return plus_idx[t], ph
            u1, q1 = sig32(j, t1)
            t2, ph2 = sig32(j, c)
            u2, q2 = gcm[b][t2]
            lhs = {}
            lhs[u1] = lhs.get(u1, 0) + ph1 * q1
            lhs[u2] = lhs.get(u2, 0) - ph2 * q2
            lhs = {k: v for k, v in lhs.items() if v != 0}
            rhs = {}
            for bb in range(10):
                if Ej[bb, b] != 0:
                    u3, q3 = gcm[bb][c]
                    rhs[u3] = rhs.get(u3, 0) + Ej[bb, b] * q3
            rhs = {k: v for k, v in rhs.items() if v != 0}
            if lhs != rhs:
                okgi = False
check("gamma intertwining [Sig_j, gamma_b] = gamma_{E_j e_b}", okgi)

print("== solve k2 = a2 a3: S-block identity, V-block confirm ==")
# E3: rho(A(x,y)) psi + (1/32) B(x,y) psi + (1/8) B(psi,y) x
#     + k2 * gamma(W(x,psi)) y = 0   for all basis (x, psi, y)
def gammaW_y(x, psi, y):
    """gamma( W(x, psi) ) applied to y (plus sector) -> minus vector."""
    out = zeros(16, 1)
    for a in range(10):
        w = Wbas[a][x][psi]
        if w == 0:
            continue
        tgt, ph = gcm[a][plus_idx[y]]
        out[mpos[tgt]] += w * ph
    return out
def E3_parts(x, psi, y):
    base = zeros(16, 1)
    Axy = Atable[(x, y)]
    for j in range(45):
        t, ph = cm[j][psi]
        base[t] += Axy[j] * ph
    base += Rational(1, 32) * B[x, y] * eye(16)[:, psi]
    ex = zeros(16, 1); ex[x] = Rational(1, 8) * B[psi, y]
    base += ex
    return base, gammaW_y(x, psi, y)
k2 = None
for x in range(16):
    for psi in range(16):
        for y in range(16):
            base, gv = E3_parts(x, psi, y)
            for r in range(16):
                if gv[r] != 0:
                    k2 = -base[r] / gv[r]
                    break
            if k2 is not None:
                break
        if k2 is not None:
            break
    if k2 is not None:
        break
check("k2 solved from the S-block identity", k2 is not None)
print(f"    k2 = a2 a3 = {k2}")
okE3 = True
for x in range(16):
    for psi in range(16):
        for y in range(16):
            base, gv = E3_parts(x, psi, y)
            if not is_zero(base + k2 * gv):
                okE3 = False
check("S-block identity holds on all 4096 triples", okE3)
# E4: k2 * W(x, gamma(v) y) = A_mat(x, y) v - (1/16) B(x, y) v
okE4 = True
for x in range(16):
    for y in range(16):
        Axy = Atable[(x, y)]
        Amat = sum((Axy[j] * E[p] for j, p in enumerate(pairs)),
                   zeros(10))
        for vb in range(10):
            # gamma(vb) y in minus sector:
            tgt, ph = gcm[vb][plus_idx[y]]
            mid = mpos[tgt]
            lhs = zeros(10, 1)
            for a in range(10):
                lhs[a] = k2 * ph * Wbas[a][x][mid]
            rhs = Amat[:, vb] - Rational(1, 16) * B[x, y] * \
                eye(10)[:, vb]
            if not eq(lhs, rhs):
                okE4 = False
check("V-block identity holds with the SAME k2 (2560 evals)", okE4)

print("== assemble the 27 action; verify all brackets ==")
# Sparse-dict operator engine on the 27 (idx 0..15 = S, 16..25 = V,
# 26 = singlet). The dense version of these checks was correct but
# slow; the sparse form is the same mathematics, exact throughout.
al, be, ga = Rational(-1, 3), Rational(2, 3), Rational(-4, 3)
a1, a2, a3, a4 = Rational(1), Rational(1), k2, Rational(-1, 8)
Ejcols = []
for j, p in enumerate(pairs):
    Ej = E[p]
    Ejcols.append([[ (bb, Ej[bb, b]) for bb in range(10)
                     if Ej[bb, b] != 0] for b in range(10)])
def dclean(d):
    return {k: v for k, v in d.items() if v != 0}
def dsub(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) - v
    return dclean(out)
def opY_s(j, d):
    out = {}
    for r, val in d.items():
        if r < 16:
            t, ph = cm[j][r]
            out[t] = out.get(t, 0) + val * ph
        elif r < 26:
            for bb, ev in Ejcols[j][r - 16]:
                out[16 + bb] = out.get(16 + bb, 0) + val * ev
    return dclean(out)
def opZ6_s(d):
    out = {}
    for r, val in d.items():
        sc = al if r < 16 else (be if r < 26 else ga)
        out[r] = val * sc
    return dclean(out)
def opTx_s(x, d):
    out = {}
    s = d.get(26, 0)
    if s != 0:
        out[x] = out.get(x, 0) + a1 * s
    for t, val in d.items():
        if t < 16:
            for a in range(10):
                w = Wbas[a][x][t]
                if w != 0:
                    out[16 + a] = out.get(16 + a, 0) + a2 * val * w
    return dclean(out)
def opTy_s(y, d):
    out = {}
    sout = 0
    for r, val in d.items():
        if 16 <= r < 26:
            tgt, ph = gcm[r - 16][plus_idx[y]]
            tp = mpos[tgt]
            out[tp] = out.get(tp, 0) + a3 * val * ph
        elif r < 16:
            sout += val * B[r, y]
    if sout != 0:
        out[26] = a4 * sout
    return dclean(out)
def densify(d):
    v = zeros(27, 1)
    for k, val in d.items():
        v[k] = val
    return v
BASIS27 = [eye(27)[:, i] for i in range(27)]
BASISD = [{i: Rational(1)} for i in range(27)]
def brk(opA, opB, d):
    return dsub(opA(opB(d)), opB(opA(d)))
ok = True
for x in range(16):
    for I in range(27):
        if brk(opZ6_s, lambda w: opTx_s(x, w), BASISD[I]) != \
           opTx_s(x, BASISD[I]):
            ok = False
        if brk(opZ6_s, lambda w: opTy_s(x, w), BASISD[I]) != \
           dclean({k: -v for k, v in opTy_s(x, BASISD[I]).items()}):
            ok = False
check("[Z6, T_x] = T_x, [Z6, T_y] = -T_y (all x, y, basis)", ok)
ok = True
for j in range(45):
    for x in range(16):
        tx, phx = cm[j][x]
        ty, phy = cmp_[j][x]
        for I in range(27):
            lhs = brk(lambda w: opY_s(j, w), lambda w: opTx_s(x, w),
                      BASISD[I])
            rhs = dclean({k: phx * v for k, v in
                          opTx_s(tx, BASISD[I]).items()})
            if lhs != rhs:
                ok = False
            lhs = brk(lambda w: opY_s(j, w), lambda w: opTy_s(x, w),
                      BASISD[I])
            rhs = dclean({k: phy * v for k, v in
                          opTy_s(ty, BASISD[I]).items()})
            if lhs != rhs:
                ok = False
check("[Y_j, T_x] = T_{rho_j x}, [Y_j, T_y] = T_{rho'_j y}", ok)
# [T_x, T_y] against the seed bracket as a block matrix per pair
ok = True
for x in range(16):
    for y in range(16):
        Axy = Atable[(x, y)]
        RA = sum((Axy[j] * rho[j] for j in range(45)), zeros(16)) \
            + lam * B[x, y] * al * eye(16)
        Amat = sum((Axy[j] * E[p] for j, p in enumerate(pairs)),
                   zeros(10)) + lam * B[x, y] * be * eye(10)
        R27 = Matrix.vstack(
            Matrix.hstack(RA, zeros(16, 11)),
            Matrix.hstack(zeros(10, 16), Amat, zeros(10, 1)),
            Matrix.hstack(zeros(1, 26),
                          Matrix([[lam * B[x, y] * ga]])))
        for I in range(27):
            lhs = brk(lambda w: opTx_s(x, w), lambda w: opTy_s(y, w),
                      BASISD[I])
            rhs = {k: R27[k, I] for k in range(27) if R27[k, I] != 0}
            if lhs != rhs:
                ok = False
check("[T_x, T_y] = action of [x, y]_seed (all 256 pairs, basis)", ok)
ok = True
for x in range(16):
    for xp in range(x+1, 16):
        for I in range(27):
            if brk(lambda w: opTx_s(x, w), lambda w: opTx_s(xp, w),
                   BASISD[I]) != {}:
                ok = False
            if brk(lambda w: opTy_s(x, w), lambda w: opTy_s(xp, w),
                   BASISD[I]) != {}:
                ok = False
check("[T_x, T_x'] = 0 and [T_y, T_y'] = 0 (all pairs, basis)", ok)

print("== cubic invariant ==")
# Tensor method (the per-triple Chat loop was too slow and timed
# out at 1700 s; this is the same mathematics, factored): build
# the symmetric trilinear tensors T1, T2 of the two ansatz forms
# ONCE from their values on basis triples (only the piece
# patterns (S,S,V) and (V,V,1) can be nonzero), then the
# infinitesimal action of each generator is a sparse tensor
# contraction delta_g T, and invariance is delta_g T = 0.
def split(vec):
    return (vec.extract(list(range(16)), [0]),
            vec.extract(list(range(16, 26)), [0]), vec[26])
def Chat(U1, U2, U3, d1, d2):
    p1, v1, s1 = split(U1)
    p2, v2, s2 = split(U2)
    p3, v3, s3 = split(U3)
    def Wof(pa, pb):
        out = zeros(10, 1)
        for a in range(10):
            tot = 0
            for r in range(16):
                if pa[r] == 0:
                    continue
                for t in range(16):
                    if pb[t] != 0 and Wbas[a][r][t] != 0:
                        tot += pa[r] * pb[t] * Wbas[a][r][t]
            out[a] = tot
        return out
    t1 = (Wof(p1, p2).dot(v3) + Wof(p1, p3).dot(v2)
          + Wof(p2, p3).dot(v1)) * Rational(1, 3)
    t2 = (s1 * v2.dot(v3) + s2 * v1.dot(v3)
          + s3 * v1.dot(v2)) * Rational(1, 3)
    return d1 * t1 + d2 * t2
def build_tensor(d1, d2):
    T = {}
    for I in range(16):
        for J in range(I, 16):
            for K in range(16, 26):
                val = Chat(BASIS27[I], BASIS27[J], BASIS27[K],
                           d1, d2)
                if val != 0:
                    T[(I, J, K)] = val
    for I in range(16, 26):
        for J in range(I, 26):
            val = Chat(BASIS27[I], BASIS27[J], BASIS27[26],
                       d1, d2)
            if val != 0:
                T[(I, J, 26)] = val
    return T
T1 = build_tensor(1, 0)
T2 = build_tensor(0, 1)
check("both ansatz tensors nonzero",
      len(T1) > 0 and len(T2) > 0)
print(f"    tensor sizes: |T1| = {len(T1)}, |T2| = {len(T2)}")
GENS = []
for j in range(45):
    GENS.append((f"Y{j}", lambda w, j=j: opY_s(j, w)))
GENS.append(("Z6", opZ6_s))
for x in range(16):
    GENS.append((f"Tx{x}", lambda w, x=x: opTx_s(x, w)))
    GENS.append((f"Ty{x}", lambda w, x=x: opTy_s(x, w)))
COLS = {}
for name, g in GENS:
    COLS[name] = [g(BASISD[I]) for I in range(27)]
# ROWS[name][M] = {I: G[M, I]}: the transpose contraction.
# delta_g T (I,J,K) = sum_M G[M,I] T[M,J,K] + (two more slots):
# the action contracts T's slot index with g's ROW index.
# (First implementation contracted with the column index --
# the transposed action; the direct Chat cross-check caught
# it: at (0,7,16) under Y0 the tensor delta gave 2i/3 where
# the direct evaluation gives 0.)
ROWS = {}
for name, g in GENS:
    rows = [dict() for _ in range(27)]
    for I in range(27):
        for M, gv in COLS[name][I].items():
            rows[M][I] = gv
    ROWS[name] = rows
def delta(T, row):
    out = {}
    for (I, J, K), val in T.items():
        for (slot, others) in ((I, (J, K)), (J, (I, K)),
                               (K, (I, J))):
            for L, gv in row[slot].items():
                key = tuple(sorted((L,) + others))
                out[key] = out.get(key, 0) + val * gv
    return {k: v for k, v in out.items() if v != 0}
# verify each form is separately k-invariant (46 gens): this is
# the equivariance content, checked at tensor level
okk = True
for name, g in GENS[:46]:
    if delta(T1, ROWS[name]) != {} or delta(T2, ROWS[name]) != {}:
        okk = False
check("each ansatz form separately invariant under k (46 gens)",
      okk)
# solve the ratio from a p-generator's delta
sol = None
for name, g in GENS[46:]:
    dT1 = delta(T1, ROWS[name])
    dT2 = delta(T2, ROWS[name])
    keys = set(dT1) | set(dT2)
    for key in sorted(keys):
        u1, u2 = dT1.get(key, 0), dT2.get(key, 0)
        if u1 != 0 or u2 != 0:
            sol = (name, key, u1, u2)
            break
    if sol is not None:
        break
check("a p-generator sees the cubic ansatz", sol is not None)
name0, key0, u1, u2 = sol
d1, d2 = u2, -u1
print(f"    first non-invariant direction: {name0} at {key0}: "
      f"(u1, u2) = ({u1}, {u2}); ratio (d1, d2) = ({d1}, {d2})")
TC = {}
for key in set(T1) | set(T2):
    v = d1 * T1.get(key, 0) + d2 * T2.get(key, 0)
    if v != 0:
        TC[key] = v
okc = True
for name, g in GENS:
    if delta(TC, ROWS[name]) != {}:
        okc = False
check("combined cubic invariant under all 78 generators "
      "(tensor deltas)", okc)
astar = next((a, r, t) for a in range(10) for r in range(16)
             for t in range(16) if Wbas[a][r][t] != 0)
Xt1 = BASIS27[astar[1]] + BASIS27[astar[2]] \
    + BASIS27[16 + astar[0]] + BASIS27[26]
Xt2 = BASIS27[16] + BASIS27[26]
Xt3 = BASIS27[astar[1]] + BASIS27[astar[2]] + BASIS27[16 + astar[0]]
check("the invariant cubic is not identically zero",
      Chat(Xt1, Xt1, Xt1, d1, d2) != 0
      or Chat(Xt2, Xt2, Xt2, d1, d2) != 0
      or Chat(Xt3, Xt3, Xt3, d1, d2) != 0)

print("== irreducibility: cyclicity from each piece ==")
GENS_S = []
for j in range(45):
    GENS_S.append(lambda w, j=j: opY_s(j, w))
GENS_S.append(opZ6_s)
for x in range(16):
    GENS_S.append(lambda w, x=x: opTx_s(x, w))
    GENS_S.append(lambda w, x=x: opTy_s(x, w))
def cyclic(seedI):
    span = Matrix.hstack(BASIS27[seedI])
    frontier = [BASISD[seedI]]
    while span.rank() < 27 and frontier:
        new = []
        for w in frontier:
            for g in GENS_S:
                gw = g(w)
                if not gw:
                    continue
                dv = densify(gw)
                trial = Matrix.hstack(span, dv)
                if trial.rank() > span.rank():
                    span = trial
                    new.append(gw)
        frontier = new
    return span.rank()
check("basis vector of S is cyclic (span = 27)", cyclic(0) == 27)
check("basis vector of V is cyclic (span = 27)", cyclic(16) == 27)
check("the singlet basis vector is cyclic (span = 27)",
      cyclic(26) == 27)

print()
print(f"charges (Z6): S {al}, V {be}, 1 {ga};  Z' = -3 Z6 gives "
      f"({-3*al}, {-3*be}, {-3*ga})")
print(f"gauge: (a1, a2, a3, a4) = (1, 1, {k2}, -1/8); "
      f"a1 a4 = -1/8, a2 a3 = {k2}")
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
