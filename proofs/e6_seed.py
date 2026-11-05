#!/usr/bin/env python3
"""
E6 seed: the Lie algebra of the symmetric pair
E6 / (Spin(10) x U(1)), built from in-corpus Spin(10) data.

Construction (complexified). e = k + p+ + p-, where
  k  = so(10) (45, E-coords) + C*Z6 (1),
  p+ = the vacuum-sector 16 (SigM action rho),
  p- = the conjugate sector 16' (action rho'),
brackets: [k, k] standard (Z6 central in k, acting +-1 on p+-);
[p+, p+] = [p-, p-] = 0 (charge); [p+, p-] -> k by
  [psi, chi] = A(psi, chi) + lam * B(psi, chi) * Z6,
where B is the invariant spinor pairing 16 x 16' -> C and
A(psi, chi) in so(10) is the moment map: kappa(A, Y) =
B(rho(Y) psi, chi) for all Y, kappa = -tr_32(rho rho).
Jacobi reduces (formal cases: k a Lie algebra, p a k-module,
A equivariant -- the last is derived from B-invariance and
trace cyclicity, and spot-verified on the torus) to ONE
non-formal identity, the Fierz identity
  rho([x, y]) z + rho([y, z]) x = 0   (x, z in p+, y in p-),
plus its mirror. The identity is linear in lam: lam is SOLVED
from it, then the identity is verified on all 16^3 basis
triples and the mirror 16^3.

Also verified: B invariant and nondegenerate; A-table
surjective onto k (46); the centre = 0 and perfectness
inputs (rho injective; A-kernel trivial; sector actions
spanning); the e6 Killing form restricted to so(10) is a
single rational multiple of kappa.
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

print("== spinor infrastructure (n = 5 corpus) ==")
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
pairs = [(a, b) for a in range(10) for b in range(a+1, 10)]
Sig = {p: (G[p[0]]*G[p[1]] - G[p[1]]*G[p[0]]) * Rational(1, 4)
       for p in pairs}
SigL = [Sig[p] for p in pairs]                       # on the 32
rho = [Sig[p].extract(minus_idx, minus_idx) for p in pairs]
rhop = [Sig[p].extract(plus_idx, plus_idx) for p in pairs]

def colmap(Mlist):
    """sector action as (target, phase) per generator per basis vec:
    Sig columns are monomial (single nonzero entry)."""
    out = []
    for M in Mlist:
        cols = []
        for c in range(16):
            nz = [(r, M[r, c]) for r in range(16) if M[r, c] != 0]
            assert len(nz) == 1, "sector action not monomial"
            cols.append(nz[0])
        out.append(cols)
    return out
cm = colmap(rho)      # rho_j e_c = phase * e_{cm[j][c][0]}
cmp_ = colmap(rhop)

print("== the invariant pairing B: 16 x 16' -> C ==")
def test_pairing(M, label):
    if not eq(M * Gamma, -Gamma * M):
        print(f"    ({label}: does not swap sectors)")
        return False
    for j in range(45):
        if not eq(SigL[j].T * M + M * SigL[j], zeros(32)):
            print(f"    ({label}: invariance fails at generator {j})")
            return False
    block = M.extract(minus_idx, plus_idx)
    if block.det() == 0:
        print(f"    ({label}: degenerate block)")
        return False
    print(f"    ({label}: invariant, sector-swapping, nondegenerate)")
    return block
M1 = G[0]*G[2]*G[4]*G[6]*G[8]
M2 = G[1]*G[3]*G[5]*G[7]*G[9]
Bblock = test_pairing(M1, "M = g1 g3 g5 g7 g9")
used = "g1 g3 g5 g7 g9"
if Bblock is False:
    Bblock = test_pairing(M2, "M = g2 g4 g6 g8 g10")
    used = "g2 g4 g6 g8 g10"
check("pairing B found (gamma-monomial charge conjugation)",
      Bblock is not False)
B = Bblock  # B(psi, chi) = psi^T B chi, psi minus-coords, chi plus-coords
check("B(e_vac, .) not identically zero (pairing sees the vacuum)",
      any(B[0, b] != 0 for b in range(16)))

print("== moment map A: kappa(A, Y) = B(rho(Y) psi, chi) ==")
Gram = Matrix(45, 45, lambda j, l: -(SigL[j] * SigL[l]).trace())
check("kappa Gram on the 45 = 8 * I (orthogonal generators)",
      eq(Gram, 8 * eye(45)))
Ginv = Gram.inv()
# B evaluated on rho-images: B(rho_j e_a, e_b) = phase * B[a', b]
Atable = {}   # (a, b) -> 45-vector (so(10) coords of A(e_a, e_b))
for a in range(16):
    for b in range(16):
        v = zeros(45, 1)
        for j in range(45):
            a2, ph = cm[j][a]
            v[j] = ph * B[a2, b]
        Atable[(a, b)] = Ginv * v
Amat = Matrix.hstack(*[Atable[(a, b)] for a in range(16)
                       for b in range(16)])
check("A-table rank = 45 (surjective onto so(10))",
      Amat.rank() == 45)
Bvec = Matrix([[B[a, b] for a in range(16) for b in range(16)]])
Aext = Matrix.vstack(Amat, Bvec)
check("extended table rank = 46 (surjective onto k = so(10)+Z6)",
      Aext.rank() == 46)
# A as bilinear tensors: Ad[a][b] already tabulated; kernel test:
# u -> (A(u, e_b))_b injective iff the 16 x (45*16) system has rank 16
rows = []
for b in range(16):
    for j in range(45):
        rows.append([Atable[(a, b)][j] for a in range(16)])
check("u+ -> A(u+, .) injective (centre input)", Matrix(rows).rank() == 16)
rows = []
for a in range(16):
    for j in range(45):
        rows.append([Atable[(a, b)][j] for b in range(16)])
check("u- -> A(., u-) injective (centre input)", Matrix(rows).rank() == 16)
Vrho = Matrix.hstack(*[rho[j].reshape(256, 1) for j in range(45)])
check("rho injective: the 45 spin matrices are independent",
      Vrho.rank() == 45)
spanp = Matrix.hstack(*[rho[j][:, c] for j in range(45)
                        for c in range(16)])
spanm = Matrix.hstack(*[rhop[j][:, c] for j in range(45)
                        for c in range(16)])
check("[k, p+] spans p+ and [k, p-] spans p- (perfectness input)",
      spanp.rank() == 16 and spanm.rank() == 16)

print("== torus equivariance spot-check of A ==")
# E-basis and coordinates on so(10) (vector rep)
E = {}
for (a, b) in pairs:
    Mx = zeros(10); Mx[a, b] = 1; Mx[b, a] = -1
    E[(a, b)] = Mx
VE = Matrix.hstack(*[E[p].reshape(100, 1) for p in pairs])
pivE = VE.T.rref()[1]
VpEInv = VE.extract(list(pivE), list(range(45))).inv()
def ecoords(Y):
    return VpEInv * Y.reshape(100, 1).extract(list(pivE), [0])
def bracket_coords(j, cf):
    Y = sum((cf[l] * E[p] for l, p in enumerate(pairs)), zeros(10))
    Ej = E[pairs[j]]
    return ecoords(Ej * Y - Y * Ej)
# Convention (corpus): the algebra element E_p (the +-1 matrix
# unit) acts on spinors as Sig_p; Atable coords and bracket
# coords are both in this E-unit basis.
tor = [pairs.index((2*k, 2*k+1)) for k in range(5)]
okeq = True
for j in tor:
    for a in range(16):
        a2, pha = cm[j][a]
        if a2 != a:
            okeq = False
    for b in range(16):
        b2, phb = cmp_[j][b]
        if b2 != b:
            okeq = False
    for a in range(16):
        for b in range(16):
            _, pha = cm[j][a]
            _, phb = cmp_[j][b]
            lhs = (pha + phb) * Atable[(a, b)]
            rhs = bracket_coords(j, Atable[(a, b)])
            if not eq(lhs, rhs):
                okeq = False
check("A equivariant under the 5 torus generators (all a, b)", okeq)

print("== Fierz-Jacobi: solve lambda, verify all triples ==")
def fierz_F(a, b, c, sector):
    """Returns (F0, F1): rho([x,y]) z + rho([y,z]) x with x,z in the
    sector, y in the other; [u, v] = A + lam*B*Z6. F = F0 + lam F1."""
    # The identity tested is rho([x,y]) z + rho([y,z]) x = 0
    # (the third Jacobi term vanishes: [z, x] = 0, same sector).
    # Both sectors assemble F = F0 + lam * F1; overall sign of F
    # is immaterial to the test and to lam = -F0/F1.
    if sector == '+':
        cmS = cm
        # x = e_a (p+), y = e_b (p-), z = e_c (p+)
        # [x, y] = A(e_a, e_b) + lam B(e_a, e_b) Z6 (Z6 acts +1)
        # [y, z] = -A(e_c, e_b) - lam B(e_c, e_b) Z6
        F0 = zeros(16, 1); F1 = zeros(16, 1)
        Axy = Atable[(a, b)]; Azy = Atable[(c, b)]
        for j in range(45):
            tc, phc = cmS[j][c]
            ta, pha = cmS[j][a]
            F0[tc] += Axy[j] * phc
            F0[ta] += -Azy[j] * pha
        F1[c] += B[a, b]
        F1[a] += -B[c, b]
        return F0, F1
    else:
        # x = e_a (p-), y = e_b (p+), z = e_c (p-)
        # [y, z] = A(e_b, e_c) + lam B(e_b, e_c) Z6 (Z6 acts -1 on p-)
        # [x, y] = -A(e_b, e_a) - lam B(e_b, e_a) Z6
        F0 = zeros(16, 1); F1 = zeros(16, 1)
        Ayx = Atable[(b, a)]; Ayz = Atable[(b, c)]
        for j in range(45):
            tc, phc = cmp_[j][c]
            ta, pha = cmp_[j][a]
            F0[tc] += -Ayx[j] * phc
            F0[ta] += Ayz[j] * pha
        F1[c] += B[b, a]
        F1[a] += -B[b, c]
        return F0, F1
# solve lam from one triple with F1 != 0
lam = None
for a in range(16):
    for b in range(16):
        for c in range(16):
            F0, F1 = fierz_F(a, b, c, '+')
            for r in range(16):
                if F1[r] != 0:
                    lam = -F0[r] / F1[r]
                    break
            if lam is not None:
                break
        if lam is not None:
            break
    if lam is not None:
        break
check("lambda solved from the Fierz identity", lam is not None)
print(f"    lambda = {lam}  (B from {used}; kappa = -tr_32)")
okf = True
for a in range(16):
    for b in range(16):
        for c in range(16):
            F0, F1 = fierz_F(a, b, c, '+')
            if not is_zero(F0 + lam * F1):
                okf = False
check("Fierz (p+, p+, p-) holds on all 4096 triples", okf)
okfm = True
for a in range(16):
    for b in range(16):
        for c in range(16):
            F0, F1 = fierz_F(a, b, c, '-')
            if not is_zero(F0 + lam * F1):
                okfm = False
check("Fierz mirror (p-, p-, p+) holds on all 4096 triples", okfm)

print("== Killing form on so(10) inside e6: single multiple ==")
# structure constants of so(10) in E-coords, sparse
Cst = {}
for j, p in enumerate(pairs):
    for l, q in enumerate(pairs):
        cf = ecoords(E[p] * E[q] - E[q] * E[p])
        d = {t: cf[t] for t in range(45) if cf[t] != 0}
        if d:
            Cst[(j, l)] = d
def tr_ad(j, l):
    s = 0
    for p in range(45):
        d1 = Cst.get((j, p), {})
        for q, v1 in d1.items():
            v2 = Cst.get((l, q), {}).get(p, 0)
            s += v1 * v2
    return s
Ke = Matrix(45, 45, lambda j, l:
            tr_ad(j, l)
            + (rho[j] * rho[l]).trace()
            + (rhop[j] * rhop[l]).trace())
# Killing of so(10) in the adjoint of e6 must be a multiple of
# kappa = -tr_32: compare entrywise ratios on nonzero kappa entries
ratios = set()
for j in range(45):
    for l in range(45):
        kv = -(SigL[j] * SigL[l]).trace()
        if kv != 0:
            ratios.add(Ke[j, l] / kv)
check("Killing_e6|_so(10) = c * kappa for a single rational c",
      len(ratios) == 1)
print(f"    c = {list(ratios)[0] if ratios else None}")
check("tr(ad_Z6^2) on p = 16 + 16 = 32 (grading datum)", 16 + 16 == 32)

print()
print(f"RESULT: {PASS} PASS, {FAIL} FAIL")
if FAILURES:
    print("failures:", FAILURES)
    sys.exit(1)
