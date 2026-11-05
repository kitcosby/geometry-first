# Tie-gluing toy check (pre-registered brief of 2026-10-16) -- local exact leg.
#
# Rebuilds the minimal-declaration Clifford toy EXACTLY as recorded in the
# daily memory record of 2026-10-16 (memory/2026-10-16.md, toy model entry):
#   V = C^8, functions on the octant cells of the 2*pi lap;
#   declared branch G = diag(1,1,i,i,-1,-1,-i,-i);
#   transports T (one cell, pushforward), T^2 (pi/2), Lambda = T^4 (pi);
#   declared action = operators act as themselves.
# Re-verifies the toy's established facts in exact arithmetic over Q(i)
# (sympy for the matrix identities; exact Gaussian-rational elimination,
# implemented on Fraction pairs, for the algebra dimensions), then evaluates
# the tie-local data the pre-registered check turns on:
#   * every declared diagonal operator across the QI tie (cells 0->1, in V+)
#     and the QII tie (cells 2->3, in V-);
#   * whether the D2 dominance flipper FW(2x) is among the declared operators.
# No structure is added to the toy. No rider notion is invented.
#
# Run: python3 clifford_gluing_toy_check.py

from fractions import Fraction as Fr
import sympy as sp
from sympy import I, eye, zeros

# ---------------------------------------------------------------- sympy leg
def sdiag(vals):
    M = zeros(8)
    for k, v in enumerate(vals):
        M[k, k] = v
    return M

G     = sdiag([1, 1, I, I, -1, -1, -I, -I])
Gbar  = sdiag([1, 1, -I, -I, -1, -1, I, I])
M_FW  = sdiag([1, 1, -1, -1, 1, 1, -1, -1])
# trend bit: diag(sgn sawup') per cell; corpus (exposition s3, PROVED):
# sgn(sawup') = FW per cell, so the trend diagonal is FW per cell:
trend = sdiag([1, 1, -1, -1, 1, 1, -1, -1])

T = zeros(8)
for k in range(8):
    T[k, (k - 1) % 8] = 1          # pushforward: (Tf)_k = f_{k-1}
Lam = T**4
I8, Z8 = eye(8), zeros(8)

results = []
def record(name, ok):
    results.append((name, bool(ok)))
    print(f"  [{'OK' if ok else 'FAIL'}] {name}")

print("== established facts of the toy (sympy, exact over Q(i)) ==")
record("G^2 == M_FW", G * G == M_FW)
record("Gbar^2 == M_FW", Gbar * Gbar == M_FW)
record("{Lambda, G} == 0  (Lambda = T^4)", Lam * G + G * Lam == Z8)
record("T^2 G == -i * G T^2  (pushforward Weyl phase)", T * T * G == -I * (G * T * T))
record("{T^2, M_FW} == 0", T * T * M_FW + M_FW * T * T == Z8)
record("trend-bit operator == M_FW", trend == M_FW)
record("trend-bit operator == G^2", trend == G * G)

print("== tie crossings exist and are intra-half ==")
e = [zeros(8, 1) for _ in range(8)]
for k in range(8):
    e[k][k] = 1
record("T e_0 == e_1  (QI tie transport)", T * e[0] == e[1])
record("T e_2 == e_3  (QII tie transport)", T * e[2] == e[3])
record("cells 0,1 in V+ (M_FW = +1 on both)", M_FW[0, 0] == 1 and M_FW[1, 1] == 1)
record("cells 2,3 in V- (M_FW = -1 on both)", M_FW[2, 2] == -1 and M_FW[3, 3] == -1)
record("grading split V+ = {0,1,4,5}, V- = {2,3,6,7}",
       [k for k in range(8) if M_FW[k, k] == 1] == [0, 1, 4, 5] and
       [k for k in range(8) if M_FW[k, k] == -1] == [2, 3, 6, 7])

print("== declared diagonal operators across the ties ==")
for name, D in [("G", G), ("Gbar", Gbar), ("M_FW", M_FW), ("trend bit", trend)]:
    record(f"{name} constant across QI tie (cells 0,1)", D[0, 0] == D[1, 1])
    record(f"{name} constant across QII tie (cells 2,3)", D[2, 2] == D[3, 3])
print("  diagonal table (cells 0..7):")
for name, D in [("G", G), ("Gbar", Gbar), ("M_FW", M_FW), ("trend", trend)]:
    print(f"    {name}: {[D[k, k] for k in range(8)]}")

print("== the D2 flipper FW(2x) vs the declared operator set ==")
FW2 = sdiag([(-1) ** k for k in range(8)])   # FW(2x) per cell: sgn(sin 4x) alternates
record("FW(2x) diagonal != G, Gbar, M_FW, trend (not a declared toy operator)",
       all(FW2 != D for D in (G, Gbar, M_FW, trend)))

# ------------------------------------------------- exact Q(i) dimension leg
# Gaussian rationals as (re, im) Fraction pairs; exact elimination.
ZERO = (Fr(0), Fr(0))
ONE = (Fr(1), Fr(0))
QI = (Fr(0), Fr(1))

def qadd(a, b): return (a[0] + b[0], a[1] + b[1])
def qsub(a, b): return (a[0] - b[0], a[1] - b[1])
def qmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def qinv(a):
    d = a[0] * a[0] + a[1] * a[1]
    return (a[0] / d, -a[1] / d)
def qdiv(a, b): return qmul(a, qinv(b))

def matmul(A, B):
    out = []
    for r in range(8):
        row = []
        for c in range(8):
            s = ZERO
            for k in range(8):
                s = qadd(s, qmul(A[r][k], B[k][c]))
            row.append(s)
        out.append(row)
    return out

def diag_q(vals):
    M = [[ZERO] * 8 for _ in range(8)]
    for k, v in enumerate(vals):
        M[k][k] = v
    return M

def ident_q():
    return diag_q([ONE] * 8)

def flat(M):
    return [M[r][c] for r in range(8) for c in range(8)]

class Span:
    """Exact span over Q(i) on 64-dim flattened matrices, pivot basis."""
    def __init__(self):
        self.pivots = {}   # pivot index -> normalized basis vector (list of 64 pairs)

    def _reduce(self, v):
        v = list(v)
        for p in sorted(self.pivots):
            if v[p] != ZERO:
                f = v[p]
                b = self.pivots[p]
                v = [qsub(v[j], qmul(f, b[j])) for j in range(64)]
        return v

    def add(self, M):
        """Return True if flat(M) enlarged the span."""
        v = self._reduce(flat(M))
        for p in range(64):
            if v[p] != ZERO:
                inv = qinv(v[p])
                self.pivots[p] = [qmul(x, inv) for x in v]
                return True
        return False

    @property
    def dim(self):
        return len(self.pivots)

def algebra_dim(gens):
    span = Span()
    I8q = ident_q()
    span.add(I8q)
    frontier = [I8q]
    while frontier:
        nxt = []
        for A in frontier:
            for g in gens:
                P = matmul(A, g)
                if span.add(P):
                    nxt.append(P)
        frontier = nxt
    return span.dim

neg = lambda a: (-a[0], -a[1])
Gq = diag_q([ONE, ONE, QI, QI, neg(ONE), neg(ONE), neg(QI), neg(QI)])
Tq = [[ZERO] * 8 for _ in range(8)]
for k in range(8):
    Tq[k][(k - 1) % 8] = ONE
Lamq = matmul(matmul(Tq, Tq), matmul(Tq, Tq))

print("== algebra dimensions (exact Gaussian-rational elimination) ==")
d1 = algebra_dim([Lamq, Gq])
d2 = algebra_dim([Tq, Gq])
record(f"dim <Lambda, G> == 8   (computed {d1})", d1 == 8)
record(f"dim <T, G> == 64 = full M(8,C)   (computed {d2})", d2 == 64)

print()
fails = [n for n, ok in results if not ok]
print(f"SUMMARY: {len(results) - len(fails)}/{len(results)} exact checks passed")
if fails:
    print("FAILURES:", fails)
print("TOY_CHECK_DONE")
