# Tie-gluing toy check -- Sage leg (toe), exact over Q(i).
# Same computation as clifford_gluing_toy_check.py (local leg).
# Targeted imports only (from sage.all import * is broken on toe);
# ends with a top-level print (stdin EOF quirk).

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField
from sage.matrix.constructor import matrix, identity_matrix, diagonal_matrix
from sage.modules.free_module_element import vector

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
K = NumberField(x**2 + 1, 'i')
i = K.gen()
one = K(1)

G = diagonal_matrix(K, [1, 1, i, i, -1, -1, -i, -i])
Gbar = diagonal_matrix(K, [1, 1, -i, -i, -1, -1, i, i])
M_FW = diagonal_matrix(K, [1, 1, -1, -1, 1, 1, -1, -1])
trend = diagonal_matrix(K, [1, 1, -1, -1, 1, 1, -1, -1])

rows = []
for k in range(8):
    r = [K(0)] * 8
    r[(k - 1) % 8] = one
    rows.append(r)
T = matrix(K, rows)
Lam = T**4
I8 = identity_matrix(K, 8)
Z8 = matrix(K, 8, 8)

results = []
def record(name, ok):
    results.append((name, bool(ok)))
    print(("  [OK] " if ok else "  [FAIL] ") + name)

print("== established facts of the toy (Sage, exact over Q(i)) ==")
record("G^2 == M_FW", G * G == M_FW)
record("Gbar^2 == M_FW", Gbar * Gbar == M_FW)
record("{Lambda, G} == 0  (Lambda = T^4)", Lam * G + G * Lam == Z8)
record("T^2 G == -i * G T^2  (pushforward Weyl phase)", T * T * G == (-i) * (G * T * T))
record("{T^2, M_FW} == 0", T * T * M_FW + M_FW * T * T == Z8)
record("trend-bit operator == M_FW", trend == M_FW)
record("trend-bit operator == G^2", trend == G * G)

print("== tie crossings exist and are intra-half ==")
def unit(k):
    v = [K(0)] * 8
    v[k] = one
    return vector(K, v)
record("T e_0 == e_1  (QI tie transport)", T * unit(0) == unit(1))
record("T e_2 == e_3  (QII tie transport)", T * unit(2) == unit(3))
record("cells 0,1 in V+ (M_FW = +1 on both)", M_FW[0, 0] == 1 and M_FW[1, 1] == 1)
record("cells 2,3 in V- (M_FW = -1 on both)", M_FW[2, 2] == -1 and M_FW[3, 3] == -1)
record("grading split V+ = {0,1,4,5}, V- = {2,3,6,7}",
       [k for k in range(8) if M_FW[k, k] == 1] == [0, 1, 4, 5] and
       [k for k in range(8) if M_FW[k, k] == -1] == [2, 3, 6, 7])

print("== declared diagonal operators across the ties ==")
for name, D in [("G", G), ("Gbar", Gbar), ("M_FW", M_FW), ("trend bit", trend)]:
    record(name + " constant across QI tie (cells 0,1)", D[0, 0] == D[1, 1])
    record(name + " constant across QII tie (cells 2,3)", D[2, 2] == D[3, 3])

print("== the D2 flipper FW(2x) vs the declared operator set ==")
FW2 = diagonal_matrix(K, [(-1) ** k for k in range(8)])
record("FW(2x) diagonal != G, Gbar, M_FW, trend (not a declared toy operator)",
       all(FW2 != D for D in (G, Gbar, M_FW, trend)))

print("== algebra dimensions (Sage, exact over Q(i)) ==")
def rank_of(vecs):
    if not vecs:
        return 0
    return matrix(K, vecs).rank()

def algebra_dim(gens):
    basis = [vector(K, I8.list())]
    r = 1
    frontier = [I8]
    while frontier:
        nxt = []
        for A in frontier:
            for g in gens:
                P = A * g
                v = vector(K, P.list())
                if rank_of(basis + [v]) > r:
                    basis.append(v)
                    r += 1
                    nxt.append(P)
        frontier = nxt
    return r

d1 = algebra_dim([Lam, G])
d2 = algebra_dim([T, G])
record("dim <Lambda, G> == 8   (computed %d)" % d1, d1 == 8)
record("dim <T, G> == 64 = full M(8,C)   (computed %d)" % d2, d2 == 64)

fails = [n for n, ok in results if not ok]
print("SUMMARY: %d/%d exact checks passed" % (len(results) - len(fails), len(results)))
if fails:
    print("FAILURES: " + str(fails))
print("TOY_CHECK_SAGE_DONE")
