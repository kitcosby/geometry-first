# Monomial (Pauli-string) engine for the Clifford constructions.
# A monomial is (f, c, mu): f = bitmask over m qubits, c in Z/4,
# mu = tuple of m values in Z/4. It represents the operator
#   e_x |-> i^(c + sum_j mu_j * bit_j(x)) * e_{x ^ f}
# All Clifford monomials built from gamma products have mu_j in {0, 2}.
# Arithmetic is exact integer arithmetic mod 4.

def mmul(M1, M2):
    """Product M1 * M2 (M2 applied first)."""
    f1, c1, mu1 = M1
    f2, c2, mu2 = M2
    m = len(mu1)
    add = 0
    for j in range(m):
        if (f2 >> j) & 1:
            add += mu1[j]
    c = (c1 + c2 + add) % 4
    mu = tuple((mu2[j] + mu1[j] * (1 - 2 * ((f2 >> j) & 1))) % 4 for j in range(m))
    return (f1 ^ f2, c, mu)

def mneg(M):
    return (M[0], (M[1] + 2) % 4, M[2])

def mscale_i(M, k=1):
    return (M[0], (M[1] + k) % 4, M[2])

def ident(m):
    return (0, 0, (0,) * m)

def is_scalar(M):
    return M[0] == 0 and all(v % 4 == 0 for v in M[2])

def trace_is_zero(M):
    f, c, mu = M
    if f != 0:
        return True
    # trace = i^c * prod_j (1 + i^mu_j); zero iff some mu_j == 2
    return any(mu[j] % 4 == 2 for j in range(len(mu)))

def mtranspose(M):
    f, c, mu = M
    m = len(mu)
    add = sum(mu[j] for j in range(m) if (f >> j) & 1) % 4
    mu2 = tuple((mu[j] * (1 - 2 * ((f >> j) & 1))) % 4 for j in range(m))
    return (f, (c + add) % 4, mu2)

def madjoint(M):
    f, c, mu = mtranspose(M)
    return (f, (-c) % 4, tuple((-v) % 4 for v in mu))

def gamma(a, m):
    """Gamma matrix a (1-based) of Cl(2m) in the Pauli recursion:
    gamma_{2j-1} = I^(j-1) x sigma1 x sigma3^(m-j),
    gamma_{2j}   = I^(j-1) x sigma2 x sigma3^(m-j)."""
    j = (a + 1) // 2
    bit = j - 1
    f = 1 << bit
    if a % 2 == 1:
        c = 0
        mu = tuple(2 if q >= j else 0 for q in range(m))
    else:
        c = 1
        mu = tuple(2 if q >= bit else 0 for q in range(m))
    return (f, c, mu)

def mprod_list(mons):
    out = None
    for M in mons:
        out = M if out is None else mmul(out, M)
    return out

def so_bracket(a, b, c, d):
    """[e_ab, e_cd] for e_xy = E_xy - E_yx (antisymmetric extension).
    Returns dict {(x,y): coeff} with x<y, coeff in {+1,-1} (or empty)."""
    def E(x, y):
        if x == y:
            return None
        if x < y:
            return ((x, y), 1)
        return ((y, x), -1)
    terms = {}
    def add(term, sgn):
        if term is None:
            return
        pair, s = term
        terms[pair] = terms.get(pair, 0) + sgn * s
    if b == c:
        add(E(a, d), +1)
    if a == d:
        add(E(b, c), +1)
    if b == d:
        add(E(a, c), -1)
    if a == c:
        add(E(b, d), -1)
    return dict((k, v) for k, v in terms.items() if v != 0)
