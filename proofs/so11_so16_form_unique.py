import sys, time
sys.path.insert(0, '.')
from pauli_monomial import gamma as mgamma, mmul, mprod_list, mscale_i
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF

VALS = ((1, 0), (0, 1), (-1, 0), (0, -1))

def bxor(a, b):
    return a + b - 2 * (a & b)

def dense_modp(M, p, ri):
    f, c, mu = M
    m = len(mu)
    d = 2 ** m
    rows = [[0] * d for _ in range(d)]
    for x in range(d):
        q = (c + sum(mu[j] for j in range(m) if (x >> j) & 1)) % 4
        re, im = VALS[q]
        rows[bxor(x, f)][x] = (re + ri * im) % p
    Fp = GF(p)
    return matrix(Fp, d, d, [rows[i][j] for i in range(d) for j in range(d)])

def solve_for_n_p(n, p):
    even = (n % 2 == 0)
    m = n // 2 if even else (n + 1) // 2
    d = 2 ** m
    Fp = GF(p)
    ri = int(Fp(-1).sqrt())
    gam = dict((a, mgamma(a, m)) for a in range(1, n + 1))
    # chirality / volume monomial
    w = mprod_list([gam[a] for a in range(1, n + 1)])
    w2 = mmul(w, w)
    G = w if w2[1] % 4 == 0 else mscale_i(w)
    Gm = dense_modp(G, p, ri)
    Idm = matrix(Fp, d, d, 1)
    plus = (Gm - Idm).right_kernel().basis()
    minus = (Gm + Idm).right_kernel().basis()
    dp = len(plus)
    assert dp + len(minus) == d
    cols = [list(v) for v in plus] + [list(v) for v in minus]
    P = matrix(Fp, d, d, [[cols[k][i] for k in range(d)] for i in range(d)])
    Pinv = P.inverse()
    # generator blocks: T_ab = 2 Sigma_ab for adjacent pairs (a, a+1)
    blocks = []
    for a in range(1, n):
        Tm = dense_modp(mmul(gam[a], gam[a + 1]), p, ri)
        B = Pinv * Tm * P
        blocks.append(B.submatrix(0, 0, dp, dp))
    out = {}
    for sym in (True, False):
        basis = []
        if sym:
            for i in range(dp):
                for j in range(i, dp):
                    E = matrix(Fp, dp, dp)
                    E[i, j] = 1
                    if i != j:
                        E[j, i] = 1
                    basis.append(E)
        else:
            for i in range(dp):
                for j in range(i + 1, dp):
                    E = matrix(Fp, dp, dp)
                    E[i, j] = 1
                    E[j, i] = p - 1
                    basis.append(E)
        nu = len(basis)
        eqs = []
        for B in blocks:
            Bt = B.transpose()
            grow = []
            for F in basis:
                L = Bt * F + F * B
                if sym:
                    grow.append([L[i, j] for i in range(dp) for j in range(i, dp)])
                else:
                    grow.append([L[i, j] for i in range(dp) for j in range(i + 1, dp)])
            eqs.append(grow)
        rows_list = [[g[k][e] for k in range(nu)]
                     for g in eqs for e in range(len(g[0]))]
        A = matrix(Fp, len(rows_list), nu, [v for row in rows_list for v in row])
        out['sym' if sym else 'antisym'] = nu - A.rank()
    return dp, out

def main():
    for n in (11, 12):
        for p in (1000033, 104729):
            t0 = time.time()
            dp, out = solve_for_n_p(n, p)
            print('UNIQUE n=%d p=%d blockdim=%d sym_dim=%d antisym_dim=%d elapsed=%.1fs'
                  % (n, p, dp, out['sym'], out['antisym'], time.time() - t0), flush=True)
    print('UNIQUE DONE', flush=True)

main()
