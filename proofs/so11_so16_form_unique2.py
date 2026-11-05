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
def solve(n, sign, sym, p):
    even = (n % 2 == 0)
    m = n // 2 if even else (n + 1) // 2
    d = 2 ** m
    Fp = GF(p)
    ri = int(Fp(-1).sqrt())
    gam = dict((a, mgamma(a, m)) for a in range(1, n + 1))
    w = mprod_list([gam[a] for a in range(1, n + 1)])
    w2 = mmul(w, w)
    G = w if w2[1] % 4 == 0 else mscale_i(w)
    Gm = dense_modp(G, p, ri)
    Idm = matrix(Fp, d, d, 1)
    blk = (Gm - Idm).right_kernel().basis() if sign > 0 else (Gm + Idm).right_kernel().basis()
    dp = len(blk)
    oth = (Gm + Idm).right_kernel().basis() if sign > 0 else (Gm - Idm).right_kernel().basis()
    cols = [list(v) for v in blk] + [list(v) for v in oth]
    P = matrix(Fp, d, d, [[cols[k][i] for k in range(d)] for i in range(d)])
    Pinv = P.inverse()
    blocks = []
    for a in range(1, n):
        Tm = dense_modp(mmul(gam[a], gam[a + 1]), p, ri)
        B = Pinv * Tm * P
        blocks.append(B.submatrix(0, 0, dp, dp))
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
    cur = matrix(Fp, nu, nu, 1)
    k = nu
    for gi, B in enumerate(blocks):
        Bt = B.transpose()
        colsA = []
        for F in basis:
            L = Bt * F + F * B
            if sym:
                colsA.append([L[i, j] for i in range(dp) for j in range(i, dp)])
            else:
                colsA.append([L[i, j] for i in range(dp) for j in range(i + 1, dp)])
        neq = len(colsA[0])
        A = matrix(Fp, neq, nu, [[colsA[kk][e] for kk in range(nu)] for e in range(neq)])
        M = A * cur
        K = M.right_kernel().basis()
        k = len(K)
        if k == 0:
            break
        Kmat = matrix(Fp, cur.ncols(), k, [[K[j][i] for j in range(k)] for i in range(cur.ncols())])
        cur = cur * Kmat
    return dp, k
def main():
    for (n, sign) in JOBS:
        for sym in (True, False):
            for p in (1000033, 104729):
                t0 = time.time()
                dp, k = solve(n, sign, sym, p)
                print('UNIQUE2 n=%d sign=%+d type=%s p=%d blockdim=%d sol_dim=%d elapsed=%.1fs'
                      % (n, sign, 'sym' if sym else 'antisym', p, dp, k, time.time() - t0), flush=True)
    print('UNIQUE2 DONE', flush=True)
main()
