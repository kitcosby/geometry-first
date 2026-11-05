import sys, time
sys.path.insert(0, '.')
import numpy as np
from pauli_monomial import gamma as mgamma, mmul, so_bracket

def dense_np(M):
    f, c, mu = M
    m = len(mu)
    d = 2 ** m
    out = np.zeros((d, d), dtype=np.complex128)
    xs = np.arange(d)
    q = np.full(d, c, dtype=int)
    for j in range(m):
        q += mu[j] * ((xs >> j) & 1)
    vals = (1j) ** (q % 4)
    out[xs ^ f, xs] = vals
    return out

def main():
    for n in (14, 15, 16):
        t0 = time.time()
        even = (n % 2 == 0)
        m = n // 2 if even else (n + 1) // 2
        gam = dict((a, mgamma(a, m)) for a in range(1, n + 1))
        pairs = [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]
        S = dict((p, dense_np(mmul(gam[p[0]], gam[p[1]])) / 2.0) for p in pairs)
        maxres = 0.0
        cnt = 0
        for ii in range(len(pairs)):
            for jj in range(ii + 1, len(pairs)):
                p, q = pairs[ii], pairs[jj]
                comm = S[p] @ S[q] - S[q] @ S[p]
                br = so_bracket(p[0], p[1], q[0], q[1])
                tgt = np.zeros_like(comm)
                for (xy, sgn) in br.items():
                    tgt += sgn * S[xy]
                res = float(np.max(np.abs(comm - tgt))) if comm.size else 0.0
                if res > maxres:
                    maxres = res
                cnt += 1
        print('NUMERIC n=%d pairs=%d max_residual=%.3e elapsed=%.1fs'
              % (n, cnt, maxres, time.time() - t0), flush=True)
    print('NUMERIC ALL DONE', flush=True)

main()
