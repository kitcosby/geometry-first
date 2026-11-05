import sys, time, random
sys.path.insert(0, '.')
from pauli_monomial import gamma as mgamma, mmul, mneg, mscale_i, mprod_list, so_bracket
from sage.rings.number_field.number_field import CyclotomicField
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.rational_field import QQ

Qi = CyclotomicField(4)
iv = Qi.gen()

def bxor(a, b):
    return a + b - 2 * (a & b)

def mdense(M):
    f, c, mu = M
    m = len(mu)
    d = 2 ** m
    rows = [[Qi(0)] * d for _ in range(d)]
    for x in range(d):
        q = c + sum(mu[j] for j in range(m) if (x >> j) & 1)
        rows[bxor(x, f)][x] = iv ** (q % 4)
    return matrix(Qi, d, d, [rows[i][j] for i in range(d) for j in range(d)])

def build(n):
    even = (n % 2 == 0)
    m = n // 2 if even else (n + 1) // 2
    gam = dict((a, mgamma(a, m)) for a in range(1, n + 1))
    pairs = [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]
    Mm = dict((p, mmul(gam[p[0]], gam[p[1]])) for p in pairs)
    return m, gam, pairs, Mm

def dense_sigmas(n):
    m, gam, pairs, Mm = build(n)
    S = dict((p, mdense(Mm[p]) / 2) for p in pairs)
    return m, gam, pairs, Mm, S

def rho_of_bracket(S, br):
    d = list(S.values())[0].nrows()
    out = matrix(Qi, d, d)
    for (xy, sgn) in br.items():
        out += sgn * S[xy]
    return out

def hom_check(n, pairlist=None):
    m, gam, pairs, Mm, S = dense_sigmas(n)
    idx = dict((p, k) for k, p in enumerate(pairs))
    todo = []
    if pairlist is None:
        for ii in range(len(pairs)):
            for jj in range(ii + 1, len(pairs)):
                todo.append((pairs[ii], pairs[jj]))
    else:
        todo = pairlist
    fails = 0
    for (p, q) in todo:
        comm = S[p] * S[q] - S[q] * S[p]
        br = so_bracket(p[0], p[1], q[0], q[1])
        if comm != rho_of_bracket(S, br):
            fails += 1
            print('  DENSE HOM FAIL', p, q, flush=True)
    return len(todo), fails

def gamma_dense(n):
    m, gam, pairs, Mm = build(n)
    Gd = dict((a, mdense(gam[a])) for a in range(1, n + 1))
    return m, Gd

def main():
    t00 = time.time()
    if MODE == 'control10':
        # dense gamma anticommutation, all pairs, n=10
        m, Gd = gamma_dense(10)
        d = 2 ** m
        Id = matrix(Qi, d, d, 1)
        ok = True
        for a in range(1, 11):
            for b in range(1, 11):
                v = Gd[a] * Gd[b] + Gd[b] * Gd[a]
                if v != (2 * Id if a == b else matrix(Qi, d, d)):
                    ok = False
        print('DENSE n=10 gamma anticommute all pairs:', ok, flush=True)
        cnt, fails = hom_check(10)
        print('DENSE HOM n=10 pairs=%d fails=%d' % (cnt, fails), flush=True)
    elif MODE == 'hom':
        cnt, fails = hom_check(N)
        print('DENSE HOM n=%d pairs=%d fails=%d' % (N, cnt, fails), flush=True)
    elif MODE == 'sample':
        m, gam, pairs, Mm, S = dense_sigmas(N)
        random.seed(20261016 + N)
        allpairs = [(pairs[i], pairs[j]) for i in range(len(pairs))
                    for j in range(i + 1, len(pairs))]
        todo = random.sample(allpairs, 50)
        fails = 0
        for (p, q) in todo:
            comm = S[p] * S[q] - S[q] * S[p]
            br = so_bracket(p[0], p[1], q[0], q[1])
            if comm != rho_of_bracket(S, br):
                fails += 1
                print('  SAMPLE FAIL', p, q, flush=True)
        print('DENSE SAMPLE n=%d pairs=50 fails=%d' % (N, fails), flush=True)
    elif MODE == 'ranks':
        for n in range(11, 17):
            t0 = time.time()
            even = (n % 2 == 0)
            m = n // 2 if even else (n + 1) // 2
            d = 2 ** m
            gam = dict((a, mgamma(a, m)) for a in range(1, n + 1))
            pairs = [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]
            # vector Casimir over QQ
            C2 = matrix(QQ, n, n)
            for (a, b) in pairs:
                E = matrix(QQ, n, n)
                E[a - 1, b - 1] = 1
                E[b - 1, a - 1] = -1
                C2 += E * E
            vec_ok = (C2 == -(n - 1) * matrix(QQ, n, n, 1))
            # mod-p rank certificate for span of T_ab = 2 Sigma_ab
            res = {}
            for p in (1000033, 104729):
                Fp = GF(p)
                r = Fp(-1).sqrt()
                ri = int(r)
                flat = [0] * (len(pairs) * d * d)
                for pi_, (a, b) in enumerate(pairs):
                    Mm = mmul(gam[a], gam[b])
                    f, c, mu = Mm
                    base = pi_ * d * d
                    for x in range(d):
                        q = (c + sum(mu[j] for j in range(m) if (x >> j) & 1)) % 4
                        re, im = ((1, 0), (0, 1), (-1, 0), (0, -1))[q]
                        flat[base + bxor(x, f) * d + x] = (re + ri * im) % p
                Mflat = matrix(Fp, len(pairs), d * d, flat)
                res[p] = Mflat.rank()
            # chirality/volume split ranks mod p
            w = mprod_list([gam[a] for a in range(1, n + 1)])
            from pauli_monomial import mscale_i as _ms
            w2 = mmul(w, w)
            G = w if w2[1] % 4 == 0 else _ms(w)
            p = 1000033
            Fp = GF(p)
            ri = int(Fp(-1).sqrt())
            f, c, mu = G
            Gp = [[0] * d for _ in range(d)]
            for x in range(d):
                q = (c + sum(mu[j] for j in range(m) if (x >> j) & 1)) % 4
                re, im = ((1, 0), (0, 1), (-1, 0), (0, -1))[q]
                Gp[bxor(x, f)][x] = (re + ri * im) % p
            Gm = matrix(Fp, d, d, [Gp[i][j] for i in range(d) for j in range(d)])
            Idm = matrix(Fp, d, d, 1)
            rk_plus = (Gm - Idm).rank()
            rk_minus = (Gm + Idm).rank()
            print('RANKS n=%d vecCasimir=-(%d)I ok=%s spanrank mod1000033=%d mod104729=%d '
                  'target=%d dim(ker G-I)=%d dim(ker G+I)=%d elapsed=%.1fs'
                  % (n, n - 1, vec_ok, res[1000033], res[104729], len(pairs),
                     d - rk_plus, d - rk_minus, time.time() - t0), flush=True)
    print('DENSE %s DONE elapsed=%.1fs' % (MODE, time.time() - t00), flush=True)

main()
