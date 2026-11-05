import time
from sage.rings.qqbar import QQbar
from sage.matrix.constructor import matrix

K = QQbar
s3 = K(3).sqrt()
s = K(1) / 2
c = s3 / 2
# operators at x = 30 deg, from the definitions (sin, cos > 0 in QI)
srx = 1 / s + c / s
sxp = 1 / srx
cxp = 1 / c + s / c
crx = 1 / cxp
urx = srx - crx
uxp = cxp - sxp
sawup = 1 / urx
sawdown = 1 / uxp
FW = 1 / urx + 1 / uxp
assert FW == 1
W = [srx, sxp, cxp, crx, urx, sawup, uxp, sawdown, FW]
print('weights ready: srx=%s cxp=%s urx=%s uxp=%s FW=%s' % (srx, cxp, urx, uxp, FW), flush=True)

class Span:
    def __init__(self, L):
        self.L = L
        self.rows = []
        self.piv = []
    def try_add(self, v):
        v = list(v)
        for r, p in zip(self.rows, self.piv):
            f = v[p]
            if f != 0:
                v = [a - f * b for a, b in zip(v, r)]
        p = -1
        for i, a in enumerate(v):
            if a != 0:
                p = i
                break
        if p < 0:
            return False
        f = v[p]
        v = [a / f for a in v]
        for i2 in range(len(self.rows)):
            f2 = self.rows[i2][p]
            if f2 != 0:
                self.rows[i2] = [a - f2 * b for a, b in zip(self.rows[i2], v)]
        self.rows.append(v)
        self.piv.append(p)
        return True
    def dim(self):
        return len(self.rows)

def run_path(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pidx = dict((p, k) for k, p in enumerate(pairs))
    def uncoords(v):
        M = matrix(K, n, n)
        for k, (i, j) in enumerate(pairs):
            M[i, j] = v[k]
            M[j, i] = -v[k]
        return M
    B = Span(len(pairs))
    for i in range(1, n):
        f = W[(i - 1) % 9]
        v = [K(0)] * len(pairs)
        v[pidx[(i - 1, i)]] = f
        B.try_add(v)
    span0 = B.dim()
    changed = True
    while changed:
        changed = False
        mats = [uncoords(B.rows[r]) for r in range(B.dim())]
        L = len(mats)
        for a in range(L):
            for b in range(a + 1, L):
                comm = mats[a] * mats[b] - mats[b] * mats[a]
                v = [comm[i, j] for (i, j) in pairs]
                if B.try_add(v):
                    changed = True
    return span0, B.dim()

def main():
    for n in range(11, 17):
        t0 = time.time()
        span0, d = run_path(n)
        print('PATH n=%d gen-span=%d closure=%d target=%d match=%s elapsed=%.1fs'
              % (n, span0, d, n * (n - 1) // 2, d == n * (n - 1) // 2, time.time() - t0), flush=True)
    print('PATH ALL DONE', flush=True)

main()
