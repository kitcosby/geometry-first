import time
from sage.rings.qqbar import QQbar
from sage.matrix.constructor import matrix
K = QQbar
IU = K(-1).sqrt()
s2 = K(2).sqrt(); s3 = K(3).sqrt(); s6 = s2 * s3
TRIG = {
 15: ((s6 - s2) / 4, (s6 + s2) / 4), 30: (K(1) / 2, s3 / 2),
 45: (s2 / 2, s2 / 2), 60: (s3 / 2, K(1) / 2),
 75: ((s6 + s2) / 4, (s6 - s2) / 4),
 105: ((s6 + s2) / 4, -(s6 - s2) / 4), 120: (s3 / 2, -K(1) / 2),
 135: (s2 / 2, -s2 / 2), 150: (K(1) / 2, -s3 / 2),
 165: ((s6 - s2) / 4, -(s6 + s2) / 4)}
DEGS = [15, 30, 45, 60, 75, 105, 120, 135, 150, 165]
def ops(deg):
    s, c = TRIG[deg]
    absc = c if deg < 90 else -c
    srx = 1 / s + c / s
    sxp = 1 / srx
    cxp = 1 / absc + s / c
    crx = 1 / cxp
    urx = srx - crx
    uxp = cxp - sxp
    FW = 1 / urx + 1 / uxp
    return (srx, sxp, cxp, crx, urx, uxp, FW, 2 * s * c, c * c - s * s)
def build3(deg):
    srx, sxp, cxp, crx, urx, uxp, FW, sin2, cos2 = ops(deg)
    M = [[K(0) for _ in range(3)] for _ in range(3)]
    diag = [urx, uxp, FW]
    for i in range(3):
        M[i][i] = diag[i]
    for i in range(1, 3):
        M[i - 1][i] = srx if i % 2 == 1 else crx
        M[i][i - 1] = cxp if i % 2 == 1 else sxp
    M[0][2] = cos2; M[2][0] = sin2
    Mx = matrix(K, 3, 3, [M[i][j] for i in range(3) for j in range(3)])
    tr = Mx.trace()
    Mx = Mx - (tr / 3) * matrix(K, 3, 3, [K(1) if i == j else K(0) for i in range(3) for j in range(3)])
    return Mx
def reim(z):
    r = (z + z.conjugate()) / 2
    i2 = (z - z.conjugate()) / (2 * IU)
    return r, i2
def rvec(M):
    v = []
    for i in range(3):
        for j in range(3):
            v.append(reim(M[i, j])[0])
    for i in range(3):
        for j in range(3):
            v.append(reim(M[i, j])[1])
    return v
class RSpan:
    def __init__(self, L):
        self.L = L; self.piv = []; self.rows = []
    def try_add(self, v):
        v = list(v)
        for (p, row) in zip(self.piv, self.rows):
            if v[p] != 0:
                f = v[p]
                v = [a - f * b for a, b in zip(v, row)]
        p = -1
        for idx in range(self.L):
            if v[idx] != 0:
                p = idx; break
        if p < 0:
            return False
        pv = v[p]
        v = [a / pv for a in v]
        self.piv.append(p); self.rows.append(v)
        return True
    def dim(self):
        return len(self.rows)
def comm(A, B):
    return A * B - B * A
def real_closure(gens):
    S = RSpan(18)
    basis = []
    for g in gens:
        if S.try_add(rvec(g)):
            basis.append(g)
    changed = True
    while changed:
        changed = False
        cur = list(basis)
        for a in range(len(cur)):
            for b in range(a + 1, len(cur)):
                C = comm(cur[a], cur[b])
                if S.try_add(rvec(C)):
                    basis.append(C); changed = True
    return basis, S.dim()
def main():
    t0 = time.time()
    samples = [build3(d) for d in DEGS]
    gval = {d: (K(1) if d < 90 else IU) for d in DEGS}
    jf = [(-gval[d]) * samples[k] for k, d in enumerate(DEGS)]
    basisT, dimT = real_closure(samples + jf)
    print('TEST ARM (samples + J-partners): real closure dim:', dimT, flush=True)
    basisR, dimR = real_closure(samples)
    print('real samples alone: real closure dim:', dimR, flush=True)
    anti = []
    sym = []
    import itertools
    B = basisR
    rowsA = []
    for a in range(9):
        pass
    def coeff_solve(want_anti):
        eqs = []
        for i in range(3):
            for j in range(3):
                row = []
                for X in B:
                    row.append(X[i, j] + X[j, i] if want_anti else X[i, j] - X[j, i])
                eqs.append(row)
        S = RSpan(8)
        sols = []
        free = []
        basis_vecs = []
        Mrows = []
        for r in eqs:
            Mrows.append(list(r))
        piv_cols = []
        rr = 0
        MM = [list(r) for r in Mrows]
        for c in range(8):
            piv = -1
            for i in range(rr, len(MM)):
                if MM[i][c] != 0:
                    piv = i; break
            if piv < 0:
                continue
            MM[rr], MM[piv] = MM[piv], MM[rr]
            pv = MM[rr][c]
            MM[rr] = [x / pv for x in MM[rr]]
            for i in range(len(MM)):
                if i != rr and MM[i][c] != 0:
                    f = MM[i][c]
                    MM[i] = [x - f * y for x, y in zip(MM[i], MM[rr])]
            piv_cols.append(c); rr += 1
        free = [c for c in range(8) if c not in piv_cols]
        out = []
        for fc in free:
            coef = [K(0)] * 8
            coef[fc] = K(1)
            for i, pc in enumerate(piv_cols):
                coef[pc] = -MM[i][fc]
            X = sum((coef[k] * B[k] for k in range(8)), matrix(K, 3, 3, [K(0)] * 9))
            out.append(X)
        return out
    anti = coeff_solve(True)
    sym = coeff_solve(False)
    print('CONTROL: antisymmetric (compact) basis dim:', len(anti), ' symmetric dim:', len(sym), flush=True)
    wick = anti + [IU * X for X in sym]
    basisW, dimW = real_closure(wick)
    print('CONTROL (standard Wick on full basis): real closure dim:', dimW, flush=True)
    if dimW == 8:
        A = [[6 * (basisW[i] * basisW[j]).trace() for j in range(8)] for i in range(8)]
        allreal = all(reim(A[i][j])[1] == 0 for i in range(8) for j in range(8))
        print('Killing Gram entries all real:', allreal, flush=True)
        A = [[reim(A[i][j])[0] for j in range(8)] for i in range(8)]
        pos = neg = zer = 0
        d = 8
        for k in range(d):
            if A[k][k] == 0:
                sw = -1
                for i in range(k + 1, d):
                    if A[i][i] != 0:
                        sw = i; break
                if sw >= 0:
                    for r in range(d):
                        A[r][k], A[r][sw] = A[r][sw], A[r][k]
                    for cc in range(d):
                        A[k][cc], A[sw][cc] = A[sw][cc], A[k][cc]
                else:
                    add = -1
                    for i in range(k + 1, d):
                        if A[i][k] != 0:
                            add = i; break
                    if add < 0:
                        zer += 1; continue
                    for cc in range(d):
                        A[k][cc] = A[k][cc] + A[add][cc]
                    for r2 in range(d):
                        A[r2][k] = A[r2][k] + A[r2][add]
            piv = A[k][k]
            if piv > 0:
                pos += 1
            elif piv < 0:
                neg += 1
            else:
                zer += 1
            if piv != 0:
                for i in range(k + 1, d):
                    f = A[i][k]
                    if f != 0:
                        f = f / piv
                        for j in range(k, d):
                            A[i][j] = A[i][j] - f * A[k][j]
        print('CONTROL Killing signature (pos, neg, zero):', (pos, neg, zer), flush=True)
    print('elapsed %.1fs' % (time.time() - t0), flush=True)
    print('RUN DONE', flush=True)
main()
