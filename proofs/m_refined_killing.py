from sage.rings.qqbar import QQbar, AA
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace
s2v = QQbar(2).sqrt(); s3v = QQbar(3).sqrt()
def sample(tt):
    b = {}
    b['sx'] = 2*tt/(1+tt**2); b['cx'] = (1-tt**2)/(1+tt**2)
    b['s2'] = 4*tt*(1-tt**2)/(1+tt**2)**2; b['c2'] = (1-6*tt**2+tt**4)/(1+tt**2)**2
    b['srx'] = 1/tt; b['sxp'] = tt; b['cxp'] = (1+tt)/(1-tt); b['crx'] = (1-tt)/(1+tt)
    b['urx'] = b['srx']-b['crx']; b['uxp'] = b['cxp']-b['sxp']
    b['sawup'] = 1/b['urx']; b['sawdown'] = 1/b['uxp']
    m = matrix(QQbar, [[b['sx'], b['cxp'], -b['sxp'], b['uxp']],
                       [b['srx'], b['cx'], b['c2'], b['sawup']],
                       [-b['crx'], b['s2'], 0, 0],
                       [b['urx'], b['sawdown'], 0, 0]])
    return m - (m.trace()/4)*matrix.identity(QQbar, 4)
V = VectorSpace(QQbar, 16)
ts = [2-s3v, s2v-1, 1/s3v, QQbar(1)/2, QQbar(3)/4]
G = [sample(tt) for tt in ts]
W = V.span([V(m.list()) for m in G]); changed = True
while changed:
    changed = False
    cur = list(G)
    for i in range(len(cur)):
        for j in range(i+1, len(cur)):
            comm = cur[i]*cur[j] - cur[j]*cur[i]
            if V(comm.list()) not in W:
                G.append(comm); W = V.span([V(m.list()) for m in G]); changed = True
print('dim:', W.dimension())
bas = W.basis()
Bm = [matrix(QQbar, 4, 4, list(v)) for v in bas]
n = len(Bm)
Bmat = matrix(QQbar, 16, n, lambda r, c: bas[c][r])
ads = []
for i in range(n):
    Ai = matrix(QQbar, n, n)
    for a in range(n):
        comm = Bm[i]*Bm[a] - Bm[a]*Bm[i]
        y = Bmat.solve_right(V(comm.list()))
        for bb in range(n): Ai[bb, a] = y[bb]
    ads.append(Ai)
K = matrix(QQbar, n, n, lambda i, j: (ads[i]*ads[j]).trace())
evs = K.eigenvalues()
print('eigenvalue count:', len(evs))
allreal = all(e.imag() == 0 for e in evs)
POS = sum(1 for e in evs if e.real() > 0); NEG = sum(1 for e in evs if e.real() < 0); ZER = sum(1 for e in evs if e.real() == 0)
print('all eigenvalues real:', allreal)
print('SIGNATURE (pos, neg, zero):', (POS, NEG, ZER))
