from sage.rings.qqbar import AA
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace
s2v = AA(2).sqrt(); s3v = AA(3).sqrt()
def sample(tt):
    b = {}
    b['sx'] = 2*tt/(1+tt**2); b['cx'] = (1-tt**2)/(1+tt**2)
    b['s2'] = 4*tt*(1-tt**2)/(1+tt**2)**2; b['c2'] = (1-6*tt**2+tt**4)/(1+tt**2)**2
    b['srx'] = 1/tt; b['sxp'] = tt; b['cxp'] = (1+tt)/(1-tt); b['crx'] = (1-tt)/(1+tt)
    b['urx'] = b['srx']-b['crx']; b['uxp'] = b['cxp']-b['sxp']
    b['sawup'] = 1/b['urx']; b['sawdown'] = 1/b['uxp']
    m = matrix(AA, [[b['sx'], b['cxp'], -b['sxp'], b['uxp']],
                    [b['srx'], b['cx'], b['c2'], b['sawup']],
                    [-b['crx'], b['s2'], 0, 0],
                    [b['urx'], b['sawdown'], 0, 0]])
    return m - (m.trace()/4)*matrix.identity(AA, 4)
V = VectorSpace(AA, 16)
ts = [2-s3v, s2v-1, 1/s3v, AA(1)/2, AA(3)/4]
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
Bm = [matrix(AA, 4, 4, list(v)) for v in bas]
n = len(Bm)
Bmat = matrix(AA, 16, n, lambda r, c: bas[c][r])
ads = []
for i in range(n):
    Ai = matrix(AA, n, n)
    for a in range(n):
        comm = Bm[i]*Bm[a] - Bm[a]*Bm[i]
        y = Bmat.solve_right(V(comm.list()))
        for bb in range(n): Ai[bb, a] = y[bb]
    ads.append(Ai)
K = matrix(AA, n, n, lambda i, j: (ads[i]*ads[j]).trace())
ev = K.eigenvalues()
pos = sum(1 for e in ev if e > 0); neg = sum(1 for e in ev if e < 0); zer = sum(1 for e in ev if e == 0)
print('Killing signature (pos, neg, zero):', (pos, neg, zer))
print('Killing eigenvalues:', sorted([str(e) for e in ev]))
