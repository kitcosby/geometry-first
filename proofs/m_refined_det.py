from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace

# ---- symbolic det over QQ(t) ----
R = PolynomialRing(QQ, 't'); t = R.gen(); F = R.fraction_field(); u = F(t)
def rat(num, den=1): return F(num)/F(den)
sx = rat(2*u, 1+u**2); cx = rat(1-u**2, 1+u**2)
s2 = rat(4*u*(1-u**2), (1+u**2)**2); c2 = rat(1-6*u**2+u**4, (1+u**2)**2)
SRX, SXP = rat(1, u), u
CXP, CRX = rat(1+u, 1-u), rat(1-u, 1+u)
URX, UXP = SRX-CRX, CXP-SXP
SAWUP, SAWDOWN = 1/URX, 1/UXP
M = matrix(F, [[sx, CXP, -SXP, UXP],[SRX, cx, c2, SAWUP],[-CRX, s2, 0, 0],[URX, SAWDOWN, 0, 0]])
d = M.det()
print('det M_refined(t) numerator factor:', d.numerator().factor())
print('det M_refined(t) denominator factor:', d.denominator().factor())
A = (-SXP)*SAWUP - UXP*c2
B = (-CRX)*SAWDOWN - s2*URX
print('A numerator:', A.numerator().factor(), '| A denom:', A.denominator().factor())
print('B numerator:', B.numerator().factor(), '| B denom:', B.denominator().factor())
print('A*B == det:', A*B == d)

# ---- exact samples over QQbar, commutator closure ----
s2v = QQbar(2).sqrt(); s3v = QQbar(3).sqrt()
def sample(tt):
    b = {}
    b['sx'] = 2*tt/(1+tt**2); b['cx'] = (1-tt**2)/(1+tt**2)
    b['s2'] = 4*tt*(1-tt**2)/(1+tt**2)**2; b['c2'] = (1-6*tt**2+tt**4)/(1+tt**2)**2
    b['srx'] = 1/tt; b['sxp'] = tt; b['cxp'] = (1+tt)/(1-tt); b['crx'] = (1-tt)/(1+tt)
    b['urx'] = b['srx']-b['crx']; b['uxp'] = b['cxp']-b['sxp']
    b['sawup'] = 1/b['urx']; b['sawdown'] = 1/b['uxp']
    return matrix(QQbar, [[b['sx'], b['cxp'], -b['sxp'], b['uxp']],
                          [b['srx'], b['cx'], b['c2'], b['sawup']],
                          [-b['crx'], b['s2'], 0, 0],
                          [b['urx'], b['sawdown'], 0, 0]])
V = VectorSpace(QQbar, 16)
def flat(m): return V(list(m))
def closure_dim(mats):
    gens = [flat(m) for m in mats]
    W = V.span(gens); mats_g = list(mats)
    changed = True
    while changed:
        changed = False
        cur = list(mats_g)
        for i in range(len(cur)):
            for j in range(i+1, len(cur)):
                comm = cur[i]*cur[j] - cur[j]*cur[i]
                v = flat(comm)
                if v not in W:
                    mats_g.append(comm); W = V.span([flat(m) for m in mats_g]); changed = True
    # verify closure fully
    ok = all(flat(a*b - b*a) in W for a in mats_g for b in mats_g)
    return W.dimension(), ok
ts = [2-s3v, s2v-1, 1/s3v, QQbar(1)/2]
Ms = [sample(tt) for tt in ts]
print('sample traces:', [m.trace() for m in Ms])
r, ok = closure_dim(Ms)
print('M_refined: closure dim =', r, '| closed =', ok)
M0 = [m - (m.trace()/4)*matrix.identity(QQbar, 4) for m in Ms]
r0, ok0 = closure_dim(M0)
print('traceless-corrected: closure dim =', r0, '| closed =', ok0)
