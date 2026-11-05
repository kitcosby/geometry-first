from sage.rings.qqbar import QQbar
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
    return matrix(QQbar, [[b['sx'], b['cxp'], -b['sxp'], b['uxp']],
                          [b['srx'], b['cx'], b['c2'], b['sawup']],
                          [-b['crx'], b['s2'], 0, 0],
                          [b['urx'], b['sawdown'], 0, 0]])
V = VectorSpace(QQbar, 16)
def closure(mats):
    W = V.span([V(m.list()) for m in mats]); G = list(mats); changed = True
    while changed:
        changed = False
        cur = list(G)
        for i in range(len(cur)):
            for j in range(i+1, len(cur)):
                comm = cur[i]*cur[j] - cur[j]*cur[i]
                if V(comm.list()) not in W:
                    G.append(comm); W = V.span([V(m.list()) for m in G]); changed = True
    ok = all(V((a*b-b*a).list()) in W for a in G for b in G)
    return W.dimension(), ok, G
ts = [2-s3v, s2v-1, 1/s3v, QQbar(1)/2, QQbar(3)/4]
Ms = [sample(tt) for tt in ts]
r, ok, G = closure(Ms)
print('M_refined: closure dim =', r, '| closed =', ok, '| traceless all gens:', all(m.trace()==0 for m in G))
M0 = [m - (m.trace()/4)*matrix.identity(QQbar, 4) for m in Ms]
r0, ok0, G0 = closure(M0)
print('traceless-corrected: closure dim =', r0, '| closed =', ok0)
