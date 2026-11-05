import sympy as sp
t = sp.symbols('t')
def forms(alpha, beta):
    srx = sp.cancel((alpha*(1+t**2) + (1-t**2))/(2*t))
    cxp = sp.cancel((beta*(1+t**2) + 2*t)/(1-t**2))
    sxp, crx = sp.cancel(1/srx), sp.cancel(1/cxp)
    urx, uxp = sp.cancel(srx-crx), sp.cancel(cxp-sxp)
    FW = 1 if alpha*beta*(1 if True else 0) is not None else None
    return srx,sxp,cxp,crx,urx,uxp
FWs = {'QI':1,'QII':-1,'QIII':1,'QIV':-1}
S = lambda e: sp.factor(sp.cancel(e))
for name,(a,b) in {'QI':(1,1),'QII':(1,-1),'QIII':(-1,-1),'QIV':(-1,1)}.items():
    srx,sxp,cxp,crx,urx,uxp = forms(a,b)
    FW = FWs[name]
    M1 = sp.Matrix([[urx, srx],[-cxp, FW]])
    M2 = sp.Matrix([[uxp, sxp],[-crx, FW]])
    C = sp.simplify(M1*M2 - M2*M1)
    print(f'--- {name} (FW={FW}) ---')
    print('  tr M1 =', S(sp.trace(M1)), '| det M1 =', S(M1.det()))
    print('  tr M2 =', S(sp.trace(M2)), '| det M2 =', S(M2.det()))
    print('  [M1,M2] =', [[S(C[i,j]) for j in range(2)] for i in range(2)])
    print('  tr [M1,M2] =', S(sp.trace(C)), '| det [M1,M2] =', S(C.det()))
    P = sp.simplify(M1*M2)
    print('  M1*M2 =', [[S(P[i,j]) for j in range(2)] for i in range(2)])
