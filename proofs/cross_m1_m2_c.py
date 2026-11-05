import sympy as sp
t = sp.symbols('t')
def forms(alpha, beta):
    srx = sp.cancel((alpha*(1+t**2) + (1-t**2))/(2*t))
    cxp = sp.cancel((beta*(1+t**2) + 2*t)/(1-t**2))
    sxp, crx = sp.cancel(1/srx), sp.cancel(1/cxp)
    return srx,sxp,cxp,crx,sp.cancel(srx-crx),sp.cancel(cxp-sxp)
FWs = {'QI':1,'QII':-1,'QIII':1,'QIV':-1}
for name,(a,b) in {'QI':(1,1),'QII':(1,-1),'QIII':(-1,-1),'QIV':(-1,1)}.items():
    srx,sxp,cxp,crx,urx,uxp = forms(a,b)
    FW = FWs[name]
    M1 = sp.Matrix([[urx, srx],[-cxp, FW]]); M2 = sp.Matrix([[uxp, sxp],[-crx, FW]])
    print(f'{name}: tr(M1M2) =', sp.simplify(sp.trace(M1*M2)),
          '| urx*uxp - srx*crx - sxp*cxp =', sp.simplify(urx*uxp - srx*crx - sxp*cxp))
