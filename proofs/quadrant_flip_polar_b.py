import sympy as sp
t = sp.symbols('t')
def q(alpha, beta, tt):
    srx = sp.cancel((alpha*(1+tt**2) + (1-tt**2))/(2*tt))
    cxp = sp.cancel((beta*(1+tt**2) + 2*tt)/(1-tt**2))
    sxp, crx = sp.cancel(1/srx), sp.cancel(1/cxp)
    return dict(srx=srx, sxp=sxp, cxp=cxp, crx=crx, urx=sp.cancel(srx-crx), uxp=sp.cancel(cxp-sxp),
                sawup=sp.cancel(1/(srx-crx)), sawdown=sp.cancel(1/(cxp-sxp)))
QI = q(1, 1, t); QIII = q(-1, -1, t)
sub = lambda e: sp.cancel(e.subs(t, -1/t))
print('Pi-shift involution t -> -1/t maps each QI formula to the SAME operator QIII formula:')
for k in ['srx','sxp','cxp','crx','urx','uxp','sawup','sawdown']:
    print(f'  {k}: {sp.simplify(sub(QI[k]) - QIII[k]) == 0}')
M1_QIII = sp.Matrix([[QIII['urx'], QIII['srx']],[-QIII['cxp'], 1]])
M2_QI = sp.Matrix([[QI['uxp'], QI['sxp']],[-QI['crx'], 1]])
D = sp.diag(1, -1)
print('Block flip: M1_QIII == diag(1,-1) * M2_QI * diag(1,-1):',
      sp.simplify(M1_QIII - D*M2_QI*D) == sp.zeros(2))
