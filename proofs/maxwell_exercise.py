import sympy as sp
t = sp.symbols('t', positive=True)
# QI: d/dx = ((1+t^2)/2) d/dt, t = tan(x/2).
D = lambda f: sp.cancel(((1+t**2)/2) * sp.diff(f, t))
srx, sxp = 1/t, t
cxp, crx = (1+t)/(1-t), (1-t)/(1+t)
urx, uxp = sp.cancel(srx-crx), sp.cancel(cxp-sxp)
print('Riccati system (QI):')
print('  srx\' = -(1+srx^2)/2 :', sp.simplify(D(srx) + (1+srx**2)/2) == 0)
print('  crx\' = -(1+crx^2)/2 :', sp.simplify(D(crx) + (1+crx**2)/2) == 0)
print('  sxp\' = +(1+sxp^2)/2 :', sp.simplify(D(sxp) - (1+sxp**2)/2) == 0)
print('  cxp\' = +(1+cxp^2)/2 :', sp.simplify(D(cxp) - (1+cxp**2)/2) == 0)
print('UNA log-derivatives (note: the diagonal cross sums from before):')
print('  urx\'/urx = -(srx+crx)/2 :', sp.simplify(D(urx)/urx + (srx+crx)/2) == 0)
print('  uxp\'/uxp = +(cxp+sxp)/2 :', sp.simplify(D(uxp)/uxp - (cxp+sxp)/2) == 0)
# Linearization: y = 2u'/u for the minus class gives u'' + u/4 = 0.
x = sp.symbols('x', real=True)
u = sp.Function('u')
print('Oscillator check: sin(x/2), cos(x/2) solve u\'\' + u/4 = 0:',
      sp.simplify(sp.diff(sp.sin(x/2), x, 2) + sp.sin(x/2)/4) == 0,
      sp.simplify(sp.diff(sp.cos(x/2), x, 2) + sp.cos(x/2)/4) == 0)
print('srx = 2 (d/dx log sin(x/2)) in QI:',
      sp.simplify(2*sp.diff(sp.log(sp.sin(x/2)), x) - 1/sp.tan(x/2)) == 0)
# Curl operator = so(3) cross matrix (uxp placement): eigenvalues 0, +/- i*theta (helicity).
a,b,c = sp.symbols('a b c')
Y = sp.Matrix([[0,-c,b],[c,0,-a],[-b,a,0]])
lam = sp.symbols('lam')
print('charpoly of curl matrix:', sp.factor(Y.charpoly(lam).as_expr()))
