import sympy as sp
t = sp.symbols('t', real=True)
# QI forms.
sxp = t; cxp = (1+t)/(1-t); uxp = sp.cancel(cxp - sxp)
th2 = sp.factor(sxp**2 + cxp**2 + uxp**2)
print('uxp =', uxp)
print('theta^2 = sxp^2+cxp^2+uxp^2 =', th2)
# Try: is theta^2 a perfect square of a rational function?
print('theta^2 = 2*(sxp^2 - sxp*cxp + cxp^2):', sp.simplify(th2 - 2*(sxp**2 - sxp*cxp + cxp**2)) == 0)
# Exact at 30 deg: t = tan(15deg) = 2 - sqrt(3).
tv = 2 - sp.sqrt(3)
vals = {k: sp.nsimplify(v.subs(t, tv), [sp.sqrt(3)]) for k, v in [('sxp', sxp), ('cxp', cxp), ('uxp', uxp)]}
print('at 30deg:', vals)
th2v = sp.nsimplify(th2.subs(t, tv), [sp.sqrt(3)])
print('theta^2 at 30deg =', th2v, ' theta =', sp.nsimplify(sp.sqrt(th2v), [sp.sqrt(3)]), '=', float(sp.sqrt(th2v)))
