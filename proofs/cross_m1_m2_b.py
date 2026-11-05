import sympy as sp
t = sp.symbols('t')
srx, sxp = 1/t, t
cxp, crx = (1+t)/(1-t), (1-t)/(1+t)
urx, uxp = sp.cancel(srx-crx), sp.cancel(cxp-sxp)
FW = 1
M1 = sp.Matrix([[urx, srx],[-cxp, FW]])
M2 = sp.Matrix([[uxp, sxp],[-crx, FW]])
A = sp.simplify(M1*M2 + M2*M1)
print('QI anticommutator {M1,M2} tr =', sp.factor(sp.trace(A)), ' (2 tr(M1M2))')
C = sp.simplify(M1*M2 - M2*M1)
print('[M1,M2]_(1,2) =', sp.factor(C[0,1]), ' = -(urx+uxp) sin^2(x)? sin^2 in t: 4t^2/(1+t^2)^2 ->',
      sp.simplify(C[0,1] + (urx+uxp)*(4*t**2/(1+t**2)**2)) == 0)
print('[M1,M2]_(1,1) =', sp.factor(C[0,0]))
# Exact values at 30 deg.
tv = 2 - sp.sqrt(3)
for nm, M in [('M1', M1), ('M2', M2), ('[M1,M2]', C)]:
    print(nm, 'at 30deg: tr =', sp.nsimplify(sp.trace(M).subs(t, tv), [sp.sqrt(3)]),
          '| det =', sp.nsimplify(M.det().subs(t, tv), [sp.sqrt(3)]))
K = sp.kronecker_product(M1, M2)
print('Kronecker: det(M1 x M2) = det(M1)^2 det(M2)^2 :', sp.simplify(K.det() - (M1.det()**2)*(M2.det()**2)) == 0,
      '| tr =', sp.factor(sp.trace(K)), '(= trM1 trM2)')
print('eigvals M1 at 30deg:', [sp.nsimplify(e, [sp.sqrt(3)]) for e in (M1.subs(t, tv)).eigenvals()])
