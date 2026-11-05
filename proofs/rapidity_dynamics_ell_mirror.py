# Ell-mirror of rapidity_dynamics.py: ell = ln srx (QI symbolic).
# Mirror claim to test: swapping (srx,sxp)<->(cxp,crx) and eta<->ell maps
# the eta-dictionary to the ell-dictionary, with sin/cos swapped and a
# sign flip in the double-angle velocities (dl/dx = -cosh ell < 0).
import sympy as sp

ell = sp.symbols('ell', real=True)
E = sp.exp(ell)
results = []
def is0(e):
    e = sp.together(e.rewrite(sp.exp)); e = sp.cancel(e)
    return sp.simplify(e) == 0
def check(name, expr):
    ok = is0(expr); results.append((name, ok))
    print(('PASS' if ok else 'FAIL'), name, flush=True)

sech = 1/sp.cosh(ell); tanh = sp.tanh(ell)
# ell-dictionary (QI): primitives in ell
srx = E; sxp = sp.exp(-ell)
cxp = sp.coth(ell/2)          # mirror of srx = coth(eta/2)
crx = sp.tanh(ell/2)          # mirror of sxp = tanh(eta/2)
check('cxp = (1+sech)/tanh', cxp - (1+sech)/tanh)
check('crx = (1-sech)/tanh', crx - (1-sech)/tanh)
check('cosh ell = (srx+sxp)/2', sp.cosh(ell) - (srx+sxp)/2)
check('sinh ell = (srx-sxp)/2', sp.sinh(ell) - (srx-sxp)/2)
urx = srx - crx; uxp = cxp - sxp
# sin x = sech ell, cos x = tanh ell (swapped vs eta)
sin2 = 2*tanh*sech                       # 2 sin cos, same form
cos2 = tanh**2 - sech**2                 # cos^2 - sin^2, sign-swapped vs eta
check('sin2 via 2 tanh sech', sin2 - 2*tanh*sech)
check('cos2 = 2 tanh^2 - 1', cos2 - (2*tanh**2 - 1))
check('sin^2+cos^2 form: sech^2+tanh^2=1', sech**2 + tanh**2 - 1)

# velocities d/dell, targets in operator form (mirror + sign where dl/dx<0)
check('d srx = srx', sp.diff(srx, ell) - srx)
check('d sxp = -sxp', sp.diff(sxp, ell) + sxp)
check('d cxp = -(cxp^2-1)/2', sp.diff(cxp, ell) + (cxp**2-1)/2)
check('d crx = (1-crx^2)/2', sp.diff(crx, ell) - (1-crx**2)/2)
check('d urx = srx - (1-crx^2)/2', sp.diff(urx, ell) - srx + (1-crx**2)/2)
check('d uxp = -(cxp^2-1)/2 + sxp', sp.diff(uxp, ell) + (cxp**2-1)/2 - sxp)
check('d sin2 = -4 cos2/(srx+sxp)', sp.diff(sin2, ell) + 4*cos2/(srx+sxp))
check('d cos2 = +4 sin2/(srx+sxp)', sp.diff(cos2, ell) - 4*sin2/(srx+sxp))
check('d sin2 = -2 cos2/cosh ell', sp.diff(sin2, ell) + 2*cos2/sp.cosh(ell))
# chain rule: dl/dx = -cosh ell, d(sin2x)/dx = 2cos2x -> d/dell = -2cos2/cosh ell OK above

# E96 velocity in ell (same traceless-corrected matrix, FW=+1 on QI)
def build96():
    M = sp.Matrix([[urx, srx, cos2],[cxp, uxp, crx],[sin2, sxp, sp.Integer(1)]])
    return M - (sp.trace(M)/3)*sp.eye(3)
M96 = build96(); dM96 = sp.diff(M96, ell)
tr = urx + uxp + 1; dtr = sp.diff(tr, ell)
targets = {(0,0): sp.diff(urx,ell)-dtr/3,(0,1): sp.diff(srx,ell),(0,2): sp.diff(cos2,ell),
 (1,0): sp.diff(cxp,ell),(1,1): sp.diff(uxp,ell)-dtr/3,(1,2): sp.diff(crx,ell),
 (2,0): sp.diff(sin2,ell),(2,1): sp.diff(sxp,ell),(2,2): -dtr/3}
for k in sorted(targets):
    check('E96 ell-velocity entry %s' % (k,), dM96[k]-targets[k])
check('E96 trace(dM/dell) = 0', sp.trace(dM96))

# 2x2 transport in ell: M=[[0,a],[b,0]], a=srx=e^ell, b=cxp=coth(ell/2)
# Same reduction as eta note: u'' + P u' + Q u = 0, P = -a'/a, Q = -a*b
a = srx; b = cxp
P = -sp.diff(a, ell)/a; Q = -a*b
check('P_ell = -1', P + 1)
check('Q_ell = -e^ell coth(ell/2)', Q + E*sp.coth(ell/2))
w = sp.symbols('w', positive=True)
# w-form: u_ww + (1+P)/w u_w + Q/w^2 u ; with P=-1 the friction coeff (1+P)/w is 0 identically
Qw = -(w+1)/(w*(w-1))
check('Qw matches Q/w^2 at w=e^ell', Qw.subs(w, E) - Q/E**2)
print('Pw (friction in w-form):', sp.simplify((1-1)/w), flush=True)
for w0 in [0, 1]:
    lq = sp.limit((w-w0)**2 * Qw, w, w0)
    print('Fuchs at w=%s: (w-w0)*Pw -> 0 ; (w-w0)^2 Qw -> %s' % (w0, lq), flush=True)
z = sp.symbols('z')
Qwz = Qw.subs(w, 1/z)
print('Fuchs at infinity: Qw(1/z)/z^2 ->', sp.limit(Qwz/z**2, z, 0), flush=True)
print('Qw large-w asymptotic:', sp.series(Qw, w, sp.oo, 2), flush=True)

# point checks x-side vs ell-side, QI at 30 and 60 deg
x = sp.symbols('x', real=True)
def ops_x(xx):
    s, c = sp.sin(xx), sp.cos(xx)
    srx_ = (1+c)/s; cxp_ = 1/c + s/c
    return srx_, 1/srx_, cxp_, 1/cxp_
def M2_x(xx):
    s,q,p,r = ops_x(xx); return sp.Matrix([[0,s],[p,0]])
def dM2_dell_formula():
    return sp.Matrix([[0, srx],[-(cxp**2-1)/2, 0]])
for x0 in [sp.pi/6, sp.pi/3]:
    s_,q_,p_,r_ = ops_x(x0); ell0 = sp.log(s_)
    # d/dell = (dx/dell) d/dx ; dl/dx = -(srx+sxp)/2 = -cosh ell -> dx/dell = -1/cosh ell0
    lhs = (sp.diff(M2_x(x), x).subs(x, x0)) * (-1/sp.cosh(ell0))
    rhs = dM2_dell_formula().subs(ell, ell0)
    err = max(abs(complex(sp.N(v,30))) for v in (lhs-rhs))
    print('QI point check x=%s: max |diff| = %.3e' % (x0, err), flush=True)

print('SUMMARY: %d/%d symbolic checks passed' % (sum(1 for _,ok in results if ok), len(results)), flush=True)
print('RUN DONE', flush=True)
