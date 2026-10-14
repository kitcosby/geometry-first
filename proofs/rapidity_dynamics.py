# Rapidity dynamics of the operator families (QI symbolic + point checks).
# eta = ln cxp (inverse Gudermannian). All identities verified by
# rewriting to exp(eta) rational functions.
import sympy as sp

eta = sp.symbols('eta', real=True)
E = sp.exp(eta)
results = []

def is0(e):
    e = sp.together(e.rewrite(sp.exp))
    e = sp.cancel(e)
    return sp.simplify(e) == 0

def check(name, expr):
    ok = is0(expr)
    results.append((name, ok))
    print(('PASS' if ok else 'FAIL'), name, flush=True)

# ---- Part A: the dictionary in eta (QI) ----
sech = 1 / sp.cosh(eta)
tanh = sp.tanh(eta)
srx = sp.coth(eta / 2)          # = (1 + cos x)/sin x
sxp = sp.tanh(eta / 2)          # = (1 - cos x)/sin x
cxp = E
crx = sp.exp(-eta)
check('srx = (1+sech)/tanh', srx - (1 + sech) / tanh)
check('sxp = (1-sech)/tanh', sxp - (1 - sech) / tanh)
check('cosh eta = (cxp+crx)/2', sp.cosh(eta) - (cxp + crx) / 2)
check('sinh eta = (cxp-crx)/2', sp.sinh(eta) - (cxp - crx) / 2)
urx = srx - crx
uxp = cxp - sxp
sin2 = 2 * tanh * sech                      # sin 2x
cos2 = sech**2 - tanh**2                    # cos 2x
check('sin2 via 2 sin cos', sin2 - 2 * tanh * (sech))
check('cos2 = 1 - 2 tanh^2', cos2 - (1 - 2 * tanh**2))

# ---- Part B: velocities d/deta, targets in operator form ----
check('d srx = -(srx^2-1)/2', sp.diff(srx, eta) + (srx**2 - 1) / 2)
check('d sxp = (1-sxp^2)/2', sp.diff(sxp, eta) - (1 - sxp**2) / 2)
check('d cxp = cxp', sp.diff(cxp, eta) - cxp)
check('d crx = -crx', sp.diff(crx, eta) + crx)
check('d urx = -(srx^2-1)/2 + crx', sp.diff(urx, eta) + (srx**2 - 1) / 2 - crx)
check('d uxp = cxp - (1-sxp^2)/2', sp.diff(uxp, eta) - cxp + (1 - sxp**2) / 2)
check('d sin2 = 4 cos2/(cxp+crx)', sp.diff(sin2, eta) - 4 * cos2 / (cxp + crx))
check('d cos2 = -4 sin2/(cxp+crx)', sp.diff(cos2, eta) + 4 * sin2 / (cxp + crx))
# chain-rule consistency: d/deta = (dx/deta) d/dx with deta/dx = cosh eta,
# and d(sin2x)/dx = 2 cos2x  =>  d sin2/deta = 2 cos2 / cosh eta
check('d sin2 = 2 cos2/cosh eta', sp.diff(sin2, eta) - 2 * cos2 / sp.cosh(eta))

# ---- Part C: E96 family in eta, traceless-corrected, and its velocity ----
def build96():
    M = sp.Matrix([[urx, srx, cos2],
                   [cxp, uxp, crx],
                   [sin2, sxp, sp.Integer(1)]])
    return M - (sp.trace(M) / 3) * sp.eye(3)

M96 = build96()
dM96 = sp.diff(M96, eta)
tr = urx + uxp + 1
dtr = sp.diff(tr, eta)
targets = {
    (0, 0): sp.diff(urx, eta) - dtr / 3,
    (0, 1): sp.diff(srx, eta),
    (0, 2): sp.diff(cos2, eta),
    (1, 0): sp.diff(cxp, eta),
    (1, 1): sp.diff(uxp, eta) - dtr / 3,
    (1, 2): sp.diff(crx, eta),
    (2, 0): sp.diff(sin2, eta),
    (2, 1): sp.diff(sxp, eta),
    (2, 2): -dtr / 3,
}
for k in sorted(targets):
    check('E96 velocity entry %s' % (k,), dM96[k] - targets[k])
check('E96 trace(dM/deta) = 0', sp.trace(dM96))

# ---- Part D: 2x2 transport ODE, w = e^eta form, singular points ----
a = srx                     # coth(eta/2)
b = cxp                     # e^eta
P = -sp.diff(a, eta) / a    # ODE u'' + P u' + Q u = 0 for the (1,1) channel
Q = -a * b
check('P = 1/sinh eta', P - 1 / sp.sinh(eta))
w = sp.symbols('w', positive=True)
Pw = (1 + 2 * w / (w**2 - 1)) / w        # coeff of u_w after w-substitution
Qw = -(w + 1) / (w * (w - 1))            # coeff of u after w-substitution
check('Pw matches (1+P)/w at w=e^eta', Pw.subs(w, E) - (1 + P) / E)
check('Qw matches Q/w^2 at w=e^eta', Qw.subs(w, E) - Q / E**2)
print('denominator Pw:', sp.factor(sp.denom(sp.together(Pw))), flush=True)
print('denominator Qw:', sp.factor(sp.denom(sp.together(Qw))), flush=True)
for w0 in [0, 1, -1]:
    lp = sp.limit((w - w0) * Pw, w, w0)
    lq = sp.limit((w - w0)**2 * Qw, w, w0)
    print('Fuchs at w=%s: (w-w0)Pw -> %s ; (w-w0)^2 Qw -> %s'
          % (w0, lp, lq), flush=True)
z = sp.symbols('z')
Pwz = Pw.subs(w, 1 / z)
Qwz = Qw.subs(w, 1 / z)
lp_inf = sp.limit(2 - Pwz / z, z, 0)
lq_inf = sp.limit(Qwz / z**2, z, 0)
print('Fuchs at infinity: 2 - Pw(1/z)/z -> %s ; Qw(1/z)/z^2 -> %s'
      % (lp_inf, lq_inf), flush=True)
# numerator of Pw: where the friction coefficient vanishes
numPw = sp.factor(sp.together(Pw).as_numer_denom()[0])
print('numerator Pw:', numPw, flush=True)
check('sinh eta = -1 at w = sqrt2-1', (w - 1 / w) / 2 + 1)

# ---- Part E: point cross-checks, x-side vs eta-side ----
x = sp.symbols('x', real=True)
def ops_x_QI(xx):
    s, c = sp.sin(xx), sp.cos(xx)
    srx_ = (1 + c) / s
    cxp_ = 1 / c + s / c
    return srx_, 1 / srx_, cxp_, 1 / cxp_

def M2_x(xx):
    s, q, p, r = ops_x_QI(xx)
    return sp.Matrix([[0, s], [p, 0]])

def dM2_deta_formula():
    return sp.Matrix([[0, -(srx**2 - 1) / 2], [cxp, 0]])

for x0 in [sp.pi / 6, sp.pi / 3]:
    s_, q_, p_, r_ = ops_x_QI(x0)
    eta0 = sp.log(p_)
    lhs = (sp.diff(M2_x(x), x).subs(x, x0)) / (1 / sp.cos(x0))
    rhs = dM2_deta_formula().subs(eta, eta0)
    err = max(abs(complex(sp.N(v, 30))) for v in (lhs - rhs))
    print('QI point check x=%s: max |diff| = %.3e' % (x0, err), flush=True)

# QII signed dictionary at x0 = 2pi/3 (M2 only)
x0 = 2 * sp.pi / 3
s0, c0 = sp.sin(x0), sp.cos(x0)
cxp0 = abs(1 / c0) + s0 / c0
eta0 = sp.log(cxp0)
print('QII: sinh(eta0) - tan(x0) =', sp.simplify(sp.sinh(eta0) - sp.tan(x0)),
      flush=True)
print('QII: tanh(eta0) + sin(x0) =', sp.simplify(sp.tanh(eta0) + s0),
      ' (signed dictionary: tanh eta = -sin x when cos x < 0)', flush=True)
srx0 = abs(1 / s0) + c0 / s0
srx_eta_QII = -1 / sp.tanh(eta) + 1 / sp.sinh(eta)   # eta < 0 branch form
print('QII: srx(x0) - srx_eta(eta0) =',
      sp.simplify(srx0 - srx_eta_QII.subs(eta, eta0)), flush=True)
dM2dx = sp.diff(M2_x(x), x)
# M2_x uses QI branch formulas; at QII the (0,1) entry differs, so build the
# QII sample matrix directly from branch-correct values as functions of x:
srx_x = abs(1 / sp.sin(x)) + sp.cos(x) / sp.sin(x)
cxp_x = abs(1 / sp.cos(x)) + sp.sin(x) / sp.cos(x)
M2g = sp.Matrix([[0, srx_x], [cxp_x, 0]])
lhs = (sp.diff(M2g, x).subs(x, x0)) / abs(1 / c0)
dsrx_eta = sp.diff(srx_eta_QII, eta)
rhs = sp.Matrix([[0, dsrx_eta.subs(eta, eta0)], [cxp0, 0]])
err = max(abs(complex(sp.N(v, 30))) for v in (lhs - rhs))
print('QII point check x=2pi/3: max |diff| = %.3e' % err, flush=True)

print('SUMMARY: %d/%d symbolic checks passed'
      % (sum(1 for _, ok in results if ok), len(results)), flush=True)
print('RUN DONE', flush=True)
