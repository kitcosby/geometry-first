import sympy as sp
w, t = sp.symbols('w t')
I = sp.I
z = w**2  # zeta = e^{ix}, w = e^{ix/2}
# srx = i (z + z^-1 + 2 alpha)/(z - z^-1); alpha = sgn(sin x) = sgn(Im zeta).
for alpha, label in [(1, 'alpha=+1'), (-1, 'alpha=-1')]:
    expr = I*(z + 1/z + 2*alpha)/(z - 1/z)
    expr = sp.cancel(expr)
    # target: alpha=+1 -> i (w+w^-1)/(w-w^-1) = cot(x/2); alpha=-1 -> i (w-w^-1)/(w+w^-1) = -tan(x/2)
    tgt = I*(w + 1/w)/(w - 1/w) if alpha == 1 else I*(w - 1/w)/(w + 1/w)
    print(f'{label}: srx(zeta) == target half-angle form: {sp.simplify(expr - sp.cancel(tgt)) == 0}   srx = {expr}')
    num = sp.expand(z + 1/z + 2*alpha)
    sq = (w + 1/w)**2 if alpha == 1 else (w - 1/w)**2
    print(f'   numerator z+z^-1{ "+" if alpha==1 else "-" }2 is a perfect square in w: {sp.expand(sq - (w**2 + 2*alpha + w**-2)) == 0}')
# pi-shift: x -> x+pi sends zeta -> -zeta, w -> i w, and t = tan(x/2) -> -1/t (Mobius involution).
# Check entry maps on QI forms under t -> -1/t vs QIII forms.
def q(alpha, beta, tt):
    srx = sp.cancel((alpha*(1+tt**2) + (1-tt**2))/(2*tt))
    cxp = sp.cancel((beta*(1+tt**2) + 2*tt)/(1-tt**2))
    sxp, crx = sp.cancel(1/srx), sp.cancel(1/cxp)
    return dict(srx=srx, sxp=sxp, cxp=cxp, crx=crx, urx=sp.cancel(srx-crx), uxp=sp.cancel(cxp-sxp))
QI = q(1, 1, t); QIII = q(-1, -1, t)
sub = lambda e: sp.cancel(e.subs(t, -1/t))
print('Under t -> -1/t (the pi-shift):')
print('  srx_QI -> srx_QIII form (-t):', sp.simplify(sub(QI['srx']) - QIII['srx']) == 0)
print('  urx_QI -> uxp_QIII:', sp.simplify(sub(QI['urx']) - QIII['uxp']) == 0)
print('  uxp_QI -> urx_QIII:', sp.simplify(sub(QI['uxp']) - QIII['urx']) == 0)
print('  sxp_QI -> crx_QIII:', sp.simplify(sub(QI['sxp']) - QIII['crx']) == 0)
print('  cxp_QI -> cxp_QIII (pi-periodic, stays):', sp.simplify(sub(QI['cxp']) - QIII['cxp']) == 0)
