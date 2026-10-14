# Fix-up verifications for rapidity_dynamics.py
import sympy as sp

w = sp.symbols('w')
# (a) friction zero: sinh at w = sqrt2-1 is exactly -1 (substitution, not identity)
val = ((w - 1 / w) / 2).subs(w, sp.sqrt(2) - 1)
print('sinh at w=sqrt2-1:', sp.simplify(val), '(expect -1)')

# (b) QII dictionary at x0 = 2pi/3 by exact rational arithmetic in w0 = e^eta0
w0 = 2 - sp.sqrt(3)                      # cxp(120 deg) = |sec| + tan = 2 - sqrt3
sinh0 = (w0 - 1 / w0) / 2
cosh0 = (w0 + 1 / w0) / 2
tanh0 = sp.simplify(sinh0 / cosh0)
print('sinh(eta0) =', sp.simplify(sinh0), '(expect -sqrt3 = tan 120)')
print('cosh(eta0) =', sp.simplify(cosh0), '(expect 2 = |sec 120|)')
print('tanh(eta0) =', tanh0, '(expect -sqrt3/2 = -sin 120)')
srx_eta = -1 / tanh0 + 1 / sinh0
print('srx_eta(eta0) =', sp.simplify(srx_eta), '(expect sqrt3/3 = srx 120)')

# (c) tie-point values for the screen: M2 and dM2/deta at eta* = ln(1+sqrt2)
a = 1 + sp.sqrt(2)
M2star = sp.Matrix([[0, a], [a, 0]])
dM2star = sp.Matrix([[0, -a], [a, 0]])
print('M2(eta*) eigenvalues:', M2star.eigenvals())
print('dM2/deta(eta*) eigenvalues:', dM2star.eigenvals())
# and the velocity (0,1) entry formula at eta*: -(srx^2-1)/2 with srx = 1+sqrt2
print('-(srx^2-1)/2 at tie:', sp.simplify(-(a**2 - 1) / 2), '(expect -(1+sqrt2))')
print('RUN DONE')
