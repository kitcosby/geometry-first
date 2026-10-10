from sage.all import *
# E26: Try to prove |U|^2 = urx^2 + uxp^2 is constant (or find formula).
# Work in QI (0<x<pi/2) where abs simplifies.
# Definitions:
# srx = csc+cot = (1+cos)/sin
# sxp = 1/srx
# cxp = sec+tan = (1+sin)/cos
# crx = 1/cxp
# urx = srx - crx, uxp = cxp - sxp.

x = var('x')
# QI: sin>0, cos>0 (for 0<x<pi/2)
srx = (1+cos(x))/sin(x)
sxp = 1/srx
cxp = (1+sin(x))/cos(x)
crx = 1/cxp
urx = srx - crx
uxp = cxp - sxp

mod2 = urx**2 + uxp**2
mod2_simpl = mod2.simplify_full()
print("urx^2 + uxp^2 simplified:")
print(mod2_simpl)
print()

# Try trig_simplify
mod2_trig = mod2.simplify_trig()
print("trig simplified:")
print(mod2_trig)
print()

# Is it constant? Take derivative.
dmod2 = diff(mod2, x).simplify_full()
print("d/dx (urx^2+uxp^2):")
print(dmod2)
print()
print("If derivative is zero, modulus is constant. If not, it's x-dependent.")

# Numerical check at non-midpoint (30 deg) vs midpoint (22.5)
import math
def num_mod2(deg):
    r=math.radians(deg)
    s=math.sin(r); c=math.cos(r)
    srx_=(1+c)/s; cxp_=(1+s)/c
    urx_=srx_-1/cxp_; uxp_=cxp_-1/srx_
    return urx_**2+uxp_**2

print()
print(f"|U|^2 at 22.5deg: {num_mod2(22.5):.6f}")
print(f"|U|^2 at 30deg:   {num_mod2(30):.6f}")
print(f"|U|^2 at 10deg:   {num_mod2(10):.6f}")
if abs(num_mod2(30)-num_mod2(22.5))>1e-6:
    print("NOT CONSTANT: varies with x. The 8-midpoint match was symmetry, not identity.")
else:
    print("Constant at tested points.")

with open('/root/modulus.txt','w') as f:
    f.write("checked\n")
print("WROTE /root/modulus.txt")
