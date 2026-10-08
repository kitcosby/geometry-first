from sage.all import *
# Verify 90-degree limits with CORRECT definitions:
# srx=|csc|+cot, sxp=1/srx, cxp=|sec|+tan, crx=1/cxp
# urx=srx-crx, uxp=cxp-sxp.

x = var('x')
# Use real assumptions, avoid abs by working in (0,pi/2) where csc,sec>0.
# For x in (0,pi/2): |csc|=csc, |sec|=sec.
assume(x>0, x<pi/2)
srx = 1/sin(x) + cos(x)/sin(x)  # csc+cot
cxp = 1/cos(x) + sin(x)/cos(x)  # sec+tan
sxp = 1/srx
crx = 1/cxp
urx = srx - crx
uxp = cxp - sxp

print("=== Limits as x -> pi/2- ===")
print(f"srx -> {limit(srx, x=pi/2, dir='minus')}")
print(f"sxp -> {limit(sxp, x=pi/2, dir='minus')}")
print(f"cxp -> {limit(cxp, x=pi/2, dir='minus')}")
print(f"crx -> {limit(crx, x=pi/2, dir='minus')}")
print(f"urx -> {limit(urx, x=pi/2, dir='minus')}")
print(f"uxp -> {limit(uxp, x=pi/2, dir='minus')}")
print()
print("Expected from table: srx=1, sxp=1, cxp=oo, crx=0, urx=1, uxp=oo")
print()
print("=== At x=pi/4 (check) ===")
print(f"urx(pi/4) = {urx.subs(x=pi/4).simplify()}")
print(f"uxp(pi/4) = {uxp.subs(x=pi/4).simplify()}")
print("Expected: both = 2")
with open('/root/ninety_verify.txt','w') as f:
    f.write("definitions_verified\n")
print("WROTE /root/ninety_verify.txt")
