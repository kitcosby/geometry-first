from sage.all import *
# Least Phase Action: delta ∫ (urx+uxp) dx = 0.
# If L = urx(x)+uxp(x) depends only on x (not x'), Euler-Lagrange gives dL/dx=0.
# Compute d(urx+uxp)/dx and find stationary points.

x = var('x')
assume(x>0, x<pi/2)
srx = 1/sin(x) + cos(x)/sin(x)
cxp = 1/cos(x) + sin(x)/cos(x)
sxp = 1/srx
crx = 1/cxp
urx = srx - crx
uxp = cxp - sxp
L = urx + uxp

dL = diff(L, x).simplify_trig()
print(f"d(urx+uxp)/dx = {dL}")
print()
# Find zeros in (0,pi/2).
# Solve dL=0.
sols = solve(dL==0, x)
print(f"Stationary points in (0,pi/2): {sols}")
print()
# But this is 1D. For Snell's law we need a PATH variational problem.
# Snell's law: n1 sin(theta1) = n2 sin(theta2).
# Fermat: minimize T = n1*L1/c + n2*L2/c over interface crossing point.
# Our analog: minimize Phase = ∫ (urx+uxp) d? over ?.
# The '?' (path parameter, geometry) is not specified.
# Without the geometric setup (what is varied? what is the path?),
# δ∫(urx+uxp)dx=0 is just dL/dx=0, which gives critical points, not Snell's law.
print("CONCLUSION:")
print("δ∫(urx+uxp)dx=0 with L=L(x) gives dL/dx=0 (critical points of L).")
print("This is NOT Snell's law; Snell's law needs a 2D path variational problem")
print("(vary the interface crossing point, like Fermat).")
print("To test against refraction, we need:")
print("  (a) the geometric setup (what path? what is varied?),")
print("  (b) how (urx+uxp) couples to the path (is it the 'refractive index'?).")
print("The guess L=urx+uxp is the Lagrangian; the variational setup is still OPEN.")

with open('/root/phase_action.txt','w') as f:
    f.write("Lagrangian_guess_OK variational_setup_OPEN\n")
print("WROTE /root/phase_action.txt")
