from sage.all import *
# Second law: when is dS/dx >= 0?
# dS/dx = (srx*|csc x| - crx*|sec x|) * S^2.
# Since S^2 >= 0, need: srx*|csc x| - crx*|sec x| >= 0.
# With: srx=|csc|+cot, crx=1/(|sec|+tan).

x = var('x')
# Work in (0,pi/2) where all trig positive, |.| drops.
# But we need the sign chart across all x. Let's do it piecewise.
# Actually, let's use the exact Q7 formula and analyze sign.

# srx = csc+cot (for x in (0,pi), sin>0)
# crx = 1/(sec+tan)
# For x in (0,pi/2): sec>0, tan>0.
# srx*|csc| = (csc+cot)*csc = csc^2 + csc*cot
# crx*|sec| = [1/(sec+tan)]*sec = sec/(sec+tan)

# Define f(x) = srx*csc - crx*sec (for x in (0,pi/2), no abs needed)
assume(x>0, x<pi/2)
srx = 1/sin(x) + cos(x)/sin(x)
crx = 1/(1/cos(x) + sin(x)/cos(x))
f = srx*(1/sin(x)) - crx*(1/cos(x))
f_simp = f.simplify_trig()
print(f"f(x) = srx*csc - crx*sec = {f_simp}")
print()
# Find where f>=0 in (0,pi/2).
# Solve f(x)=0.
sols = solve(f_simp==0, x)
print(f"Zeros in (0,pi/2): {sols}")
# Sample sign.
for deg in [15,30,45,60,75]:
    r = RDF(deg*pi/180)
    val = float(f_simp.subs(x=r))
    print(f"  x={deg}deg: f={val:.4f} {'(dS/dx>=0)' if val>=0 else '(dS/dx<0)'}")
print()
print("If f<0 somewhere, entropy DECREASES there (second law violated in x).")
print("But x is phase, not time. The physical second law needs t, not x.")
print("This is the GAP: dS/dx is PROVED, but dS/dt needs the x->t map (Gamma).")

with open('/root/second_law.txt','w') as f_out:
    f_out.write("second_law_x_analyzed\n")
print("WROTE /root/second_law.txt")
