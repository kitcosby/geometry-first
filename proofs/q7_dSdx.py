from sage.all import *
import numpy as np

# Q7: Symbolic derivation of dS/dx where S = sawup = 1/urx.
# Claim (from RUN144 lemma): dS/dx = (srx*|csc x| - crx*|sec x|) * S^2.
# Verify symbolically per-octant (fix signs), then numerically.

print("=== Symbolic per-octant verification ===")
# On an octant, sigma_s = sign(sin x), sigma_c = sign(cos x) are constant.
# srx = sigma_s*(1+sigma_s*cos x)/sin x
# crx = 1/cxp, cxp = sigma_c*(1+sigma_c*sin x)/cos x
# We verify: d/dx [1/(srx - crx)] = (srx*|csc| - crx*|sec|) / (srx-crx)^2
# using l' = -|csc|, eta' = |sec|.

x = var('x')
for (ss, sc, name) in [(1,1,'Q1'), (1,-1,'Q2'), (-1,-1,'Q3'), (-1,1,'Q4')]:
    # define with fixed signs (no abs, using sigma explicitly)
    # |csc| = 1/|sin|; on octant, |sin| = ss*sin (since ss=sign(sin))
    # srx = |csc| + cot = 1/(ss*sin) + cos/sin
    srx = 1/(ss*sin(x)) + cos(x)/sin(x)
    cxp = 1/(sc*cos(x)) + sin(x)/cos(x)
    crx = 1/cxp
    urx = srx - crx
    S = 1/urx
    dS = diff(S, x)
    # claimed form
    abs_csc = 1/(ss*sin(x))  # = |csc| on this octant
    abs_sec = 1/(sc*cos(x))  # = |sec| on this octant
    claimed = (srx*abs_csc - crx*abs_sec) * S**2
    # simplify difference
    diff_expr = (dS - claimed).simplify_full()
    # check numerically at a test point in the octant
    test_x = {'Q1':0.7, 'Q2':2.0, 'Q3':3.8, 'Q4':5.5}[name]
    v1 = float(dS.subs(x=test_x))
    v2 = float(claimed.subs(x=test_x))
    match = abs(v1-v2) < 1e-10
    print(f"{name}: dS/dx={v1:.6f}, claimed={v2:.6f}, match={match}")

print()
print("=== Numerical global check (all octants) ===")
def ops_num(x):
    s=np.sin(x); c=np.cos(x)
    srx=abs(1/s)+c/s; cxp=abs(1/c)+s/c
    crx=1/cxp; urx=srx-crx
    S=1/urx
    # claimed coefficient
    coef = srx*abs(1/s) - crx*abs(1/c)
    return S, coef

xs = np.linspace(0.01, 4*np.pi-0.01, 5000)
k = np.round(xs/(np.pi/2)); valid = np.abs(xs-k*np.pi/2) > 1e-3
xs = xs[valid]
S = np.array([ops_num(x)[0] for x in xs])
coef = np.array([ops_num(x)[1] for x in xs])
dS_num = np.gradient(S, xs[1]-xs[0])
dS_claim = coef * S**2
err = np.abs(dS_num - dS_claim)
# exclude near-singularities where gradient is noisy
mask = np.abs(S) < 10
print(f"max err (non-singular): {err[mask].max():.2e}")
print(f"mean err (non-singular): {err[mask].mean():.2e}")
print(f"FORMULA CONFIRMED: {err[mask].max() < 1e-4}")

# Analyze the coefficient c(x) = srx*|csc| - crx*|sec|
# Is it expressible via S? Check correlation.
print()
print("=== Coefficient analysis ===")
print("c(x) range:", coef[mask].min(), "to", coef[mask].max())
# c(x) vs S: if Riccati with constant coeffs, c would be const. It's not.
print("c(x) std:", coef[mask].std(), "(non-constant => not simple Riccati)")

with open('/root/q7_dSdx.txt','w') as f:
    f.write(f"maxerr {err[mask].max():.2e}\nconfirmed {err[mask].max()<1e-4}\n")
print("WROTE /root/q7_dSdx.txt")
