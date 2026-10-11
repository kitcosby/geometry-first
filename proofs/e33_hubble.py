from sage.all import *
# E33: Hubble tension calibration.
# H0 late (SH0ES) ~73, H0 early (Planck) ~67. Ratio = 73/67 = 1.0896.
# Hypothesis: H(x) proportional to (d(uxp)/dx)/uxp (like standard Hubble).
# Then H(x2)/H(x1) = [uxp'(x2)/uxp(x2)] / [uxp'(x1)/uxp(x1)].
# Solve for x_now given x_cmb (early) and ratio=1.0896.

import math

def uxp_and_deriv(deg):
    r = math.radians(deg)
    s = math.sin(r); c = math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None,None
    srx = (1+c)/s; cxp = (1+s)/c
    u = cxp - 1/srx
    # numerical derivative
    h = 1e-5
    def u_(d):
        rr=math.radians(d); ss=math.sin(rr); cc=math.cos(rr)
        srx_=(1+cc)/ss; cxp_=(1+ss)/cc
        return cxp_ - 1/srx_
    dup = (u_(deg+h)-u_(deg-h))/(2*h)  # per degree; convert later if needed
    return u, dup

def H_ratio(x1_deg, x2_deg):
    # H ~ (d(uxp)/dx)/uxp. Ratio is dimensionless, h cancels.
    u1, du1 = uxp_and_deriv(x1_deg)
    u2, du2 = uxp_and_deriv(x2_deg)
    if u1 is None or u2 is None: return None
    H1 = du1/u1
    H2 = du2/u2
    return H2/H1

target = 73/67
print(f"Target H ratio (late/early): {target:.4f}")
print()
print("Assuming x_cmb (early) is small, solve for x_now:")
print("x_cmb | x_now solving H(x_now)/H(x_cmb)=1.0896")
print()

for x_cmb in [1, 5, 10, 15]:
    # Find x_now > x_cmb such that ratio matches
    # Scan x_now from x_cmb+1 to 89
    best = None
    best_err = 1e9
    for x_now in [x_cmb + i*0.5 for i in range(1, int((89-x_cmb)*2))]:
        r = H_ratio(x_cmb, x_now)
        if r is None: continue
        err = abs(r - target)
        if err < best_err:
            best_err = err
            best = x_now
    if best:
        print(f"{x_cmb:5.1f} | {best:5.1f} (err {best_err:.4f})")
    else:
        print(f"{x_cmb:5.1f} | no solution")

print()
print("INTERPRETATION:")
print("If CMB at x~5deg, we are at x~? (see above).")
print("This is RELATIVE calibration (no absolute units).")
print("Assumes H ~ uxp'/uxp and x_cmb is 'early'.")
print("The x<->t map is still needed for absolute time.")

with open('/root/hubble.txt','w') as f:
    f.write(f"target_ratio={target}\n")
print("WROTE /root/hubble.txt")
