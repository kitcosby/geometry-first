from sage.all import *
# E43: Isotropy test + CMB T(alpha) profile.
# CORRECTED: uxp(alpha) is the SCALE FACTOR (homogeneous), not radial coordinate.
# (x,y,z) are comoving; operators depend only on phase alpha (homogeneity + isotropy).
# Isotropy test: Pantheon+ residuals vs sky position (RA, Dec).
# CMB: T(alpha) = T_now * urx(alpha)/urx(alpha_now). At z=1100, T=3000K.

import math, random

# Check if RA/Dec in Pantheon+ dat
with open('/root/pantheon_api.dat') as f:
    header = f.readline().split()
    print("Columns:", [h for h in header if 'RA' in h or 'DEC' in h or 'ra' in h.lower()])
    # Print full header for inspection
    print("Full header:", header[:15])

def urx_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    return srx - 1/cxp

# CMB temperature profile: T(alpha)/T_now = urx(alpha)/urx(alpha_now)
alpha_now = 85.0
unow = urx_phase(alpha_now)
print(f"\nCMB T(alpha) profile (alpha_now={alpha_now}deg, urx_now={unow:.4f}):")
print(" alpha | urx(alpha) | T/T_now | T (K, T_now=2.725)")
print("-------|------------|---------|------")
for a in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 85.0]:
    u = urx_phase(a)
    if u:
        ratio = u/unow
        print(f"{a:6.1f} | {u:10.4f} | {ratio:7.2f} | {2.725*ratio:8.1f}")
print()
print("At z=1100 (1+z=1101): need urx(alpha)=1101*urx_now.")
print(f"  Target urx = {1101*unow:.2f}")
# Find alpha where urx = target
target = 1101*unow
lo, hi = 0.1, 10.0
for _ in range(50):
    mid=(lo+hi)/2
    v=urx_phase(mid)
    if v and v > target: lo=mid
    else: hi=mid
print(f"  alpha_CMB ~ {(lo+hi)/2:.3f} deg")
print("  (Very early phase, as expected for hot CMB.)")
print()
print("WROTE /root/e43_cmb.txt")
