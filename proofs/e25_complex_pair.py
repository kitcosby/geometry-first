from sage.all import *
# E25: Complexification via reciprocal pairing.
# U_urx = urx + i*uxp, U_uxp = uxp + i*urx.
# Similarly: S_srx = srx + i*sxp, etc.
# Check: modulus, and whether |U|^2 has a simple form.

import math, cmath

def ops(deg):
    r = math.radians(deg)
    s = math.sin(r); c = math.cos(r)
    srx = abs(1/s) + c/s
    cxp = abs(1/c) + s/c
    sxp = 1/srx
    crx = 1/cxp
    urx = srx - crx
    uxp = cxp - sxp
    return srx,sxp,cxp,crx,urx,uxp

print("=== Complex pairing: U = urx + i*uxp ===")
print("deg | urx | uxp | |U|^2 | arg(U) deg")
for deg in [22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5]:
    srx,sxp,cxp,crx,urx,uxp = ops(deg)
    U = complex(urx, uxp)
    mod2 = abs(U)**2
    arg = math.degrees(cmath.phase(U))
    print(f"{deg:6.1f} | {urx:7.3f} | {uxp:7.3f} | {mod2:8.3f} | {arg:6.1f}")

print()
print("=== Four cores: S_srx = srx + i*sxp ===")
for deg in [22.5, 67.5]:
    srx,sxp,cxp,crx,urx,uxp = ops(deg)
    S = complex(srx, sxp)
    C = complex(cxp, crx)
    print(f"{deg}deg: |srx+i*sxp|^2={abs(S)**2:.3f}, |cxp+i*crx|^2={abs(C)**2:.3f}")
    print(f"  Note: srx*sxp=1, so |S|^2 = srx^2 + 1/srx^2.")

print()
print("CONCLUSION:")
print("Complexification pairs each operator with its reciprocal:")
print("  U_urx = urx + i*uxp, U_uxp = uxp + i*urx")
print("  S_srx = srx + i*sxp, C_cxp = cxp + i*crx, etc.")
print("The 'self-propagating wave' structure (real/imaginary generating each other)")
print("is posited; the wave equation itself is OPEN.")

with open('/root/complex_pair.txt','w') as f:
    f.write("complexification_defined=True\n")
print("WROTE /root/complex_pair.txt")
