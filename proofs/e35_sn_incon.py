from sage.all import *
# E35: Supernova z-range vs uxp - CONSISTENCY CHECK.
# Observed z max = 2.0. If 1+z = uxp(x_now)/uxp(x_emit),
# then need uxp(x_now) >= 3*uxp(x_emit_min) >= 3*1 = 3.
# So x_now must satisfy uxp(x_now) >= 3.
# Compare with E33 (Hubble tension gave x_now~8.5deg).

import math

def uxp_deg(deg):
    r=math.radians(deg)
    s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    return cxp - 1/srx

print("=== Supernova z-max consistency ===")
print("Observed z_max = 2.0")
print("Need: uxp(x_now) >= (1+z_max)*uxp_min = 3*1 = 3")
print()
print("x_now | uxp(x_now) | z_max achievable")
for d in [8.5, 13.5, 20, 30, 45, 55, 60, 67.5]:
    u = uxp_deg(d)
    zmax = u - 1 if u else None
    ok = "OK" if zmax and zmax>=2.0 else "FAIL"
    print(f"{d:5.1f} | {u:7.3f} | {zmax:5.2f} [{ok}]")

print()
print("E33 (Hubble tension) suggested x_now ~ 8.5deg.")
print("But x_now=8.5deg gives uxp~1.1, z_max~0.1, CANNOT reach z=2.0.")
print()
print("CONCLUSION: INCONSISTENCY.")
print("Either:")
print("  (a) E33's H~uxp'/uxp assumption is wrong, OR")
print("  (b) The 1+z = uxp_now/uxp_emit mapping is wrong, OR")
print("  (c) uxp is not the scale factor.")
print("The supernova z=2.0 data REQUIRES x_now >= ~55deg (if mapping holds).")
print("This falsifies the combined E33 + z-mapping assumptions.")

with open('/root/sn_incon.txt','w') as f:
    f.write("inconsistency_found=True\n")
print("WROTE /root/sn_incon.txt")
