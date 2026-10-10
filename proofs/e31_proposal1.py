from sage.all import *
# E31: Test proposal one: T(x) = cos(phi)*E0 + sin(phi)*E4, where phi=arg(U).
# Question: does phi(x) select E0 (cos-dominated) when FW>0, E4 when FW<0?
# Compute phi(x), cos(phi), sin(phi), FW sign across all quadrants.

import math

def ops_full(deg):
    # Full definitions with abs, valid in all quadrants
    r = math.radians(deg)
    s = math.sin(r); c = math.cos(r)
    if abs(s) < 1e-12 or abs(c) < 1e-12:
        return None
    srx = abs(1/s) + c/s
    cxp = abs(1/c) + s/c
    sxp = 1/srx if abs(srx)>1e-12 else float('inf')
    crx = 1/cxp if abs(cxp)>1e-12 else float('inf')
    urx = srx - crx
    uxp = cxp - sxp
    if abs(urx)<1e-12 or abs(uxp)<1e-12:
        return None
    fw = 1/urx + 1/uxp
    # phi = atan2(uxp, urx), in degrees, normalized to [0,360)
    phi = math.degrees(math.atan2(uxp, urx)) % 360
    return {'urx':urx,'uxp':uxp,'fw':fw,'phi':phi,
            'cos_phi':math.cos(math.radians(phi)),
            'sin_phi':math.sin(math.radians(phi))}

print("=== Proposal One Test: phi(x) vs FW sign ===")
print("deg | phi | cos(phi) | sin(phi) | FW | E0-dom? | FW>0?")
print("If proposal works: cos>sin (E0-dom) iff FW>0.")
print()
mismatch = 0
total = 0
for deg in range(5, 360, 10):
    o = ops_full(deg)
    if o is None: continue
    total += 1
    e0_dom = o['cos_phi'] > o['sin_phi']
    # Note: cos_phi, sin_phi can be negative; "E0-dom" means |cos|>|sin| and cos>0?
    # Simpler: check if cos_phi > 0 and cos_phi > abs(sin_phi)?
    # Actually proposal was T=cos*E0+sin*E4. "Select E0" means cos dominates.
    # Let's use: cos_phi > sin_phi (as real numbers).
    fw_pos = o['fw'] > 0
    match = (e0_dom == fw_pos)
    if not match: mismatch += 1
    flag = "OK" if match else "MISMATCH"
    if not match or deg % 45 == 5:
        print(f"{deg:3d} | {o['phi']:6.1f} | {o['cos_phi']:7.3f} | {o['sin_phi']:7.3f} | {'+' if fw_pos else '-'} | {str(e0_dom):5s} | {str(fw_pos):5s} [{flag}]")

print()
print(f"Mismatches: {mismatch}/{total}")
if mismatch == 0:
    print("PROPOSAL CONFIRMED: phi selects E0 iff FW>0.")
else:
    print("PROPOSAL FALSIFIED (or needs refinement): phi does not cleanly select by FW.")
    print("The coupling via phi(x) does not implement the E0/FW>0 rule.")

with open('/root/proposal1.txt','w') as f:
    f.write(f"mismatch={mismatch}\ntotal={total}\n")
print("WROTE /root/proposal1.txt")
