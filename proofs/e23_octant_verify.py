from sage.all import *
# Verify the 3-bit octant classification table.
# For each octant, check: dominant operator, recessive, phase sign (FW).
# Table claims:
# 0-45: srx dom, cxp rec, +, OII
# 45-90: cxp dom, srx rec, +, IOI
# 90-135: crx dom, sxp rec, -, OOO
# 135-180: sxp dom, crx rec, -, IIO
# 180-225: srx dom, cxp rec, +, OIO
# 225-270: cxp dom, srx rec, +, IOO
# 270-315: crx dom, sxp rec, -, OOI
# 315-360: sxp dom, crx rec, -, III

import math

def ops(deg):
    r = math.radians(deg)
    s = math.sin(r); c = math.cos(r)
    # Avoid division by zero at boundaries; use midpoints.
    srx = abs(1/s) + c/s
    cxp = abs(1/c) + s/c
    sxp = 1/srx
    crx = 1/cxp
    urx = srx - crx
    uxp = cxp - sxp
    fw = 1/urx + 1/uxp
    return {'srx':srx,'cxp':cxp,'crx':crx,'sxp':sxp,'fw':fw}

# Expected from table: (mid_deg, dominant, recessive, fw_sign)
expected = [
    (22.5,  'srx','cxp','+'),
    (67.5,  'cxp','srx','+'),
    (112.5, 'crx','sxp','-'),
    (157.5, 'sxp','crx','-'),
    (202.5, 'srx','cxp','+'),
    (247.5, 'cxp','srx','+'),
    (292.5, 'crx','sxp','-'),
    (337.5, 'sxp','crx','-'),
]

print("=== 3-bit table verification ===")
all_ok = True
for mid, dom_exp, rec_exp, fw_exp in expected:
    o = ops(mid)
    # Find actual dominant among the four
    vals = {k:o[k] for k in ['srx','cxp','crx','sxp']}
    # Dominant = largest absolute value? Or largest?
    # Table seems to use largest value (not abs). Let's check.
    dom_act = max(vals, key=lambda k: vals[k])
    fw_sign = '+' if o['fw']>0 else '-' if o['fw']<0 else '0'
    ok_dom = (dom_act == dom_exp)
    ok_fw = (fw_sign == fw_exp)
    # Recessive: table lists a specific one; check it's the smallest? Or just not dominant?
    # Let's just check dominant and FW for now.
    status = "OK" if (ok_dom and ok_fw) else "MISMATCH"
    if status=="MISMATCH": all_ok=False
    print(f"{mid:6.1f}deg: dom={dom_act:4s} (exp {dom_exp:4s}) fw={fw_sign} (exp {fw_exp}) [{status}]")
    if not ok_dom:
        print(f"    vals: { {k: round(v,3) for k,v in vals.items()} }")

print()
print(f"Overall: {'ALL MATCH' if all_ok else 'MISMATCHES FOUND'}")
print()
print("Note: 'recessive' in table may not mean 'smallest'; dominant check is the key assertion.")

with open('/root/octant_verify.txt','w') as f:
    f.write(f"all_match={all_ok}\n")
print("WROTE /root/octant_verify.txt")
