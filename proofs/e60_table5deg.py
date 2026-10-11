from sage.all import *
# E60: Clean 5-degree increment operator table.
# Columns: Degree, srx, sxp, cxp, crx, uxp, urx, sawx(1/uxp), sawr(1/urx), Flatwave
# Excludes exact boundaries (0,90,180,270,360) where operators undefined.

import math

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-12 or abs(c)<1e-12:
        return None
    srx=abs(1/s)+c/s; sxp=1/srx
    cxp=abs(1/c)+s/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    sawx=1/uxp if uxp!=0 else None
    sawr=1/urx if urx!=0 else None
    # Flatwave = 1/urx + 1/uxp = sgn(sin2x) in theory; compute directly
    fw=(1/urx+1/uxp) if (urx!=0 and uxp!=0) else None
    return srx,sxp,cxp,crx,uxp,urx,sawx,sawr,fw

print("Generating 5-degree table...")
rows=[]
for deg in range(5, 360, 5):
    if deg%90==0: continue  # skip 90,180,270
    v=ops(deg)
    if v is None: continue
    rows.append((deg,)+v)

# Write CSV
with open('/root/operator_table_5deg.csv','w') as f:
    f.write("Degree,srx,sxp,cxp,crx,uxp,urx,sawx,sawr,Flatwave\n")
    for r in rows:
        f.write(f"{r[0]},{r[1]:.6f},{r[2]:.6f},{r[3]:.6f},{r[4]:.6f},{r[5]:.6f},{r[6]:.6f},{r[7]:.6f},{r[8]:.6f},{r[9]:.6f}\n")
print(f"Wrote {len(rows)} rows to /root/operator_table_5deg.csv")
# Preview first 5 and last 5
print("\nPreview (first 3):")
for r in rows[:3]:
    print(f"  {r[0]}: srx={r[1]:.4f}, uxp={r[5]:.4f}, urx={r[6]:.4f}, FW={r[9]:.4f}")
print("Preview (last 3):")
for r in rows[-3:]:
    print(f"  {r[0]}: srx={r[1]:.4f}, uxp={r[5]:.4f}, urx={r[6]:.4f}, FW={r[9]:.4f}")
print("\nWROTE /root/operator_table_5deg.csv")
