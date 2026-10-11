from sage.all import *
# E36: Finer calibration using full (z, mu) shape.
# For candidate x_now, compute x_emit for each z via 1+z=uxp_now/uxp_emit.
# Then check correlation between observed mu and (x_now - x_emit).
# Best x_now maximizes the monotonic relationship.

import math, csv

data = []
with open('/root/Updated_Model_Dataset.csv') as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            z = float(row['Redshift (z)'])
            mu = float(row['Observed Distance Modulus (w/ Uncertainty)'])
            data.append((z, mu))
        except: pass

def uxp_deg(deg):
    r=math.radians(deg)
    s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    return cxp - 1/srx

def find_x_emit(z, x_now, uxp_now):
    target = uxp_now/(1+z)
    # uxp is increasing in QI, so invert by scanning
    best=None; best_err=1e9
    for d in [0.5 + i*0.5 for i in range(int((x_now-0.5)*2))]:
        u = uxp_deg(d)
        if u is None: continue
        err=abs(u-target)
        if err<best_err: best_err=err; best=d
    return best

def correlation(x_now):
    uxp_now = uxp_deg(x_now)
    if uxp_now is None: return -1
    xs=[]; mus=[]
    for z, mu in data:
        xe = find_x_emit(z, x_now, uxp_now)
        if xe is None: continue
        # Distance proxy: (x_now - xe) -- further back = larger
        # mu should increase with (x_now - xe)
        xs.append(x_now - xe)
        mus.append(mu)
    if len(xs)<10: return -1
    # Pearson correlation
    n=len(xs)
    mx=sum(xs)/n; mm=sum(mus)/n
    num=sum((x-mx)*(m-mm) for x,m in zip(xs,mus))
    den=math.sqrt(sum((x-mx)**2 for x in xs)*sum((m-mm)**2 for m in mus))
    if den==0: return -1
    return num/den

print("=== Finer calibration: correlation(mu, x_now - x_emit) ===")
print("x_now | correlation")
best=None; best_c=-1
for x_now in [60, 62, 65, 67.5, 70, 75, 80, 85]:
    c = correlation(x_now)
    print(f"{x_now:5.1f} | {c:6.4f}")
    if c>best_c: best_c=c; best=x_now

print()
print(f"Best x_now by correlation: {best} deg (r={best_c:.4f})")
print()
print("This uses the FULL (z,mu) shape, not just z_max.")
print("Assumes mu increases monotonically with (x_now - x_emit).")

with open('/root/fine_cal.txt','w') as f:
    f.write(f"best_x_now={best}\ncorrelation={best_c}\n")
print("WROTE /root/fine_cal.txt")
