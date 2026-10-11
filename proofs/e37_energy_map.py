from sage.all import *
# E37: Energy-based redshift mapping from photon-to-space mechanism.
# 1+z = urx(x_emit)/urx(x_now)  [energy ratio, not scale ratio]
# urx decreases in QI: inf -> 2 -> 1 as x: 0 -> 45 -> 90 deg.
# For each z, find x_emit. Check consistency and correlation.

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

def urx_deg(deg):
    r=math.radians(deg)
    s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    crx=1/cxp
    return srx - crx

def find_x_emit_energy(z, x_now, urx_now):
    # Solve urx(x_emit) = (1+z)*urx_now for x_emit in (0, x_now)
    target = (1+z)*urx_now
    best=None; best_err=1e9
    # urx decreases with x, so x_emit < x_now. Scan from 0.1 to x_now.
    for d in [0.1 + i*0.2 for i in range(int((x_now-0.1)*5))]:
        u = urx_deg(d)
        if u is None: continue
        err=abs(u-target)
        if err<best_err: best_err=err; best=d
    # Check if best is reasonable (err small relative to target)
    if best_err > 0.1*target:
        return None
    return best

def correlation_energy(x_now):
    urx_now = urx_deg(x_now)
    if urx_now is None: return -1, 0
    xs=[]; mus=[]; count=0
    for z, mu in data:
        xe = find_x_emit_energy(z, x_now, urx_now)
        if xe is None: continue
        count+=1
        xs.append(x_now - xe)
        mus.append(mu)
    if len(xs)<10: return -1, count
    n=len(xs)
    mx=sum(xs)/n; mm=sum(mus)/n
    num=sum((x-mx)*(m-mm) for x,m in zip(xs,mus))
    den=math.sqrt(sum((x-mx)**2 for x in xs)*sum((m-mm)**2 for m in mus))
    if den==0: return -1, count
    return num/den, count

print("=== Energy-based mapping: 1+z = urx_emit/urx_now ===")
print("x_now | urx_now | correlation | N_fit")
best=None; best_c=-1
for x_now in [8.5, 13.5, 20, 30, 45, 60, 67.5, 75]:
    urx_now = urx_deg(x_now)
    c, n = correlation_energy(x_now)
    urx_str = f"{urx_now:.3f}" if urx_now else "None"
    print(f"{x_now:5.1f} | {urx_str:7s} | {c:6.4f} | {n:3d}")
    if c>best_c: best_c=c; best=x_now

print()
print(f"Best x_now (energy mapping): {best} deg (r={best_c:.4f})")
print()
print("Compare with scale-based (E36): x~75deg.")
print("If energy-based gives different x, the mapping matters.")

with open('/root/energy_map.txt','w') as f:
    f.write(f"best_x_now={best}\n")
print("WROTE /root/energy_map.txt")
