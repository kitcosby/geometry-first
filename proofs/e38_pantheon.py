from sage.all import *
# E38: Re-run calibration with Pantheon+ (1701 SNe).
# File: /tmp/pantheon_api.dat, columns include zHD, MU_SH0ES.
# Test both mappings: scale-based (1+z=uxp0/uxpe) and energy-based (1+z=urxe/urx0).

import math

# Parse Pantheon+ .dat (space-separated, header row)
data = []
with open('/root/pantheon_api.dat') as f:
    header = f.readline().split()
    try:
        iz = header.index('zHD')
        imu = header.index('MU_SH0ES')
    except ValueError:
        # Fallback: try zCMB
        iz = header.index('zCMB')
        imu = header.index('MU_SH0ES')
    print(f"Using columns: z={header[iz]}, mu={header[imu]}")
    for line in f:
        parts = line.split()
        if len(parts) <= max(iz, imu): continue
        try:
            z = float(parts[iz])
            mu = float(parts[imu])
            # Filter: z>0.01 (avoid peculiar velocities), z<2.5, mu reasonable
            if 0.01 < z < 2.5 and 30 < mu < 50:
                data.append((z, mu))
        except: pass

print(f"Loaded {len(data)} Pantheon+ supernovae.")
print(f"z range: {min(d[0] for d in data):.4f} to {max(d[0] for d in data):.4f}")
print()

def uxp_deg(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    return (1+s)/c - 1/((1+c)/s)

def urx_deg(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    return srx - 1/cxp

def find_x_for_target(target_fn, target_val, x_now, fn):
    # Find x_emit in (0, x_now) such that fn(x_emit)=target_val
    best=None; best_err=1e9
    for d in [0.2 + i*0.2 for i in range(int((x_now-0.2)*5))]:
        v = fn(d)
        if v is None: continue
        err=abs(v-target_val)
        if err<best_err: best_err=err; best=d
    # Require reasonable match
    if best_err > 0.15*abs(target_val):
        return None
    return best

def corr_for_xnow(x_now, mode):
    # mode: 'scale' (1+z=uxp0/uxpe) or 'energy' (1+z=urxe/urx0)
    if mode=='scale':
        base_now = uxp_deg(x_now)
        fn = uxp_deg
        # target: uxp_emit = uxp_now/(1+z)
        def get_target(z): return base_now/(1+z)
    else:
        base_now = urx_deg(x_now)
        fn = urx_deg
        def get_target(z): return base_now*(1+z)
    if base_now is None: return -1, 0
    xs=[]; mus=[]
    for z, mu in data:
        xe = find_x_for_target(None, get_target(z), x_now, fn)
        if xe is None: continue
        xs.append(x_now - xe); mus.append(mu)
    if len(xs)<50: return -1, len(xs)
    n=len(xs); mx=sum(xs)/n; mm=sum(mus)/n
    num=sum((x-mx)*(m-mm) for x,m in zip(xs,mus))
    den=math.sqrt(sum((x-mx)**2 for x in xs)*sum((m-mm)**2 for m in mus))
    return (num/den if den!=0 else -1), len(xs)

print("=== Scale-based mapping (1+z = uxp0/uxpe) ===")
print("x_now | corr | N")
best_s=None; best_cs=-1
for x_now in [60, 65, 70, 75, 80, 85]:
    c, n = corr_for_xnow(x_now, 'scale')
    print(f"{x_now:5.1f} | {c:6.4f} | {n:4d}")
    if c>best_cs: best_cs=c; best_s=x_now
print(f"Best (scale): {best_s} deg")
print()

print("=== Energy-based mapping (1+z = urxe/urx0) ===")
print("x_now | corr | N")
best_e=None; best_ce=-1
for x_now in [45, 60, 67.5, 75, 80, 85]:
    c, n = corr_for_xnow(x_now, 'energy')
    print(f"{x_now:5.1f} | {c:6.4f} | {n:4d}")
    if c>best_ce: best_ce=c; best_e=x_now
print(f"Best (energy): {best_e} deg")
print()

print(f"CONCLUSION: Pantheon+ (N={len(data)}) calibration.")
print(f"  Scale-based best: x~{best_s}deg, Energy-based best: x~{best_e}deg.")

with open('/root/pantheon_cal.txt','w') as f:
    f.write(f"best_scale={best_s}\nbest_energy={best_e}\nN={len(data)}\n")
print("WROTE /root/pantheon_cal.txt")
