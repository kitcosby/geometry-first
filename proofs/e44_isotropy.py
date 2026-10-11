from sage.all import *
# E44: Isotropy test on Pantheon+.
# If homogeneous+isotropic, mu residuals (from best urx model) should NOT
# correlate with sky position (RA, Dec).
# Test: split sky into hemispheres, compare mean residual.

import math, random

data = []  # (z, mu, err, ra, dec)
with open('/root/pantheon_api.dat') as f:
    header = f.readline().split()
    iz=header.index('zHD'); imu=header.index('MU_SH0ES')
    ierr=header.index('MU_SH0ES_ERR_DIAG')
    ira=header.index('RA'); idec=header.index('DEC')
    for line in f:
        parts=line.split()
        try:
            z=float(parts[iz]); mu=float(parts[imu]); err=float(parts[ierr])
            ra=float(parts[ira]); dec=float(parts[idec])
            if 0.01<z<2.5 and 30<mu<50 and err>0:
                data.append((z,mu,err,ra,dec))
        except: pass
print(f"Loaded {len(data)} SNe with RA/Dec.")

def urx_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    return (1+c)/s - 1/((1+s)/c)

def find_xe(z, x_now):
    # 1+z = urx(xe)/urx(x_now)
    unow=urx_phase(x_now)
    target=(1+z)*unow
    best=None; be=1e9
    for d in [0.5+i*0.25 for i in range(int((x_now-0.5)*4))]:
        v=urx_phase(d)
        if v is None: continue
        e=abs(v-target)
        if e<be: be=e; best=d
    return best if be<0.1*target else None

# Best model from E39: urx, x_now=85, mu = alpha*(x_now-xe)+beta, alpha=0.259, beta=35.143
x_now=85.0; alpha=0.259; beta=35.143
resids=[]  # (ra, dec, resid, err)
for z,mu,err,ra,dec in data:
    xe=find_xe(z,x_now)
    if xe is None: continue
    mu_pred=alpha*(x_now-xe)+beta
    resids.append((ra,dec,(mu-mu_pred)/err,err))

print(f"Computed {len(resids)} residuals.")
# Split by RA hemisphere: RA<180 vs RA>=180
r1=[r for ra,dec,r,e in resids if ra<180]
r2=[r for ra,dec,r,e in resids if ra>=180]
m1=sum(r1)/len(r1); m2=sum(r2)/len(r2)
# Standard error of mean
import statistics
s1=statistics.stdev(r1)/math.sqrt(len(r1)) if len(r1)>1 else 0
s2=statistics.stdev(r2)/math.sqrt(len(r2)) if len(r2)>1 else 0
print(f"\nRA<180:  N={len(r1)}, mean resid={m1:.4f} +/- {s1:.4f}")
print(f"RA>=180: N={len(r2)}, mean resid={m2:.4f} +/- {s2:.4f}")
print(f"Difference: {abs(m1-m2):.4f} (significance: {abs(m1-m2)/math.sqrt(s1*s1+s2*s2):.2f} sigma)")
print()
# Split by Dec: North vs South
n1=[r for ra,dec,r,e in resids if dec>=0]
n2=[r for ra,dec,r,e in resids if dec<0]
mn1=sum(n1)/len(n1); mn2=sum(n2)/len(n2)
sn1=statistics.stdev(n1)/math.sqrt(len(n1)) if len(n1)>1 else 0
sn2=statistics.stdev(n2)/math.sqrt(len(n2)) if len(n2)>1 else 0
print(f"Dec>=0 (N): N={len(n1)}, mean resid={mn1:.4f} +/- {sn1:.4f}")
print(f"Dec<0 (S):  N={len(n2)}, mean resid={mn2:.4f} +/- {sn2:.4f}")
print(f"Difference: {abs(mn1-mn2):.4f} (significance: {abs(mn1-mn2)/math.sqrt(sn1*sn1+sn2*sn2):.2f} sigma)")
print()
if abs(m1-m2)/math.sqrt(s1*s1+s2*s2) < 3 and abs(mn1-mn2)/math.sqrt(sn1*sn1+sn2*sn2) < 3:
    print("CONCLUSION: No significant anisotropy (<3 sigma). Isotropy HOLDS.")
else:
    print("CONCLUSION: Significant anisotropy detected. Isotropy FALSIFIED.")
print("WROTE /root/e44_iso.txt")
