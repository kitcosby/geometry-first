from sage.all import *
# E40: Test quasi-linear saw readouts.
# User hypothesis: the correct Gamma_theta uses saw coordinates (quasi-linear).
# sawr = 1/urx = lambda in QI (transfer fraction, 0->1).
# Test: mu = alpha*(sawr_now - sawr_emit) + beta, and log version.

import math, random

data = []
with open('/root/pantheon_api.dat') as f:
    header = f.readline().split()
    iz = header.index('zHD')
    imu = header.index('MU_SH0ES')
    ierr = header.index('MU_SH0ES_ERR_DIAG')
    for idx, line in enumerate(f):
        parts = line.split()
        try:
            z = float(parts[iz]); mu = float(parts[imu]); err = float(parts[ierr])
            if 0.01 < z < 2.5 and 30 < mu < 50 and err > 0:
                data.append((z, mu, err, idx))
        except: pass

print(f"Loaded {len(data)} SNe.")
N = len(data)
random.seed(42)
indices = list(range(N)); random.shuffle(indices)
n_fit = N // 2
fit_idx = set(indices[:n_fit])
fit_data = [d for i,d in enumerate(data) if i in fit_idx]
val_data = [d for i,d in enumerate(data) if i not in fit_idx]
print(f"Fit: {len(fit_data)}, Validation: {len(val_data)}")
print()

def urx_deg(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    srx=(1+c)/s; cxp=(1+s)/c
    return srx - 1/cxp

def find_x_for_urx_ratio(target_ratio, x_now):
    unow = urx_deg(x_now)
    if unow is None: return None
    target = target_ratio * unow
    best=None; best_err=1e9
    for d in [0.5 + i*0.25 for i in range(int((x_now-0.5)*4))]:
        v = urx_deg(d)
        if v is None: continue
        err=abs(v-target)
        if err<best_err: best_err=err; best=d
    if best_err > 0.1*target: return None
    return best

def model_sawr_lin(x_now, dataset):
    # sawr = 1/urx. mu = alpha*(sawr_now - sawr_emit) + beta.
    unow = urx_deg(x_now)
    if unow is None: return None
    saw_now = 1/unow
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_urx_ratio(1+z, x_now)
        if xe is None: return None
        ue = urx_deg(xe)
        if ue is None or ue==0: return None
        preds.append((saw_now - 1/ue, mu, err))
    return preds

def model_sawr_log(x_now, dataset):
    # mu = alpha*ln(sawr_now/sawr_emit) + beta.
    unow = urx_deg(x_now)
    if unow is None: return None
    saw_now = 1/unow
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_urx_ratio(1+z, x_now)
        if xe is None: return None
        ue = urx_deg(xe)
        if ue is None or ue==0: return None
        saw_emit = 1/ue
        if saw_emit<=0 or saw_now<=0: return None
        preds.append((math.log(saw_now/saw_emit), mu, err))
    return preds

def fit_alpha_beta(preds):
    Sx=Sy=Sxx=Sxy=Sw=0
    for x, mu, err in preds:
        w=1/(err*err)
        Sx+=w*x; Sy+=w*mu; Sxx+=w*x*x; Sxy+=w*x*mu; Sw+=w
    denom = Sw*Sxx - Sx*Sx
    if abs(denom)<1e-12: return None
    alpha = (Sw*Sxy - Sx*Sy)/denom
    beta = (Sy - alpha*Sx)/Sw
    chi2 = sum(((mu-(alpha*x+beta))/err)**2 for x,mu,err in preds)
    return alpha, beta, chi2

models = {
    'sawr_lin': model_sawr_lin,
    'sawr_log': model_sawr_log,
}

print("=== FITTING (quasi-linear saws) ===")
print("Model      | x_now | alpha  | beta   | chi2/dof")
print("-----------|-------|--------|--------|----------")
results={}
for name, mfn in models.items():
    best=None
    for x_now in [60, 65, 70, 75, 80, 85]:
        preds = mfn(x_now, fit_data)
        if preds is None: continue
        fb = fit_alpha_beta(preds)
        if fb is None: continue
        alpha, beta, chi2 = fb
        dof = len(preds)-3
        if best is None or chi2 < best[3]:
            best = (x_now, alpha, beta, chi2, dof, len(preds))
    if best:
        x_now, alpha, beta, chi2, dof, n = best
        print(f"{name:10s} | {x_now:5.1f} | {alpha:6.3f} | {beta:6.3f} | {chi2:7.1f}/{dof}={chi2/dof:.3f}")
        results[name]=best
    else:
        print(f"{name:10s} | FAILED")
print()

print("=== VALIDATION ===")
print("Model      | x_now | chi2/dof (val)")
print("-----------|-------|----------------")
for name, mfn in models.items():
    if name not in results: continue
    x_now, alpha, beta, chi2_fit, dof_fit, n_fit = results[name]
    preds = mfn(x_now, val_data)
    if preds is None:
        print(f"{name:10s} | FAILED"); continue
    chi2_val = sum(((mu-(alpha*x+beta))/err)**2 for x,mu,err in preds)
    print(f"{name:10s} | {x_now:5.1f} | {chi2_val:7.1f}/{len(preds)}={chi2_val/len(preds):.3f}")
print()
print("WROTE /root/e40_results.txt")
with open('/root/e40_results.txt','w') as f:
    for name, b in results.items():
        f.write(f"{name}: x_now={b[0]}, chi2/dof={b[3]/b[4]:.4f}\n")
