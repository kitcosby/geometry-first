from sage.all import *
# E39: S1-compliant cosmological fit.
# - Four readouts: urx (log), 1/urx (linear), uxp (log), sawx (linear)
# - Chi-squared with diagonal errors (declared uncertainty model)
# - Held-out validation: fit on 50%, validate on 50%
# - Same fitting record for all models (fair comparison)

import math, random

# Load Pantheon+ data
data = []  # (z, mu, mu_err, idx)
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

# Split into fit (50%) and validation (50%), deterministic
random.seed(42)
indices = list(range(N))
random.shuffle(indices)
n_fit = N // 2
fit_idx = set(indices[:n_fit])
fit_data = [d for i,d in enumerate(data) if i in fit_idx]
val_data = [d for i,d in enumerate(data) if i not in fit_idx]
print(f"Fit: {len(fit_data)}, Validation: {len(val_data)}")
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

def find_x_for_urx_ratio(target_ratio, x_now):
    # Find x_emit such that urx(x_emit)/urx(x_now) = target_ratio
    # target_ratio = 1+z > 1, urx decreases with x, so x_emit < x_now
    unow = urx_deg(x_now)
    if unow is None: return None
    target = target_ratio * unow
    best=None; best_err=1e9
    # x_emit in (0.5, x_now)
    for d in [0.5 + i*0.25 for i in range(int((x_now-0.5)*4))]:
        v = urx_deg(d)
        if v is None: continue
        err=abs(v-target)
        if err<best_err: best_err=err; best=d
    if best_err > 0.1*target: return None
    return best

def find_x_for_uxp_ratio(target_ratio, x_now):
    # Find x_emit such that uxp(x_now)/uxp(x_emit) = target_ratio
    # target_ratio = 1+z > 1, uxp increases with x, so x_emit < x_now
    unow = uxp_deg(x_now)
    if unow is None: return None
    target = unow / target_ratio
    best=None; best_err=1e9
    for d in [0.5 + i*0.25 for i in range(int((x_now-0.5)*4))]:
        v = uxp_deg(d)
        if v is None: continue
        err=abs(v-target)
        if err<best_err: best_err=err; best=d
    if best_err > 0.1*target: return None
    return best

# Four readout models. Each: given x_now, compute mu_pred for each SN.
# Then fit (alpha, beta) by linear least squares, compute chi2.

def model_urx_log(x_now, dataset):
    # 1+z = urx_emit/urx_now. mu = alpha*(x_now - x_emit) + beta.
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_urx_ratio(1+z, x_now)
        if xe is None: return None
        preds.append((x_now - xe, mu, err))
    return preds

def model_invurx_lin(x_now, dataset):
    # mu = alpha*(1/urx_now - 1/urx_emit) + beta. (entropy difference)
    unow = urx_deg(x_now)
    if unow is None: return None
    sinv_now = 1/unow
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_urx_ratio(1+z, x_now)
        if xe is None: return None
        ue = urx_deg(xe)
        if ue is None: return None
        preds.append((sinv_now - 1/ue, mu, err))
    return preds

def model_uxp_log(x_now, dataset):
    # 1+z = uxp_now/uxp_emit. mu = alpha*(x_now - x_emit) + beta.
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_uxp_ratio(1+z, x_now)
        if xe is None: return None
        preds.append((x_now - xe, mu, err))
    return preds

def model_sawx_lin(x_now, dataset):
    # sawx = 1/uxp. mu = alpha*(sawx_emit - sawx_now) + beta.
    # sawx decreases with x (uxp increases), so sawx_emit > sawx_now.
    unow = uxp_deg(x_now)
    if unow is None: return None
    saw_now = 1/unow
    preds=[]
    for z, mu, err, idx in dataset:
        xe = find_x_for_uxp_ratio(1+z, x_now)
        if xe is None: return None
        ue = uxp_deg(xe)
        if ue is None: return None
        preds.append((1/ue - saw_now, mu, err))
    return preds

def fit_alpha_beta(preds):
    # Weighted linear least squares: mu = alpha*x + beta.
    # Returns (alpha, beta, chi2).
    # chi2 = sum((mu - (alpha*x+beta))/err)^2
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
    'urx_log': model_urx_log,
    'invurx_lin': model_invurx_lin,
    'uxp_log': model_uxp_log,
    'sawx_lin': model_sawx_lin,
}

print("=== FITTING (50% of data) ===")
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
        dof = len(preds)-3  # x_now, alpha, beta
        if best is None or chi2 < best[3]:
            best = (x_now, alpha, beta, chi2, dof, len(preds))
    if best:
        x_now, alpha, beta, chi2, dof, n = best
        print(f"{name:10s} | {x_now:5.1f} | {alpha:6.3f} | {beta:6.3f} | {chi2:7.1f}/{dof}={chi2/dof:.3f}")
        results[name]=best
    else:
        print(f"{name:10s} | FAILED")
print()

print("=== VALIDATION (held-out 50%) ===")
print("Model      | x_now | alpha  | beta   | chi2/dof (val)")
print("-----------|-------|--------|--------|---------------")
for name, mfn in models.items():
    if name not in results: continue
    x_now, alpha, beta, chi2_fit, dof_fit, n_fit = results[name]
    preds = mfn(x_now, val_data)
    if preds is None:
        print(f"{name:10s} | FAILED")
        continue
    chi2_val = sum(((mu-(alpha*x+beta))/err)**2 for x,mu,err in preds)
    dof_val = len(preds)  # parameters fixed from fit
    print(f"{name:10s} | {x_now:5.1f} | {alpha:6.3f} | {beta:6.3f} | {chi2_val:7.1f}/{dof_val}={chi2_val/dof_val:.3f}")
print()
print("Uncertainty model: diagonal (MU_SH0ES_ERR_DIAG). Full covariance available for refinement.")
print("WROTE /root/e39_results.txt")
with open('/root/e39_results.txt','w') as f:
    for name, b in results.items():
        f.write(f"{name}: x_now={b[0]}, alpha={b[1]:.4f}, beta={b[2]:.4f}, chi2_fit/dof={b[3]/b[4]:.4f}\n")
