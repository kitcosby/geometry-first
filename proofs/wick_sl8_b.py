import numpy as np
rt3 = np.sqrt(3)
a = np.array([2+rt3, 2-rt3, rt3, rt3/3, 2+2*rt3/3, 3/4-rt3/4, -2+2*rt3, 1/4+rt3/4])
def E(i,j):
    M = np.zeros((8,8), dtype=complex); M[i,j] = 1; return M
W = []
for i in range(8):
    Xi = a[i]*E(i, (i+1) % 8)
    W.append((Xi - Xi.T)/2 + 1j*(Xi + Xi.T)/2)
# Claim: U W_i U^dagger = (a_i/sqrt2)(E_{i,i+1} - E_{i+1,i}), U = diag(e^{-i k pi/4}).
ph = np.exp(1j*np.pi/4*np.arange(8))
U = np.diag(ph)
ok = True
for i in range(8):
    M = U @ W[i] @ U.conj().T
    target = (a[i]/np.sqrt(2))*(np.eye(8)[ :, [(i+1)%8]]*0 + E(i,(i+1)%8).real - E((i+1)%8, i).real)
    ok &= np.allclose(M, target) and np.allclose(M.imag, 0)
print('U W_i U^dag = (a_i/sqrt2)(E_{i,i+1}-E_{i+1,i}) real antisym, all i:', ok)
# So generated algebra = U^dag so(8) U, dim 28. Trace form on a basis of it: negative definite (compact).
def real_vec(M): return np.concatenate([M.real.flatten(), M.imag.flatten()])
def rank_of(mats): return np.linalg.matrix_rank(np.array([real_vec(m) for m in mats]).T, tol=1e-7) if mats else 0
G = list(W); r = rank_of(G)
for _ in range(14):
    cur = list(G); added = False
    for i in range(len(cur)):
        for j in range(i+1, len(cur)):
            comm = cur[i]@cur[j] - cur[j]@cur[i]
            if rank_of(G+[comm]) > r: G.append(comm); r = rank_of(G); added = True
    if not added: break
print('dim =', r)
basis = G[:28]
T = np.array([[np.trace(basis[i]@basis[j]).real for j in range(28)] for i in range(28)])
ev = np.linalg.eigvalsh(T)
print('tr-form: max eigenvalue =', round(ev.max(),3), ' => negative definite (compact):', ev.max() < 0)
