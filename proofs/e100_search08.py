from sage.all import *
# E100: Random search for (0,8) Killing from trig operators.
# Try random 3x3 traceless matrices from operator pool.
# Test rank 8, closure, Killing. Look for (0,8).
import random

K = QQbar
sqrt2 = K(2).sqrt()
sqrt3 = K(3).sqrt()
sqrt6 = sqrt2*sqrt3
trig = {
    15: ((sqrt6-sqrt2)/4, (sqrt6+sqrt2)/4),
    30: (K(1)/2, sqrt3/2),
    45: (sqrt2/2, sqrt2/2),
    60: (sqrt3/2, K(1)/2),
}
def ops_at(deg):
    s,c = trig[deg]
    srx = (1+c)/s; sxp = 1/srx
    cxp = (1+s)/c; crx = 1/cxp
    urx = srx-crx; uxp = cxp-sxp
    sawup = 1/urx; sawdown = 1/uxp
    FW = sawup + sawdown
    # Return dict.
    return {
        'srx': srx, 'sxp': sxp, 'cxp': cxp, 'crx': crx,
        'urx': urx, 'uxp': uxp, 'sawup': sawup, 'sawdown': sawdown,
        'FW': FW, 'zero': K(0), 'one': K(1),
    }

OP_NAMES = ['srx','sxp','cxp','crx','urx','uxp','sawup','sawdown','FW','zero']

def random_matrix(deg, seed):
    # Seed is a list of 9 op names for positions (0,0)...(2,2).
    # Make traceless by adjusting (2,2).
    ops = ops_at(deg)
    vals = [ops[name] for name in seed[:8]]  # first 8 positions
    # Positions: (0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)
    # (2,2) = -(sum of diag)
    a00,a01,a02,a10,a11,a12,a20,a21 = vals
    a22 = -(a00 + a11)
    M = matrix(K, 3, 3, [
        [a00,a01,a02],
        [a10,a11,a12],
        [a20,a21,a22],
    ])
    return M

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

def test_seed(seed):
    # Build 4 matrices, test rank and closure.
    G = [random_matrix(p, seed) for p in [15,30,45,60]]
    # Quick rank check.
    if mrank(G) < 4:
        return None  # not even rank 4
    for i in [1,2,3]:
        G.append(G[0]*G[i]-G[i]*G[0])
    r=mrank(G); cnt=0
    while r<8 and cnt<20:
        cnt+=1; found=False
        for i in range(len(G)):
            for j in range(i+1,len(G)):
                comm=G[i]*G[j]-G[j]*G[i]
                if mrank(G+[comm])>r:
                    G.append(comm); r=mrank(G); found=True; break
            if found: break
        if not found: break
    if r!=8:
        return None
    # Closure.
    G8=[]
    for g in G:
        if len(G8)>=8: break
        if mrank(G8+[g])>mrank(G8):
            G8.append(g)
    n=3
    Gflat=[[G8[k][a,b] for a in range(n) for b in range(n)] for k in range(8)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    if len(piv)!=8:
        return None
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(8)] for a in range(8)])
    try:
        Ninv=Nmat.inverse()
    except:
        return None
    for i in range(8):
        for j in range(i+1,8):
            comm=G8[i]*G8[j]-G8[j]*G8[i]
            dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(8)])
            cf=Ninv*dv
            rc=sum(cf[k]*G8[k] for k in range(8))
            if not (rc-comm).is_zero():
                return None  # no closure
    # Killing.
    f={}
    for i in range(8):
        for j in range(8):
            comm=G8[i]*G8[j]-G8[j]*G8[i]
            dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(8)])
            fij=Ninv*dv
            for k in range(8):
                if fij[k]!=0: f[(i,j,k)]=fij[k]
    Kmat=matrix(K, 8, 8)
    for i in range(8):
        for j in range(i,8):
            s=K(0)
            for k in range(8):
                for l in range(8):
                    s+=f.get((i,k,l),K(0))*f.get((j,l,k),K(0))
            Kmat[i,j]=s; Kmat[j,i]=s
    ev=Kmat.eigenvalues()
    pos=sum(1 for e in ev if e>0); neg=sum(1 for e in ev if e<0); z=sum(1 for e in ev if e==0)
    return (pos,neg,z)

print("E100: Random search for (0,8).", flush=True)
random.seed(42)
found=False
for trial in range(100):
    seed = [random.choice(OP_NAMES) for _ in range(8)]
    # Skip if all zeros on diagonal? Just try.
    result = test_seed(seed)
    if result is not None:
        pos,neg,z = result
        print(f"  Trial {trial}: closed, Killing ({pos},{neg},{z}), seed={seed}", flush=True)
        if pos==0 and neg==8 and z==0:
            print(f"  *** FOUND (0,8)! seed={seed}", flush=True)
            found=True
            break
    if trial%20==19:
        print(f"  ... {trial+1}/100", flush=True)
if not found:
    print("  No (0,8) found in 100 trials.", flush=True)
with open('/root/e100_done.txt','w') as fh:
    fh.write(f"found={found}\n")
