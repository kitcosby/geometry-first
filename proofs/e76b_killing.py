from sage.all import *
# E76b: Killing signature for (lam,mu) != (0,0).
# Is it still sl(3,R)?

import math
import numpy as np

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=(1+c)/s; sxp=1/srx
    cxp=(1+s)/c; crx=1/cxp
    urx=srx-crx; uxp=cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def Ta(deg, lam, mu):
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    return np.array([[urx, srx, lam],[cxp, -urx, uxp],[mu, -uxp, 0]], dtype=float)

def killing_signature(lam, mu):
    G=[]
    for p in [15,30,45,60]:
        G.append(Ta(p, lam, mu))
    G.append(G[0]@G[1]-G[1]@G[0])
    G.append(G[0]@G[2]-G[2]@G[0])
    G.append(G[0]@G[3]-G[3]@G[0])
    A=G[0]@G[1]-G[1]@G[0]
    G.append(A@G[2]-G[2]@A)
    n=8
    # Structure constants via least squares.
    M=np.array([g.flatten() for g in G]).T
    f=np.zeros((n,n,n))
    for i in range(n):
        for j in range(n):
            comm=(G[i]@G[j]-G[j]@G[i]).flatten()
            coeffs=np.linalg.lstsq(M, comm, rcond=None)[0]
            f[i,j,:]=coeffs
    # Killing: K_{ij} = sum_{k,l} f_{ikl} f_{jlk}
    Kmat=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            s=0
            for k in range(n):
                for l in range(n):
                    s+=f[i,k,l]*f[j,l,k]
            Kmat[i,j]=s
    ev=np.linalg.eigvalsh(Kmat)
    pos=np.sum(ev>1e-6); neg=np.sum(ev<-1e-6); zero=np.sum(np.abs(ev)<=1e-6)
    return pos, neg, zero, ev

print("Killing signature for various (lam, mu):")
for lam, mu in [(0,0), (1,0), (0,1), (1,1), (2,3), (-1,0.5)]:
    pos, neg, zero, ev = killing_signature(lam, mu)
    print(f"  ({lam},{mu}): ({pos},{neg},{zero})")
    if not (pos==5 and neg==3 and zero==0):
        print(f"    Eigenvalues: {ev}")
        print(f"    NOT sl(3,R)!")
print()
print("If all (5,3,0), the corners don't change the Lie algebra type.")
