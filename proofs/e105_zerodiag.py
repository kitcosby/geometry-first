from sage.all import *
# E105: Zero-diagonal, all 6 off-diagonals with core operators.
# M = [[0, srx, urx],[cxp, 0, sxp],[uxp, crx, 0]]
# Test rank and closure.

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
def ops_exact(deg):
    s,c = trig[deg]
    srx = (1+c)/s; sxp = 1/srx
    cxp = (1+s)/c; crx = 1/cxp
    urx = srx-crx; uxp = cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

def M_zero_diag(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    M = matrix(K, 3, 3, [
        [0,   srx, urx],
        [cxp, 0,   sxp],
        [uxp, crx, 0  ],
    ])
    return M

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print("E105: Zero-diagonal, 6 off-diagonals.", flush=True)
G = [M_zero_diag(p) for p in [15,30,45,60]]
print(f"Trace (should be 0): {[g.trace() for g in G]}", flush=True)
r0 = mrank(G)
print(f"Initial rank (4 samples): {r0}", flush=True)
# Add commutators.
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
r1 = mrank(G)
print(f"After seed commutators: {r1}", flush=True)
r=r1; cnt=0
while r<8 and cnt<20:
    cnt+=1; found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True; break
        if found: break
    if not found: break
print(f"Final rank: {r}", flush=True)
if r==8:
    print("Rank 8 achieved.", flush=True)
    # Closure check.
    n=3
    G8=[]
    for g in G:
        if len(G8)>=8: break
        if mrank(G8+[g])>mrank(G8):
            G8.append(g)
    Gflat=[[G8[k][a,b] for a in range(n) for b in range(n)] for k in range(8)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    if len(piv)==8:
        Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(8)] for a in range(8)])
        Ninv=Nmat.inverse()
        closed=True
        for a in range(8):
            for b in range(a+1,8):
                comm=G8[a]*G8[b]-G8[b]*G8[a]
                dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(8)])
                cf=Ninv*dv
                rc=sum(cf[k]*G8[k] for k in range(8))
                if not (rc-comm).is_zero():
                    closed=False; break
            if not closed: break
        print(f"Closure: {'YES' if closed else 'NO'}", flush=True)
        if closed:
            # Killing.
            f={}
            for a in range(8):
                for b in range(8):
                    comm=G8[a]*G8[b]-G8[b]*G8[a]
                    dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(8)])
                    fij=Ninv*dv
                    for k in range(8):
                        if fij[k]!=0: f[(a,b,k)]=fij[k]
            Kmat=matrix(K, 8, 8)
            for a in range(8):
                for b in range(a,8):
                    s=K(0)
                    for k in range(8):
                        for l in range(8):
                            s+=f.get((a,k,l),K(0))*f.get((b,l,k),K(0))
                    Kmat[a,b]=s; Kmat[b,a]=s
            ev=Kmat.eigenvalues()
            pos=sum(1 for e in ev if e>0); neg=sum(1 for e in ev if e<0); zer=sum(1 for e in ev if e==0)
            print(f"Killing signature: ({pos},{neg},{zer})", flush=True)
            if (pos,neg,zer)==(5,3,0):
                print("RESULT: sl(3,R) PROVED. Fourth construction, zero-diagonal.", flush=True)
else:
    print(f"RESULT: Rank {r}, not 8.", flush=True)
with open('/root/e105_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
