from sage.all import *
# E99: Block construction with w = -v (antisymmetric off-diagonal).
# X4 = [[X3, v], [-v^T, 0]]. Test rank, closure, Killing.
# Hypothesis: may give a different real form (more compact).

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

def X4_wneg(deg):
    # w = -v.
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    A3 = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    S2 = 1/(urx+uxp); C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    X4 = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            X4[a,b]=A3[a,b]
    # v = [S2, C2, 0], w = -v.
    X4[0,3]=S2; X4[1,3]=C2; X4[2,3]=0
    X4[3,0]=-S2; X4[3,1]=-C2; X4[3,2]=0
    X4[3,3]=0
    return X4

print("E99: w=-v block.", flush=True)
G=[X4_wneg(p) for p in [15,30,45,60]]
# Check trace.
for g in G:
    if g.trace()!=0:
        print(f"  NOT TRACELESS!", flush=True)
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()
print(f"  Starting rank: {mrank(G)}", flush=True)
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
A01=G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)
r=mrank(G); cnt=0
while r<15 and cnt<50:
    cnt+=1; found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True; break
        if found: break
    if not found: break
print(f"  Final rank: {r} ({len(G)} mats)", flush=True)
if r==15:
    G15=[]
    for g in G:
        if len(G15)>=15: break
        if mrank(G15+[g])>mrank(G15):
            G15.append(g)
    n=4
    Gflat=[[G15[k][a,b] for a in range(n) for b in range(n)] for k in range(15)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(15)] for a in range(15)])
    try:
        Ninv=Nmat.inverse()
        closed=True
        for i in range(15):
            for j in range(i+1,15):
                comm=G15[i]*G15[j]-G15[j]*G15[i]
                dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
                cf=Ninv*dv
                rc=sum(cf[k]*G15[k] for k in range(15))
                if not (rc-comm).is_zero():
                    closed=False; break
            if not closed: break
        print(f"  Closure: {'YES' if closed else 'NO'}", flush=True)
        if closed:
            f={}
            for i in range(15):
                for j in range(15):
                    comm=G15[i]*G15[j]-G15[j]*G15[i]
                    dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
                    fij=Ninv*dv
                    for k in range(15):
                        if fij[k]!=0: f[(i,j,k)]=fij[k]
            Kmat=matrix(K, 15, 15)
            for i in range(15):
                for j in range(i,15):
                    s=K(0)
                    for k in range(15):
                        for l in range(15):
                            s+=f.get((i,k,l),K(0))*f.get((j,l,k),K(0))
                    Kmat[i,j]=s; Kmat[j,i]=s
            ev=Kmat.eigenvalues()
            pos=sum(1 for e in ev if e>0); neg=sum(1 for e in ev if e<0); z=sum(1 for e in ev if e==0)
            print(f"  Killing: ({pos},{neg},{z})", flush=True)
            if pos==0 and neg==15:
                print("  RESULT: (0,15) = su(4). Compact!", flush=True)
            elif pos==9 and neg==6:
                print("  RESULT: (9,6) = sl(4,R). Same as before.", flush=True)
            else:
                print(f"  RESULT: Different real form.", flush=True)
    except Exception as e:
        print(f"  Error: {e}", flush=True)
else:
    print(f"  RESULT: rank {r}, not 15.", flush=True)
with open('/root/e99_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
