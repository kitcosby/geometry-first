from sage.all import *
# E98: FW matrix with sawup/sawdown on diagonal.
# M = [[sawup, srx, 0],[cxp, sawdown, crx],[0, sxp, FW]]
# where sawup=1/urx, sawdown=1/uxp, FW=sawup+sawdown.
# Traceless-corrected. Test rank+closure.

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
    sawup = 1/urx; sawdown = 1/uxp
    FW = sawup + sawdown
    return srx,sxp,cxp,crx,urx,uxp,sawup,sawdown,FW

def M_saw(deg):
    srx,sxp,cxp,crx,urx,uxp,sawup,sawdown,FW = ops_exact(deg)
    M = matrix(K, 3, 3, [
        [sawup, srx, 0],
        [cxp, sawdown, crx],
        [0, sxp, FW],
    ])
    tr = M.trace()
    return M - (tr/3)*identity_matrix(K, 3)

print("E98: Sawup/Sawdown diagonal matrix.", flush=True)
for deg in [15,30,45,60]:
    M = M_saw(deg)
    print(f"  {deg}deg: trace={M.trace()}", flush=True)

G = [M_saw(p) for p in [15,30,45,60]]
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print(f"\nStarting rank: {mrank(G)}", flush=True)
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
print(f"After seed: {mrank(G)} ({len(G)})", flush=True)
r=mrank(G); cnt=0
while r<8 and cnt<30:
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
    G8=[]
    for g in G:
        if len(G8)>=8: break
        if mrank(G8+[g])>mrank(G8):
            G8.append(g)
    n=3
    Gflat=[[G8[k][a,b] for a in range(n) for b in range(n)] for k in range(8)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(8)] for a in range(8)])
    try:
        Ninv=Nmat.inverse()
        closed=True
        for i in range(8):
            for j in range(i+1,8):
                comm=G8[i]*G8[j]-G8[j]*G8[i]
                dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(8)])
                cf=Ninv*dv
                rc=sum(cf[k]*G8[k] for k in range(8))
                if not (rc-comm).is_zero():
                    closed=False; break
            if not closed: break
        print(f"Closure: {'YES' if closed else 'NO'}", flush=True)
        if closed:
            print("RESULT: Saw matrix gives sl(3,R). PROVED.", flush=True)
        else:
            print("RESULT: Rank 8 but NO closure.", flush=True)
    except Exception as e:
        print(f"Error: {e}", flush=True)
else:
    print(f"RESULT: Rank {r}.", flush=True)
with open('/root/e98_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
