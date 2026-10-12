from sage.all import *
# E96: Traceless-corrected version of Kit's FW matrix.
# M_traceless = M - (tr(M)/3)*I.
# Test if it closes.

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

def M_kit_traceless(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    S2 = 1/(urx+uxp)
    C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    sin2x = 4*S2; cos2x = C2
    FW = 1/urx + 1/uxp
    M = matrix(K, 3, 3, [
        [urx, srx, cos2x],
        [cxp, uxp, crx],
        [sin2x, sxp, FW],
    ])
    # Make traceless.
    tr = M.trace()
    M_tl = M - (tr/3)*identity_matrix(K, 3)
    return M_tl

print("E96: Traceless Kit FW matrix.", flush=True)
for deg in [15,30,45,60]:
    M = M_kit_traceless(deg)
    print(f"  {deg}deg: trace = {M.trace()} (should be 0)", flush=True)

G = [M_kit_traceless(p) for p in [15,30,45,60]]
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print(f"\nStarting rank: {mrank(G)}", flush=True)
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
print(f"After seed: {mrank(G)} ({len(G)} mats)", flush=True)
r = mrank(G)
cnt=0
while r<8 and cnt<30:
    cnt+=1; found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True; break
        if found: break
    if not found:
        print(f"Stalled at {r}.", flush=True); break
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
            print("RESULT: Traceless FW matrix gives sl(3,R). PROVED.", flush=True)
        else:
            print("RESULT: Rank 8 but NO closure. Not a Lie algebra.", flush=True)
    except Exception as e:
        print(f"Error: {e}", flush=True)
else:
    print(f"RESULT: Rank {r}, not 8.", flush=True)
with open('/root/e96_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
