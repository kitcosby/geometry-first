from sage.all import *
# E97: Kit's FW matrix variants.
# A: sin(2x),cos(2x) -> 0.
# B: sin(2x),cos(2x) -> sin(x),cos(x).
# Both traceless-corrected. Test rank+closure.

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
    return srx,sxp,cxp,crx,urx,uxp,s,c

def M_variant(deg, mode):
    # mode: 'zero' or 'sincos'.
    srx,sxp,cxp,crx,urx,uxp,s,c = ops_exact(deg)
    FW = 1/urx + 1/uxp
    if mode == 'zero':
        a13 = K(0); a31 = K(0)
    elif mode == 'sincos':
        a13 = c; a31 = s  # cos(x), sin(x)
    M = matrix(K, 3, 3, [
        [urx, srx, a13],
        [cxp, uxp, crx],
        [a31, sxp, FW],
    ])
    tr = M.trace()
    return M - (tr/3)*identity_matrix(K, 3)

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

def test_mode(mode, label):
    print(f"\n=== {label} ===", flush=True)
    G = [M_variant(p, mode) for p in [15,30,45,60]]
    for i in [1,2,3]:
        G.append(G[0]*G[i]-G[i]*G[0])
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
        if not found: break
    print(f"  Rank: {r}", flush=True)
    if r!=8:
        print(f"  RESULT: rank {r}, not sl(3).", flush=True)
        return r, False
    # Closure.
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
    except:
        print(f"  RESULT: singular.", flush=True)
        return r, False
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
    print(f"  Closure: {'YES' if closed else 'NO'}", flush=True)
    if closed:
        print(f"  RESULT: {label} gives sl(3,R). PROVED.", flush=True)
    else:
        print(f"  RESULT: rank 8 but NO closure.", flush=True)
    return r, closed

print("E97: FW matrix variants.", flush=True)
r1,c1 = test_mode('zero', 'A: zeros')
r2,c2 = test_mode('sincos', 'B: sin(x),cos(x)')
with open('/root/e97_done.txt','w') as fh:
    fh.write(f"zero: rank={r1} closed={c1}\nsincos: rank={r2} closed={c2}\n")
