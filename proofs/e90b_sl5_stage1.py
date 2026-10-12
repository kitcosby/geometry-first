from sage.all import *
# E90b: Exact sl(5,R), staged. Stage 1: rank 24 + closure. (Killing later.)
# With progress logging.

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

def X5_exact(deg):
    # Build via X4 then extend.
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    # X4.
    A3 = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    S2 = 1/(urx+uxp); C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    X4 = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            X4[a,b]=A3[a,b]
    X4[0,3]=S2; X4[1,3]=C2; X4[3,0]=S2; X4[3,1]=C2
    # X5: 4x.
    S4 = 8*S2*C2; C4 = C2*C2 - 16*S2*S2
    X5 = matrix(K, 5, 5)
    for a in range(4):
        for b in range(4):
            X5[a,b]=X4[a,b]
    X5[0,4]=S4; X5[1,4]=C4; X5[4,0]=S4; X5[4,1]=C4
    return X5

print("E90b: Stage 1 (rank + closure).", flush=True)
G = [X5_exact(p) for p in [15,30,45,60]]
print(f"  4 base.", flush=True)
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
print(f"  +3 commutators = {len(G)}.", flush=True)
A01=G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)
print(f"  +1 nested = {len(G)}.", flush=True)

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

r=mrank(G)
print(f"  Rank: {r}.", flush=True)
# Grow.
cnt=0
while r<24 and cnt<50:
    cnt+=1
    found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True
                print(f"  Added, rank {r} (total {len(G)}).", flush=True)
                break
        if found: break
    if not found:
        print("  No more independent.", flush=True)
        break
print(f"Final: {len(G)} mats, rank {r}.", flush=True)
# Select 24.
G24=[]
for g in G:
    if len(G24)>=24: break
    if mrank(G24+[g])>mrank(G24):
        G24.append(g)
print(f"Selected {len(G24)}.", flush=True)
if len(G24)==24:
    print("PROVED: rank 24.", flush=True)
    # Closure.
    n=5
    Gflat=[[G24[k][a,b] for a in range(n) for b in range(n)] for k in range(24)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(24)] for a in range(24)])
    Ninv=Nmat.inverse()
    print("  Checking 276 commutators...", flush=True)
    closed=True
    for i in range(24):
        for j in range(i+1,24):
            comm=G24[i]*G24[j]-G24[j]*G24[i]
            dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(24)])
            cf=Ninv*dv
            rc=sum(cf[k]*G24[k] for k in range(24))
            if not (rc-comm).is_zero():
                print(f"  FAIL [{i},{j}].", flush=True); closed=False; break
        if not closed: break
        if i%6==5:
            print(f"  ... {i+1}/24 done.", flush=True)
    if closed:
        print("PROVED: closure. 24-dim Lie algebra.", flush=True)
print("WROTE /root/e90b_done.txt", flush=True)
with open('/root/e90b_done.txt','w') as fh:
    fh.write(f"rank={len(G24)}\n")
