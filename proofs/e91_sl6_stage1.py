from sage.all import *
# E91: Exact sl(6,R), Stage 1: rank 35 + closure (595 commutators).
# X6(x) = [[X5(x), v5],[w5^T,0]], v5=w5=[S8_op, C8_op, 0,0,0].
# S8=sin(8x), C8=cos(8x) via iterated operator double-angle.

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

def X6_exact(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    # X3.
    A3 = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    # 2x.
    S2 = 1/(urx+uxp); C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    # X4.
    X4 = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            X4[a,b]=A3[a,b]
    X4[0,3]=S2; X4[1,3]=C2; X4[3,0]=S2; X4[3,1]=C2
    # 4x.
    S4 = 8*S2*C2; C4 = C2*C2 - 16*S2*S2
    # X5.
    X5 = matrix(K, 5, 5)
    for a in range(4):
        for b in range(4):
            X5[a,b]=X4[a,b]
    X5[0,4]=S4; X5[1,4]=C4; X5[4,0]=S4; X5[4,1]=C4
    # 8x.
    S8 = 2*S4*C4; C8 = C4*C4 - S4*S4
    # X6.
    X6 = matrix(K, 6, 6)
    for a in range(5):
        for b in range(5):
            X6[a,b]=X5[a,b]
    X6[0,5]=S8; X6[1,5]=C8; X6[5,0]=S8; X6[5,1]=C8
    return X6

print("E91: sl(6) Stage 1.", flush=True)
G = [X6_exact(p) for p in [15,30,45,60]]
print(f"  4 base.", flush=True)
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
A01=G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)
print(f"  Total {len(G)}.", flush=True)

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

r=mrank(G)
print(f"  Rank {r}.", flush=True)
cnt=0
while r<35 and cnt<80:
    cnt+=1
    found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True
                if r%5==0:
                    print(f"  Rank {r} ({len(G)} mats).", flush=True)
                break
        if found: break
    if not found:
        print("  Stalled.", flush=True); break
print(f"Final rank {r}, {len(G)} mats.", flush=True)

G35=[]
for g in G:
    if len(G35)>=35: break
    if mrank(G35+[g])>mrank(G35):
        G35.append(g)
print(f"Selected {len(G35)}.", flush=True)
if len(G35)==35:
    print("PROVED: rank 35.", flush=True)
    n=6
    Gflat=[[G35[k][a,b] for a in range(n) for b in range(n)] for k in range(35)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(35)] for a in range(35)])
    Ninv=Nmat.inverse()
    print("  595 commutators...", flush=True)
    closed=True
    for i in range(35):
        for j in range(i+1,35):
            comm=G35[i]*G35[j]-G35[j]*G35[i]
            dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(35)])
            cf=Ninv*dv
            rc=sum(cf[k]*G35[k] for k in range(35))
            if not (rc-comm).is_zero():
                print(f"  FAIL [{i},{j}].", flush=True); closed=False; break
        if not closed: break
        if i%7==6:
            print(f"  ... {i+1}/35.", flush=True)
    if closed:
        print("PROVED: closure. 35-dim Lie algebra.", flush=True)
with open('/root/e91_done.txt','w') as fh:
    fh.write(f"rank={len(G35)}\n")
print("WROTE /root/e91_done.txt", flush=True)
