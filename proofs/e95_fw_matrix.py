from sage.all import *
# E95: Test Kit's FW matrix.
# M = [[urx, srx, C2],[cxp, uxp, crx],[S2, sxp, FW]]
# where C2=cos(2x), S2=sin(2x), FW=1/urx+1/uxp.
# Check: trace, rank, closure, what algebra does it generate?

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

def M_kit(deg):
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    # Double-angle via Lemma (PROVED).
    S2 = 1/(urx+uxp)  # sin(2x)/4
    C2 = (srx+crx-sxp-cxp)/(urx+uxp)  # cos(2x)
    sin2x = 4*S2
    cos2x = C2
    FW = 1/urx + 1/uxp
    M = matrix(K, 3, 3, [
        [urx, srx, cos2x],
        [cxp, uxp, crx],
        [sin2x, sxp, FW],
    ])
    return M

print("E95: Kit's FW matrix.", flush=True)
for deg in [15,30,45,60]:
    M = M_kit(deg)
    tr = M.trace()
    print(f"  {deg}deg: trace = {complex(tr):.4f} {'(traceless)' if tr==0 else '(NOT traceless)'}", flush=True)

# Generate Lie algebra via commutators.
G = [M_kit(p) for p in [15,30,45,60]]
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print(f"\nStarting rank: {mrank(G)} (4 mats)", flush=True)
# Seed.
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
print(f"After seed: {mrank(G)} ({len(G)} mats)", flush=True)
# Grow.
r = mrank(G)
cnt = 0
while r < 8 and cnt < 30:
    cnt += 1
    found = False
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            comm = G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm]) > r:
                G.append(comm); r = mrank(G); found = True
                break
        if found: break
    if not found:
        print(f"Stalled at rank {r}.", flush=True)
        break
print(f"Final rank: {r} ({len(G)} mats)", flush=True)

if r == 8:
    # Check closure and Killing.
    G8 = []
    for g in G:
        if len(G8) >= 8: break
        if mrank(G8+[g]) > mrank(G8):
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
        print(f"Closure: {'YES (28)' if closed else 'NO'}", flush=True)
        if closed:
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
            print(f"Killing: ({pos},{neg},{z})", flush=True)
            if pos==5 and neg==3:
                print("RESULT: sl(3,R) from Kit's FW matrix. PROVED.", flush=True)
            else:
                print(f"RESULT: 8-dim Lie algebra, not sl(3,R).", flush=True)
    except Exception as e:
        print(f"Singular: {e}", flush=True)
elif r < 8:
    print(f"RESULT: Only rank {r}, not sl(3,R).", flush=True)
else:
    print(f"RESULT: Rank {r} > 8, unexpected.", flush=True)

with open('/root/e95_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
