from sage.all import *
# E107: 2x2 with srx, cxp. Does it give sl(2,R)?
# M = [[0, srx],[cxp, 0]]. Test rank and closure.

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
    srx = (1+c)/s
    cxp = (1+s)/c
    return srx, cxp

def M_sl2(deg):
    srx, cxp = ops_exact(deg)
    M = matrix(K, 2, 2, [[0, srx],[cxp, 0]])
    return M

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print("E107: 2x2 [[0,srx],[cxp,0]] for sl(2,R).", flush=True)
G = [M_sl2(p) for p in [15,30,45,60]]
r0 = mrank(G)
print(f"Initial rank (4 samples): {r0} (sl2 needs 3)", flush=True)
# Commutators.
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
r1 = mrank(G)
print(f"After seed: {r1}", flush=True)
r=r1; cnt=0
while r<3 and cnt<10:
    cnt+=1; found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>r:
                G.append(comm); r=mrank(G); found=True; break
        if found: break
    if not found: break
print(f"Final rank: {r}", flush=True)
if r==3:
    print("Rank 3 achieved.", flush=True)
    # Closure: sl(2) has 3-dim, check 3 commutators.
    n=2; G3=[]
    for g in G:
        if len(G3)>=3: break
        if mrank(G3+[g])>mrank(G3):
            G3.append(g)
    Gflat=[[G3[k][a,b] for a in range(n) for b in range(n)] for k in range(3)]
    Mmat=matrix(K, Gflat); piv=Mmat.pivots()
    if len(piv)==3:
        Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(3)] for a in range(3)])
        Ninv=Nmat.inverse()
        closed=True
        for a in range(3):
            for b in range(a+1,3):
                comm=G3[a]*G3[b]-G3[b]*G3[a]
                dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(3)])
                cf=Ninv*dv
                rc=sum(cf[k]*G3[k] for k in range(3))
                if not (rc-comm).is_zero():
                    closed=False; break
            if not closed: break
        print(f"Closure: {'YES' if closed else 'NO'}", flush=True)
        if closed:
            # Killing for sl2: should be (2,1,0).
            f={}
            for a in range(3):
                for b in range(3):
                    comm=G3[a]*G3[b]-G3[b]*G3[a]
                    dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(3)])
                    fij=Ninv*dv
                    for k in range(3):
                        if fij[k]!=0: f[(a,b,k)]=fij[k]
            Kmat=matrix(K, 3, 3)
            for a in range(3):
                for b in range(a,3):
                    s=K(0)
                    for k in range(3):
                        for l in range(3):
                            s+=f.get((a,k,l),K(0))*f.get((b,l,k),K(0))
                    Kmat[a,b]=s; Kmat[b,a]=s
            ev=Kmat.eigenvalues()
            pos=sum(1 for e in ev if e>0); neg=sum(1 for e in ev if e<0)
            print(f"Killing: ({pos},{neg})", flush=True)
            if (pos,neg)==(2,1):
                print("RESULT: sl(2,R) PROVED from 2 functions.", flush=True)
else:
    print(f"RESULT: Rank {r}, not 3.", flush=True)
with open('/root/e107_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
