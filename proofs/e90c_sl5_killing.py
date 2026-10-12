from sage.all import *
# E90c: Killing signature for exact sl(5). Target (14,10,0).
# Rebuilds the 24 generators, computes 24x24 Killing exactly.

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
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    A3 = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    S2 = 1/(urx+uxp); C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    X4 = matrix(K, 4, 4)
    for a in range(3):
        for b in range(3):
            X4[a,b]=A3[a,b]
    X4[0,3]=S2; X4[1,3]=C2; X4[3,0]=S2; X4[3,1]=C2
    S4 = 8*S2*C2; C4 = C2*C2 - 16*S2*S2
    X5 = matrix(K, 5, 5)
    for a in range(4):
        for b in range(4):
            X5[a,b]=X4[a,b]
    X5[0,4]=S4; X5[1,4]=C4; X5[4,0]=S4; X5[4,1]=C4
    return X5

print("E90c: Killing for sl(5).", flush=True)
G = [X5_exact(p) for p in [15,30,45,60]]
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
A01=G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

while mrank(G)<24:
    found=False
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            comm=G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm])>mrank(G):
                G.append(comm); found=True; break
        if found: break
    if not found: break

G24=[]
for g in G:
    if len(G24)>=24: break
    if mrank(G24+[g])>mrank(G24):
        G24.append(g)
print(f"Got {len(G24)} generators.", flush=True)

# Structure constants.
n=5
Gflat=[[G24[k][a,b] for a in range(n) for b in range(n)] for k in range(24)]
Mmat=matrix(K, Gflat); piv=Mmat.pivots()
Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(24)] for a in range(24)])
Ninv=Nmat.inverse()
print("Computing f_ijk...", flush=True)
f={}
for i in range(24):
    for j in range(24):
        comm=G24[i]*G24[j]-G24[j]*G24[i]
        dv=vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(24)])
        fij=Ninv*dv
        for k in range(24):
            if fij[k]!=0:
                f[(i,j,k)]=fij[k]
    if i%6==5:
        print(f"  {i+1}/24.", flush=True)
print(f"Non-zero: {len(f)}.", flush=True)

print("Killing matrix (24x24)...", flush=True)
Kmat=matrix(K, 24, 24)
for i in range(24):
    for j in range(i,24):
        s=K(0)
        for k in range(24):
            for l in range(24):
                s+=f.get((i,k,l),K(0))*f.get((j,l,k),K(0))
        Kmat[i,j]=s; Kmat[j,i]=s
    if i%6==5:
        print(f"  row {i+1}/24.", flush=True)

print("Eigenvalues (exact)...", flush=True)
ev=Kmat.eigenvalues()
pos=sum(1 for e in ev if e>0)
neg=sum(1 for e in ev if e<0)
z=sum(1 for e in ev if e==0)
print(f"Signature: ({pos},{neg},{z}).", flush=True)
if pos==14 and neg==10 and z==0:
    print("PROVED: (14,10). This is sl(5,R).", flush=True)
else:
    print(f"Expected (14,10,0).", flush=True)
with open('/root/e90c_done.txt','w') as fh:
    fh.write(f"sig=({pos},{neg},{z})\n")
print("WROTE /root/e90c_done.txt", flush=True)
