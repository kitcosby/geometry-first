from sage.all import *
# E93: sl(4,R) with non-zero v[2]. Test if zeros are necessary.
# X4(x) with v = [S2, C2, a], a != 0. Try a=1.
# Check: rank 15, closure (105), Killing (9,6).

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

def X4_exact_a(deg, a):
    # a = v[2], non-zero.
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    A3 = matrix(K, 3, 3, [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]])
    S2 = 1/(urx+uxp); C2 = (srx+crx-sxp-cxp)/(urx+uxp)
    X4 = matrix(K, 4, 4)
    for i in range(3):
        for j in range(3):
            X4[i,j] = A3[i,j]
    # v = [S2, C2, a].
    X4[0,3] = S2; X4[1,3] = C2; X4[2,3] = a
    X4[3,0] = S2; X4[3,1] = C2; X4[3,2] = a
    # (4,4) = 0 (trace).
    X4[3,3] = 0
    return X4

def test_a(a_val, label):
    print(f"\n=== Testing a = {label} ===", flush=True)
    G = [X4_exact_a(p, a_val) for p in [15,30,45,60]]
    # Check traceless.
    for g in G:
        if g.trace() != 0:
            print(f"  NOT TRACELESS. Fail.", flush=True)
            return False
    print(f"  Traceless OK.", flush=True)
    # Grow to 15.
    def mrank(mats):
        if not mats: return 0
        n = mats[0].nrows()
        M = matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
        return M.rank()
    # Seed commutators.
    for i in [1,2,3]:
        G.append(G[0]*G[i]-G[i]*G[0])
    A01 = G[0]*G[1]-G[1]*G[0]
    # Grow.
    r = mrank(G)
    cnt = 0
    while r < 15 and cnt < 50:
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
            print(f"  Stalled at rank {r}.", flush=True)
            break
    print(f"  Rank {r} ({len(G)} mats).", flush=True)
    if r < 15:
        print(f"  FAIL: rank < 15.", flush=True)
        return False
    # Select 15.
    G15 = []
    for g in G:
        if len(G15) >= 15: break
        if mrank(G15+[g]) > mrank(G15):
            G15.append(g)
    # Closure.
    n = 4
    Gflat = [[G15[k][a,b] for a in range(n) for b in range(n)] for k in range(15)]
    Mmat = matrix(K, Gflat); piv = Mmat.pivots()
    Nmat = matrix(K, [[Gflat[k][piv[a]] for k in range(15)] for a in range(15)])
    try:
        Ninv = Nmat.inverse()
    except:
        print(f"  FAIL: singular.", flush=True)
        return False
    closed = True
    for i in range(15):
        for j in range(i+1, 15):
            comm = G15[i]*G15[j]-G15[j]*G15[i]
            dv = vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
            cf = Ninv*dv
            rc = sum(cf[k]*G15[k] for k in range(15))
            if not (rc-comm).is_zero():
                closed = False; break
        if not closed: break
    if not closed:
        print(f"  FAIL: not closed.", flush=True)
        return False
    print(f"  Closure OK (105).", flush=True)
    # Killing.
    f = {}
    for i in range(15):
        for j in range(15):
            comm = G15[i]*G15[j]-G15[j]*G15[i]
            dv = vector(K, [comm[piv[a]//n, piv[a]%n] for a in range(15)])
            fij = Ninv*dv
            for k in range(15):
                if fij[k] != 0:
                    f[(i,j,k)] = fij[k]
    Kmat = matrix(K, 15, 15)
    for i in range(15):
        for j in range(i, 15):
            s = K(0)
            for k in range(15):
                for l in range(15):
                    s += f.get((i,k,l), K(0)) * f.get((j,l,k), K(0))
            Kmat[i,j] = s; Kmat[j,i] = s
    ev = Kmat.eigenvalues()
    pos = sum(1 for e in ev if e > 0); neg = sum(1 for e in ev if e < 0); z = sum(1 for e in ev if e == 0)
    print(f"  Killing: ({pos},{neg},{z}).", flush=True)
    if pos == 9 and neg == 6 and z == 0:
        print(f"  SUCCESS: sl(4,R) with a={label}.", flush=True)
        return True
    else:
        print(f"  FAIL: not (9,6,0).", flush=True)
        return False

print("E93: Non-zero v[2] for sl(4).", flush=True)
# Try a=1.
test_a(K(1), "1")
# Try a=srx at 45deg? Actually a must be constant across x? No, a can be a function.
# For now, constant a=1. If works, try operator-based.
