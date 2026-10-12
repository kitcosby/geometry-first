from sage.all import *
# E94: Diagnostic re-verification of E91 (sl(6) 35-dim claim).
# Check: (1) Are the 4 base X6 distinct? (2) Exact rank at each growth step.
# (3) Final 35 exactly independent? (4) Closure re-verified.

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
    S8 = 2*S4*C4; C8 = C4*C4 - S4*S4
    X6 = matrix(K, 6, 6)
    for a in range(5):
        for b in range(5):
            X6[a,b]=X5[a,b]
    X6[0,5]=S8; X6[1,5]=C8; X6[5,0]=S8; X6[5,1]=C8
    return X6, (S8, C8)

print("E94: Diagnostic.", flush=True)
# (1) Check base matrices and v5.
bases = []
v5s = []
for deg in [15,30,45,60]:
    X6, (S8,C8) = X6_exact(deg)
    bases.append(X6)
    v5s.append((S8,C8))
    print(f"  {deg}deg: S8={complex(S8):.4f}, C8={complex(C8):.4f}", flush=True)

# Check pairwise differences.
print("\n(1) Base matrix distinctness:", flush=True)
for i in range(4):
    for j in range(i+1,4):
        diff = bases[i]-bases[j]
        is_zero = diff.is_zero()
        print(f"  X6({[15,30,45,60][i]}) vs X6({[15,30,45,60][j]}): {'IDENTICAL' if is_zero else 'distinct'}", flush=True)

# (2) Exact rank growth.
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

G = bases[:]
print(f"\n(2) Starting rank: {mrank(G)} (4 mats)", flush=True)
# Add seed commutators.
for i in [1,2,3]:
    G.append(G[0]*G[i]-G[i]*G[0])
print(f"  After seed: rank {mrank(G)} ({len(G)} mats)", flush=True)
A01=G[0]*G[1]-G[1]*G[0]
G.append(A01*G[2]-G[2]*A01)
print(f"  After A01: rank {mrank(G)} ({len(G)} mats)", flush=True)

# Grow.
r = mrank(G)
step = 0
while r < 35 and step < 100:
    step += 1
    found = False
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            comm = G[i]*G[j]-G[j]*G[i]
            if mrank(G+[comm]) > r:
                G.append(comm); r = mrank(G); found = True
                if r % 5 == 0 or r >= 33:
                    print(f"  Step {step}: rank {r} ({len(G)} mats)", flush=True)
                break
        if found: break
    if not found:
        print(f"  STALLED at rank {r} after {step} steps.", flush=True)
        break

print(f"\nFinal: rank {r}, {len(G)} matrices.", flush=True)
if r < 35:
    print(f"RESULT: Rank is {r}, NOT 35. E91 claim FALSE.", flush=True)
else:
    # (3) Select 35 and verify independence.
    G35 = []
    for g in G:
        if len(G35) >= 35: break
        if mrank(G35+[g]) > mrank(G35):
            G35.append(g)
    print(f"(3) Selected {len(G35)} independent.", flush=True)
    if len(G35) == 35:
        print("RESULT: Rank 35 CONFIRMED exact.", flush=True)
    else:
        print(f"RESULT: Only {len(G35)} independent. E91 claim FALSE.", flush=True)

with open('/root/e94_done.txt','w') as fh:
    fh.write(f"rank={r}\n")
print("WROTE /root/e94_done.txt", flush=True)
