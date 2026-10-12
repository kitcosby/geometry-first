from sage.all import *
# E104: Parameter theorem exploration.
# Test: for T_a skeleton, which diagonal functions (d1,d2,d3) give sl(3,R)?
# Vary d1,d2,d3 as multiples of urx, or constants, etc.
# Find the condition.

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

def M_param(deg, d1_fn, d2_fn, d3_fn):
    # d1_fn etc are functions of (srx,sxp,cxp,crx,urx,uxp) returning diagonal.
    srx,sxp,cxp,crx,urx,uxp = ops_exact(deg)
    d1 = d1_fn(srx,sxp,cxp,crx,urx,uxp)
    d2 = d2_fn(srx,sxp,cxp,crx,urx,uxp)
    d3 = d3_fn(srx,sxp,cxp,crx,urx,uxp)
    # Traceless check: d1+d2+d3 should be 0.
    # T_a skeleton: [[d1, srx, 0],[cxp, d2, uxp],[0, sxp, d3]]? 
    # Wait, original T_a: [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]].
    # Let's use: [[d1, srx, 0],[cxp, d2, crx],[0, sxp, d3]]? 
    # Simpler: keep off-diagonals as srx,cxp,sxp,crx, vary diagonal.
    M = matrix(K, 3, 3, [
        [d1, srx, 0],
        [cxp, d2, crx],
        [0, sxp, d3],
    ])
    return M

def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

def test_diagonal(name, d1_fn, d2_fn, d3_fn):
    G = [M_param(p, d1_fn, d2_fn, d3_fn) for p in [15,30,45,60]]
    # Check traceless.
    tr_ok = all((g.trace()==0) for g in G)
    r0 = mrank(G)
    # Add commutators.
    for i in [1,2,3]:
        G.append(G[0]*G[i]-G[i]*G[0])
    r1 = mrank(G)
    r=r1; cnt=0
    while r<8 and cnt<20:
        cnt+=1; found=False
        for i in range(len(G)):
            for j in range(i+1,len(G)):
                comm=G[i]*G[j]-G[j]*G[i]
                if mrank(G+[comm])>r:
                    G.append(comm); r=mrank(G); found=True; break
            if found: break
        if not found: break
    # Closure check (simplified: just report rank).
    print(f"{name}: traceless={tr_ok}, rank4={r0}, final={r}", flush=True)
    return r==8

print("E104: Diagonal parameter exploration.", flush=True)
# Case 1: Original (urx, -urx, 0).
test_diagonal("urx,-urx,0",
    lambda srx,sxp,cxp,crx,urx,uxp: urx,
    lambda srx,sxp,cxp,crx,urx,uxp: -urx,
    lambda srx,sxp,cxp,crx,urx,uxp: K(0))
# Case 2: (1,-1,0) constant.
test_diagonal("1,-1,0",
    lambda srx,sxp,cxp,crx,urx,uxp: K(1),
    lambda srx,sxp,cxp,crx,urx,uxp: K(-1),
    lambda srx,sxp,cxp,crx,urx,uxp: K(0))
# Case 3: (srx, cxp, -(srx+cxp)).
test_diagonal("srx,cxp,-(srx+cxp)",
    lambda srx,sxp,cxp,crx,urx,uxp: srx,
    lambda srx,sxp,cxp,crx,urx,uxp: cxp,
    lambda srx,sxp,cxp,crx,urx,uxp: -(srx+cxp))
# Case 4: (0,0,0) - should FAIL (zero diagonal).
test_diagonal("0,0,0",
    lambda srx,sxp,cxp,crx,urx,uxp: K(0),
    lambda srx,sxp,cxp,crx,urx,uxp: K(0),
    lambda srx,sxp,cxp,crx,urx,uxp: K(0))
# Case 5: (urx, uxp, -(urx+uxp)).
test_diagonal("urx,uxp,-(urx+uxp)",
    lambda srx,sxp,cxp,crx,urx,uxp: urx,
    lambda srx,sxp,cxp,crx,urx,uxp: uxp,
    lambda srx,sxp,cxp,crx,urx,uxp: -(urx+uxp))
print("E104 done.", flush=True)
with open('/root/e104_done.txt','w') as fh:
    fh.write("done\n")
