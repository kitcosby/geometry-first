from sage.all import *
# E101: Wick rotation sl(3,R) -> su(3). Get (0,8).
# Take Chevalley basis: 3 antisymmetric (keep), 5 symmetric (multiply by i).
# Verify closure and Killing (0,8).

K = QQbar
i = K(-1).sqrt()
print(f"i**2 = {i**2}", flush=True)

# Chevalley basis for sl(3,R) (integer entries).
# Compact (antisymmetric):
A1 = matrix(K, 3, 3, [[0,1,0],[-1,0,0],[0,0,0]])
A2 = matrix(K, 3, 3, [[0,0,0],[0,0,1],[0,-1,0]])
A3 = matrix(K, 3, 3, [[0,0,1],[0,0,0],[-1,0,0]])
# Non-compact (symmetric traceless):
H1 = matrix(K, 3, 3, [[1,0,0],[0,-1,0],[0,0,0]])
H2 = matrix(K, 3, 3, [[0,0,0],[0,1,0],[0,0,-1]])
S1 = matrix(K, 3, 3, [[0,1,0],[1,0,0],[0,0,0]])
S2 = matrix(K, 3, 3, [[0,0,0],[0,0,1],[0,1,0]])
S3 = matrix(K, 3, 3, [[0,0,1],[0,0,0],[1,0,0]])

# Wick: multiply symmetric by i.
G = [A1, A2, A3, i*H1, i*H2, i*S1, i*S2, i*S3]
print(f"8 generators (3 real antisym + 5 imag sym).", flush=True)
# Check anti-Hermitian: X.conjugate_transpose() == -X?
# For QQbar, conjugate is .conjugate()? Might not exist. Check: (i*H1).H?
# Sage: M.H is conjugate transpose? Actually M.conjugate_transpose().
# For exact, let's just verify closure and Killing.
def mrank(mats):
    if not mats: return 0
    n=mats[0].nrows()
    M=matrix(K, [[m[a,b] for a in range(n) for b in range(n)] for m in mats])
    return M.rank()

print(f"Rank: {mrank(G)} (should be 8)", flush=True)
# Closure: check all 28 commutators in span.
n=3
Gflat=[[G[k][a,b] for a in range(n) for b in range(n)] for k in range(8)]
Mmat=matrix(K, Gflat); piv=Mmat.pivots()
print(f"Pivots: {len(piv)}", flush=True)
if len(piv)==8:
    Nmat=matrix(K, [[Gflat[k][piv[a]] for k in range(8)] for a in range(8)])
    Ninv=Nmat.inverse()
    closed=True
    for a in range(8):
        for b in range(a+1,8):
            comm=G[a]*G[b]-G[b]*G[a]
            dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(8)])
            cf=Ninv*dv
            rc=sum(cf[k]*G[k] for k in range(8))
            if not (rc-comm).is_zero():
                closed=False; break
        if not closed: break
    print(f"Closure: {'YES' if closed else 'NO'}", flush=True)
    if closed:
        # Killing. Note: for complex matrices, Killing is complex bilinear.
        # For su(3) as REAL algebra, we need to treat as real 8-dim.
        # But our matrices are complex; the REAL span is what matters.
        # Actually, the 8 matrices are linearly independent OVER QQbar (complex).
        # As a REAL Lie algebra, su(3) has real dimension 8.
        # Our G has 8 elements, but they're complex. The REAL span would be 16-dim.
        # Hmm. This is the issue: Wick rotation gives COMPLEX matrices,
        # but su(3) as a REAL algebra needs 8 REAL dimensions.
        # The 8 complex matrices A1,A2,A3,iH1,... are linearly independent over C,
        # but as a REAL vector space, {A1,A2,A3,iH1,...} spans... 
        # Actually, A1 is real, iH1 is imaginary. Together they span a REAL space?
        # A1 (real) and iH1 (imag) are independent over R? Yes, because one is real, one is imag.
        # So the REAL span of G is 8-dim (3 real + 5 imag = 8 real).
        # Good. Now Killing on this REAL space.
        # For Killing, use the complex bilinear form restricted? Actually for real Lie algebra,
        # Killing is real bilinear. Our structure constants are... complex?
        # [A1, iH1] = i[A1,H1]. [A1,H1] is real (bracket of real matrices).
        # So structure constants are either real or imag.
        # This is getting messy. Let's just compute the Killing as complex and check signature?
        # For su(3), the Killing form (as complex bilinear on sl(3,C)) restricted to su(3) (real form)
        # should be negative definite REAL.
        # Simpler: compute using real and imaginary parts separately?
        # Actually, let's just compute the 8x8 complex Killing and see.
        f={}
        for a in range(8):
            for b in range(8):
                comm=G[a]*G[b]-G[b]*G[a]
                dv=vector(K, [comm[piv[x]//n, piv[x]%n] for x in range(8)])
                fij=Ninv*dv
                for k in range(8):
                    if fij[k]!=0: f[(a,b,k)]=fij[k]
        Kmat=matrix(K, 8, 8)
        for a in range(8):
            for b in range(a,8):
                s=K(0)
                for k in range(8):
                    for l in range(8):
                        s+=f.get((a,k,l),K(0))*f.get((b,l,k),K(0))
                Kmat[a,b]=s; Kmat[b,a]=s
        # Eigenvalues may be complex. For su(3), they should be real negative.
        ev=Kmat.eigenvalues()
        # Count real negative.
        real_neg=sum(1 for e in ev if e.imag()==0 and e.real()<0)
        print(f"Eigenvalues: {ev}", flush=True)
        print(f"Real negative: {real_neg}/8", flush=True)
        if real_neg==8:
            print("RESULT: (0,8) su(3). PROVED via Wick.", flush=True)
        else:
            print("RESULT: Not (0,8).", flush=True)
with open('/root/e101_done.txt','w') as fh:
    fh.write("wick done\n")
