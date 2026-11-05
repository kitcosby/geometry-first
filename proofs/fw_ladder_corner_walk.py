# Corner-walk parity check for proofs/fw_ladder_induction.md (2026-10-16).
# Exact rational arithmetic (fractions.Fraction), desktop, no Sage needed.
# A1..A7 are the coefficient matrices of the Ladder-FW family in the
# basis {srx, crx, cxp, sxp, FW, cos 2x, sin 2x}, read off
# fw_ladder.py's build(): diag (urx = srx - crx, uxp = cxp - sxp, FW),
# path edges alternating (srx,cxp) odd / (crx,sxp) even, corners
# (1,n) = cos 2x, (n,1) = sin 2x, traceless correction -(tr/n) I with
# tr = urx + uxp + FW distributed by coefficient.
# Output (verified 2026-10-16 for n = 5..10):
#   [A6,A7] = E11 - Enn at every n.
#   Even n: [A6,A3] = E_{1,n-1} - E_{2,n} (the pair), [A6,A4] = 0;
#           [A7,A1] = -E_{n-1,1} + E_{n,1} + E_{n,2}, [A7,A2] = -E_{n,1}.
#   Odd n:  [A6,A3] = -E_{2,n} and [A6,A4] = +E_{1,n-1} (roots separate);
#           [A7,A1] = E_{n,1} + E_{n,2}, [A7,A2] = -E_{n-1,1} - E_{n,1}
#           (roots separate modulo the corner multiple E_{n,1} = A7).
from fractions import Fraction as F

def mat(n):
    return [[F(0)]*n for _ in range(n)]
def addE(A,i,j,v=F(1)): A[i][j]+=v

def build_A(n):
    A=[mat(n) for _ in range(7)]
    addE(A[0],0,0); addE(A[1],0,0,F(-1))          # urx = srx - crx
    addE(A[2],1,1); addE(A[3],1,1,F(-1))          # uxp = cxp - sxp
    if n>=3: addE(A[4],2,2)                        # FW
    for i in range(1,n):
        if i%2==1:
            addE(A[0],i-1,i); addE(A[2],i,i-1)
        else:
            addE(A[1],i-1,i); addE(A[3],i,i-1)
    addE(A[5],0,n-1); addE(A[6],n-1,0)             # corners
    for idx,c in [(0,F(-1,n)),(1,F(1,n)),(2,F(-1,n)),(3,F(1,n)),(4,F(-1,n))]:
        for i in range(n): A[idx][i][i]+=c         # -(tr/n) I by coefficient
    return A

def mul(X,Y):
    n=len(X); return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def br(X,Y):
    P,Q=mul(X,Y),mul(Y,X); return [[P[i][j]-Q[i][j] for j in range(len(X))] for i in range(len(X))]
def support(X):
    return {(i+1,j+1):X[i][j] for i in range(len(X)) for j in range(len(X)) if X[i][j]!=0}

for n in range(5,11):
    A=build_A(n)
    print(f"=== n={n}")
    print(" [A6,A7] =", support(br(A[5],A[6])))
    print(" [A6,A3] =", support(br(A[5],A[2])), " [A6,A4] =", support(br(A[5],A[3])))
    print(" [A7,A1] =", support(br(A[6],A[0])), " [A7,A2] =", support(br(A[6],A[1])))
    print(" [A6,A1] =", support(br(A[5],A[0])), " [A6,A2] =", support(br(A[5],A[1])), " [A6,A5] =", support(br(A[5],A[4])))
