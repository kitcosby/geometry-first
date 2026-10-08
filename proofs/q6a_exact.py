from sage.all import *
# T3 (Q6a exact): Find X in intertwiner space with X^2 = +I.
# Strategy: get exact basis, try simple rational combinations, verify exactly.
# If simple search fails, attempt Groebner on reduced system.

print("=== Exact intertwiner basis ===")
I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
R = E[0:4]; J = E[4:8]
n=16
rows=[]
for mu in range(4):
    K = R[mu].transpose().tensor_product(I16) - I16.tensor_product(J[mu])
    rows.append(K)
M = block_matrix([[r] for r in rows])
ns = M.right_kernel()
basis = ns.basis()
print(f"Basis size: {len(basis)}")
# Convert vec to matrix
def vec_to_mat(v):
    return matrix(QQ, 16, 16, list(v))
B = [vec_to_mat(v) for v in basis]

print()
print("=== Trying single basis elements ===")
found=None
for i,Bi in enumerate(B):
    sq = Bi*Bi
    # Check if scalar multiple of I
    # sq = c*I  <=>  sq - (sq[0,0])*I == 0 and sq[0,0]!=0
    c = sq[0,0]
    if c!=0 and sq == c*I16:
        print(f"B[{i}]^2 = {c}*I")
        if c>0:
            # normalize: (Bi/sqrt(c))^2 = I. But sqrt(c) may be irrational.
            # Check if c is a rational square.
            cq = QQ(c)
            s = cq.sqrt(all=True)
            # s is list of roots (may be empty if not a square in QQ)
            print(f"  c={cq}, sqrt in QQ? {len(s)>0}")
            if len(s)>0:
                Xcand = Bi / s[0]
                if Xcand*Xcand == I16:
                    print(f"  PROVED: (B[{i}]/{s[0]})^2 = I")
                    found = Xcand
                    break
        elif c<0:
            print(f"  (negative: would give X^2=-I, not our target)")

if found is None:
    print("No single basis element works. Trying pairs...")
    # Try B_i + B_j for small i,j
    import itertools
    for i,j in itertools.combinations(range(len(B)),2):
        Xcand = B[i]+B[j]
        sq = Xcand*Xcand
        c = sq[0,0]
        if c!=0 and sq == c*I16 and c>0:
            cq=QQ(c); s=cq.sqrt(all=True)
            if len(s)>0:
                Xn = Xcand/s[0]
                if Xn*Xn==I16:
                    print(f"PROVED: ((B[{i}]+B[{j}])/{s[0]})^2 = I")
                    found=Xn
                    break
        if found: break

if found is None:
    print("Simple search failed. Attempting Groebner (may be heavy)...")
    # Parametrize X = sum c_i B_i. Impose (X^2 - I)[k,k]=0 for k=0..15 (diagonal, 16 eqns).
    # This is necessary but not sufficient; if no solution, then no involution.
    # If solution, verify full X^2=I.
    pr = PolynomialRing(QQ, 16, 'c')
    c = pr.gens()
    Xsym = sum(c[i]*B[i] for i in range(16))
    Xsq = Xsym*Xsym
    eqns = [(Xsq[k,k]-1) for k in range(16)]
    # Add a few off-diagonal to constrain
    for k in range(15):
        eqns.append(Xsq[k,k+1])
    print(f"System: {len(eqns)} equations, 16 vars, degree 2")
    try:
        I = pr.ideal(eqns)
        gb = I.groebner_basis()
        print(f"Groebner basis computed, {len(gb)} elements")
        if gb == [1] or (len(gb)==1 and gb[0]==1):
            print("INCONSISTENT: no solution to diagonal+offdiag (INCOMPLETE for full)")
        else:
            print("System may have solutions; full verification needed (INCOMPLETE)")
    except Exception as ex:
        print(f"Groebner failed: {ex} (INCOMPLETE)")
else:
    print()
    print("T3 PROVED: explicit involution found.")
    with open('/root/t3_q6a_exact.txt','w') as f:
        f.write("proved True\n")
