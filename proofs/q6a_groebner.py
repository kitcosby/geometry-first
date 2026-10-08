from sage.all import *
# T3 (Q6a exact, attempt 2): Find X in intertwiner space with X^2=+I via Groebner.
# X = sum c_i B_i. Equations: (X^2 - I)[i,j] = 0. 16 vars, 256 quadrics.

print("=== Exact intertwiner basis ===")
I2 = matrix(QQ, [[1,0],[0,1]])
Xm = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [Xm.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(Xm), Z.tensor_product(Y)]
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
def vec_to_mat(v):
    return matrix(QQ, 16, 16, list(v))
B = [vec_to_mat(v) for v in basis]
print(f"Basis: {len(B)} matrices")

print()
print("=== Polynomial system: X^2 = I ===")
pr = PolynomialRing(QQ, 16, 'c', order='degrevlex')
c = pr.gens()
Xsym = sum(c[i]*B[i] for i in range(16))
Xsq = Xsym*Xsym
# Use all 256 equations (X^2 - I == 0)
eqns = [(Xsq - I16).list()[k] for k in range(256)]
# Remove zero polynomials
eqns = [e for e in eqns if e != 0]
print(f"Nonzero equations: {len(eqns)}")
# Deduplicate (some may be identical)
seen=set(); uniq=[]
for e in eqns:
    s=str(e)
    if s not in seen:
        seen.add(s); uniq.append(e)
print(f"Unique equations: {len(uniq)}")

print()
print("=== Groebner basis (this may take a while) ===")
try:
    I = pr.ideal(uniq)
    gb = I.groebner_basis()
    print(f"Groebner basis: {len(gb)} elements")
    if len(gb)==1 and gb[0]==1:
        print("RESULT: Ideal is [1] => NO involution exists. (Unexpected! Q6 found one numerically.)")
        print("T3: FALSIFIED (no X^2=+I)")
    else:
        print("Ideal is proper => solutions may exist.")
        # Try to find one solution via triangular decomposition or saturation
        # For now, try sage's solve with a subset, or use variety over QQbar (may be slow)
        print("Attempting to find an explicit solution...")
        # Use elimination to get univariate polynomial, then solve
        # This is heuristic; may not complete.
        try:
            # Try variety over algebraic closure (can be very slow; limit)
            sols = I.variety(ring=QQbar)
            print(f"Found {len(sols)} solutions over QQbar")
            if sols:
                s0=sols[0]
                cc=[s0[pr.gen(i)] for i in range(16)]
                Xsol=sum(cc[i]*B[i] for i in range(16))
                # Verify X^2=I over QQbar
                if Xsol*Xsol == identity_matrix(QQbar,16):
                    print("T3 PROVED: explicit involution over QQbar, X^2=I verified.")
                    with open('/root/t3_q6a_exact.txt','w') as f:
                        f.write("proved True\n")
                else:
                    print("Solution does not verify (numerical issue). INCOMPLETE.")
        except Exception as ex:
            print(f"Variety computation failed: {ex}")
            print("T3: INCOMPLETE (Groebner done, explicit solution not extracted)")
    with open('/root/t3_groebner.txt','w') as f:
        f.write(f"gb_len {len(gb)}\nis_one {len(gb)==1 and gb[0]==1}\n")
except Exception as ex:
    print(f"Groebner failed: {ex}")
    print("T3: INCOMPLETE")
