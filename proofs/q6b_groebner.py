from sage.all import *
# T4 (Q6b exact): Prove NO X in intertwiner space with X^2 = -I.
# X = sum c_i B_i. Equations: (X^2 + I)[i,j] = 0. If Groebner = [1], PROVED.

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
print("=== Polynomial system: X^2 = -I ===")
pr = PolynomialRing(QQ, 16, 'c', order='degrevlex')
c = pr.gens()
Xsym = sum(c[i]*B[i] for i in range(16))
Xsq = Xsym*Xsym
eqns = [(Xsq + I16).list()[k] for k in range(256)]
eqns = [e for e in eqns if e != 0]
seen=set(); uniq=[]
for e in eqns:
    s=str(e)
    if s not in seen:
        seen.add(s); uniq.append(e)
print(f"Unique nonzero equations: {len(uniq)}")

print()
print("=== Groebner basis ===")
try:
    I = pr.ideal(uniq)
    gb = I.groebner_basis()
    print(f"Groebner basis: {len(gb)} elements")
    if len(gb)==1 and gb[0]==1:
        print("T4 PROVED: Ideal is [1] => NO X with X^2=-I exists in intertwiner space.")
        print("No invariant complex structure via intertwiners. (Matches Q6 numerical.)")
        with open('/root/t4_q6b_exact.txt','w') as f:
            f.write("proved True\nnonexistence True\n")
    else:
        print("Ideal is proper => solutions MAY exist (unexpected).")
        print("T4: INCOMPLETE (or FALSIFIED if explicit solution found)")
        with open('/root/t4_q6b_exact.txt','w') as f:
            f.write("proved False\n")
    with open('/root/t4_groebner.txt','w') as f:
        f.write(f"gb_len {len(gb)}\nis_one {len(gb)==1 and gb[0]==1}\n")
except Exception as ex:
    print(f"Groebner failed: {ex}")
    print("T4: INCOMPLETE")
