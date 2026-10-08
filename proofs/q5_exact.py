from sage.all import *
# T2 (Q5 exact): Sector intertwiner dimension = 16.
# PROVE: dim{ X (16x16) : X*R_mu = J_mu*X for mu=0..3 } = 16,
# where R_mu = real-sector generators, J_mu = imaginary-sector generators.
# Method: exact linear algebra over QQ.

print("=== Building generators (same as T1) ===")
I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
R = E[0:4]  # real sector
J = E[4:8]  # imaginary sector

print("=== Exact intertwiner computation ===")
# X*R_mu = J_mu*X  <=>  (R_mu^T ⊗ I - I ⊗ J_mu) vec(X) = 0
n=16
rows=[]
for mu in range(4):
    K = R[mu].transpose().tensor_product(I16) - I16.tensor_product(J[mu])
    rows.append(K)
M = block_matrix([[r] for r in rows])
print(f"System: {M.nrows()} equations, {M.ncols()} unknowns")
ns = M.right_kernel()
d = ns.dimension()
print(f"Intertwiner dimension (exact over QQ): {d}")
if d==16:
    print("PROVED: intertwiner space is 16-dimensional (= M4(R)).")
else:
    print(f"NOT 16: dim={d}")

with open('/root/t2_q5_exact.txt','w') as f:
    f.write(f"dim {d}\nproved {d==16}\n")
print("WROTE /root/t2_q5_exact.txt")
