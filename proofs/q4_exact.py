from sage.all import *
# T1 (Q4 exact): Commutant of Cl(4,4) 16-dim representation.
# PROVE: only scalars commute with all 8 generators (dim = 1).
# Method: exact linear algebra over QQ.

print("=== Building Cl(4,4) generators (exact, over QQ) ===")
# Use tensor products of 2x2 matrices to get 16x16.
# Signature (4,4): E_i^2 = +I for i<4, -I for i>=4 (adjust to match Q2's +,-,+,- pattern per sector).
I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])   # X^2=+I
Z = matrix(QQ, [[1,0],[0,-1]])  # Z^2=+I
Y = matrix(QQ, [[0,-1],[1,0]])  # Y^2=-I
W = matrix(QQ, [[0,1],[-1,0]])  # W^2=-I

# Cl(2,2) 4x4 generators (signature +,-,+,-): e0^2=+, e1^2=-, e2^2=+, e3^2=-
# Graded tensor: e0 = X⊗I, e1 = Y⊗I, e2 = Z⊗X, e3 = Z⊗Y. Check: 
# (Z⊗X)^2 = Z^2⊗X^2 = I⊗I = +I. Good.
# (Z⊗Y)^2 = Z^2⊗Y^2 = I⊗(-I) = -I. Good.
# Anticommutators: {X⊗I, Y⊗I} = {X,Y}⊗I = 0. Good.
# {X⊗I, Z⊗X} = XZ⊗X + ZX⊗X = {X,Z}⊗X = 0. Good.
# {Y⊗I, Z⊗X} = YZ⊗X + ZY⊗X = {Y,Z}⊗X = 0? {Y,Z} = YZ+ZY = [[0,-1],[1,0]]*[[1,0],[0,-1]] + [[1,0],[0,-1]]*[[0,-1],[1,0]] = [[0,1],[1,0]] + [[0,-1],[-1,0]] = 0. Good.
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]

print("Verifying Cl(2,2) relations:")
ok = True
for i in range(4):
    sq = e[i]*e[i]
    expected = I2.tensor_product(I2) if i in [0,2] else -I2.tensor_product(I2)
    if sq != expected:
        print(f"  e[{i}]^2 WRONG"); ok=False
for i in range(4):
    for j in range(i+1,4):
        if e[i]*e[j] + e[j]*e[i] != zero_matrix(QQ,4):
            print(f"  {{e[{i}],e[{j}]}} != 0"); ok=False
print("  Cl(2,2) OK" if ok else "  Cl(2,2) FAILED")

# Cl(4,4) 16x16: E_mu = e_mu ⊗ I4 (mu=0..3), E_{4+nu} = Omega ⊗ e_nu (nu=0..3)
# Omega = volume element of first Cl(2,2) = e0*e1*e2*e3 (should square to +I or -I; for the twist we need Omega^2=+I and {Omega, e_mu}=0).
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
print(f"Omega^2 = +I? {Omega*Omega == I4}")
# Check Omega anticommutes with each e_mu (required for graded tensor)
for mu in range(4):
    if Omega*e[mu] + e[mu]*Omega != zero_matrix(QQ,4):
        print(f"  {{Omega, e[{mu}]}} != 0; using e[0]*e[1]*e[2]*e[3] may need sign fix")

E = []
for mu in range(4):
    E.append(e[mu].tensor_product(I4))
for nu in range(4):
    E.append(Omega.tensor_product(e[nu]))

I16 = identity_matrix(QQ,16)
print("Verifying Cl(4,4) relations:")
ok=True
for i in range(8):
    sq = E[i]*E[i]
    # Signature: first 4 from e ( +,-,+,- ), last 4 from Omega⊗e_nu.
    # (Omega⊗e_nu)^2 = Omega^2 ⊗ e_nu^2.
    if sq != I16 and sq != -I16:
        print(f"  E[{i}]^2 not ±I"); ok=False
for i in range(8):
    for j in range(i+1,8):
        if E[i]*E[j] + E[j]*E[i] != zero_matrix(QQ,16):
            print(f"  {{E[{i}],E[{j}]}} != 0"); ok=False
print("  Cl(4,4) anticommutators OK" if ok else "  FAILED")
sqs = [1 if E[i]*E[i]==I16 else -1 for i in range(8)]
print(f"  Signature (E_i^2): {sqs}")

print()
print("=== Exact commutant computation ===")
# Solve [X, E_i] = 0 for 16x16 X (256 unknowns).
# Vectorize: [X,E_i]=0  <=>  (I⊗E_i - E_i^T⊗I) vec(X) = 0.
n=16
rows=[]
for i in range(8):
    # Kronecker: (I⊗E_i - E_i^T⊗I)
    K = I16.tensor_product(E[i]) - E[i].transpose().tensor_product(I16)
    rows.append(K)
M = block_matrix([[r] for r in rows])  # (8*256) x 256
print(f"System: {M.nrows()} equations, {M.ncols()} unknowns")
ns = M.right_kernel()
print(f"Commutant dimension (exact over QQ): {ns.dimension()}")
if ns.dimension()==1:
    print("PROVED: commutant is 1-dimensional (scalars only). Cl(4,4) 16 is irreducible.")
else:
    print(f"NOT 1: dim={ns.dimension()}")

with open('/root/t1_q4_exact.txt','w') as f:
    f.write(f"dim {ns.dimension()}\nproved {ns.dimension()==1}\n")
print("WROTE /root/t1_q4_exact.txt")
