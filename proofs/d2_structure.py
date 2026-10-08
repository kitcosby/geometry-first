from sage.all import *
# D2: Structure of the 16-dim intertwiner space V.
# Q1: Is V closed under matrix multiplication? (Is it an algebra?)
# Q2: What is the Gram matrix (trace inner product)?
# Q3: Explicit V ~= M4(R) as vector spaces.

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
print(f"dim V = {len(B)}")

def is_intertwiner(X):
    # Check X*R_mu == J_mu*X for all mu
    for mu in range(4):
        if X*R[mu] != J[mu]*X:
            return False
    return True

print()
print("=== Q1: Closure under multiplication ===")
# Check if B_i * B_j is in V for all i,j
# First: is identity in V? (Quick negative test)
print(f"I16 in V? {is_intertwiner(I16)} (expect False: R != J)")
closed=True
fail=[]
for i in range(16):
    for j in range(16):
        if not is_intertwiner(B[i]*B[j]):
            closed=False
            fail.append((i,j))
            break
    if not closed:
        break
print(f"V closed under multiplication? {closed}")
if not closed:
    print(f"First failure: B[{fail[0][0]}]*B[{fail[0][1]}] not in V")
    print("=> V is NOT an algebra under 16x16 matrix multiplication.")
    print("   (M4(R) isomorphism is vector-space only, not algebra.)")

print()
print("=== Q2: Gram matrix G[i,j] = tr(B_i^T B_j) ===")
G = matrix(QQ, 16, 16, [[(B[i].transpose()*B[j]).trace() for j in range(16)] for i in range(16)])
print(f"Rank of Gram: {G.rank()} (expect 16 if basis is independent)")
# Is Gram a scalar multiple of identity? (orthonormal up to scale?)
# Check G / G[0,0]
c = G[0,0]
print(f"G[0,0] = {c}")
# Check if G is diagonal
offdiag = sum(1 for i in range(16) for j in range(16) if i!=j and G[i,j]!=0)
print(f"Nonzero off-diagonal entries: {offdiag}")
if offdiag==0:
    diag_vals = set(G[i,i] for i in range(16))
    print(f"G is diagonal. Distinct diagonal values: {diag_vals}")

print()
print("=== Q3: What IS V? Summary ===")
print("V = Hom_{Cl(2,2)}(R, J), dim 16.")
print("As vector space: V ~= M4(R) (via R=4S, J=4S, S=4-dim irred).")
print(f"As algebra under 16x16 mult: {'CLOSED' if closed else 'NOT CLOSED'}.")
print("Natural bilinear form: <X,Y> = tr(X^T Y).")
# Check: does V contain any nonzero idempotent? Any X with X^2=X?
# (Just for structure; not required)
print()
print("=== Additional: trace of basis elements ===")
traces = [B[i].trace() for i in range(16)]
print(f"tr(B_i): {traces[:8]}...")

with open('/root/d2_structure.txt','w') as f:
    f.write(f"closed {closed}\n")
    f.write(f"gram_rank {G.rank()}\n")
    f.write(f"offdiag {offdiag}\n")
print()
print("WROTE /root/d2_structure.txt")
