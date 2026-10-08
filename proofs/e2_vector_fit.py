from sage.all import *
# Fix the vector intertwiner: directly fit [S_mn, E_rho] coefficients.

I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]

print("=== E_mu^2 (metric signature) ===")
for mu in range(8):
    M2 = E[mu]*E[mu]
    # check scalar
    d = M2[0,0] if all(M2[i,j]==0 for i in range(16) for j in range(16) if i!=j) else "nonscalar"
    print(f"E_{mu}^2 = {d}")

print()
print("=== Fitting [S_mn, E_rho] = sum_sigma C^{mn}_{rho,sigma} E_sigma ===")
S = {}
for mu in range(8):
    for nu in range(mu+1,8):
        S[(mu,nu)] = (E[mu]*E[nu] - E[nu]*E[mu])/2

# For each (mn, rho), express [S_mn, E_rho] in the E_sigma basis.
# Since E_sigma are linearly independent (check), solve for coefficients.
# Build the coordinate map: flatten E_sigma, solve linear system.
import itertools
vecs = [M.list() for M in E]  # 8 vectors in QQ^256
B = matrix(QQ, vecs).transpose()  # 256x8
# Check rank
print(f"rank of {{E_sigma}}: {B.rank()}")

# Fit one example to see the pattern
mn = (0,1); rho = 0
comm = S[mn]*E[rho] - E[rho]*S[mn]
v = vector(QQ, comm.list())
# Solve B*c = v
try:
    c = B.solve_right(v)
    print(f"[S_01, E_0] = {list(c)} (coeffs of E_0..E_7)")
except Exception as ex:
    print(f"solve failed: {ex}")

# Now fit ALL and tabulate the structure constants
print()
print("=== Structure constants C^{mn}_{rho -> sigma} ===")
# C[(mn,rho)] = coefficient vector
C = {}
for mn in S:
    for rho in range(8):
        comm = S[mn]*E[rho] - E[rho]*S[mn]
        v = vector(QQ, comm.list())
        c = B.solve_right(v)
        C[(mn,rho)] = c

# Check: is it of the form a*delta_{nu,rho} E_mu + b*delta_{mu,rho} E_nu?
# Print a few
for mn in [(0,1),(0,4),(4,5)]:
    mu,nu = mn
    print(f"S_{mn}:")
    for rho in range(8):
        c = C[(mn,rho)]
        nz = [(s, c[s]) for s in range(8) if c[s]!=0]
        print(f"  rho={rho}: {nz}")

with open('/root/vector_fit.txt','w') as f:
    for mn in sorted(S.keys()):
        for rho in range(8):
            c = C[(mn,rho)]
            nz = [(s,int(c[s])) for s in range(8) if c[s]!=0]
            f.write(f"{mn} {rho} {nz}\n")
print("WROTE /root/vector_fit.txt")
