from sage.all import *
# D3: Natural substructures of the Cl(4,4) 16.
# Q1: Spin(4,4)-invariant bilinear forms (symmetric? skew? dimensions?)
# Q2: Volume element action (scalar? Schur check)
# Q3: Any other natural structure.

print("=== Building Cl(4,4) generators ===")
I2 = matrix(QQ, [[1,0],[0,1]])
Xm = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [Xm.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(Xm), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
n=16

print()
print("=== Q2: Volume element (check first, quick) ===")
vol = E[0]
for i in range(1,8):
    vol = vol*E[i]
print(f"vol^2 = +I? {vol*vol == I16}")
print(f"vol^2 = -I? {vol*vol == -I16}")
# Schur: vol is central (8 gens, even), so vol = scalar. Which?
# Check if vol is scalar matrix
v00 = vol[0,0]
is_scalar = (vol == v00*I16)
print(f"vol is scalar? {is_scalar}, value: {v00}")
if is_scalar:
    print("=> No Z2-grading from volume element (Schur). 16 stays irreducible.")

print()
print("=== Q1: Spin(4,4)-invariant bilinear forms ===")
print("Spin Lie algebra spanned by S_mn = [E_m, E_n].")
print("Want M with S_mn^T M + M S_mn = 0 for all m<n.")
# Build the linear system for M (16x16 = 256 unknowns).
# For symmetric M: 136 unknowns. For skew: 120.
# We'll do general M first (256), then decompose.
S = {}
for m in range(8):
    for nn in range(m+1,8):
        S[(m,nn)] = E[m]*E[nn] - E[nn]*E[m]
print(f"Spin generators: {len(S)}")

# Vectorize: (S^T M + M S = 0) <=> (I (x) S^T + S^T (x) I) vec(M) = 0
# vec(S^T M) = (I (x) S^T) vec(M). vec(M S) = (S^T (x) I) vec(M).
# So: (I Kron S^T + S^T Kron I) vec(M) = 0.
rows=[]
for key in S:
    Smn = S[key]
    K = I16.tensor_product(Smn.transpose()) + Smn.transpose().tensor_product(I16)
    rows.append(K)
Msys = block_matrix([[r] for r in rows])
print(f"System: {Msys.nrows()} eqns, {Msys.ncols()} unknowns")
ns = Msys.right_kernel()
print(f"Dim of invariant bilinear forms (general): {ns.dimension()}")

# Now separate symmetric and skew.
# Symmetric: M^T = M. Skew: M^T = -M.
# We can intersect with the symmetric/skew subspaces.
# Symmetric subspace: vec(M) with M[i,j]=M[j,i]. Basis: E_{ij}+E_{ji}.
# Instead: compute nullspace, then project.
if ns.dimension() > 0:
    basis_ns = ns.basis()
    def vec_to_mat(v):
        return matrix(QQ,16,16,list(v))
    Mb = [vec_to_mat(v) for v in basis_ns]
    nsym=0; nskew=0
    for Mi in Mb:
        if Mi.transpose() == Mi:
            nsym+=1
        elif Mi.transpose() == -Mi:
            nskew+=1
    print(f"Symmetric invariant forms: {nsym}, Skew: {nskew}")
    # Note: basis may not align with sym/skew; this counts those that are purely sym/skew.
    # For exact dims, we'd intersect subspaces. Let's do it properly.
    print()
    print("Exact symmetric/skew dimensions via subspace intersection:")
    # Symmetric constraint: M - M^T = 0. This is linear in vec(M).
    # Build matrix for (M - M^T). vec(M - M^T) = (I - P) vec(M) where P is transpose permutation.
    # P vec(M) = vec(M^T). P is the commutation matrix.
    # Simpler: directly parametrize symmetric M (136 vars) and impose invariance.
    # Symmetric basis: for i<=j, S_{ij} = E_{ij}+E_{ji} (with E_{ii}).
    sym_basis=[]
    for i in range(16):
        for j in range(i,16):
            Mij = zero_matrix(QQ,16)
            Mij[i,j]=1
            if i!=j:
                Mij[j,i]=1
            sym_basis.append(Mij)
    # Express invariance as linear system in 136 coeffs.
    # For each Smn, Smn^T M + M Smn = 0. M = sum c_k S_k.
    # This is linear in c. Build matrix.
    nsk=len(sym_basis)
    # vec(Smn^T M + M Smn) = sum c_k vec(Smn^T S_k + S_k Smn)
    # Stack over all (m,n).
    erows=[]
    for key in S:
        Smn=S[key]
        # For each k, compute vec(Smn^T S_k + S_k Smn), as column k.
        cols=[]
        for Sk in sym_basis:
            V = Smn.transpose()*Sk + Sk*Smn
            cols.append(vector(QQ, list(V)))
        erows.append(block_matrix([[c.column()] for c in cols]).transpose())
        # Actually we want: sum_k c_k * vec(...) = 0. So matrix with columns vec(...).
    # Stack all (m,n) vertically
    Msym = block_matrix([[r] for r in erows])
    print(f"Symmetric system: {Msym.nrows()} eqns, {Msym.ncols()} unknowns")
    ns_sym = Msym.right_kernel()
    print(f"Dim symmetric invariant: {ns_sym.dimension()}")

    # Skew: 120 vars, A_{ij} = E_{ij}-E_{ji} for i<j.
    skew_basis=[]
    for i in range(16):
        for j in range(i+1,16):
            Aij = zero_matrix(QQ,16)
            Aij[i,j]=1; Aij[j,i]=-1
            skew_basis.append(Aij)
    erows2=[]
    for key in S:
        Smn=S[key]
        cols=[]
        for Ak in skew_basis:
            V = Smn.transpose()*Ak + Ak*Smn
            cols.append(vector(QQ, list(V)))
        erows2.append(block_matrix([[c.column()] for c in cols]).transpose())
    Mskew = block_matrix([[r] for r in erows2])
    print(f"Skew system: {Mskew.nrows()} eqns, {Mskew.ncols()} unknowns")
    ns_skew = Mskew.right_kernel()
    print(f"Dim skew invariant: {ns_skew.dimension()}")

with open('/root/d3_forms.txt','w') as f:
    f.write(f"general {ns.dimension()}\n")
    try:
        f.write(f"sym {ns_sym.dimension()}\nskew {ns_skew.dimension()}\n")
    except:
        pass
    f.write(f"vol_scalar {is_scalar}\nvol_val {v00}\n")
print()
print("WROTE /root/d3_forms.txt")
