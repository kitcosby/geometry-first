from sage.all import *
# D3b: Fixed version. Get invariant forms explicitly.
print("=== Generators ===")
I2 = matrix(QQ, [[1,0],[0,1]])
Xm = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [Xm.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(Xm), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)

print("=== Volume element ===")
vol = prod(E)
print(f"vol^2=+I? {vol*vol==I16}, vol^2=-I? {vol*vol==-I16}")
v00=vol[0,0]
print(f"vol scalar? {vol==v00*I16}, val={v00}")

print()
print("=== Spin-invariant forms (fixed) ===")
Slist=[]
for m in range(8):
    for nn in range(m+1,8):
        Slist.append(E[m]*E[nn]-E[nn]*E[m])
print(f"{len(Slist)} spin generators")
# General M: vec(S^T M + M S)=0.
rows=[]
for Smn in Slist:
    K = I16.tensor_product(Smn.transpose()) + Smn.transpose().tensor_product(I16)
    rows.append(K)
Msys = block_matrix([[r] for r in rows])
ns = Msys.right_kernel()
print(f"General dim: {ns.dimension()}")
Mb=[matrix(QQ,16,16,list(v)) for v in ns.basis()]
for i,Mi in enumerate(Mb):
    sym = (Mi.transpose()==Mi)
    sk = (Mi.transpose()==-Mi)
    print(f"  M[{i}]: symmetric={sym}, skew={sk}, rank={Mi.rank()}, det_zero={Mi.det()==0}")

with open('/root/d3b.txt','w') as f:
    f.write(f"gendim {ns.dimension()}\nvol_scalar {vol==v00*I16}\nvol_val {v00}\n")
print("WROTE /root/d3b.txt")
