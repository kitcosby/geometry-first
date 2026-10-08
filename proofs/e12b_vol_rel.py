from sage.all import *
# Fix: define volW, check W-W' volume relation.

I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)

volW = E[0]*E[1]*E[4]*E[5]
volWp = E[2]*E[3]*E[6]*E[7]
vol_full = E[0]*E[1]*E[2]*E[3]*E[4]*E[5]*E[6]*E[7]

print(f"[volW,volWp]=0? {(volW*volWp-volWp*volW).is_zero()}")
print(f"volW*volWp = +vol_full? {(volW*volWp - vol_full).is_zero()}")
print(f"volW*volWp = -vol_full? {(volW*volWp + vol_full).is_zero()}")
print(f"vol_full^2 = {(vol_full*vol_full)[0,0]}")
print(f"tr(vol_full) = {vol_full.trace()}")
# vol_full charpoly?
cp = vol_full.charpoly()
print(f"vol_full charpoly: {cp.factor()}")
with open('/root/vol_rel.txt','w') as f:
    f.write("done\n")
print("WROTE /root/vol_rel.txt")
