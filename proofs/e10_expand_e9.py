from sage.all import *
# Expand E9: explicit U, and the 16 sectors from 4 commuting involutions.

I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
S = lambda mu,nu: (E[mu]*E[nu] - E[nu]*E[mu])/2

T = S(0,4)*S(1,5) - S(0,5)*S(1,4)
U = T/2
volW = E[0]*E[1]*E[4]*E[5]
H1 = S(0,1); H2 = S(4,5)

print("=== Explicit U in Cl(2,2) basis ===")
Widx = [0,1,4,5]
basis = [I16]
bas_names = ["I"]
for mu in Widx:
    basis.append(E[mu]); bas_names.append(f"E{mu}")
for i,mu in enumerate(Widx):
    for nu in Widx[i+1:]:
        basis.append(E[mu]*E[nu]); bas_names.append(f"E{mu}E{nu}")
for i,mu in enumerate(Widx):
    for j,nu in enumerate(Widx):
        for k,rh in enumerate(Widx):
            if i<j<k:
                basis.append(E[mu]*E[nu]*E[rh]); bas_names.append(f"E{mu}E{nu}E{rh}")
basis.append(E[0]*E[1]*E[4]*E[5]); bas_names.append("volW")
B = matrix(QQ, [M.list() for M in basis]).transpose()
v = vector(QQ, U.list())
c = B.solve_right(v)
print("U = sum c_i * basis_i:")
for i,ci in enumerate(c):
    if ci != 0:
        print(f"  {ci} * {bas_names[i]}")
print()

print("=== 16 sectors from (volW, H1, H2, U) ===")
# P(e1,e2,e3,e4) = [(1+e1*volW)/2][(1+e2*H1)/2][(1+e3*H2)/2][(1+e4*U)/2]
import itertools
sectors = {}
for e1,e2,e3,e4 in itertools.product([1,-1], repeat=4):
    P = ((I16 + e1*volW)/2) * ((I16 + e2*H1)/2) * ((I16 + e3*H2)/2) * ((I16 + e4*U)/2)
    tr = P.trace()
    # tr = dimension of sector (since P is a projector)
    sectors[(e1,e2,e3,e4)] = tr
    if tr != 0:
        print(f"  ({e1:+d},{e2:+d},{e3:+d},{e4:+d}): dim {int(tr)}")
print()
nonzero = sum(1 for d in sectors.values() if d!=0)
print(f"Nonzero sectors: {nonzero}/16")
total = sum(sectors.values())
print(f"Total dim: {int(total)} (expect 16)")

print()
print("=== Are the 4 involutions independent? ===")
# If all 16 sectors nonzero and dim 1 each, they're independent and maximal.
# If fewer, there are relations.
print(f"If 16 sectors of dim 1: independent, maximal.")
print(f"If fewer: relations exist.")
# Check pairwise products: is U = volW*H1*H2 or similar?
print()
print("Check relations:")
print(f"U == volW*H1*H2? {(U - volW*H1*H2).is_zero()}")
print(f"U == volW*H1? {(U - volW*H1).is_zero()}")
print(f"U == H1*H2? {(U - H1*H2).is_zero()}")
# U might be a product of the others, or independent.

with open('/root/e9_expand.txt','w') as f:
    f.write(f"nonzero_sectors {nonzero}\n")
    f.write(f"U_is_product {((U-volW*H1*H2).is_zero() or (U-volW*H1).is_zero() or (U-H1*H2).is_zero())}\n")
print("WROTE /root/e9_expand.txt")
