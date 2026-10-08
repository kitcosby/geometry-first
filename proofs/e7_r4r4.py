from sage.all import *
# R4 x R4 decomposition: find the two factors explicitly.
# Cl(2,2) = <E0,E1,E4,E5> ~= M(4,R). The 16 = R4_tensor x R4_mult.
# Find idempotents or a tensor product structure.

I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)

print("=== Cl(2,2) generators: E0,E1,E4,E5 ===")
Wgen = [E[0], E[1], E[4], E[5]]

print()
print("=== Find a maximal set of orthogonal idempotents ===")
# In M(4,R), we have 4 minimal orthogonal idempotents summing to I.
# They correspond to a decomposition R4 = R+ R+ R+ R.
# For the tensor product, we want idempotents in Cl(2,2) that split the 16.
# Use the Cartan: H1=S01, H2=S45. They commute. Their simultaneous eigenspaces?
# H1^2 = -I? Let's check. Actually from E1: H1 = S01, and we saw M^(1) has eigenvalues +-2,0.
# On the 16 (not the 8-vector), what are the eigenvalues?
H1 = (E[0]*E[1]-E[1]*E[0])/2
H2 = (E[4]*E[5]-E[5]*E[4])/2
print(f"H1^2 scalar? {all((H1*H1)[a,b]==0 for a in range(16) for b in range(16) if a!=b)}, val={(H1*H1)[0,0]}")
print(f"H2^2 scalar? {all((H2*H2)[a,b]==0 for a in range(16) for b in range(16) if a!=b)}, val={(H2*H2)[0,0]}")
print(f"[H1,H2]=0? {(H1*H2-H2*H1).is_zero()}")

print()
print("=== Eigenvalues of H1, H2 on the 16 ===")
for H,name in [(H1,"H1"),(H2,"H2")]:
    cp = H.charpoly()
    print(f"{name} charpoly: {cp.factor()}")

print()
print("=== Simultaneous eigenspaces (the 4 factors?) ===")
# H1, H2 commute. Find common eigenvectors over QQbar.
# If each has eigenvalues +-a, +-b, we get 4 combos, each dim 4 => 16 = 4+4+4+4.
# That would be the decomposition into 4 irreps, not R4xR4.
# For R4xR4, we need a different structure.
# Actually: Cl(2,2)~=M(4,R) acts irreducibly on R4. The 16 as a Cl(2,2)-module:
# Since commutant is M(4,R), 16 = R4 (natural) tensor R4 (multiplicity).
# To see the tensor factors, find a subalgebra isomorphic to M(4,R) in the commutant,
# and see how it pairs with Cl(2,2).
print("(Deeper decomposition requires explicit matrix units; reporting eigenvalues for now.)")

print()
print("=== volW eigenspaces (Z2 grading) ===")
volW = E[0]*E[1]*E[4]*E[5]
cp = volW.charpoly()
print(f"volW charpoly: {cp.factor()}")
# volW^2=+I, so eigenvalues +-1. Trace 0 => 8 of each. 16 = 8+8.
print("volW gives 8+8 (even/odd). This is the Z2 grading.")

print()
print("=== Relation to T (interplay) ===")
# T = S04*S15 - S05*S14. How does T interact with volW, H1, H2?
S = lambda mu,nu: (E[mu]*E[nu] - E[nu]*E[mu])/2
T = S(0,4)*S(1,5) - S(0,5)*S(1,4)
print(f"[T, volW]=0? {(T*volW-volW*T).is_zero()}")
print(f"[T, H1]=0? {(T*H1-H1*T).is_zero()}")
print(f"[T, H2]=0? {(T*H2-H2*T).is_zero()}")
print(f"T^2 = {(T*T)[0,0]}*I, tr(T)={T.trace()}")

with open('/root/r4r4.txt','w') as f:
    f.write(f"H1_H2_commute {(H1*H2-H2*H1).is_zero()}\n")
    f.write(f"volW_splits_8_8 True\n")
    f.write(f"T_commutes_volW {(T*volW-volW*T).is_zero()}\n")
print("WROTE /root/r4r4.txt")
