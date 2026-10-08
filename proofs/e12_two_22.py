from sage.all import *
# The second (2,2): W' = span{E2,E3,E6,E7} (P- complement).
# How do W and W' interact?

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

print("=== W = span{E0,E1,E4,E5} (P+), W' = span{E2,E3,E6,E7} (P-) ===")
W = [0,1,4,5]
Wp = [2,3,6,7]
print("W signature:", [(mu, int((E[mu]*E[mu])[0,0])) for mu in W])
print("W' signature:", [(mu, int((E[mu]*E[mu])[0,0])) for mu in Wp])

print()
print("=== Spin(2,2)' generators (from W') ===")
S22p = [S(mu,nu) for i,mu in enumerate(Wp) for nu in Wp[i+1:]]
print(f"{len(S22p)} generators")
H1p = S(2,3); H2p = S(6,7)
print(f"H1'=S23, H2'=S67. [H1',H2']=0? {(H1p*H2p-H2p*H1p).is_zero()}")
print(f"H1'^2 scalar? {all((H1p*H1p)[a,b]==0 for a in range(16) for b in range(16) if a!=b)}, val={(H1p*H1p)[0,0]}")

print()
print("=== volW' ===")
volWp = E[2]*E[3]*E[6]*E[7]
print(f"volW'^2 scalar? {all((volWp*volWp)[a,b]==0 for a in range(16) for b in range(16) if a!=b)}, val={(volWp*volWp)[0,0]}")
print(f"tr(volW')={volWp.trace()}")

print()
print("=== How do W and W' interact? ===")
print("Do [Spin(2,2), Spin(2,2)'] commute?")
# Take H1=S01 (from W) and H1'=S23 (from W'). Do they commute?
H1 = S(0,1)
print(f"[H1,H1']=[S01,S23]=0? {(H1*H1p-H1p*H1).is_zero()}")
print(f"[H1,H2']=0? {(H1*H2p-H2p*H1).is_zero()}")
H2 = S(4,5)
print(f"[H2,H1']=0? {(H2*H1p-H1p*H2).is_zero()}")
print(f"[H2,H2']=0? {(H2*H2p-H2p*H2).is_zero()}")
print()
print("Do volW and volW' commute?")
print(f"[volW,volW']=0? {(volW*volWp-volWp*volW).is_zero()}")
print()
print("What is volW * volWp?")
# volW*volWp = (E0E1E4E5)(E2E3E6E7) = E0E1E2E3E4E5E6E7 (up to sign).
# This is the FULL volume element of Cl(4,4)!
vol_full = E[0]*E[1]*E[2]*E[3]*E[4]*E[5]*E[6]*E[7]
print(f"volW*volWp = +/- vol_full? {(volW*volWp - vol_full).is_zero() or (volW*volWp + vol_full).is_zero()}")
print(f"vol_full^2 = {(vol_full*vol_full)[0,0]} (scalar? {all((vol_full*vol_full)[a,b]==0 for a in range(16) for b in range(16) if a!=b)})")

with open('/root/two_22.txt','w') as f:
    f.write(f"W_Wp_commute_H {(H1*H1p-H1p*H1).is_zero()}\n")
    f.write(f"volW_volWp_commute {(volW*volWp-volWp*volW).is_zero()}\n")
print("WROTE /root/two_22.txt")
