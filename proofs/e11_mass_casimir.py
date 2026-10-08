from sage.all import *
# Can mass be derived? Look for Casimir invariants distinguishing the 4 sectors.
# If Spin(2,2) Casimirs differ across sectors, that's a geometric mass splitting.

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

print("=== Spin(2,2) generators ===")
Widx = [0,1,4,5]
S22 = [S(mu,nu) for i,mu in enumerate(Widx) for nu in Widx[i+1:]]
print(f"{len(S22)} generators")

print()
print("=== Quadratic Casimir C2 = sum_{a} eta^{aa} S_a^2 ===")
# For so(2,2), need the Killing form. Simpler: sum over all S_mn S^mn with metric.
# Actually, for the spinor rep, the Casimir acts as scalar (irrep).
# But our 16 is NOT irreducible under Spin(2,2) (it's 4 copies or 4 sectors).
# Compute C2 = sum_{mu<nu} sgn * S_mu_nu^2, where sgn from Killing.
# For so(p,q), Killing on S_mn: proportional to (eta_mm eta_nn).
# Let's just compute sum S_mn^2 and sum eta_mm eta_nn S_mn^2.
eta = {0:1, 1:-1, 4:1, 5:-1}
C2_naive = sum(M*M for M in S22)
print(f"sum S^2: scalar? {all(C2_naive[a,b]==0 for a in range(16) for b in range(16) if a!=b)}")
if all(C2_naive[a,b]==0 for a in range(16) for b in range(16) if a!=b):
    print(f"  = {C2_naive[0,0]}*I")
# With Killing signs:
C2_kill = sum(QQ(eta[mu]*eta[nu]) * S(mu,nu)*S(mu,nu) for i,mu in enumerate(Widx) for nu in Widx[i+1:])
print(f"Killing-weighted: scalar? {all(C2_kill[a,b]==0 for a in range(16) for b in range(16) if a!=b)}")
if all(C2_kill[a,b]==0 for a in range(16) for b in range(16) if a!=b):
    print(f"  = {C2_kill[0,0]}*I")

print()
print("=== Project onto 4 sectors, check Casimir per sector ===")
volW = E[0]*E[1]*E[4]*E[5]
H1 = S(0,1); H2 = S(4,5)
# Sectors from E10: (+,+,+), (+,-,-), (-,+,-), (-,-,+) in (volW,H1,H2), each dim 4.
import itertools
for e1,e2,e3 in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]:
    P = ((I16+e1*volW)/2)*((I16+e2*H1)/2)*((I16+e3*H2)/2)
    # Restrict C2 to this sector: P*C2*P. Is it scalar on the sector?
    # Compute trace(P*C2)/trace(P) = average eigenvalue.
    trP = P.trace()
    trPC = (P*C2_kill).trace()
    avg = trPC/trP if trP!=0 else None
    print(f"Sector volW={e1:+d},H1={e2:+d},H2={e3:+d}: dim={int(trP)}, <C2>={avg}")
print()
print("If <C2> differs across sectors: geometric mass splitting PROVED.")
print("If same: no mass from Casimir; mass needs external scale (wall).")

with open('/root/mass_casimir.txt','w') as f:
    f.write("casimir_check_done\n")
print("WROTE /root/mass_casimir.txt")
