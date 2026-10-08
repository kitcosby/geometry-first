from sage.all import *
# Hunt for the natural bridge: Spin(4,4) generators + spinor-vector intertwiner.
# Goal: find the algebraic object that makes operator<->physics mapping natural.

print("=== Building Cl(4,4) (same as T1) ===")
I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
print(f"8 generators E_mu, 16x16 over QQ")

print()
print("=== 28 Spin(4,4) generators S_mn = [E_m,E_n]/2 ===")
S = {}
for mu in range(8):
    for nu in range(mu+1,8):
        S[(mu,nu)] = (E[mu]*E[nu] - E[nu]*E[mu])/2
print(f"Built {len(S)} generators")

# Find which S_mn are diagonalizable over QQ, and their spectra
print()
print("=== Spectra (looking for Cartan-like / coordinate-like) ===")
diag_count = 0
for key, M in S.items():
    # characteristic polynomial
    cp = M.charpoly()
    # check if splits over QQ with distinct linear factors (diagonalizable)
    # just report the trace and whether it's nilpotent-ish
    tr = M.trace()
    # S_mn are antisymmetric in the Clifford sense; check M^2 scalar?
    M2 = M*M
    # is M^2 a scalar matrix?
    is_scalar = all(M2[i,j]==0 for i in range(16) for j in range(16) if i!=j)
    if is_scalar:
        d = M2[0,0]
        print(f"S_{key}: M^2 = {d}*I (scalar!), tr={tr}")
        diag_count += 1
print(f"{diag_count}/28 have scalar square")

print()
print("=== Maximal commuting subset (Cartan hunt) ===")
# Greedy: find largest set of mutually commuting S_mn
keys = list(S.keys())
# Build commutation graph
import itertools
commute = {}
for a,b in itertools.combinations(keys,2):
    c = S[a]*S[b] - S[b]*S[a]
    commute[(a,b)] = c.is_zero()

# Greedy maximal commuting set
chosen = [keys[0]]
for k in keys[1:]:
    if all(commute.get(tuple(sorted((k,c))), False) or commute.get(tuple(sorted((c,k))), False) for c in chosen):
        # check actual commutation (handle ordering)
        ok = True
        for c in chosen:
            cc = S[k]*S[c] - S[c]*S[k]
            if not cc.is_zero():
                ok=False; break
        if ok:
            chosen.append(k)
print(f"Greedy maximal commuting set size: {len(chosen)}")
print(f"Keys: {chosen}")

print()
print("=== Intertwiner: 16x16 -> 8 (spinor bilinear to vector) ===")
# The vector rep of SO(4,4) is 8-dim. We want explicit maps v: 16x16 -> 8
# such that v(S_mn . M) = (vector action) . v(M).
# Simpler: find 8 linearly independent 16x16 matrices V_a such that
# [S_mn, V_a] = sum_b (vector_rep_mn)_ab V_b.
# This is the statement that {V_a} spans a submodule isomorphic to the vector.
# We can find it via: the E_mu themselves transform as vectors!
# Check: [S_mn, E_rho] = ?
# For Clifford: [E_mu E_nu, E_rho] gives the vector action.
print("Testing if {E_rho} spans the vector representation:")
print("(i.e., [S_mn, E_rho] = eta_nu_rho E_mu - eta_mu_rho E_nu)")
# Metric: first 4 are +1, last 4 are -1 (for Cl(4,4))
eta = [1]*4 + [-1]*4
ok_count=0; tot=0
for (mu,nu), Smn in S.items():
    for rho in range(8):
        tot+=1
        lhs = Smn*E[rho] - E[rho]*Smn
        # RHS: eta_nu_rho E_mu - eta_mu_rho E_nu, but S_mn = [E_mu,E_nu]/2
        # [E_mu E_nu, E_rho] = E_mu{E_nu,E_rho} - {E_mu,E_rho}E_nu
        # = E_mu(2 eta_nu_rho) - (2 eta_mu_rho)E_nu = 2(eta_nu_rho E_mu - eta_mu_rho E_nu)
        # So [S_mn, E_rho] = eta_nu_rho E_mu - eta_mu_rho E_nu.
        rhs = (eta[nu] if nu==rho else 0)*E[mu] - (eta[mu] if mu==rho else 0)*E[nu]
        # Actually eta is diagonal: eta_nu_rho = eta[nu]*delta_nu_rho
        rhs = (eta[nu] if nu==rho else 0)*E[mu] - (eta[mu] if mu==rho else 0)*E[nu]
        if (lhs - rhs).is_zero():
            ok_count+=1
print(f"Vector intertwiner via {{E_rho}}: {ok_count}/{tot} correct")
if ok_count==tot:
    print("PROVED: {E_0..E_7} is an explicit SO(4,4)-vector intertwiner (16x16 matrices -> 8-dim vector).")
    print("This IS the natural bridge: spinor bilinears map to vectors via the E_mu themselves.")

with open('/root/bridge_hunt.txt','w') as f:
    f.write(f"scalar_square {diag_count}/28\ncartan_greedy {len(chosen)}\nvector_intertwiner {ok_count}/{tot}\n")
print("WROTE /root/bridge_hunt.txt")
