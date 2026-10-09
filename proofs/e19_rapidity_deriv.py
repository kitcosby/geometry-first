from sage.all import *
# Build exp(phi * S01) and test derivatives.
# S01 is the boost generator (E0 timelike, E1 spacelike).

print("=== Rebuilding Cl(4,4) E_mu ===")
I2 = matrix(QQ, [[1,0],[0,1]])
X = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [X.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(X), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
print("E_mu built (8 generators, 16x16)")

# S01 = [E0,E1]/2. E0 is timelike (+), E1 is spacelike (-). So S01 is a BOOST.
S01 = (E[0]*E[1] - E[1]*E[0])/2
print(f"S01 built. S01^2 = scalar? Checking...")
S01_sq = S01*S01
# Check if scalar
is_scalar = all(S01_sq[i,j]==0 for i in range(16) for j in range(16) if i!=j)
print(f"  S01^2 is scalar: {is_scalar}, value: {S01_sq[0,0] if is_scalar else 'N/A'}")

# For a boost, S^2 = +I (or +scalar). Then exp(phi*S) = cosh(phi)*I + sinh(phi)*S.
# Let's verify and construct symbolically.
phi = var('phi')
if is_scalar:
    s = S01_sq[0,0]  # S01^2 = s*I
    print(f"  S01^2 = {s}*I")
    # exp(phi*S01) = cosh(phi*sqrt(s))*I + sinh(phi*sqrt(s))/sqrt(s) * S01
    # If s=1: exp = cosh(phi)*I + sinh(phi)*S01.
    # If s=-1: exp = cos(phi)*I + sin(phi)*S01 (rotation).
    if s == 1:
        print("  S01 is a BOOST (S^2=+I). exp(phi*S01) = cosh(phi)*I + sinh(phi)*S01.")
    elif s == -1:
        print("  S01 is a ROTATION (S^2=-I). exp(phi*S01) = cos(phi)*I + sin(phi)*S01.")
    else:
        print(f"  S01^2 = {s}*I, unusual.")

print()
print("=== Derivative w.r.t. rapidity ===")
print("d/dphi exp(phi*S01) = S01 * exp(phi*S01).")
print("This is a derivative w.r.t. the GROUP PARAMETER phi.")
print("It gives back the Lie algebra generator S01, not a spacetime derivative.")
print()
print("=== Does this 'allow derivatives'? ===")
print("NO for spacetime derivatives P_mu = -i d/dx^mu.")
print("The rapidity phi is a parameter, not a coordinate.")
print("S01 is a matrix (Lie algebra), not a differential operator.")
print("exp(phi*S01) is a matrix function of phi, not a field on spacetime.")
print()
print("To get P_mu (translations), we need:")
print("  (a) a manifold with coordinates x^mu, AND")
print("  (b) an action of translations on that manifold.")
print("The (2,2) gives (b) the Lorentz action, but not (a) the manifold.")
print("The E8 lattice did not supply (a) either; its absence does not create it.")
print()
print("CONCLUSION: The rapidity picture clarifies the BOOSTS, but does not")
print("by itself supply spacetime derivatives. P_mu remains OPEN (E18).")

with open('/root/rapidity_deriv.txt','w') as f:
    f.write("rapidity_gives_boosts_not_derivatives\n")
print("WROTE /root/rapidity_deriv.txt")
