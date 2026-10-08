from sage.all import *
# T4b: Extract explicit solutions to X^2=-I and check if any are REAL.
print("=== Rebuilding system ===")
I2 = matrix(QQ, [[1,0],[0,1]])
Xm = matrix(QQ, [[0,1],[1,0]])
Z = matrix(QQ, [[1,0],[0,-1]])
Y = matrix(QQ, [[0,-1],[1,0]])
e = [Xm.tensor_product(I2), Y.tensor_product(I2), Z.tensor_product(Xm), Z.tensor_product(Y)]
I4 = identity_matrix(QQ,4)
Omega = e[0]*e[1]*e[2]*e[3]
E = [e[mu].tensor_product(I4) for mu in range(4)] + [Omega.tensor_product(e[nu]) for nu in range(4)]
I16 = identity_matrix(QQ,16)
R = E[0:4]; J = E[4:8]
rows=[]
for mu in range(4):
    K = R[mu].transpose().tensor_product(I16) - I16.tensor_product(J[mu])
    rows.append(K)
M = block_matrix([[r] for r in rows])
ns = M.right_kernel()
basis = ns.basis()
def vec_to_mat(v):
    return matrix(QQ, 16, 16, list(v))
B = [vec_to_mat(v) for v in basis]

pr = PolynomialRing(QQ, 16, 'c', order='degrevlex')
c = pr.gens()
Xsym = sum(c[i]*B[i] for i in range(16))
Xsq = Xsym*Xsym
eqns = [(Xsq + I16).list()[k] for k in range(256)]
eqns = [e for e in eqns if e != 0]
seen=set(); uniq=[]
for e in eqns:
    s=str(e)
    if s not in seen:
        seen.add(s); uniq.append(e)
I = pr.ideal(uniq)
print("Computing variety over QQbar...")
sols = I.variety(ring=QQbar)
print(f"Found {len(sols)} solutions over QQbar")
nreal=0
for idx,s0 in enumerate(sols[:10]):  # check first 10
    cc=[s0[pr.gen(i)] for i in range(16)]
    # Check if all coefficients are real (imag part zero)
    is_real = all(abs(complex(ci).imag) < 1e-10 for ci in cc)
    if is_real:
        nreal+=1
        print(f"Solution {idx}: REAL coefficients")
        # Verify X^2=-I with real matrix
        Xreal = sum(float(complex(cc[i]).real)*B[i] for i in range(16))
        # Check X^2 = -I numerically
        err = (Xreal*Xreal + I16).norm()
        print(f"  ||X^2+I|| = {err:.2e}")
    else:
        print(f"Solution {idx}: complex (non-real) coefficients")
print()
print(f"Real solutions among first 10: {nreal}")
if nreal>0:
    print("T4 FALSIFIED: real complex structure exists in intertwiner space.")
else:
    print("No real solutions found among checked. T4 likely holds (no real J), but not yet PROVED.")
    print("(Would need real algebraic geometry to prove no real points.)")
with open('/root/t4_solutions.txt','w') as f:
    f.write(f"nsols {len(sols)}\nnreal_first10 {nreal}\n")
