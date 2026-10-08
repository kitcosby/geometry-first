from sage.all import *
# Follow the algebra: how do E_mu act on the 4 sectors?
# Sectors from (H1,H2): (++,--,+-,-+), each dim 4.
# Do E_mu preserve sectors, or swap them?

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

H1 = S(0,1); H2 = S(4,5)
volW = E[0]*E[1]*E[4]*E[5]

print("=== 4 sectors (H1,H2) ===")
# Projectors
Ppp = ((I16+H1)/2)*((I16+H2)/2)
Pmm = ((I16-H1)/2)*((I16-H2)/2)
Ppm = ((I16+H1)/2)*((I16-H2)/2)
Pmp = ((I16-H1)/2)*((I16+H2)/2)
print(f"dim(++,--,+-,-+): {int(Ppp.trace())},{int(Pmm.trace())},{int(Ppm.trace())},{int(Pmp.trace())}")

print()
print("=== How does E0 act? ===")
# E0 * Ppp * E0: does it stay in (++) or go elsewhere?
# Compute E0*Ppp: take a vector in (++), apply E0, see which sector it's in.
# Simpler: compute P_q * E0 * P_p for all p,q. Nonzero => E0 maps p->q.
sectors = {"++":Ppp, "--":Pmm, "+-":Ppm, "-+":Pmp}
for mu in [0,1,4,5]:
    print(f"E{mu}:")
    for pn,Pp in sectors.items():
        for qn,Pq in sectors.items():
            M = Pq * E[mu] * Pp
            # Check if zero
            if not M.is_zero():
                # What's the rank? (dimension of image)
                print(f"  {pn} -> {qn}: nonzero (rank {M.rank()})")
    print()

with open('/root/sector_action.txt','w') as f:
    f.write("sector_action_computed\n")
print("WROTE /root/sector_action.txt")
