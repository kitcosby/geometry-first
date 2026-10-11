from sage.all import *
# E51: Verify multi-rank Delta_{alpha,beta} Flatwave quotient formula.
# Doc claims:
#   Delta FW = DeltaN/D_beta - N_alpha*DeltaD/(D_alpha*D_beta)
# where N=srx+cxp, D=srx*cxp-1, DeltaN=N_a-N_b, DeltaD=D_a-D_b.
# Verify symbolically with t_a=tan(x_a/2), t_b=tan(x_b/2).

print("=== Abstract verification (pure algebra) ===")
# Let N_a, D_a, N_b, D_b be symbols. Verify identity.
Na, Da, Nb, Db = var('Na Da Nb Db')
LHS = Na/Da - Nb/Db
RHS = (Na-Nb)/Db - Na*(Da-Db)/(Da*Db)
diff = (LHS - RHS).simplify_full()
print(f"LHS - RHS = {diff}")
print(f"Identity holds: {diff==0}")
print()

print("=== Concrete verification (operator definitions) ===")
# t_a, t_b in (0,1) for QI.
ta, tb = var('ta tb')
assume(ta>0); assume(ta<1); assume(tb>0); assume(tb<1)
def N(t): return 1/t + (1+t)/(1-t)  # srx+cxp
def D(t): return (1/t)*((1+t)/(1-t)) - 1  # srx*cxp-1
Na_c=N(ta); Da_c=D(ta); Nb_c=N(tb); Db_c=D(tb)
FW_a=Na_c/Da_c; FW_b=Nb_c/Db_c
Delta_FW=(FW_a-FW_b).simplify_full()
RHS_c=((Na_c-Nb_c)/Db_c - Na_c*(Da_c-Db_c)/(Da_c*Db_c)).simplify_full()
check=(Delta_FW-RHS_c).simplify_full()
print(f"Delta_FW - RHS = {check}")
# Factor to see if zero
print(f"Simplified: {check.factor()}")
print(f"Identity holds: {check==0}")
print()
print("=== Resonant gate analysis ===")
# Doc: "(cxp*srx - 1) acts as resonant gate, spikes when product -> 1."
# D = srx*cxp - 1. When does D=0? i.e., srx*cxp=1.
# srx=1/ta, cxp=(1+ta)/(1-ta). Product = (1+ta)/(ta*(1-ta)).
# Set =1: (1+ta)=ta*(1-ta) => 1+ta=ta-ta^2 => ta^2+1=0. No real solution!
# So D != 0 for real ta. The "gate" never closes (no spike) in QI.
prod=(1/ta)*((1+ta)/(1-ta))
eq=solve(prod==1, ta)
print(f"srx*cxp=1 solves to: {eq}")
print("No real solution in (0,1). D never zero in QI. No resonant spike.")
print("(The gate is always open; Flatwave is smooth on QI.)")
print()
print("WROTE /root/e51_deltarank.txt")
with open('/root/e51_deltarank.txt','w') as f:
    f.write(f"Abstract identity: {diff==0}\n")
    f.write(f"Concrete identity: {check==0}\n")
    f.write("Resonant gate: no real zero in QI.\n")
