from sage.all import *
# E30: Analyze T_a(x) = [[urx, srx, 0],[cxp, -urx, uxp],[0, -uxp, 0]]
# Check: trace, determinant, characteristic polynomial, su(3) properties.

x = var('x')
# QI definitions
srx = (1+cos(x))/sin(x)
cxp = (1+sin(x))/cos(x)
urx = srx - cos(x)/(1+sin(x))  # srx - crx
uxp = cxp - sin(x)/(1+cos(x))  # cxp - sxp

# Build matrix
T = matrix([[urx, srx, 0],
            [cxp, -urx, uxp],
            [0, -uxp, 0]])
print("T_a(x) =")
print(T)
print()

tr = T.trace().simplify_full()
print(f"Trace: {tr}")
print(f"  Traceless: {tr==0}")
print()

det = T.determinant().simplify_full()
print(f"Determinant: {det}")
print()

cp = T.characteristic_polynomial()
print(f"Characteristic polynomial (in l):")
print(cp)
print()

# su(3) generators are 3x3 traceless. Check if T is anti-Hermitian (for su(3))
# or Hermitian (for Gell-Mann convention).
# T is real, so Hermitian = symmetric, anti-Hermitian = antisymmetric.
# T_12=srx, T_21=cxp. Symmetric needs srx=cxp. Antisymmetric needs cxp=-srx.
# Neither holds generally.
print("Symmetry check:")
print(f"  T12=srx, T21=cxp. Symmetric (srx=cxp)? Only at x=45deg.")
print(f"  Antisymmetric (cxp=-srx)? No (both positive in QI).")
print()
print("CONCLUSION:")
print("T_a(x) is 3x3 traceless (trace=0).")
print("It is NOT (anti-)Hermitian in general, so not a standard su(3) generator.")
print("The sparsity pattern [*,*,0; *,*,*; 0,*,0] is notable.")
print("Determinant and char poly computed above.")

with open('/root/ta.txt','w') as f:
    f.write("ta_analyzed=True\n")
print("WROTE /root/ta.txt")
