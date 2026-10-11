from sage.all import *
# E48: uxp-based q-derivative.
# Idea: use q = uxp(x_now) as the deformation parameter.
# D_q f(x) = (f(q*x) - f(x)) / ((q-1)*x).
# At x_now=85deg, uxp~22, so q~22. Test on f(x)=x^2, f(x)=srx(x).

import math

def uxp_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9 or abs(c)<1e-9: return None
    return (1+s)/c - 1/((1+c)/s)

def srx_phase(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    if abs(s)<1e-9: return None
    return abs(1/s)+c/s

x_now=85.0
q=uxp_phase(x_now)
print(f"x_now={x_now}deg, q=uxp(x_now)={q:.4f}")
print()

def Dq(f, x, q, h_scale=1.0):
    # D_q f(x) = (f(q*x)-f(x))/((q-1)*x). x in degrees, scale q*x also in degrees.
    # Need q*x within domain (0,90) for srx.
    qx=q*x
    fx=f(x); fqx=f(qx) if qx<89.9 else None
    if fx is None or fqx is None: return None
    if abs((q-1)*x)<1e-12: return None
    return (fqx-fx)/((q-1)*x)

# Test f(x)=x^2 (x in radians for derivative check)
# D_q x^2 = (q^2 x^2 - x^2)/((q-1)x) = (q+1)x. As q->1, ->2x (correct).
print("Test f(x)=x (x in degrees):")
print("  D_q x should be 1 (since (qx - x)/((q-1)x) = 1).")
for x in [10, 30, 50]:
    print(f"  x={x}: D_q={Dq(lambda d: d, x, q)}")
print()
print("Test f(x)=x^2:")
print("  D_q x^2 = (q+1)*x.")
for x in [10, 30]:
    dq=Dq(lambda d: d*d, x, q)
    expected=(q+1)*x
    print(f"  x={x}: D_q={dq:.4f}, expected={expected:.4f}, match={abs(dq-expected)<1e-6}")
print()
# Test on srx (our operator). q*x may exceed 90deg, so use small x.
print("Test f=srx (operator):")
print("  Note: q*x must stay <90deg, so x < 90/q.")
xmax=90/q
print(f"  q={q:.2f}, xmax={xmax:.2f}deg")
for x in [1, 2, 3]:
    if x>=xmax: continue
    dq=Dq(srx_phase, x, q)
    # Ordinary derivative for comparison (numerical)
    h=0.01
    d_ord=(srx_phase(x+h)-srx_phase(x-h))/(2*h)
    print(f"  x={x}deg: D_q(srx)={dq:.4f}, d/dx(srx)={d_ord:.4f}")
print()
print("Interpretation:")
print("  D_uxp uses q=uxp(x_now)~22. The 'lattice spacing' is large (q>>1).")
print("  As x_now->0, uxp->1, q->1, D_q -> d/dx (continuous limit).")
print("  As x_now->90deg, uxp->inf, q->inf, highly discrete.")
print("WROTE /root/e48_uxp_q.txt")
