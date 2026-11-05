"""Entropy-selector self-consistency check.

Kit's claim (2026-10-10, flagged by him as untested -- "I haven't checked this
with a calculator"): applying the entropy gradient to the 'other' (non-selected)
time direction produces an absurd, reciprocal result -- "reciprocal under one so
only a negative feedback" -- i.e. the entropy selector is self-validating
because the wrong choice refutes itself.

Reading used here (stated so Kit can correct it):
  selected direction : U = urx + i*uxp, entropy S  = 1/urx
  "other" direction  : U_tilde = uxp + i*urx, entropy S_t = 1/uxp
  ("reciprocal under one" = 1/uxp; "negative feedback" = dS_t/dx <= 0)

Operator definitions (per proofs/q7_dSdx.py):
  srx = |csc x| + cot x,  sxp = |csc x| - cot x
  cxp = |sec x| + tan x,  crx = |sec x| - tan x
  urx = srx - crx,        uxp = cxp - sxp,   S = 1/urx

Result (symbolic, sympy):
  (1) sxp*srx = 1 and crx*cxp = 1                (reciprocal pairs)
  (2) uxp(x) = urx(pi/2 - x)                    (mirror identity, exact)
  (3) dS_t/dx(x) = -dS/dx(pi/2 - x)             (chain rule on (2))
  With E16 (dS/dx >= 0 on QI, proved): for x in (0, pi/2),
  dS_t/dx(x) = -(positive) <= 0.  The wrong choice yields exactly the
  negative of the selected direction's gradient at the mirror point.
Numeric corroboration (numpy, abs-correct branches):
  QI   : min dS/dx  = +0.51 > 0 ; max dS_t/dx = -0.51 < 0
  QIII : min dS/dx  = +0.51 > 0 ; max dS_t/dx = -0.51 < 0
  |uxp(x) - urx(pi/2-x)| < 4e-13 on both quadrants.

Scope: the math (dS_t/dx = -S'(pi/2-x) <= 0) is PROVED given E16.
The physical reading -- that entropy-decrease "refutes" the choice --
is Kit's postulate (ASSERTED), not established here.
"""
from sympy import *
import numpy as np

x = symbols('x', real=True)

# --- Symbolic part: QI branch (sin>0, cos>0), abs dropped ---
srx = 1/sin(x) + cos(x)/sin(x)
sxp = 1/sin(x) - cos(x)/sin(x)
cxp = 1/cos(x) + sin(x)/cos(x)
crx = 1/cos(x) - sin(x)/cos(x)
urx = srx - crx
uxp = cxp - sxp
S = 1/urx
St = 1/uxp

r1a = simplify(sxp*srx - 1)
r1b = simplify(crx*cxp - 1)
y = pi/2 - x
r2 = simplify(uxp - (srx.subs(x, y) - crx.subs(x, y)))
r3 = simplify(diff(St, x) + diff(S, x).subs(x, y))
print("symbolic: sxp*srx-1 =", r1a)
print("symbolic: crx*cxp-1 =", r1b)
print("symbolic: uxp(x)-urx(pi/2-x) =", r2)
print("symbolic: dSt/dx + dS/dy|(y=pi/2-x) =", r3)
sym_ok = (r1a == 0 and r1b == 0 and r2 == 0 and r3 == 0)
print("SYMBOLIC IDENTITIES EXACT:", sym_ok)

# --- Numeric part: abs-correct branches on QI and QIII ---
def ops(xv):
    s = np.sin(xv); c = np.cos(xv)
    srx_ = np.abs(1/s) + c/s; sxp_ = np.abs(1/s) - c/s
    cxp_ = np.abs(1/c) + s/c; crx_ = np.abs(1/c) - s/c
    return (srx_ - crx_), (cxp_ - sxp_)

num_ok = True
for name, a, b in [('QI', 0.02, np.pi/2 - 0.02),
                   ('QIII', np.pi + 0.02, 3*np.pi/2 - 0.02)]:
    xs = np.linspace(a, b, 2000)
    urx_, uxp_ = ops(xs)
    urx_m, _ = ops(np.pi/2 - xs)
    S_ = 1/urx_; St_ = 1/uxp_
    h = xs[1] - xs[0]
    dS = np.gradient(S_, h); dSt = np.gradient(St_, h)
    mirror_err = np.abs(uxp_ - urx_m).max()
    ok = dS.min() > 0 and dSt.max() < 0 and mirror_err < 1e-9
    num_ok = num_ok and ok
    print(f"{name}: min dS/dx={dS.min():.6f} (>0), "
          f"max dSt/dx={dSt.max():.6f} (<0), "
          f"mirror err={mirror_err:.2e} -> {'OK' if ok else 'FAIL'}")

print("NUMERIC CORROBORATION:", num_ok)
print("OVERALL:", "PASS" if (sym_ok and num_ok) else "FAIL")
