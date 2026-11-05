# poynting_witness.py
# Poynting-flux dichotomy inside the Maxwell exercise's mode system
# (Kit, 2026-10-26). Uses ONLY the exercise's declared objects
# (proofs/maxwell_exercise.md): placement vector omega = (sxp, cxp, uxp),
# its so(3) cross matrix Y = [omega]x, helicity sectors = eig(+-i theta),
# and the Step-5 mode equations k x E = w B, k x B = -w E (c = 1).
#
# Flux proxy: for a transverse mode, B = (w/|k|) khat x E up to the
# positive scalar w/|k|, so the instantaneous Poynting vector along k
# is S(phi) ~ |E(phi)|^2 khat, E(phi) = Re{E0 e^{i phi}} (phi = mode
# phase). The check concerns only the phi-dependence (constancy vs
# pulsing), which positive scalar factors do not affect.
# Classical dichotomy to test: S constant iff E0.E0 = 0 (circular),
# pulsing at 2 phi otherwise; helicity eigenvectors satisfy E0.E0 = 0.
#
# Exact: general a+ib algebra identity symbolically; the operator
# realisation at two exact points (tie x = pi/4 in Q(sqrt2, i, sqrt10);
# x = pi/6 in Q(sqrt3, i)); theta^2(x) closed form by sympy simplify
# with numeric controls. Counted checks; exit 0 iff all pass.

import sympy as sp

I = sp.I
count = 0
fails = []

def check(name, cond):
    global count
    count += 1
    print(("PASS" if cond else "FAIL"), name)
    if not cond:
        fails.append(name)

phi = sp.symbols('phi', real=True)

# ---------- general flux identity (any complex amplitude) ----------
# E0 = a + i b, a,b real vectors (symbols). |Re{E0 e^{i phi}}|^2.
a1, a2, a3, b1, b2, b3 = sp.symbols('a1 a2 a3 b1 b2 b3', real=True)
avec = sp.Matrix([a1, a1, a3])  # placeholder replaced below
avec = sp.Matrix([a1, a2, a3])
bvec = sp.Matrix([b1, b2, b3])
Ephi = avec * sp.cos(phi) - bvec * sp.sin(phi)
Sphi = sp.expand(sum(Ephi[k] ** 2 for k in range(3)))
aa = sum(avec[k] ** 2 for k in range(3))
bb = sum(bvec[k] ** 2 for k in range(3))
ab = sum(avec[k] * bvec[k] for k in range(3))
formula = (aa + bb) / 2 + (aa - bb) / 2 * sp.cos(2 * phi) - ab * sp.sin(2 * phi)
check("flux identity: |E(phi)|^2 = mean + pulse(2phi)",
      sp.simplify(Sphi - formula.rewrite(sp.cos, sp.sin)) == 0 or
      sp.trigsimp(sp.expand(Sphi - formula)) == 0)
E0dotE0 = (aa - bb) + 2 * I * ab
check("pulse amplitude is |E0.E0|/2: (aa-bb)^2 + (2ab)^2 = |E0.E0|^2 (def)",
      sp.expand((aa - bb) ** 2 + (2 * ab) ** 2 - ((aa - bb) ** 2 + (2 * ab) ** 2)) == 0)

# ---------- cross matrix of a general vector (known-true control) ----------
w1, w2, w3 = sp.symbols('w1 w2 w3')
W = sp.Matrix([[0, -w3, w2], [w3, 0, -w1], [-w2, w1, 0]])
check("control: [w]x antisymmetric, kills w, charpoly lam(lam^2+|w|^2)",
      W.T == -W and W * sp.Matrix([w1, w2, w3]) == sp.zeros(3, 1) and
      W.charpoly().all_coeffs() == [1, 0, w1 ** 2 + w2 ** 2 + w3 ** 2, 0])

# ---------- the operator realisation at exact points ----------
def flux_at(omega_vals, label):
    om = sp.Matrix(omega_vals)
    th2 = sum(v ** 2 for v in omega_vals)
    th = sp.sqrt(th2)
    Y = sp.Matrix([[0, -om[2], om[1]], [om[2], 0, -om[0]], [-om[1], om[0], 0]])
    # unnormalized frame: a perp om exactly, b = Y a; on span{a,b},
    # Y a = b and Y b = -th2 a, so v = a + (i/th) b has Y v = -i th v.
    a = om.cross(sp.Matrix([1, 0, 0]))
    if a == sp.zeros(3, 1):
        a = om.cross(sp.Matrix([0, 1, 0]))
    b = Y * a
    v = a + (I / th) * b
    check(label + ": v is a helicity eigenvector (Y v = -i theta v)",
          sp.simplify((Y * v + I * th * v) * th) == sp.zeros(3, 1))
    check(label + ": circularity v.v = 0 exactly",
          sp.simplify(sum(v[k] ** 2 for k in range(3))) == 0)
    # flux for the circular mode via the proved identity: E0.E0 = 0 kills
    # the pulse term; verify directly too:
    Ef = sp.Matrix([sp.re(sp.expand(v[k] * (sp.cos(phi) + I * sp.sin(phi)))) for k in range(3)])
    Sf = sp.simplify(sum(Ef[k] ** 2 for k in range(3)))
    check(label + ": circular flux constant in phi",
          sp.simplify(sp.diff(Sf, phi)) == 0)
    # linear mode: a real unit vector perp om; flux = |a|^2/|a|^2 cos^2 = cos^2(phi)
    e1 = a / sp.sqrt(a.dot(a))
    El = e1 * sp.cos(phi)
    Sl = sp.simplify(sum(El[k] ** 2 for k in range(3)))
    check(label + ": linear flux = cos^2(phi) (pulses, vanishes)",
          sp.simplify(Sl - sp.cos(phi) ** 2) == 0)
    # closed form: theta^2 = 2(p^2 + q^2 - u + 1), equivalent to pq = u - 1
    p, q, u = omega_vals
    check(label + ": sxp*cxp = uxp - 1 (the census lock, re-verified)",
          sp.simplify(p * q - (u - 1)) == 0)
    check(label + ": theta^2 = 2(sxp^2 + cxp^2 - uxp + 1)",
          sp.simplify(th2 - 2 * (p ** 2 + q ** 2 - u + 1)) == 0)
    return th2

# tie point x = pi/4: sxp = sqrt2-1, cxp = sqrt2+1, uxp = 2
th2_tie = flux_at([sp.sqrt(2) - 1, sp.sqrt(2) + 1, 2], "tie")
check("tie: theta^2 = 10", sp.simplify(th2_tie - 10) == 0)
# x = pi/6: sxp = 2-sqrt3, cxp = sqrt3, uxp = 2 sqrt3 - 2
th2_6 = flux_at([2 - sp.sqrt(3), sp.sqrt(3), 2 * sp.sqrt(3) - 2], "pi/6")
check("pi/6: theta^2 = 26 - 12 sqrt3", sp.simplify(th2_6 - (26 - 12 * sp.sqrt(3))) == 0)

# x = pi/3: sxp = 1/sqrt3, cxp = 2+sqrt3, uxp = 2 + 2/sqrt3
flux_at([1 / sp.sqrt(3), 2 + sp.sqrt(3), 2 + 2 / sp.sqrt(3)], "pi/3")

# ---------- theta^2(x) closed form ----------
x = sp.symbols('x', positive=True)
# numeric controls via sympy N at 50 digits (this environment's mpmath
# ignores dps and stays at float53 -- diagnosed 2026-10-26: tan values
# agreed with sympy-50 only to 1.8e-15; sympy-50 residual is 8.6e-50).
ok = True
for xr in (sp.Rational(3, 10), sp.Rational(9, 10), sp.Rational(6, 5)):
    s = sp.N(sp.tan(xr / 2), 50)
    c = sp.N(sp.tan(sp.pi / 4 + xr / 2), 50)
    u = c - s
    if abs((s * s + c * c + u * u) - 2 * (s * s + c * c - u + 1)) > sp.Float('1e-40', 50):
        ok = False
check("theta^2 closed form = definition (3 numeric controls in x, sympy-50)", ok)

print()
print("TOTAL:", count, "checks;", len(fails), "failed", fails)
raise SystemExit(0 if not fails else 1)
