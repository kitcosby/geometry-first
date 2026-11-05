#!/usr/bin/env python3
"""
prism_454590.py -- Kit's prism configuration, worked exactly.

Configuration (Kit, 2026-10-28): a 45-45-90 triangular prism, white
source with the beam PARALLEL TO THE HYPOTENUSE, refractive index set
to water.

Geometry: apex A = (0,1) (the 90-degree corner), B = (-1,0), C = (1,0).
The hypotenuse BC is the base (y = 0, length 2); legs AB, AC have
length sqrt(2). Beam direction d0 = (1, 0), entering leg AB. Because
the legs stand at 45 degrees to the hypotenuse, the entry incidence is
forced to exactly i1 = pi/4 -- the tie angle.

Method: (1) exact symbolic derivation for GENERAL n (mu = 1/n,
sin r1 = (sqrt2/2) mu, c = cos r1 with c^2 = 1 - mu^2/2), reduced
modulo the c-relation; (2) an independent numeric segment ray tracer
(Snell vector form, TIR by the critical test) with known-true and
known-false controls, per the symbolic-checker validation discipline.

Headline (proved below): CHANNEL A (entry low on the leg) emerges
exactly parallel to the beam for every n -- zero angular dispersion;
the prism works there as a total-internal-reflection mirror. CHANNEL B
(entry high on the leg) totally reflects off the far leg and exits
through the hypotenuse at incidence pi/4 - r1, with
sin i2 = n sin(pi/4 - r1): that channel DOES disperse (water F-C spread
~0.29 deg). Newton's classic fan (direct leg-to-leg passage) exists in
the same prism only at steep entry (> ~62 deg for water).

CORRECTION DISCLOSED: the first hand derivation of channel B claimed a
second TIR at the base and a (0,1) emergent ray. The independent
numeric tracer refuted it: the base incidence in channel B is
pi/4 - r1 (~13 deg for water), not pi/4 + r1 -- the quantity sc + mu/2
is cos(pi/4 - r1), not a sine. The symbolic identities P10-P12 of the
first pass were true algebra about a hypothetical reflected ray; the
TIR premise was false. Channel B below is the corrected treatment.
"""

import math
import sympy as sp

CHECKS = 0
FAILS = []


def check(name, ok):
    global CHECKS
    CHECKS += 1
    if not ok:
        FAILS.append(name)
        print(f"FAIL: {name}")
    else:
        print(f"ok:   {name}")


# ----------------------------------------------------------------------
# Part 1 -- exact symbolic derivation, general n
# ----------------------------------------------------------------------
mu, c = sp.symbols("mu c", positive=True)
s = sp.sqrt(2) / 2                      # sin(pi/4) = cos(pi/4)
C_REL = 1 - mu**2 / 2                   # c^2 = cos^2 r1 = 1 - sin^2 r1


def red(expr):
    """Reduce a polynomial in c modulo c^2 = 1 - mu^2/2, then simplify."""
    expr = sp.expand(expr)
    p = sp.Poly(expr, c)
    rem = sp.rem(p, sp.Poly(c**2 - C_REL, c))
    return sp.simplify(rem.as_expr())


def red_rat(expr):
    """Same for rational expressions: reduce the numerator."""
    num, den = sp.together(expr).as_numer_denom()
    return sp.simplify(red(num) / den)


# Refraction at AB (outward normal (-s, s), i1 = pi/4, sin r1 = s*mu):
# d1 = mu*d0 + (mu*s - c)*(-s, s)  (vector Snell), expanded by hand:
d1 = sp.Matrix([mu / 2 + s * c, mu / 2 - s * c])
check("P1 |d1|^2 = 1", red(d1.dot(d1) - 1) == 0)

# Hypotenuse (base, outward normal (0,-1)): incidence angle pi/4 + r1.
cos_base = s * c - mu / 2               # = cos(pi/4 + r1)
sin_base = s * c + mu / 2               # = sin(pi/4 + r1)
check("P2 base incidence: cos^2 + sin^2 = 1",
      red(cos_base**2 + sin_base**2 - 1) == 0)
# TIR at the base iff sin_base > mu; both sc, mu/2 > 0, so
# sin_base > mu  <=>  sc > mu/2  <=>  c^2 > mu^2/2:
check("P3 base TIR reduces to c^2 - mu^2/2 = 1 - mu^2 (TIR iff n > 1)",
      sp.expand(c**2 - mu**2 / 2).subs(c**2, C_REL) == 1 - mu**2)

# Channel boundary on AB: entry point (-u, 1-u); d1x + d1y = mu exactly.
check("P4 d1x + d1y = mu", sp.simplify(d1[0] + d1[1] - mu) == 0)
# tau to AC = 2u/mu ; tau to base = (1-u)/(sc - mu/2). Equal at
# u* = mu/(2 s c):
u_star = mu / (2 * s * c)
diff = 2 * u_star / mu - (1 - u_star) / (s * c - mu / 2)
check("P5 channel boundary u* = mu/(2 s c)", red_rat(diff) == 0)

# ---- Channel A (u > u*): base reflection, exit through AC ----
d2 = sp.Matrix([d1[0], -d1[1]])         # reflect across the base line
n2 = sp.Matrix([s, s])                  # outward normal of AC
check("P6 channel A exit incidence: d2.n2 = c (incidence = r1)",
      red(d2.dot(n2) - c) == 0)
n = 1 / mu
# Exit refraction: sin i2 = n sin r1 = s, so i2 = pi/4, cos i2 = s.
d3 = n * d2 - (n * c - s) * n2
check("P7 channel A emergent = (1, 0): exactly parallel, all n",
      red_rat(d3[0] - 1) == 0 and red_rat(d3[1]) == 0)

# ---- Channel B (u < u*): TIR at AC, TIR at base, exit through AC ----
cos_ac_first = d1.dot(n2)
check("P8 channel B first AC incidence: d1.n2 = s*mu (incidence = pi/2 - r1)",
      sp.simplify(cos_ac_first - s * mu) == 0)
# TIR there iff cos r1 > mu, i.e. c^2 > mu^2:
check("P9 AC TIR reduces to c^2 - mu^2 = 1 - 3mu^2/2 (TIR iff n > sqrt(3/2))",
      sp.expand(c**2 - mu**2).subs(c**2, C_REL) == 1 - 3 * mu**2 / 2)
d2b = d1 - 2 * cos_ac_first * n2        # reflect across AC
check("P10 channel B reflected ray = (sc - mu/2, -sc - mu/2)",
      sp.simplify(d2b[0] - (s * c - mu / 2)) == 0
      and sp.simplify(d2b[1] + (s * c + mu / 2)) == 0)
# At the base the incidence is pi/4 - r1 (NOT pi/4 + r1):
# sin(pi/4 - r1) = s*c - s*sin r1 = sc - mu/2 (the tangential component),
# cos(pi/4 - r1) = s*c + s*sin r1 = sc + mu/2 (the normal component).
check("P11 channel B base incidence: sin = sc - mu/2, cos = sc + mu/2",
      sp.simplify((s * c - s * (s * mu)) - (s * c - mu / 2)) == 0
      and sp.simplify((s * c + s * (s * mu)) - (s * c + mu / 2)) == 0)
# For water-range indices n sin(pi/4 - r1) < 1 (tracer, Part 2), so the
# ray refracts OUT through the hypotenuse: exit angle i2 obeys the exact
# law sin i2 = n (sc - mu/2), direction angle i2 - pi/2 (n-dependent:
# this channel disperses). No closed n-free direction exists here.

# ----------------------------------------------------------------------
# Part 2 -- independent numeric ray tracer (controls + confirmation)
# ----------------------------------------------------------------------
SQ2 = math.sqrt(2.0)
FACES = [((-1.0, 0.0), (0.0, 1.0), (-1 / SQ2, 1 / SQ2)),    # AB
         ((0.0, 1.0), (1.0, 0.0), (1 / SQ2, 1 / SQ2)),      # AC
         ((1.0, 0.0), (-1.0, 0.0), (0.0, -1.0))]            # BC (base)


def refract_or_reflect(d, nl, eta):
    cosi = -(d[0] * nl[0] + d[1] * nl[1])
    sin2t = eta * eta * (1.0 - cosi * cosi)
    if sin2t > 1.0:
        return (d[0] + 2 * cosi * nl[0], d[1] + 2 * cosi * nl[1]), True
    cost = math.sqrt(1.0 - sin2t)
    return (eta * d[0] + (eta * cosi - cost) * nl[0],
            eta * d[1] + (eta * cosi - cost) * nl[1]), False


def trace(d0, entry_u, n_in):
    """Ray parallel-ish d0 hitting AB at (-u, 1-u); returns exit angle deg."""
    p = (-entry_u - 0.5 * d0[0], 1.0 - entry_u - 0.5 * d0[1])
    d = d0
    inside = False
    for _ in range(12):
        best = None
        for (a, b, nrm) in FACES:
            ex, ey = b[0] - a[0], b[1] - a[1]
            den = d[0] * ey - d[1] * ex
            if abs(den) < 1e-14:
                continue
            t = ((a[0] - p[0]) * ey - (a[1] - p[1]) * ex) / den
            w = ((a[0] - p[0]) * d[1] - (a[1] - p[1]) * d[0]) / den
            if t > 1e-9 and -1e-9 <= w <= 1 + 1e-9:
                if best is None or t < best[0]:
                    best = (t, nrm)
        if best is None:
            return None
        t, nrm = best
        p = (p[0] + t * d[0], p[1] + t * d[1])
        nl = nrm if d[0] * nrm[0] + d[1] * nrm[1] < 0 else (-nrm[0], -nrm[1])
        eta = (1.0 / n_in) if not inside else n_in
        d, tir = refract_or_reflect(d, nl, eta)
        if not tir:
            inside = not inside
            if not inside:
                return math.degrees(math.atan2(d[1], d[0]))
    return None


def channel_b_angle(nv):
    """Closed form: exit through hypotenuse, direction angle in degrees."""
    r1 = math.asin((SQ2 / 2) / nv)
    i2 = math.asin(nv * math.sin(math.pi / 4 - r1))
    return math.degrees(i2) - 90.0


aA = trace((1.0, 0.0), 0.8, 4 / 3)
aB = trace((1.0, 0.0), 0.3, 4 / 3)
check("T1 tracer, n=4/3, u=0.8 (channel A): exit angle 0 deg",
      aA is not None and abs(aA) < 1e-9)
check("T2 tracer, n=4/3, u=0.3 (channel B): exit = asin(n sin(pi/4-r1)) - 90",
      aB is not None and abs(aB - channel_b_angle(4 / 3)) < 1e-9)
aA2 = trace((1.0, 0.0), 0.8, 1.3330)
aB2 = trace((1.0, 0.0), 0.3, 1.3330)
check("T3 tracer, n=1.3330: channel A at 0, channel B at closed form",
      abs(aA2) < 1e-9 and abs(aB2 - channel_b_angle(1.3330)) < 1e-9)
a15 = trace((1.0, 0.0), 0.8, 1.5)
check("T4 known-true instance, n=1.5, u=0.8: exit angle 0 deg",
      a15 is not None and abs(a15) < 1e-9)
# Known-false control: tilt the beam 5 degrees off parallel (i1 = 50 deg);
# the parallel-exit claim must FAIL.
d50 = (math.cos(math.radians(-5.0)), math.sin(math.radians(-5.0)))
a50 = trace(d50, 0.8, 4 / 3)
check("T5 known-false control, beam 5 deg off parallel: exit NOT at 0 deg",
      a50 is not None and abs(a50) > 0.5)
# Channel B disperses: F and C exit angles differ materially.
bF = trace((1.0, 0.0), 0.3, 1.3371)
bC = trace((1.0, 0.0), 0.3, 1.3311)
check("T6 channel B dispersion: F-C exit-angle spread > 0.1 deg",
      bF is not None and bC is not None and abs(bF - bC) > 0.1)
# Channel boundary by bisection on the tracer (n = 4/3): channel A exits
# at 0 deg, channel B at ~-72.6 deg; classify B by angle < -45 deg.
lo, hi = 0.05, 0.95
for _ in range(40):
    mid = 0.5 * (lo + hi)
    ang = trace((1.0, 0.0), mid, 4 / 3)
    if ang is not None and ang < -45.0:
        lo = mid
    else:
        hi = mid
u_num = 0.5 * (lo + hi)
r1_43 = math.asin((SQ2 / 2) / (4 / 3))
u_star_43 = (3 / 4) / (2 * (SQ2 / 2) * math.cos(r1_43))
check("T7 channel boundary: tracer bisection matches u* = mu/(2 s c)",
      abs(u_num - u_star_43) < 2e-3)

# ----------------------------------------------------------------------
# Part 3 -- water numbers and the Newton (leg-to-leg) regime
# ----------------------------------------------------------------------
print("\nWater / index table (entry i1 = 45 deg exactly):")
print(f"{'n':>8} {'r1 deg':>9} {'crit deg':>9} {'base inc':>9} {'u*':>8}")
for nval, tag in [(4 / 3, "4/3 nominal"), (1.3311, "C 656nm"),
                  (1.3330, "d 589nm"), (1.3371, "F 486nm")]:
    r1 = math.degrees(math.asin((SQ2 / 2) / nval))
    crit = math.degrees(math.asin(1 / nval))
    print(f"{nval:8.4f} {r1:9.4f} {crit:9.4f} {45 + r1:9.4f} "
          f"{(1 / nval) / (2 * (SQ2 / 2) * math.cos(math.radians(r1))):8.5f}"
          f"   {tag}")

print("\nNewton leg-to-leg regime (same prism, steep entry, apex 90 deg):")
for nval, tag in [(1.3311, "C"), (1.3330, "d"), (1.3371, "F")]:
    i_min = math.degrees(math.asin(nval * SQ2 / 2))
    dev_min = 2 * i_min - 90.0
    # direct passage needs r1 > 90 - crit:
    r1_need = 90.0 - math.degrees(math.asin(1 / nval))
    i_need = math.degrees(math.asin(nval * math.sin(math.radians(r1_need))))
    print(f"  n={nval} ({tag}): min deviation {dev_min:.3f} deg at entry "
          f"{i_min:.2f} deg; direct passage needs entry > {i_need:.2f} deg")
dF = 2 * math.degrees(math.asin(1.3371 * SQ2 / 2)) - 90.0
dC = 2 * math.degrees(math.asin(1.3311 * SQ2 / 2)) - 90.0
print(f"  F-C angular spread at min deviation: {dF - dC:.3f} deg")

r1F = math.degrees(math.asin((SQ2 / 2) / 1.3371))
r1C = math.degrees(math.asin((SQ2 / 2) / 1.3311))
print(f"\nInside the prism the fan exists: r1(F)={r1F:.4f} vs r1(C)={r1C:.4f} "
      f"deg (spread {abs(r1F - r1C):.4f} deg). Channel A cancels it exactly "
      f"at exit; channel B keeps part of it:")
for nval, tag in [(1.3311, "C"), (1.3330, "d"), (1.3371, "F")]:
    print(f"  channel B exit (n={nval}, {tag}): "
          f"{channel_b_angle(nval):.4f} deg from +x axis")
print(f"  channel B F-C spread: "
      f"{abs(channel_b_angle(1.3371) - channel_b_angle(1.3311)):.4f} deg")

print(f"\n{CHECKS - len(FAILS)}/{CHECKS} checks passed")
if FAILS:
    print("FAILURES:", FAILS)
    raise SystemExit(1)
print("ALL PASS")
