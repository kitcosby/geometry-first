# The 27 of E6, with its Cubic Invariant

Date: 2026-10-23. Status: CONSTRUCTED and PROVED exact
(23 counted checks, exit 0, ~85 s; script
proofs/e6_rep27.py), on the seed algebra of e6_seed.md.
The representation is of the complexified e6 built there.

## The space and the forced charges

V27 = S + V + 1 (16 + 10 + 1): S the vacuum-sector 16,
V = C^10, 1 a singlet. k = so(10) + Z6 acts in the obvious
way (rho, vector, trivial) with Z6 charges (al, be, ga).
The bracket [Z6, T_x] = T_x for x in p+ forces the charges:
the nonzero components of T_x must raise charge by exactly
1, and with (al, be, ga) proportional to (1, -2, 4) the
admissible differences give -3c = 1, so
  (al, be, ga) = (-1/3, 2/3, -4/3).
With Z' = -3 Z6 these are the textbook charges (1, -2, 4).

## The p-action: four maps, two determined products

By so(10)-equivariance (and charge), T_x for x in p+ has
exactly two nonzero components and T_y for y in p- two:
  T_x:  1 -> S : s |-> a1 s x;
        S -> V : psi |-> a2 W(x, psi),
          W(x, psi)_a = B(x, gamma_a psi);
  T_y:  V -> S : v |-> a3 gamma(v) y;
        S -> 1 : psi |-> a4 B(psi, y).
Every other component is zero by Hom-vanishing and charge
together (e.g. 16 (x) 10 contains 16' and 144, not 16).
W is symmetric in its spinor arguments and
B(gamma_a ., .) is symmetric on 16' (matrix identities,
verified) -- these two symmetries ARE the conditions
[T_x, T_x'] = 0 and [T_y, T_y'] = 0. W is equivariant and
[Sig_j, gamma_b] = gamma_{E_j e_b} (verified), which are
the [k, p] bracket conditions.

Only the products of the a's are determined (rescaling the
three piece-bases is gauge). They are SOLVED, not chosen:
  a1 a4 = -1/8, from [T_x, T_y] on the singlet
    (= lam * ga * B(x, y), seed lam = -3/32);
  a2 a3 = -1/16, from the S-block identity
    rho(A(x,y)) psi + (1/32) B(x,y) psi
      + (1/8) B(psi,y) x + a2a3 gamma(W(x,psi)) y = 0,
    verified on all 4096 basis triples (x, psi, y);
  and CONFIRMED independently by the V-block identity
    a2a3 W(x, gamma(v) y) = A(x,y) v - (1/16) B(x,y) v,
    verified on all 2560 basis evaluations -- the same
    constant from two different Fierz-type identities.
Gauge fixed (a1, a2, a3, a4) = (1, 1, -1/16, -1/8), the
assembled action of all 78 generators is verified bracket
by bracket on all 27 basis vectors: [Z6, T], [Y_j, T],
[T_x, T_y] against the seed bracket (all 256 pairs),
[T, T] same-side zero. Irreducibility: a basis vector from
each of the three pieces is cyclic under the 78 operators
(complete reducibility of e6-modules cited).

## The cubic invariant

The k-invariant symmetric cubics on V27 are exactly the
span of C1 = W(psi,psi).v and C2 = s (v.v): every other
piece content has nonzero Z6 charge or vanishing Hom
(charge computation in the script header; the two Hom
spaces are 1-dimensional). Each form is separately
k-invariant (verified by tensor deltas). p-invariance
imposes one condition, solved:
  C = 2 C1 - C2   (ratio (d1, d2) = (2/3, -1/3)),
i.e. C(psi, v, s) is proportional to
2 W(psi,psi).v - s (v.v), verified invariant under all 78
generators and not identically zero. Since the ansatz span
is ALL k-invariant cubics, the e6-invariant cubic is
unique up to scale: Sym^3(27*)^{e6} is 1-dimensional,
spanned by C. This is the determinant-type invariant of
the 27, in these coordinates.

## Disclosures

Four defects on the way, all in my machinery, none in the
mathematics on record: (i) the first assembled run crashed
on a sympy slicing quirk (M[:16] returns a list in this
version) and a wrong nonvanishing probe -- W(e_vac, e_vac)
= 0 necessarily, by the subset-complement pairing; W is
nonzero on other basis pairs. (ii) The first cubic
verification (per-triple evaluations) timed out at 1700 s
and was replaced by tensor contraction: the trilinear
tensors are built once (80 and 10 nonzero entries) and the
generator action is a sparse contraction. (iii) That
contraction initially used the transposed action; a direct
cross-check caught it (tensor delta 2i/3 vs direct 0 at
(0, 7, 16) under Y0) and the row-contraction form gives
the recorded ratio (2, -1) -- the (1, -1) from the buggy
pass is superseded. (iv) A probe-design lesson in the
record: basis-diagonal values C(gX, X, X) see only
W(e_r, e_r) = 0 and are blind to C1; the trilinear form
itself must be probed.

## Scope

Built: the 27 as an explicit e6-module (operator form, all
brackets verified) and its invariant cubic, for the
complexified algebra of the seed. Not built: real forms;
the Freudenthal/Jordan realisation (the 27 as J_3(O_C)
with C as determinant is the classical identification --
consistent with everything computed here, but the Jordan
product itself is not constructed); the EIII homogeneous
geometry. No FW/saw content. Local file; not pushed.
