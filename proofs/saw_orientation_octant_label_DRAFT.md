# DRAFT — Saw orientation invariant and the octant label (mod π)

Status: DRAFT, not yet in theorems/INDEX.md. Assembled 2026-10-14 from the
Algebra-first thread discussion of the same day. Every claim below is labeled;
the derivations use only identities already in the record.

## 1. The saw orientation invariant — PROVED

Definitions (Pythagorean Foundation, in record): sawdown = 1/uxp,
sawup = 1/urx, on the domain D (quadrant interiors; seams x = kπ/2 excluded).

Logarithmic-derivative identities (maxwell_exercise, PROVED):
  urx'/urx = −(srx+crx)/2,   uxp'/uxp = +(cxp+sxp)/2,
with srx, sxp, crx, cxp all positive on D.

Sign facts (Foundation, PROVED): urx·uxp = 4/|sin 2x| > 0 and
urx + uxp = 4/sin 2x, so sgn(urx) = sgn(uxp) = sgn(FW) = sgn(sin 2x).

Therefore:
  d(sawdown)/dx = −(cxp+sxp)/(2·uxp)
  d(sawup)/dx   =  (srx+crx)/(2·urx)

In the FW = −1 sectors (QII, QIV): uxp < 0, so sawdown is signed-increasing —
but sawdown is negative there, so the increase approaches the x-axis and
|sawdown| strictly decreases. Mirror: sawup is signed-decreasing there while
|sawup| strictly increases, away from the axis.

Global invariant (all quadrant interiors): |sawdown| strictly decreases
toward the axis and |sawup| strictly increases away from it. The names
"down"/"up" track magnitude-vs-axis in ALL quadrants; only the signed
monotonicity flips, with sgn(FW).
Spot value (in record): at 135°, d(sawdown)/dx = +√2/2, d(sawup)/dx = −√2/2.

## 2. Octant label completeness mod π — PROVED

Two internal bits at each x ∈ D:
  (a) sgn(FW) = sgn(sin 2x) — separates {QI, QIII} from {QII, QIV}.
  (b) saw dominance, |sawup| vs |sawdown| — equivalently sgn(p) in the UNA
      forms, |sawup|/|sawdown| = (1−p)/(1+p) — separates the halves of each
      quadrant (tie at the 45° points, both = 1/2).
Together: a complete octant label on [0, π).

Obstruction — PROVED: no function of the operators (values, derivatives,
any combination) distinguishes x from x + π; every operator is π-periodic
from its definition (|csc|, cot, |sec|, tan). QI vs QIII and QII vs QIV are
information-theoretically invisible internally. One external bit
(e.g. sgn(sin x)) completes the full quadrant, but that bit is an input to
the definitions, not an output.

## 3. Which handedness this resolves — and which it does not

1. Frame handedness of the matrix families A(x), B(x): RESOLVED internally.
   FW sits in the generator slots; "right-handed frame iff FW = +1" is a
   datum the operators supply, refined to the octant by bit (b).
2. Algebraic handedness: THERE IS NONE to resolve. The generation theorems
   are sign- and weight-blind — A(x) generates the same sl(3,R) in every
   quadrant (in record).
3. Spinor chirality (16 vs 16′, sign of Γ) and the π-flip: NOT resolved.
   No proved bridge maps the saw datum onto the Clifford side; which block
   is called "16" stays a convention. The operators carry exactly one
   orientation Z₂ (FW) and provably erase the other (π).

Precise summary: the saw characteristic is a complete internal orientation
for everything the theory actually computes — every branch formula in the
corpus needs exactly the octant-mod-π label — and it proves its own limit:
a full 2π handedness cannot be operator-internal; it needs exactly one
added bit.

## Open / not claimed
- No bridge from this datum to the Cl(4,4) spinor chirality sign (see the
  companion hunt spec for the sector-blending intertwiner C).
- Seams and 45° tie points are excluded above; boundary behaviour is not
  asserted here.
