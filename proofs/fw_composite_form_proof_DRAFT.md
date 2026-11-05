# FlatWave composite form — promotion from CHECKED to PROVED (DRAFT, 2026-10-12)

## Claim

For all x where srx, cxp are defined and srx·cxp != 1:

    (srx + cxp) / (srx·cxp − 1)  =  1/urx + 1/uxp  =  FlatWave = sgn(sin 2x)

Kit proposed the left-hand form on 2026-10-12; it was verified numerically
at 10 angles and left labeled CHECKED. The proof below uses only results
already PROVED in the corpus, so the label can be promoted.

## Proof

Definitions (Pythagorean Foundation): urx = srx − crx, uxp = cxp − sxp.

INVPAIR (already in theorems/INDEX.md, PROVED): srx·sxp = 1 and cxp·crx = 1,
i.e. sxp = 1/srx and crx = 1/cxp. (Both are the identities
|csc x|^2 − cot^2 x = 1 and |sec x|^2 − tan^2 x = 1 in disguise.)

Step 1. urx = srx − 1/cxp = (srx·cxp − 1)/cxp.
Step 2. uxp = cxp − 1/srx = (srx·cxp − 1)/srx.
Step 3. Therefore
    1/urx + 1/uxp = cxp/(srx·cxp − 1) + srx/(srx·cxp − 1)
                  = (srx + cxp)/(srx·cxp − 1).  ∎

The last equality to FlatWave, 1/urx + 1/uxp = sgn(sin 2x), is the
already-PROVED FlatWave identity from the Pythagorean Foundation.

## Why this form matters

It computes FlatWave — the composite encoding all eight octant words
(+1 on OII, IOI, OIO, IOO; −1 on OOO, IIO, OOI, III) — from only the two
operators srx and cxp, not the four (srx, crx, cxp, sxp) needed to form
urx and uxp separately. It is the 2-operator form of the FW identity,
matching the minimality theme of E106/E107 (2 independent functions).

## Scope / domain note

Requires srx·cxp != 1 (else urx = uxp = 0 and both sides are undefined)
and x away from the operator singularities (sin x, cos x, sin 2x nonzero),
the same domain discipline used in the cycle theorem. On that domain the
argument is pure algebra over the proved INVPAIR identities — no numeric
sampling is used, per the proof-over-sampling rule.

Numeric spot-check performed while drafting (sanity only, not part of the
proof): both sides agree and equal sgn(sin 2x) at 15, 30, 45, 60, 120,
200, and 300 degrees, and srx·sxp = cxp·crx = 1 at each.
