# FW Ladder — n = 11 rung: DRAFT addendum (not yet merged into fw_ladder.md)

**Date:** 2026-10-16 · **Status of this draft:** MERGED 2026-10-21 (Kit's
order) as the "Addendum 2026-10-21" of `proofs/fw_ladder.md`, after the
Sage-exact n = 11 run on toe (FW 407.7 s, DA 369.2 s; closure 120,
traceless, Killing (65, 55, 0) both variants) raised the rung to PROVED —
the merged addendum records that run and supersedes this draft's
CHECKED-grade account, which is retained here for the record. Not pushed.
The rung itself
was run at Kit's order in the Algebra-first thread (desktop leg) and is filed
in the daily memory (2026-10-16). This draft exists so the geometry-first
record does not silently stop at n = 10.

## What was run

Same construction as `proofs/fw_ladder.py`, copied verbatim, N = 11, both
variants (FW and DA). No Sage on the desktop machine, so the computation was
done in exact arithmetic over GF(p) at two primes p ≡ 1 mod 24
(p = 1000033 and p = 999983; the sample field Q(√2, √3) embeds in both).
Script: `/tmp/ladder11.py` (desktop leg).

## Results (identical at both primes, both variants)

| Variant | Sample span | Closure | Basis traceless |
|---|---|---|---|
| FW | 7 | **120** = 11² − 1 | yes |
| DA | 7 | **120** = 11² − 1 | yes |

Sample self-checks (FW = sgn(sin 2x), SINPROD) passed mod p.

## Scope — read this before citing the rung

- Identification uses the ladder's own logic (real samples + dimension
  120 + traceless basis ⇒ sl(11, R)); the mod-p independence of the
  certifying minors lifts the dimension statement to Q(√2, √3). This is
  exact-grade for **dimension and tracelessness**, the two quantities the
  identification needs.
- **The Killing signature cross-check was NOT computed in this leg.**
  Rungs n = 3..10 carry exact split signatures ((44,36,0) at n = 9,
  (54,45,0) at n = 10); the n = 11 split signature would be (65,55,0),
  and it is currently **unverified** — do not cite it as checked.
- Grade accordingly: n = 3..10 are PROVED (Sage-exact, with Killing);
  n = 11 is CHECKED at the grade above. A Sage-exact n = 11 run
  (fw_ladder.py with N = 11 on toe, which would also supply the Killing
  signature) has not been run; at the n = 10 runtimes (103–221 s per
  rung) it is a bounded job, and it is the step that would bring n = 11
  to the same grade as the rest of the ladder.
- Universal E104 and the all-n induction remain OPEN, unchanged.
