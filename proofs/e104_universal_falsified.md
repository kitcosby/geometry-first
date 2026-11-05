# E104 — the universal diagonal claim is FALSE (2026-10-26)

**Result (PROVED kill).** The universal E104 claim — that samples of the
E103 skeleton with *any* traceless diagonal pattern Lie-generate
sl(3,ℝ) (rank 8) — is false. The claim was withdrawn 2026-10-12 and
stood OPEN; it is now resolved in the negative. The five specific
patterns proved in `e100_e103_e104_proved.md` are untouched.

The skeleton, fixed by E103's decomposition
T_a(x) = urx·H + srx·E12 + cxp·E21 + uxp·(E23 − E32), is

    M_d(x) = [[d1, srx, 0], [cxp, d2, uxp], [0, −uxp, d3]],   d1+d2+d3 = 0.

## Counterexample

Take the traceless pattern built from core functions in exactly the
style of the tested pattern (srx, cxp, −srx−cxp):

    d1 = (cxp − 2·srx)/3
    d2 = (srx − 2·cxp)/3 − uxp
    d3 = (srx + cxp)/3 + uxp          (trace ≡ 0 identically)

Put λ(x) = (srx + cxp)/3. For every x,

    M_d(x) · (1,1,1)ᵀ = λ(x) · (1,1,1)ᵀ,

by three lines of algebra: row 1 gives d1 + srx = λ; row 2 gives
cxp + d2 + uxp = λ; row 3 gives −uxp + d3 = λ; and the three row
equations are consistent precisely because d1+d2+d3 = 0 forces
3λ = srx + cxp. Hence *every* sample of this pattern — any sample
set, any number of samples — lies in the stabilizer of the line
ℝ(1,1,1), a parabolic subalgebra of dimension 6 inside sl(3,ℝ).
The stabilizer is a subalgebra, so brackets stay inside it: the
generated algebra has rank ≤ 6 < 8 for this pattern. ∎

**Exact verification.** The eigenvector identity and trace were
confirmed in exact ℚ(√2,√3) arithmetic at x = π/6, π/4, π/3, π/12;
the Lie closure of the four samples has dimension exactly 6 — the
samples generate the entire line-stabilizer, so the bound is sharp.
Controls at the same four points — patterns (0,0,0) (E105),
(urx,−urx,0) (T_a itself), (1,−1,0) — all give closure dimension 8.
Scripts: `r-theory-db/ansatz_tests_2026-10-26/scripts/`
(`field.py`, `e104.py`); full test record:
`r-theory-db/ansatz_tests_2026-10-26/cluster_B_geometry_first.md` §1.

## What survives — the constant-pattern subcase (PROVED)

If the pattern is a *constant* diagonal D, then D cancels in the
sample differences M_i − M_j = O(x_i) − O(x_j). Whenever the samples'
off-diagonal coefficient vectors (srx, cxp, uxp) affinely span ℝ³,
the differences span V = span{E12, E21, E23−E32}, and Lie(V) = sl(3)
by explicit brackets: [E12,E21] = diag(1,−1,0);
[E12, E23−E32] = E13; [E21, E23−E32] = E31;
[diag(1,−1,0), E23−E32] = −E23−E32, which with E23−E32 separates
E23, E32; [E13,E31] = diag(1,0,−1) — all eight basis elements.
Hence for constant patterns, generation holds for every traceless D,
on any sample set with affine span 3 — a checkable property of the
sample set (verified exactly, rank 3, for the sample set used in the
test). Supporting CHECKED evidence: an exact grid scan of 49 constant
patterns diag(p, q, −p−q), p,q ∈ {−2,−1,−1/2,0,1/2,1,2}, at the same
samples gives closure 8 for all 49.

Caveat: the original E104 test's sample points are not recorded in
the source files, so this subcase's application to *those specific
samples* rests on the affine-span condition, which their record does
not state.

## Database record

- `gf-e104-open`: IN → IC (2026-10-26, Kit-approved), statement
  carries the counterexample.
- `gf-e104-constant-patterns`: new row, CP — the subcase above.
