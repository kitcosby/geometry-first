# Proof-archive → Pascal conversion plan

Kit's directive (2026-10-13): Pascal is the language of this
repository — new proof programs are written in Pascal, and the
existing `proofs/` archive is converted. Pascal versions are
**additions, not replacements**: the Sage/Python scripts and the
`.tex`/`.md` records stay as the proofs of record; a Pascal port
never upgrades a scope label.

Scope labels for this file: **CONVERTED-EXACT** (Pascal re-derivation
in exact arithmetic, zero-residual equality), **CONVERTED-CHECK**
(Pascal numeric recheck, CHECKED strength), **PENDING**, **NOT
PORTABLE** (needs symbolic machinery Pascal does not have here).

## Classification key

- **(a) Exactly portable** — the computation is linear algebra /
  polynomial algebra over ℚ or over a small algebraic field
  (here ℚ(√2,√3)) that Pascal can carry exactly with a rational /
  quadratic-field number type.
- **(a′) Exactly portable, engine not yet built** — exact, but needs
  a ℚ(t) rational-function engine (entries are rational functions of
  t = tan(x/2)) or a large ℚ-linear solve; portable in principle, a
  bigger build than batch 1.
- **(b) Numeric CHECK only** — as written the script is floating
  point (numpy/scipy), or its content is an identity the Pascal layer
  already re-checks numerically; no exact content is lost by a
  Double-precision port, and none is gained either.
- **(c) Not portable** — Gröbner bases, algebraic-closure root
  extraction/factorization, or general QQbar machinery. These stay
  Sage-only; Pascal can at most host the registry object.

## Inventory of `proofs/*.py` (29 scripts present)

| Script | Establishes | Class | Pascal status |
|---|---|---|---|
| q4_exact.py | T1: commutant of Cl(4,4) 16-spinor is scalars (dim 1), exact over ℚ | (a′) — 256-unknown ℚ-linear nullspace; exact but large | **CONVERTED-EXACT (batch 2)** — T1-EXACT in the gate; commutant dim 1 exact |
| q5_exact.py | T2: sector intertwiner space dim = 16, exact over ℚ | (a′) — same engine as q4 | **CONVERTED-EXACT (batch 2)** for the dimension (dim 16 exact); **finding:** the script's kernel basis, reshaped row-major by its own `vec_to_mat`, satisfies R_μᵀ·B = B·J_μᵀ, not B·R_μ = J_μ·B — the Kronecker system encodes column-stacking vec while the scripts reshape row-major, so the basis matrices are transposes of intertwiners (transposes verified exact intertwiners, 16/16). T2-EXACT object reports this and is not in the pass/fail gate |
| q6a_exact.py | T3 (attempt 1): explicit X in intertwiner space with X² = +I via rational combinations, verified exactly | (a′) — needs the q4/q5 intertwiner basis over ℚ | **CONVERTED (batch 2)** — the exact search (16 singles + 120 pairs B_i+B_j) finds NO candidate whose square is a positive rational-square multiple of I (none has a scalar square at all); this matches the Sage script's own trajectory, whose simple search also fails and falls through to Gröbner. T3's establishment therefore rests on the class-(c) scripts, exactly as before the port; the port is a verified negative for the rational-combination attempt, not a proof of T3 |
| q6a_groebner.py | T3 (attempt 2): X² = +I via Gröbner basis (16 vars, 256 quadrics) | (c) Gröbner | NOT PORTABLE (registry object only; T3 itself is (a′) via q6a_exact) |
| q6b_groebner.py | T4*: no real X with X² = −I — Gröbner side | (c) Gröbner | NOT PORTABLE |
| q6b_solutions.py | T4*: extracts the 2 algebraic-closure solutions, checks non-reality | (c) QQbar root extraction | NOT PORTABLE |
| d2_structure.py | D2CLOSE / D2GRAM: intertwiner space V not closed under mult.; Gram = 4·I₁₆, exact | (a′) — needs intertwiner basis over ℚ (shares q5 engine) | **CONVERTED-EXACT (batch 2)** — D2-EXACT in the gate; Gram = 4·I₁₆ confirmed componentwise; V not closed (first failure B[0]·B[0]); I₁₆ ∉ V |
| d3_forms.py | D3 first attempt (Spin(4,4)-invariant forms, volume element) — superseded by d3b | (a′) but superseded; convert d3b instead | PENDING (superseded) |
| d3b_forms.py | D3FORM / D3VOL: exactly 2-dim invariant symmetric forms (rank 8), 0 skew; vol² = +I, anticommutes | (a′) — ℚ-linear algebra on 16×16 generators | **CONVERTED-EXACT (batch 2)** — D3B-EXACT in the gate; form space dim 2, both forms symmetric rank 8; vol² = +I, vol not scalar |
| q7_dSdx.py | Q7: d(sawup)/dx Bernoulli identity, symbolic per-octant + numeric global check | (b) — Pascal layer already rechecks Q7/DURX numerically (RTheory.Verify) | CONVERTED-CHECK (already in layer) |
| run150b_doubleangle_product.py | SINPROD corrected global form urx·uxp = 4·\|csc(2x)\| (sympy trigsimp + numeric) | (b) — already rechecked in RTheory.Verify as SINPROD | CONVERTED-CHECK (already in layer) |
| e70_exact.py | E70: T_a at 15°/30°/45°/60° + 4 commutators: exact rank 8, all 28 commutators in span ⇒ sl(3,ℝ) | **(a)** — all entries in ℚ(√2,√3) | **CONVERTED-EXACT (batch 1)** |
| e71_exact_casimir.py | E71: exact structure constants, Killing signature (5,3), C₂ = (4/9)·I | **(a)** — signature via exact Sylvester inertia over ℚ(√2,√3) (replaces QQbar eigenvalue signs) | **CONVERTED-EXACT (batch 1)** |
| e72_exact_rest.py | E72: Cartan centralizer dim 2, so(3) triple (antisym. parts rank 3, closed), A2 roots via ad char poly | **(a)** for centralizer / so(3) / char-poly coefficients (Faddeev–LeVerrier); the Sage *factorization* of the char poly is (c) — the port instead proves the exact root-pairing form (char poly = λ²·(even polynomial), rank(ad) = 6) | **CONVERTED-EXACT (batch 1, root-pairing form)** |
| m_refined_closure.py | MREF4: M_refined closes at dim 16 (gl(4)), traceless correction at dim 15 = sl(4,ℝ) | (a′) — closure over ℚ(t)/QQbar | **CONVERTED-EXACT (batch 2)** — MREF-CLOSURE in the gate; samples lie in ℚ(√2,√3) so the QQbar computation ports exactly onto TField; dims 16/15, both closed (traceless-all-gens for the original set computed False, as the script prints it) |
| m_refined_closure_aa_attempt.py | MREF4 attempt over AA (real algebraic) — superseded working attempt | (a′)/(c) — AA machinery | PENDING (superseded attempt) |
| m_refined_det.py | MREF4: det M_refined(t) factors exactly in ℚ(t), zero at t² = 1/5 | (a′) — ℚ(t) factorization engine | **CONVERTED-EXACT (batch 2)** — MREF-DET in the gate; det = −(3t+5)(5t²−1)/((t+1)(t²+1)²) from the Pascal ℚ(t) engine, A·B = det exact as scripted; 4-sample closure re-run 16/15 |
| m_refined_killing.py | MREF4: Killing signature of the 15-dim algebra via QQbar eigenvalues | (a′) — structure constants over ℚ(t) + exact inertia | **CONVERTED-EXACT (batch 2)** — MREF-KILLING in the gate; 225/225 structure-constant reconstructions exact; signature (9,6,0) by exact Sylvester inertia over ℚ(√2,√3) (replaces the QQbar eigenvalue-sign count) |
| cross_m1_m2_a.py / _b.py / _c.py | M1M2: tr(M1·M2) = 3 identically; entry forms — sympy exact rational algebra in t = tan(x/2) per quadrant | (a′) exact over ℚ(t); Pascal layer already has the numeric scalar recheck (TM1M2ScalarProof) | CONVERTED-CHECK (already); exact ℚ(t) port PENDING |
| expmap_uxp_so3_a.py | SO3EXP: θ² identity, sympy symbolic in t | (a′) exact over ℚ(t); numeric θ² recheck already in layer (TSo3ThetaProof) | CONVERTED-CHECK (already); exact port PENDING |
| expmap_uxp_so3_b.py | SO3EXP numeric matrix at 30° via scipy expm | (b) floating point as written | CONVERTED-CHECK (already, via θ² object) |
| quadrant_flip_polar_a.py / _b.py | QFLIP: quadrant flip = perfect-square swap, sympy symbolic in w and t | (a′) exact rational/symbolic in t | **CONVERTED-EXACT (batch 2)** as a faithful port — QFLIP-EXACT object replicates both scripts as written and FAILS its conjunction on a **finding in the originals**: polar_a's three swapped-role π-shift assertions (urx_QI→uxp_QIII, uxp_QI→urx_QIII, sxp_QI→crx_QIII) are FALSE as written (the original sympy script itself prints False for exactly these three); polar_b's same-name form (8/8 + block flip) and both w-identities/perfect squares over ℚ(i)-Laurent are exact. Not in the pass/fail gate; not harmonized |
| wick_sl8_a.py / _b.py | WICKWALL numeric side: Wick-rotated cycle generators at 30°, closure dim 28 (numpy) | (b) floating point as written; Pascal wall probe TWallWickProof already rechecks closure dim 28 numerically | CONVERTED-CHECK (already in layer) |
| maxwell_exercise.py | Maxwell exercise: Riccati identities in QI (sympy) — operator identities only; E/B dictionary is ASSERTED, not a theorem | (a′) for the identities; the entry is on the not-theorems list | **CONVERTED-EXACT (batch 2, identities only; no scope change)** — MAXWELL-EXACT in the gate; 4 Riccati + 2 UNA log-derivative identities exact over ℚ(t); curl char poly λ(λ²+a²+b²+c²) verified on the exact {0,1,2}³ grid; E/B dictionary stays out |
| entropy_selector_selfconsistency.py | Entropy-selector self-consistency identities (sympy, exact given E16); the underlying claim was flagged untested by Kit | (a′) symbolic identities in t | **CONVERTED-EXACT (batch 2, identities only)** — ENTROPY-EXACT in the gate; reciprocal pairs, mirror identity, and r3 derivative identity (no chain factor, as scripted) all exact over ℚ(t); sign conclusion still needs E16 and the physical reading stays ASSERTED (not an INDEX theorem) |
| **Referenced but not present as files** | E88 / E90b / E90c / E91 / E94 (block ladder sl(4)/sl(5)/sl(6), exact over ℚ̄ — record is `proofs/sl_n_exact.tex` in the workspace root, scripts not in this repo); E105/E106/E107 scripts (GitHub-only, SHAs in INDEX) | E106/E107 already have numeric Pascal rechecks (TE106Proof, TE107Proof); ladder (a′) once scripts are vendored | PENDING |

Counts: 29 scripts present — (a) batch-1 exact: 3 (e70, e71, e72);
batch-2 exact: 12 files (q4, q5, q6a_exact, d2, d3b, the m_refined
trio, the quadrant_flip_polar pair, maxwell, entropy-selector);
already converted as CHECKED rechecks: 6 (q7, run150b, cross_m1_m2
counted once as a trio + expmap_uxp_so3_b + wick pair counted once —
see table); (a′) still pending an exact port: cross_m1_m2 exact
ℚ(t) forms and expmap_uxp_so3_a (both already CHECKED in the layer),
plus the two superseded attempts (d3_forms, m_refined_closure_aa);
(c) not portable: 3 (q6a_groebner, q6b_groebner, q6b_solutions).

## Batch 1 (this pass) — the sl(3) exact core, E70/E71/E72

New unit `src/RTheory.ExactField.pas`: exact arithmetic in the
biquadratic field ℚ(√2,√3) — every matrix entry in e70–e72 lies in
it (15° values are (√6±√2)/4). Elements are a + b√2 + c√3 + d√6
over normalized Int64 rationals; the field has an exact sign
function (recursive: ℚ(√2,√3) → ℚ(√2) → ℚ), so Sylvester inertia
gives the Killing signature exactly — no floating point anywhere.

New unit `src/RTheory.ExactLie.pas`: verifiable proof objects
- `TE70ExactProof` — the 8 generators exactly as in e70 (4 samples +
  3 commutators + 1 double commutator); exact rank = 8; all 28
  commutators exactly in the span (augmented-rank test).
- `TE71ExactProof` — structure constants via the pivot solve exactly
  as e71 (with full reconstruction check of every commutator, all 9
  entries — stronger than e71's single-commutator spot check);
  Killing matrix; signature (5,3,0) by exact inertia; C₂ computed
  from the inverse Killing form equals (4/9)·I **exactly** (zero
  residual, componentwise).
- `TE72ExactProof` — centralizer of G₀ inside the span has exact
  dimension 2 (nullity of the 9×8 commutator matrix); its two basis
  elements commute exactly; antisymmetric parts have exact rank 3
  and the chosen triple closes exactly; ad_{H₀} has exact rank 6 and
  its Faddeev–LeVerrier characteristic polynomial is exactly
  λ²·(λ⁶ + c₄λ⁴ + c₂λ² + c₀) — roots in ± pairs with exactly two
  zero roots (the A2 root pattern), proved by the exact coefficient
  form, not by QQbar factorization.

These objects are EXACT re-derivations in Pascal; the Sage scripts
remain the proofs of record named in `ProofSource`, and the registry
entries (E70-E72 etc.) keep their existing labels untouched.

## Batch 1 verification results (2026-10-13)

- VM engine (FPC 3.2.2, `~/workspace/fpc-toolchain`): **31/31 PASS,
  exit 0** — the previous 28 numeric checks plus E70-EXACT,
  E71-EXACT, E72-EXACT, all zero-residual in ℚ(√2,√3).
- Laptop engine (FPC 3.2.2, WSL `lampy` `/usr/bin/fpc`, sources
  shipped to `/root/pascal`): **31/31 PASS, exit 0**, identical
  registry totals (87 objects). This run also exercised the v4
  registry build on the laptop for the first time, closing the
  laptop leg that pascal/README.md had marked PENDING for v4.
- Engineering note, disclosed: the first version of the field unit
  used Int64 rationals and E71/E72 **overflowed loudly**
  (EIntOverflow, overflow checking on) inside the structure-constant
  / Killing computations. The rational type was rebuilt on a
  base-10^9 arbitrary-precision integer inside
  `RTheory.ExactField.pas`; the bigint core was unit-tested against
  Python's exact `fractions` arithmetic (telescoping products,
  Pell-type sign near-misses 7−5√2 < 0 and 1351−955√2 > 0, inverse
  round-trips) before the full run. No result in this batch rests
  on floating point.

## Batch 2 (2026-10-17) — the ℚ-linear family and the ℚ(t) family

New unit `src/RTheory.RatSparse.pas`: sparse exact linear algebra
over ℚ (bigint `TRational` from batch 1) — sparse-row RREF and
kernel bases for the vectorized Kronecker systems (2048×256 for
q4's commutant, 1024×256 for the intertwiner system, 7168×256 for
d3b's invariant forms; every system row starts with ≤ 2 nonzeros),
dense 16×16 ℚ helpers, and an exact rational-square test (integer
sqrt on `TBigInt`).

New unit `src/RTheory.ExactClifford.pas`: T1-EXACT (q4), T2-EXACT
(q5), T3-EXACT (q6a), D2-EXACT (d2), D3B-EXACT (d3b).

New unit `src/RTheory.RatFunc.pas`: the ℚ(t) engine — `TPoly` over
ℚ (DivMod, monic Gcd, derivative, Horner evaluation at a `TField`
point), `TRatFunc` (gcd-cancelled, monic-denominator normal form;
arithmetic, derivative, general Möbius substitution, ≤4×4
determinant), and `TLaurent` over ℚ(i) for the polar w-identities
(ratio equality by cross-multiplication).

New unit `src/RTheory.ExactQt.pas`: MREF-CLOSURE, MREF-DET,
MREF-KILLING, QFLIP-EXACT, MAXWELL-EXACT, ENTROPY-EXACT. The
m_refined sample points (2−√3, √2−1, 1/√3, 1/2, 3/4) all lie in
ℚ(√2,√3), so the Sage/QQbar sample computations port exactly onto
batch 1's `TField`; the script's QQbar eigenvalue-sign count for
the Killing form is replaced by exact Sylvester inertia
(`FSymmetricInertia`), which is basis-independent under
congruence.

### Batch 2 verification results (2026-10-17)

- VM engine (FPC 3.2.2, `~/workspace/fpc-toolchain`): **39/39 PASS,
  exit 0** — the 31 batch-1 checks plus T1-EXACT, D2-EXACT,
  D3B-EXACT, MREF-CLOSURE, MREF-DET, MREF-KILLING, MAXWELL-EXACT,
  ENTROPY-EXACT, all zero-residual. Registry unchanged: 87 objects.
- Laptop engine (FPC 3.2.2, WSL `lampy` `/usr/bin/fpc`, sources
  shipped to `/root/pascal`): **39/39 PASS, exit 0**, identical
  registry totals (87 objects).
- Three converted objects are deliberately NOT in the VerifyProofs
  pass/fail gate, because their honest verdicts are findings, and
  the gate's exit-0 semantics ("ALL CHECKS PASS") must stay
  meaningful. They are exercised by the standalone batch test
  programs:
  - **T2-EXACT (q5)** — kernel dimension exactly 16 (the theorem),
    but the scripts' basis matrices are transposes of
    intertwiners: the Sage system encodes column-stacking vec,
    the scripts reshape row-major. Disclosed in the table row
    above; the dimension claim is unaffected.
  - **T3-EXACT (q6a)** — the rational-combination search
    (singles + pairs) is a verified negative, matching the Sage
    script's own trajectory to its Gröbner fallback.
  - **QFLIP-EXACT** — faithful port; polar_a's three swapped-role
    assertions are false as written in the original (see table).

### Disclosed infrastructure repair (batch-1 unit, found by batch 2)

Batch 2's field-RREF computations (rationals with ~10²³
denominators, confirmed legitimate against a Python `fractions`
model) raised `ERationalOverflow: bigint division not exact`
inside batch 1's `RTheory.ExactField.pas`. Root cause, in
`MagDivMod`'s multiply–subtract loop: the negative-limb
correction applied `T := T + BIGBASE` only once, but
`T = U[J+I] − Borrow − (Prod mod BIGBASE)` can fall to
−2·BIGBASE+2, so a negative limb survived and was wrapped by
`UInt32(T)`, silently corrupting the division. Evidence: minimal
repro `RatDiv(11981083046996584384200, 155808076914)` raised
pre-fix and returns `665615724833143576900/8656004273` post-fix
(agrees with Python `fractions`); an exact Python replica of the
shipped algorithm underflows on 1,120/15,403 random and
5,279/15,017 crafted multi-limb divisions, and produces zero
wrong results whenever the underflow does not fire. Fix (minimal,
2026-10-17): the borrow step loops the `+BIGBASE` correction
(`while T < 0 do begin T := T + BIGBASE; Inc(Borrow); end;`,
at most two iterations; Borrow stays ≤ BIGBASE, which the
existing add-back step tolerates). `FMul` was independently
probe-verified correct (√2·√3 = √6 both orderings, √2² = 2,
√6² = 6, √2·√6 = 2√3) and is untouched. Regression after the fix:
batch-1 suite 31/31 PASS exit 0 (E71's signature (5,3,0) and
C₂ = (4/9)I unchanged — batch-1 claims are zero-residual
equalities, hence self-certifying, and its data never hit the
digit pattern); 300 public-API division cases against Python
`fractions` clean. A separate alarm raised during this batch —
that `FMul` might compute a conjugate ordering — was tested and
REFUTED by the probe above; it is recorded here so it is not
re-raised.

## Remaining after batch 2

1. Exact ℚ(t) ports still pending (both already CONVERTED-CHECK in
   the layer): cross_m1_m2 exact entry forms, expmap_uxp_so3_a.
   The ℚ(t) engine now exists, so these are ports, not builds.
2. Vendored ladder scripts (E88–E94) if/when they land in the repo.
3. Never: the three Gröbner/QQbar scripts (class (c)).
