# Method — What "Geometry First" Means

## Scope labels

Every claim in this repository carries exactly one label:

- **PROVED** — Complete mathematical proof. Either exact computer algebra (SageMath over ℚ/ℚ̄, Groebner bases, exact nullspaces) or symbolic algebra verifiable by hand. The proof script in `proofs/` reproduces it.
- **CHECKED** — Numerical verification only (floating-point, sampling). May support a PROVED theorem but is never presented as proof. Labeled explicitly.

Labels never used here: ASSERTED, OPEN, INCOMPLETE, HYPOTHESIS. If a result has one of those labels, it is not in this repository.

## Proof standards

1. **Exact over numeric.** Where Sage is used, it is exact (ℚ, not ℝDF; Groebner, not optimization; rational nullspace, not SVD). Numerical checks are supplementary, never primary.
2. **No hidden axioms.** Every proof states its inputs. If a proof uses only the operator definitions and standard algebra, it says so. Nothing enters through the back door.
3. **Failures disclosed.** T4* is included as "strong evidence" (not PROVED) because the real-algebraic emptiness proof is pending. The gap is stated, not hidden.
4. **Reproducible.** Every PROVED theorem has a script in `proofs/`. Run it; get the same result.

## What "geometry first" means

1. **Start from the operators.** srx, sxp, cxp, crx — defined, not hypothesized.
2. **Derive the structure.** Clifford algebras, intertwiners, symmetries, conic geometry — all *found* in the operators, not imposed.
3. **Mark the wall.** When derivation stops, say so precisely (see `WALL.md`).
4. **Add nothing.** No E8, no Axiom Zero, no SM, no coupling — unless and until proven.

## What this is not

- Not a physical theory (no Γ, no units, no spacetime).
- Not a unification (no fermions, no gauge groups).
- Not complete (the wall is real).

It is a **foundation**: the geometry the operators prove about themselves, cleaned of everything else. Unification — if it comes — starts here, not from imported beauty.
