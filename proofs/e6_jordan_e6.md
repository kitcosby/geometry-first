# 𝔢₆ on the Jordan side: the determinant-preserving algebra of H₃(𝕆_ℂ)

`e6_jordan_e6.py` — 20/20 exact checks over ℚ(i), exit 0. This is
the first half of the Φ bridge between the two models of the 27:
the representation-theoretic one (`e6_rep27.py`, 𝔢₆ built from the
Spin(10) spinor) and the Jordan one (`e6_jordan.py`, J = H₃(𝕆_ℂ)
with its determinant). Here the Lie algebra 𝔢₆ is realized
**intrinsically on J**, as the algebra of determinant-preserving
operators, with its full frame-relative structure; and its
𝔰𝔬(10) ⊕ 𝔲(1) part is matched quantitatively against the rep27
side. The intertwiner Φ itself (a basis-level dictionary) is the
remaining step — see the end.

## The solve

T ∈ End(J) preserves det infinitesimally iff δ_T Tdet = 0, where
Tdet is the det tensor of `e6_jordan.py` (89 nonzero entries).
The Peirce grading operator Z_J = 6L_E − 2·id (eigenvalues 4, 1, −2
on J₁, J₁⸝₂, J₀) grades End(J) by charge differences in
{0, ±3, ±6}, so the system splits into five independent sectors,
each solved by sparse exact Gaussian elimination over ℚ(i):

| charge | unknowns | equations | solution dim |
|---|---|---|---|
| 0 | 357 | 1359 | **46** |
| +3 | 176 | 864 | **16** |
| −3 | 176 | 864 | **16** |
| +6 | 10 | 90 | **0** |
| −6 | 10 | 90 | **0** |

Total **78**. Every assembled basis operator is re-verified against
δ_T det = 0 by direct tensor contraction (J2), as is the known
solution Z_J itself (J0, a regression guard).

## The algebra

- **Closure (J4a):** all 3,081 pairwise commutators are expressed
  in the basis via a pivot minor and verified as matrices — the
  solution space is a Lie algebra, computationally, in full. (The
  general reason: T ↦ δ_T is a Lie homomorphism from 𝔤𝔩(J) to
  operators on trilinear forms, δ_[A,B] = [δ_A, δ_B]; J3 verifies
  δ_[T_a,T_b] det = 0 independently for 780 commutators.)
- **Centre 0 (J4b), perfect (J4c):** [𝔤, 𝔤] spans all 78.

A 78-dimensional, centreless, perfect Lie algebra realized as the
determinant algebra of the exceptional Jordan determinant, with
the graded structure below, is 𝔢₆ in its Jordan (Freudenthal)
realization; the seed's 𝔢₆ is the same algebra in spinor
coordinates. What identifies the two *modules* (not just the
algebras) is the Casimir match in J9.

## The frame: parabolic, Levi, nilradical

The stabilizer of the **line** ℂE (E = diag(1,0,0)) has dimension
**62**: it is the maximal parabolic

𝔭 = Levi (46) ⊕ nilradical (16),

- the **Levi is exactly the charge-0 sector** (charge-0 operators
  preserve the 1-dimensional charge-4 piece ℂE automatically);
- the **nilradical is the charge-(+3) sector**: those operators
  kill e₀ for support reasons (charge 4 + 3 = 7 does not occur),
  form an abelian ideal of 𝔭 ([P⁺, P⁺] = 0 as matrices, J8a;
  [𝔭, P⁺] ⊆ P⁺ verified on all 62×16 pairs, J6b).

Consistency: the orbit of the line ℂE has dimension 78 − 62 = 16,
the complex dimension of the EIII variety E₆/(Spin(10)×U(1)).

Inside the Levi:

- **Centre 1-dimensional (J5c)**, spanned by Ẑ with Peirce block
  scalars (4, 1, −2)·s, **s = i** exactly (J5d) — the rep27 Z′
  charges, in an anti-Hermitian normalization.
- **Derived algebra 45-dimensional (J5e)**, and it is 𝔰𝔬(10)
  concretely: all 45 preserve the det₂ bilinear form
  B_Q(d₁,d₂,C) = (d₁d₂′ + d₂d₁′)/2 − bN(C, C′) on J₀
  (antisymmetry on all 45×55 basis pairs, J7).
- **[Ẑ, P±] = ±3s·P± (J8c)** — the 3-term grading
  𝔢₆ = 𝔭⁻ ⊕ Levi ⊕ 𝔭⁺ with [P⁺, P⁻] ⊆ Levi on all 256 pairs
  (J8b) and [Levi, P±] ⊆ P± on all 46×32 pairs (J8d). (The full
  parabolic does not preserve P⁻: [P⁺, P⁻] lands in the Levi.)

## The Casimir weld (J9)

One prescription, applied identically to both models: on the
derived 45, Gram_{ab} = −tr₂₇(T_a T_b), dual bases, Casimir
C = Σ_a T_a T^a. Scalars by piece:

| piece | Jordan side | rep27 side |
|---|---|---|
| singlet | J₁: 0 | 1: 0 |
| spinor 16 | J₁⸝₂: **−15/8** | S: **−15/8** |
| vector 10 | J₀: **−3/2** | V: **−3/2** |

Exact agreement. (These numbers are normalization-bound — the
DIII ladder's −45/4 used a different form on so(10); the content
here is that under *one* prescription the two independently built
𝔨-modules give identical Casimir spectra piece by piece.) The
𝔨-modules are thereby identified at the level of this complete
set of invariants; combined with the charge match (4, 1, −2) and
the bracket structure on both sides, the Φ bridge is reduced to
basis matching.

## What remains for Φ

Both 𝔢₆ actions are now explicit matrices: the rep27 one (78
operators on S ⊕ V ⊕ 1) and this one (78 operators on J). An
intertwiner is forced piece by piece: match the two 𝔰𝔬(10)
actions' Cartan weight systems (both concrete), which fixes Φ on
weight vectors up to 27 scalars; one 𝔭-bracket per piece then
fixes the scalars, and det∘Φ = c·C becomes a verification rather
than a search. Not built yet. The real forms and the EIII
homogeneous geometry also remain; note the parabolic above is
the EIII isotropy data (Levi = Spin(10)×U(1) algebra), so the
geometry layer now has its algebraic half in place.

## Disclosures

- **First run gave nullities (30, 0, 0, 0, 0).** Cause: the pair
  index PAIRL listed a repeated pattern index twice, doubling
  those coefficients in every δ equation and over-constraining
  the system. Caught because the expected dimensions failed; the
  fixed version carries the known solution Z_J as a regression
  guard (J0). Solver and checker shared the bug, so the run was
  self-consistent — the dimension expectations and the
  independent guard are what caught it.
- A scaling helper conflated rational and ℚ(i) scalars
  (TypeError mid-run); fixed uniformly.
- **My structural expectation for the stabilizer was wrong, the
  computation right:** I expected the line stabilizer to be the
  46 of Spin(10)×U(1). It is 62 — the maximal parabolic; 46 is
  its Levi factor. The orbit-dimension count (78 − 62 = 16 =
  dim_ℂ EIII) confirms the reading. The note and checks were
  rewritten to the parabolic/Levi statements.
- J8d's first form ("the parabolic preserves P±") was also wrong
  for P⁻; the correct statement is the Levi form, which passes.
