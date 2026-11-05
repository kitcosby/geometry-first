# The exceptional Jordan algebra H₃(𝕆_ℂ), built intrinsically

`e6_jordan.py` — 15/15 exact checks over ℚ(i), exit 0. This is the
Jordan-algebra layer of the E₆ stack: the 27-dimensional space of the
previous build (`e6_rep27.md`) carries a classical second life as
3×3 Hermitian matrices over the complexified octonions, with a Jordan
product and a determinant cubic. This build constructs that structure
from the octonions up, in its own coordinates, and verifies it
completely. The identification with the representation-theoretic 27
(the intertwining map Φ with det∘Φ ∝ C) is the remaining bridge and
is **not** built here — see the last section for exactly what it
requires and what this build already pins down about it.

## Construction

Octonions 𝕆_ℂ over ℚ(i) from the Fano table (lines (1,2,3), (1,4,5),
(1,7,6), (2,4,6), (2,5,7), (3,4,7), (3,6,5); eᵢ² = −1, cyclic products
+1). Conjugation is the ℂ-linear extension of e₀ ↦ e₀, eᵢ ↦ −eᵢ.

J = {3×3 Hermitian matrices over 𝕆_ℂ}, coordinates:

- diagonal d₀, d₁, d₂ ∈ ℂ (indices 0, 1, 2),
- slot A = X₀₁ (indices 3–10), slot B = X₀₂ (indices 11–18,
  X₂₀ = B̄), slot C = X₁₂ (indices 19–26, X₂₁ = C̄),

each slot a full octonion: 3 + 8 + 8 + 8 = 27. Product
X ∘ Y = (XY + YX)/2 with octonion multiplication in the entries.

## What is verified

**Octonions (S0).** Alternativity [x,x,y] = [x,y,y] = [x,y,x] = 0 on
all basis pairs — the property the Jordan identity for H₃ rests on.

**Product (S1, S2).** On all 378 basis pairs the product is again
Hermitian with complex diagonal (closure in J) and commutative, and
I = diag(1,1,1) is a two-sided identity on all basis elements.

**Jordan identity (S3, S3b, S4) — PROVED.** The identity
(X²∘Y)∘X = X²∘(Y∘X) is verified directly on all 729 basis pairs and
on all pair-sum X against all basis Y (9,477 evaluations). In operator
form it is the vanishing of [L_X, L_{X²}] on Y. Its full polarization
F̃(X₁,X₂,X₃) = Σ_cyc [L_{X₁}, L_{X₂∘X₃}] — a symmetric trilinear
operator-valued form with F̃(X,X,X) = 3[L_X, L_{X²}] — vanishes on all
3,654 basis triples with repetition (sparse exact operator algebra,
matrices as {(row, col): ℚ(i)} dicts). A symmetric multilinear form is
determined by its basis values, so [L_X, L_{X²}] = 0 for **all** X,
hence the Jordan identity for all X, Y. Together with commutativity
and the identity element, J is a Jordan algebra — verified, not cited.

**Determinant (S5).** det is defined by the trace formula
det X = (1/3)tr X³ − (1/2)(tr X²)(tr X) + (1/6)(tr X)³. The explicit
H₃ formula

det X = d₀d₁d₂ − d₀·N(C) − d₁·N(B) − d₂·N(A) + 2·Re₀((A·C)·B̄),
N(u) = u·ū,

agrees with the trace formula on the polarization-complete test set
(all basis elements, pair sums, and distinct triple sums: 3,303
elements; two cubics agreeing there agree everywhere). Exactly the
cyclic class of triple orderings matches — (A·C)·B̄, A·(C·B̄),
(B̄·A)·C — and the other three orderings fail: the computation
recovers the classical fact that the real part of an octonion triple
product is invariant under cyclic permutation and reassociation, but
not under reversal. The det tensor on basis triples has 89 nonzero
entries (1 diagonal triple, 24 norm patterns, 64 triple-product
patterns).

**Sharp map and adjoint identity (S6, S6p) — PROVED.**
X^♯ = X² − (tr X)X + σ(X)I, σ(X) = ((tr X)² − tr X²)/2. The adjoint
identity X ∘ X^♯ = det(X)·I is verified directly on all basis elements
and pair sums, and its full polarization
(1/6)Σ_perm X_{p₁} ∘ (X_{p₂} × X_{p₃}) − det̃(X₁,X₂,X₃)·I
(with × the polarized sharp map) vanishes on all 3,654 basis triples,
hence for all X.

**Peirce frame (S7).** For the idempotent E = diag(1,0,0) (checked:
E∘E = E, E^♯ = 0, det E = 0 — rank one), L_E satisfies
L_E(2L_E − I)(L_E − I) = 0 with Peirce spaces

J₁ = 1 (the d₀ line), J₁⸝₂ = 16 (slots A, B), J₀ = 10 (d₁, d₂, slot C),

dims read off the spectral projectors. Peirce product rules hold on
all 361 cross-block basis pairs: J₁ ∘ J₀ = 0, J₁⸝₂ ∘ J₁⸝₂ ⊆ J₁ ⊕ J₀,
J₀ ∘ J₁⸝₂ ⊆ J₁⸝₂, J₀ ∘ J₀ ⊆ J₀.

**Grading (S8).** The operator Z_J = 6L_E − 2·id has Peirce
eigenvalues **4, 1, −2** on (J₁, J₁⸝₂, J₀) — exactly the Z′ charges of
the (1, 16, 10) pieces of the rep27 build. And det is Z_J-invariant:
every nonzero det monomial has Peirce grades summing to 1, so its
charge sum is 6·1 − 2·3 = 0 (checked entry by entry on the tensor).

## Scope and the remaining bridge

- This is the **complex** Jordan algebra. No real form is chosen or
  built.
- E₆ preserves the determinant, not the Jordan product (the product's
  automorphism algebra is 𝔣₄-type). Nothing here claims the 𝔢₆ action
  of `e6_rep27.py` acts by Jordan derivations — it does not.
- **The Φ bridge (not built).** What remains is a linear isomorphism
  Φ: V₂₇ → J intertwining the 𝔨 = 𝔰𝔬(10) ⊕ ℂZ₆ actions with
  det_J ∘ Φ = c·C, c ≠ 0, C the rep27 cubic. This build pins two of
  its anchors: the grading forces Φ to send the charge-(4, 1, −2)
  decomposition of V₂₇ to the Peirce decomposition (J₁, J₁⸝₂, J₀);
  in particular the singlet line goes to ℂE, **not** to ℂI — and
  consistently det(E) = 0 = C on the pure singlet, while the Jordan
  identity decomposes as I = E + (I − E) across J₁ ⊕ J₀, predicting
  Φ⁻¹(I) = singlet + one specific vector-piece element. The missing
  ingredient is the dictionary between the Pauli-recursion spinor
  basis (subsets of {1,…,5}) and the octonion-pair model
  J₁⸝₂ = 𝕆_ℂ², i.e. the Spin(9) ⊂ Spin(10) octonionic model of the
  16 and the 10. That is the next layer, together with the real
  forms and the EIII homogeneous geometry.

## Disclosures

First run: 12/15, three failures — all machinery, all diagnosed
against the direct computations (which passed throughout) and fixed:

1. S4 failed at (3,3,3): the sparse `matmul` grouped the second
   factor's entries by column instead of by row, so every commutator
   was wrong. The direct Jordan checks S3/S3b passing on 10,206
   evaluations is what localized the bug to the operator machinery.
2. S6p failed at (0,19,19): the polarization sum deduplicated
   permutations via `set()`, underweighting repeated indices by
   their multiplicity (hand computation at that triple: correct value
   −I/3 − (−1/3)I = 0). Fixed by summing all 6 permutations with
   multiplicity.
3. S5's assertion "exactly one ordering matches" was wrong as a
   prediction: three match — the cyclic class above. The expectation
   was corrected by the computation, and the matching set is now
   itself a checked statement (the three non-cyclic orderings fail
   on the test set).
