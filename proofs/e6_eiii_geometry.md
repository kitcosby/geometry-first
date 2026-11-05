# The EIII homogeneous geometry

`e6_eiii_geometry.py` — **41/41 exact checks over ℚ(i), exit 0.**
The last layer of the E₆ stack: the Hermitian symmetric pair of
type EIII, realized on the Jordan-side model at the basepoint
given by the Peirce decomposition of J with respect to
E = diag(1, 0, 0). Front end and conjugation machinery rebuilt
from `e6_jordan_e6.py` / `e6_real_forms.py` (checks J0–J4c and
R1–R3a are included; Killing = 4·tr₂₇ is re-verified on all
3,081 basis pairs).

## The symmetric pair

Every basis operator of 𝔢₆ is ad(Z′)-homogeneous for the grading
operator Z′ = 6L_E − 2·id, with charges in {0, +3, −3} and counts
**(46, 16, 16)** (G1). Charges add under brackets and the ±6
brackets vanish identically (G2), giving the 3-grading

𝔢₆ = 𝔭⁻ ⊕ 𝔨 ⊕ 𝔭⁺, [𝔭±, 𝔭±] = 0, [𝔭⁺, 𝔭⁻] ⊂ 𝔨, [𝔨, 𝔭±] ⊂ 𝔭±.

𝔨 (charge 0) has 1-dimensional centre ℂZ′ and 45-dimensional
derived algebra (G3), and 𝔨 is **exactly** the algebra of
operators preserving the three Peirce pieces (J₁, J½, J₀) of the
basepoint (G4 — an independent linear solve: the
piece-preserving operators have nullity 46 and the same span as
𝔨). This is the Levi of the 62-dimensional maximal parabolic
from the Jordan build, in its second role: the parabolic is
Levi ⊕ 𝔭⁺ (the orbit picture, E₆/P), while EIII's isotropy is
the Levi alone.

## The tangent module

𝔭⁺ and 𝔭⁻ are Killing-isotropic and dually paired (G5:
Killing = 4·tr₂₇ vanishes on 𝔭⁺×𝔭⁺ and 𝔭⁻×𝔭⁻; the 𝔭⁺×𝔭⁻
pairing has rank 16). As a 𝔨-module, 𝔭⁺ carries Z′-charge +3
uniformly (G6a), is irreducible (commutant of the restricted
action = scalars, G6b), and its Casimir under the J9
prescription (Gram = −tr₂₇ on the derived 45, dual bases) is
the scalar **−15/8** (G6c) — identical to the spinor piece J½
of J under the same prescription. The holomorphic tangent space
of EIII is the 16-dimensional spinor module of 𝔰𝔬(10), its
U(1)-charge fixed by Z′.

## The real form and the Hermitian structure

σ_EIII (the conjugation of the real-forms build, F4) preserves
𝔨 and 𝔭 (G7a–c), and its 78 real fixed operators split as
**46 in 𝔨 (𝔨₀) and 32 in 𝔭 (𝔭₀)** (G7d). Ad_ψ — conjugation by
the Peirce grading involution of J — is the symmetric-pair
involution: +1 on 𝔨, −1 on 𝔭, at operator level on all 78 basis
operators (G7e). The Killing form is negative definite on 𝔨₀,
signature (0, 46), and positive definite on 𝔭₀, signature
(32, 0) (G7f, G7g): the Cartan decomposition of 𝔢₆(−14),
character 32 − 46 = −14. The quotient has real dimension 32,
complex dimension 16.

The complex structure is internal to the model: H = iZ′/3 ∈ 𝔨₀
is central in 𝔨₀ (G8d), and ad(H) preserves 𝔭₀, is real in the
fixed basis (G8a), squares to −id (G8b), and is
Killing-orthogonal, K(JX, JY) = K(X, Y) on 𝔭₀ (G8c) — the
invariant Hermitian structure of EIII, with the Killing form as
metric (positive on 𝔭₀ for this noncompact dual, the bounded
symmetric domain of type EIII).

**Compact dual (G9).** σ_compact fixes 𝔨₀ pointwise in
coordinates and fixes i·𝔭₀ elementwise, so the compact form's
Cartan decomposition is 𝔨₀ ⊕ i𝔭₀ with the *same* isotropy
algebra 𝔨₀: compact EIII, E₆^compact/(Spin(10) × U(1)), is the
compact dual of the 𝔢₆(−14) space, and both are computed here
from the single complex pair. (σ_compact's bracket
compatibility is F3 of `e6_real_forms.py`; here its span and
involution are rechecked, G9a.) The real bracket relations
[𝔨₀, 𝔨₀] ⊂ 𝔨₀, [𝔨₀, 𝔭₀] ⊂ 𝔭₀, [𝔭₀, 𝔭₀] ⊂ 𝔨₀ are verified on
all fixed-basis pairs (G10: 1,035 + 1,472 + 496 pairs).

## Scope

Lie-algebra-level homogeneous geometry, exact over ℚ(i);
signatures by exact rational congruence. Cited, not computed:
the name "EIII" in the classification of Hermitian symmetric
spaces, and the group-level isotropy (Spin(10) × U(1))/ℤ₄.
Everything else — the pair and its charges, the isotropy
characterization, the tangent module and its Casimir, both real
Cartan decompositions, the signatures, the complex structure —
is computed.

With this layer the E₆ stack in geometry-first is complete:
seed algebra, the 27, the Jordan realization, the Φ weld
between the two models, the five real forms, and the EIII
geometry — all exact, all local.

## Disclosures

- Run 1 failed at G1 and aborted in G2: the charge extraction
  read the empty bracket dict {} (the charge-0 case, [Z′, T] = 0)
  as "not homogeneous", reporting counts (0, 16, 16) with the 46
  Levi operators unclassified, and G2 then crashed on the
  unclassified entries. A defect in the check code, not the
  mathematics: treating {} as charge 0 (and guarding G2) fixed
  it, and run 2 passed 41/41 with the expected (46, 16, 16).
  The run-1 log was replaced by run 2's.
- σ_compact's full bracket verification is not repeated in this
  file; it is F3 of `e6_real_forms.py` (83/83 there). G9a
  rechecks its span preservation and coordinate involution.
