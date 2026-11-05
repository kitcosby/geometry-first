# Real forms of the E₆ stack

`e6_real_forms.py` — **83/83 exact checks over ℚ(i), exit 0.**
Built on the Jordan-side model: J = H₃(𝕆_ℂ), its determinant, and
𝔢₆ as the 78 determinant-preserving operators (front end rebuilt;
checks J0–J4c as in `e6_jordan_e6.py`, including all 3,081 pairwise
commutators expressed in the basis). All five real forms of 𝔢₆ are
constructed as fixed algebras of explicit conjugations, and the
real forms of the 27 are settled: two exist (Jordan-induced), and
for the other three forms a compatible conjugation on J is
**proved obstructed**.

## The trace-form adjoint and θ

The trace form B(X, Y) = tr(X ∘ Y) on J has **diagonal Gram
(1, 1, 1, 2 ×24)** in the coordinate basis (R1a), so the B-adjoint
of an operator is an explicit weighted transpose (weights 1, 2, ½).
θ(T) := −T^{*B} is an involutive automorphism of 𝔢₆ (R2b: θ of a
bracket = bracket of θ's, on all 3,003 basis pairs, at operator
level; R2c: span preserved). Its eigenspaces are the
structure-algebra decomposition of J, each side identified
independently of θ:

- **+1 eigenspace = Der(J), dim 52.** The derivation space is
  solved from the Leibniz identity on all 378 basis pairs (R4a,
  nullity 52); every solution is θ-fixed (R4b) and B-antisymmetric
  (R4c).
- **−1 eigenspace = {L_a : tr a = 0}, dim 26.** All 26 traceless
  multiplication operators (24 slot units + 2 diagonal differences)
  lie in the algebra (R4d), are θ-anti-fixed and B-symmetric (R4e).
  Der (52) + L (26) spans the full 78 (R4f).

Also computed: the grading operator Z′ = 6L_E − 2·id lies in the
algebra with Peirce block scalars (4, 1, −2) (R1c), θ(Z′) = −Z′
(R2d), and τ₀ commutes with θ on all 78 basis operators (R2e).

## Real forms of J

Two antilinear Jordan automorphisms of J (each verified on all
378 basis pairs; det(τX) = conj(det X) on all 3,303 test elements;
fixed bases of 27 elements, closed under the product):

- **τ₀** — coefficient conjugation. Fixed algebra H₃(𝕆):
  octonion norm signature **(8, 0)**, Jordan trace-form signature
  **(27, 0)** (the Euclidean Albert algebra).
- **τ_split** — coefficient conjugation composed with the Fano
  sign character u ↦ (−1)^bit₀(u) on octonion units (an octonion
  automorphism: the unit labels 1..7 are the nonzero vectors of
  𝔽₂³ under xor, and bit₀ is a linear functional). Fixed algebra
  H₃(𝕆_split): norm signature **(4, 4)**, Jordan trace-form
  signature **(15, 12)**.

## Real forms of 𝔢₆ — all five

The Killing form is pinned to the representation first:
**Killing(X, Y) = 4 · tr₂₇(XY)**, the scalar computed from Z′ and
the identity verified on all 3,081 basis pairs, adjoint trace from
COMM against the 27-trace (R3a).

For each conjugation σ: span preserved (express), involution in
coordinates (S·conj(S) = I), all 3,003 brackets respected in
coordinates; a real fixed basis of 78 operators extracted
(u_a = T_a + σT_a, w_a = i(T_a − σT_a), real pivot extraction —
78 elements, each σ-fixed in coordinates and as matrices, complex
rank 78); the Killing matrix on it is real, and its signature is
computed by exact rational congruence (Sylvester). Results:

| σ | definition | signature | character | form |
|---|---|---|---|---|
| σ₀ | τ₀ | (26, 52) | −26 | 𝔢₆(−26), E IV |
| σ_split | τ_split | (42, 36) | +6 | 𝔢₆(6), split, E I |
| σ_compact | τ₀ ∘ θ | (0, 78) | −78 | compact 𝔢₆ |
| σ_EIII | τ₀ ∘ θ ∘ Ad_ψ | (32, 46) | −14 | 𝔢₆(−14), E III |
| σ₂ | τ₀ ∘ θ ∘ Ad_s ∘ Ad_ψ | (40, 38) | +2 | 𝔢₆(2), E II |

Ad_ψ is conjugation by the Peirce grading involution of J
(+1 on J₁ ⊕ J₀, −1 on J₁⸝₂) and Ad_s conjugation by the split sign
automorphism; both are Jordan automorphisms preserving det (R2f:
378 basis pairs, 3,303 test elements). Naming the real form from
its Killing signature/character is the standard classification of
the five real forms of E₆ (cited); every signature above is
computed, and the five signatures are pairwise distinct, so the
five fixed algebras are pairwise non-isomorphic real forms and
the list is complete.

Screening note (exploration pass, before assembly): the diagonal
sign automorphisms Ad_D, Ad_D2 and the index permutation Ad_π of
J all induce algebra involutions with the same eigenspace split
as Ad_ψ, (46, 32) — they give E III-type conjugations again, and
in fact σ_compact ∘ Ad_D ∘ Ad_ψ = σ_compact in coordinates. The
(38, 40) involution that 𝔢₆(2) requires comes from Ad_s ∘ Ad_ψ
(split sign × Peirce grading), and from nothing else tried.

## Real forms of the 27

- **σ₀ and σ_split have compatible conjugations on J** — τ₀ and
  τ_split themselves: τ(T X) = σ(T) τ(X) verified on all
  78 × 27 basis pairs (F7), and σ(Z′) = +Z′ for both (F6). Their
  fixed spaces in J are the two real Albert algebras above.
- **σ_compact, σ_EIII, σ₂ have none — PROVED obstruction (F6, F8).**
  Each of the three satisfies σ(Z′) = −Z′. A conjugation τ on J
  compatible with σ would satisfy τ ∘ Z′ = σ(Z′) ∘ τ = −Z′ ∘ τ,
  hence map the λ-eigenspace of Z′ onto the (−λ)-eigenspace. The
  spectrum of Z′ on J is (4, 1, −2) with multiplicities (1, 16, 10)
  (R1c and the Peirce dimensions); none of (−4, −1, 2) occurs, so
  τ = 0, contradiction. For these three real forms the 27 admits
  no real structure: it remains a genuinely complex representation.

## Scope

Lie-algebra-level real forms, exact over ℚ(i); signatures by exact
congruence, no floating point. The cited step is only the naming
of a real form by its character (the fivefold classification).
Not addressed: group-level real forms, and the EIII homogeneous
geometry E₆/(Spin(10) × U(1)) as a real symmetric space — the last
named layer of the stack; its complex isotropy data (the
62-parabolic with 46 Levi) is already built.

## Disclosures

- In the exploration pass the adjoint's weight ratio was first
  written inverted (G_c/G_r for G_r/G_c). It was caught there —
  σ_compact failed its span check — and fixed before the final
  file was assembled; the Z′ sign test had passed either way
  because Z′ only joins coordinates of equal Gram weight. R1b
  (adjoint is an involution on all 78) and R2b are the standing
  guards in the final file, and the assembled file passed
  83/83 on its first full run.
- One test expectation in developing the signature routine was
  wrong (the matrix [[0,2,0],[2,0,3],[0,3,0]] has determinant 0
  and inertia (1, 1, 1), not (2, 1, 0)); the routine was right.
  The controls shipped in the file (R1d) are the verified ones.
