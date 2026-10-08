# Geometry-First Unification

**A clean break.** This repository contains ONLY proven results — theorems with complete proofs, tables with validated entries, documents describing derived structure. No hypotheses. No E8. No speculation. No "representation hosts."

## Method

**Geometry first:** Start from the native operators (srx, sxp, cxp, crx) and derive, using exact computer algebra (SageMath), how far the mathematics goes on its own. Mark the wall where derivation stops. Add nothing beyond it.

**Scope labels (every claim carries one):**
- **PROVED** — complete mathematical proof (exact Sage, symbolic algebra).
- **CHECKED** — numerical verification, no proof yet. (Included only when supporting a PROVED theorem; labeled as such.)
- Nothing else enters. No ASSERTED. No OPEN. No hypotheses.

## What is here

### Theorems (`theorems/`)
27 PROVED theorems, each with statement, proof sketch, and Sage verification script in `proofs/`. See `theorems/INDEX.md`.

**Clifford structure:**
- T1: Cl(4,4) 16-spinor is irreducible (commutant dim 1).
- T2: Sector intertwiner space is 16-dimensional.
- T3: An involution X²=+I exists in the intertwiner space.
- T4*: No REAL X²=−I in intertwiners (2 complex solutions, 0 real — strong evidence, proof of emptiness via real algebraic geometry pending).

**Operator algebra:**
- Reciprocal-partition: sawup + sawdown = FlatWave.
- One-generator collapse: every operator is a function of z=srx.
- Flip identities: crx = |sec|−tan, sxp = |csc|−cot.
- Inverse pairs: c₁·c₄=1, c₂·c₃=1.

**Symmetry:**
- V₄ action on the primitive quartet (permutation matrices).
- V₄ Klein-four closure.
- FlatWave character χ_F = (+1,−1,−1,+1).

**Möbius and conic geometry:**
- Möbius: cxp = (srx+ε)/(ε·srx−1), ε=FlatWave.
- Double-angle: sin(2x)=4/(urx·uxp); K₂/D₂=cos(2x).
- Hyperbola: D₂²−K₂²=16; circle map.
- Conic matrices: Q_conic=diag(16,4,−1), Q_hyp=diag(1,−1,−16).
- Veronese quadric and spin-1 expectations.

**Dynamics:**
- dS/dx Bernoulli exact form; d(urx)/dx exact form.
- RUN144 8-type classification.
- Harmonic evolution dF/dx.

**Intertwiner geometry (D2/D3):**
- V is a vector space, NOT an algebra (not closed under multiplication).
- Gram matrix = 4·I₁₆ (orthogonal basis).
- Two invariant symmetric bilinear forms, both rank 8; zero skew.

### Tables (`tables/`)
Validated computational tables supporting the theorems. E8-specific tables EXCLUDED (see below). Each table carries its validation status (NC/SC) from the source catalog.

### Proofs (`proofs/`)
Exact SageMath scripts. Every PROVED theorem has a corresponding script that anyone can run to verify.

### Documents (`docs/`)
- `WALL.md` — The derive-to-the-wall program: what the base gives, where it stops.
- `METHOD.md` — Scope labels, proof standards, what "geometry first" means.

## What is NOT here (deliberately excluded)

- **E8** — all E8-specific tables (30380, adjoint, simple roots), E8(−24) gateway claims, 30380 branching. E8 sits on the far side of the wall: an inspired hypothesis, not a derived consequence.
- **Axiom Zero** — the "Normal and Ideal" imaginary sector. Not used in any proof here.
- **SU(8)** — the 8 Z's proposal. ASSERTED, not PROVED. Excluded.
- **Standard Model** — 10⊕5̄⊕1, quantum numbers, fermion content. Downstream of the wall.
- **Dynamic coupling, theological interpretations** — hypotheses, not theorems.

These are not declared false. They are declared *not proven* and therefore not in this repository. If and when they are proven, they enter.

## Provenance

Built 2026-10-08 by Kit Cosby and Muse Cosby, as a clean break from the E8-augmented R Theory program. The trigonometry thread continues the manuscript audit separately. This repository is the geometry-first attempt: what the operators prove about themselves, nothing more.
