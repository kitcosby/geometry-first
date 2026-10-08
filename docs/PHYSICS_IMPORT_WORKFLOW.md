# Workflow: Physics Import as Trigonometric Identity

**Reframing (Kit, 2026-10-08):** Treat each import as a trigonometric identity. The established physics is the **known result** (the RHS). Our job is to **provide the proof** — starting from operator definitions, via valid steps, arriving at the physical law.

This is not "translation" (loose). It is **proof obligation**:
- **Given:** R Theory operators (srx, sxp, cxp, crx, urx, uxp, saws, FW) and their proven properties.
- **Prove:** [operator expression] = [known physical law].

Like sin²x+cos²x=1: you don't "translate" it. You prove it from definitions.

## Proof Template

For each target:

1. **Known result** (established physics, with source). This is the RHS — not to be proven, but to be matched.
2. **Operator statement** (conjecture in R Theory notation). This is the LHS — what we claim equals the RHS.
3. **Proof** (from operator definitions, via valid algebraic/analytic steps).
4. **Status:** PROVED (proof complete) / OPEN (no proof yet) / INCOMPLETE (partial) / FALSIFIED (operator statement refuted).

**Rules:**
- The known result is never "derived" — it is the target.
- The proof must start from operator definitions and proven theorems only.
- No physical constants (G, ħ, c) enter the operators. They live in Γ (calibration), not in the proof.
- If the proof needs an unproven lemma, mark INCOMPLETE and state the lemma.
- If the operator statement is wrong, mark FALSIFIED — do not adjust the physics to fit.

---

## Import Queue

### P1: Newtonian Gravitation
- **Known (RHS):** F = Gm₁m₂/r². Equivalently, potential Φ(r) = −GM/r satisfies ∇²Φ = 4πGρ.
- **Operator statement (LHS, conjectured):** Φ ∝ −urx, via the phase-to-radius map r(x)=R√(2/urx(x)).
- **Proof obligation:** From urx definitions, show that −urx (as a function of r) satisfies the Poisson equation (up to Γ-supplied G).
- **Blocker:** The r(x) map is definitional, not derived. The proof must either derive r(x) or explicitly assume it as a Γ input.
- **Status:** OPEN.

### P2: General Relativity
- **Known (RHS):** G_μν = 8πT_μν. Schwarzschild: ds² = −(1−2M/r)dt² + (1−2M/r)⁻¹dr² + r²dΩ².
- **Operator statement (LHS, conjectured):** TBD. Candidate: the conic chain (Q_hyp → Q_conic) mirrors causal → spatial structure; the phase-to-radius map gives r.
- **Proof obligation:** TBD (needs P1 first).
- **Status:** OPEN.

### P3: Electromagnetism
- **Known (RHS):** Maxwell's equations. Coulomb: E = Q/(4πϵ₀r²)r̂.
- **Operator statement (LHS, conjectured):** TBD. Candidate: E/B from (cxp, crx) via the E/B witnesses (cxp·crx=1 factorizes mass-shell).
- **Proof obligation:** TBD.
- **Status:** OPEN.

### P4: Quantum Mechanics
- **Known (RHS):** Dirac equation (iγ^μ∂_μ − m)ψ = 0.
- **Operator statement (LHS, conjectured):** Dirac operator = Σ E_μ ∂_μ using Cl(4,4) generators E_μ (standard construction).
- **Proof obligation:** Show the Cl(4,4) E_μ satisfy the Dirac algebra {E_μ,E_ν}=2η_μν (this is by construction — PROVED in T1 setup). Then: what does R Theory's specific 16 add beyond standard Dirac theory?
- **Status:** OPEN (the algebra is standard; the R Theory content is the question).

### P5: Thermodynamics
- **Known (RHS):** dS ≥ 0 (second law). S = k_B log Ω (Boltzmann).
- **Operator statement (LHS, conjectured):** S = sawup = 1/urx; dS/dx = (srx|csc| − crx|sec|)S² (Bernoulli, PROVED in Q7).
- **Proof obligation:** The operator dynamics (Q7) is PROVED. The missing step is the physical identification: what is "temperature" in operator terms? What is the Γ map from x (phase) to t (time)?
- **Status:** INCOMPLETE (operator side PROVED; physical identification OPEN).

---

## What "proof" means here

We are not proving the physics (it is established). We are proving the **identity**: that the operator expression equals the physical law. The proof lives entirely in the operator algebra; the physics is the fixed target.

If we cannot prove the identity, we do not have an import. We have an OPEN problem.
