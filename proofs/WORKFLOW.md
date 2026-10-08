# Workflow: From CHECKED to PROVED via Exact Sage

**Purpose:** Systematically convert numerically-CHECKED results into algebraically-PROVED results using exact SageMath (symbolic/rational arithmetic, not floating-point).

**Principle:** A numerical check (SVD, optimization, sampling) establishes CHECKED. An exact computation (rational nullspace, Groebner basis, symbolic simplification) establishes PROVED.

**Machine:** toe (laptop), WSL `lampy` distro, `/usr/bin/sage` (SageMath 9.5).
**Invocation:** `toetop-ssh-kitco.sh 'powershell -c "wsl -d lampy -u root --exec sage /root/<script>.py"'`

---

## Target List (priority order)

### T1: Q4 — Commutant dimension = 1 (Cl(4,4) irreducibility)
- **Current:** CHECKED (numerical SVD, dim 1).
- **Exact statement:** The only 16×16 real matrices X with [X, E_μ]=0 for all 8 generators E_μ are scalar multiples of I.
- **Method:** Form the linear system [X,E_μ]=0 (8×256 equations, 256 unknowns) over ℚ. Compute exact nullspace dimension via rational Gaussian elimination.
- **Script:** `q4_exact.py`
- **PROVED if:** nullspace dim = 1 exactly.

### T2: Q5 — Sector intertwiner dimension = 16
- **Current:** CHECKED (numerical SVD, dim 16).
- **Exact statement:** dim{ X : X·R_μ = J_μ·X ∀μ } = 16.
- **Method:** Linear system (8×256 equations) over ℚ. Exact nullspace.
- **Script:** `q5_exact.py`
- **PROVED if:** nullspace dim = 16 exactly.

### T3: Q6a — Involution X²=+I exists in intertwiner space
- **Current:** CHECKED (numerical optimization found candidate).
- **Exact statement:** ∃ X in intertwiner space with X² = +I₁₆.
- **Method:** Parametrize 16-dim intertwiner basis {B_i}. Let X = Σ c_i B_i. Impose X²−I=0 (256 quadratic equations in 16 unknowns c_i). Solve via Groebner basis over ℚ. 
- **Script:** `q6a_exact.py`
- **PROVED if:** Groebner basis yields explicit rational (or algebraic) solution; verify X²=I exactly.
- **Note:** May require algebraic (not rational) numbers. If Groebner is too heavy, use exact verification of the numerical candidate (rational approximation + exact check).

### T4: Q6b — No X²=−I in intertwiner space (complex structure absent)
- **Current:** NOT PROVED (optimizer failed to find; not a proof).
- **Exact statement:** ∄ X in intertwiner space with X² = −I₁₆.
- **Method:** Same parametrization. Impose X²+I=0. Compute Groebner basis. If basis = [1], no solution exists (PROVED nonexistence).
- **Script:** `q6b_exact.py`
- **PROVED if:** Groebner basis is [1] (inconsistent system).
- **Note:** This is the hardest (nonexistence). May be computationally intensive. Alternative: representation-theoretic argument (if commutant is ℝ and representation is irreducible of real type, no invariant complex structure exists).

### T5: Inverse pairs c1·c4=1, c2·c3=1
- **Current:** CHECKED (numerical).
- **Exact statement:** sqrt(srx·cxp)·(1/sqrt(srx·cxp)) = 1, etc.
- **Method:** Direct symbolic simplification (trivial).
- **Script:** `t5_exact.py`
- **PROVED if:** simplifies to 1.

### T6: SU(8) 8 Z's linearly independent over ℂ
- **Current:** ASSERTED (proposal).
- **Exact statement:** The 8 Z_{i,M} = c_i + i·c_i' are linearly independent over ℂ.
- **Method:** Requires explicit c_i' via A0.3. Then form 8×8 (or 16×8) matrix, compute rank over ℂ (or ℚ(i)).
- **Script:** `t6_exact.py` (BLOCKED on defining c_i')
- **PROVED if:** rank = 8.

---

## Execution Log

| Target | Date | Result | Status |
|--------|------|--------|--------|
| T1 (Q4) | 2026-10-07 | dim=1 exact over QQ | PROVED |
| T2 (Q5) | 2026-10-07 | dim=16 exact over QQ | PROVED |
| T3 (Q6a) | 2026-10-08 | explicit X^2=+I over QQbar, verified | PROVED |
| T4 (Q6b) | 2026-10-08 | 2 sols over QQbar, both non-real; 0 real | NO REAL J (strong) |
| D2 (V struct) | 2026-10-08 | V not closed under mult; Gram=4I (orthogonal, uniform) | VECTOR SPACE, NOT ALGEBRA |
| D3 (forms) | 2026-10-08 | 2-dim symmetric invariant forms, both rank 8; 0 skew; vol^2=+I, not scalar | NATURAL 8+8 STRUCTURE |
| T5 | | _pending_ | |
| T6 | | BLOCKED (needs c_i') | |

---

## Notes
- All scripts use `from sage.all import *` and exact arithmetic (QQ, not RDF).
- Numerical candidates (from Q6) can be verified exactly by rational reconstruction.
- If a Groebner basis computation times out, record as INCOMPLETE (do not silently downgrade).
