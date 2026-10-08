# Theorems Index — All PROVED

Each entry: statement, status, proof script in `../proofs/`.

## Clifford structure

| ID | Statement | Proof |
|----|-----------|-------|
| T1 | Cl(4,4) 16-spinor is irreducible (commutant dim 1, exact over ℚ) | `proofs/q4_exact.py` |
| T2 | Sector intertwiner space dim = 16 (exact over ℚ) | `proofs/q5_exact.py` |
| T3 | ∃ X in intertwiner space with X²=+I (explicit, over ℚ̄, verified) | `proofs/q6a_groebner.py` |
| T4* | X²=−I: 2 solutions over ℚ̄, both non-real; 0 real found | `proofs/q6b_groebner.py`, `proofs/q6b_solutions.py` |

*T4* is strong evidence (not a complete real-algebraic proof of emptiness).

## Operator algebra

| ID | Statement | Proof |
|----|-----------|-------|
| RECIP | sawup + sawdown = FlatWave (exact) | algebra |
| ONEGEN | Every operator is a function of z=srx alone | algebra |
| FLIP1 | crx = \|sec x\| − tan x (exact) | sec²−tan²=1 |
| FLIP2 | sxp = \|csc x\| − cot x (exact) | csc²−cot²=1 |
| INVPAIR | c₁·c₄=1, c₂·c₃=1 (exact) | definitions |

## Symmetry

| ID | Statement | Proof |
|----|-----------|-------|
| V4ACT | V₄ acts on (srx,sxp,cxp,crx) via ι=(0 1)(2 3), τ=(0 3)(1 2), ι_d=(0 2)(1 3) | trig + reciprocal |
| V4GRP | {id,ι,τ,ι_d} is Klein V₄ (closure verified) | permutation computation |
| CHIF | FlatWave character χ_F=(+1,−1,−1,+1) | FW definition |

## Möbius and double-angle

| ID | Statement | Proof |
|----|-----------|-------|
| MOBIUS | cxp = (srx+ε)/(ε·srx−1), ε=FlatWave | equivalent to RECIP |
| KRED | u=(urx,uxp)=K·P, K=[[1,0,0,−1],[0,−1,1,0]] | definitions |
| AB | A=srx−sxp=2cot x, B=cxp−crx=2tan x | FLIP1, FLIP2 |
| COSRAT | K₂/D₂=cos(2x), K₂=A−B, D₂=A+B | AB + trig |
| SINPROD | sin(2x)=4/(urx·uxp); urx·uxp=D₂ | algebra (Q1), numeric globally |
| HYPER | D₂²−K₂²=16 | COSRAT |
| CIRCLE | (sin2x,cos2x)=(4/D₂,K₂/D₂) | HYPER |

## Conic geometry

| ID | Statement | Proof |
|----|-----------|-------|
| QCONIC | Q_conic=diag(16,4,−1); xᵀQx=0 ⟺ 4H_R²+V_R²=1/4 | algebra |
| QHYP | Q_hyp=diag(1,−1,−16); yᵀQy=0 ⟺ Σ²−Δ²=16 | algebra |
| HYP2ELL | (Σ,Δ)↦(1/Σ,Δ/2Σ)=(H_R,V_R) maps hyperbola→ellipse; seam Σ→∞ ↦ (0,1/2) | algebra |
| VERON | ν₂(cos x,sin x) satisfies X₀X₂−X₁²=0 | algebra |
| SPIN1 | ⟨J_z⟩=cos2x, ⟨J_x⟩=sin2x; ⟨J_z⟩/2=V_R, ⟨J_x⟩/4=H_R | spin-1 representation |

## Dynamics

| ID | Statement | Proof |
|----|-----------|-------|
| Q7 | dS/dx=(srx·\|csc x\|−crx·\|sec x\|)·S², S=sawup (Bernoulli, exact per-octant) | `proofs/q7_dSdx.py` |
| DURX | d(urx)/dx=crx·\|sec x\|−srx·\|csc x\| (exact) | Q7 + S=1/urx |
| RUN144 | 8 (M,P,G) types, uniquely identified; G flips for c₁,c₄ at octant boundaries | exact lemma + census |
| HARM | F=(sin2x/4,cos2x/2), dF/dx=[[0,1],[−4,0]]·F | trig |

## Intertwiner geometry (D2/D3)

| ID | Statement | Proof |
|----|-----------|-------|
| D2CLOSE | V not closed under 16×16 multiplication (B₀·B₀∉V; I₁₆∉V) | `proofs/d2_structure.py` |
| D2GRAM | Gram G[i,j]=tr(B_iᵀB_j)=4·I₁₆ (orthogonal, uniform norm 2) | `proofs/d2_structure.py` |
| D3FORM | Exactly 2-dim Spin(4,4)-invariant symmetric forms, both rank 8; 0 skew | `proofs/d3b_forms.py` |
| D3VOL | vol=E₀···E₇: vol²=+I, not scalar, anticommutes with E_μ | `proofs/d3b_forms.py` |

**Total: 27 PROVED (+1 strong, T4*).**
