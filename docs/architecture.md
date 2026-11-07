Living document — update these diagrams when adding features.

# geometry-first — Architecture

The clean-break proven core of the R Theory program: **only** theorems with complete proofs, validated tables, and the documents that describe the derived structure. No E8, no Axiom Zero, no SU(8), no Standard Model content — everything on the far side of the wall (see `docs/WALL.md`) is excluded until it is proven. Every claim carries exactly one scope label: **PROVED** (exact SageMath over ℚ/ℚ̄, Gröbner bases, exact nullspaces, or hand-verifiable symbolic algebra) or **CHECKED** (numerical only, included solely when supporting a PROVED theorem). Labels ASSERTED / OPEN / INCOMPLETE / HYPOTHESIS never appear as content here.

The corpus: `theorems/INDEX.md` (27 PROVED pre-Lie theorems, plus a separately counted Lie campaign section), `proofs/` (243 files at HEAD `83c10989` — exact Sage/Python proof scripts and theorem write-ups, including the 2026-11 sync areas below), `tables/` (29 validated CSVs + `MANIFEST.md`), `pascal/` (a Free Pascal re-verification layer: 7 source units + the `VerifyProofs` program), and `docs/` (`METHOD.md`, `WALL.md`, `PHYSICS_IMPORT_WORKFLOW.md`).

**2026-11-07 revision:** two staged syncs (Kit's go-aheads, 2026-11-02 and 2026-11-04/05) landed the local proof pile on `main`: the 2026-11-02 sync brought INDEX status fixes, the Lie campaign section, the Q7/run150b scripts, and the current Pascal layer sources; the 2026-11-05 sync added 96 proof scripts/notes (E6 campaign, so(11)/so(16) ladder, blending, Clifford gluing, Kostant-Dirac, pure spinor, fermion miniature) plus 3 corrections, including narrowing E104's scope to the five tested diagonal patterns (its universal form was falsified 2026-10-26; the falsification is recorded in the index itself). The component architecture — index, proof scripts, tables, Pascal layer, method docs — is unchanged; what changed is corpus size and the presence of new campaign working areas in `proofs/` that are not yet indexed theorems. Counts and grounding below are updated to HEAD `83c10989`.

## Context diagram

```mermaid
flowchart LR
    KIT["Kit Cosby (author and reviewer)"]
    SAGE["SageMath on toe (laptop WSL, exact arithmetic)"]
    FPC["Free Pascal compilers (VM and laptop WSL)"]
    PHYS["Known physics results (import targets)"]
    READ["Public readers"]

    subgraph SYS["geometry-first repository"]
        CORE["Proven corpus: theorems, proofs, tables, Pascal layer, method docs"]
    end

    KIT -->|"theorem candidates and review decisions"| CORE
    CORE -->|"exact proof scripts"| SAGE
    SAGE -->|"exact verification results"| CORE
    CORE -->|"Pascal verification source"| FPC
    FPC -->|"compile and run verdicts (two engines)"| CORE
    PHYS -->|"target results to prove toward (P1 to P5)"| CORE
    CORE -->|"published proven theorems and tables"| READ
```

## Level-1 DFD

```mermaid
flowchart TD
    KIT["Kit Cosby (author and reviewer)"]
    SAGE["SageMath on toe (laptop WSL, exact arithmetic)"]
    FPC["Free Pascal compilers (VM and laptop WSL)"]
    PHYS["Known physics results (import targets)"]
    READ["Public readers"]

    D1[("D1 Theorem index (theorems/INDEX.md)")]
    D2[("D2 Proof scripts and write-ups (proofs/)")]
    D3[("D3 Validated tables and manifest (tables/)")]
    D4[("D4 Pascal verification layer (pascal/)")]
    D5[("D5 Method documents (docs/)")]

    KIT -->|"candidate results"| P1["1.0 Derive structure from the native operators"]
    P1 -->|"exact proof scripts"| SAGE
    SAGE -->|"exact results"| P2["2.0 Verify proofs exactly"]
    P2 -->|"verified proof scripts"| D2
    P2 -->|"proved statements"| P4["4.0 Record theorem with scope label"]
    P2 -->|"proof obligations restated"| D4
    D4 -->|"verification source"| FPC
    FPC -->|"two-engine verdicts"| P3["3.0 Cross-check in Pascal"]
    P3 -->|"cross-check outcome"| P4
    P4 -->|"indexed theorem entries"| D1
    P1 -->|"computed diagonals and spectra"| P5["5.0 Validate and catalog tables"]
    P5 -->|"validated tables and manifest entries"| D3
    P4 -->|"wall position and method changes"| P6["6.0 Document method and the wall"]
    P6 -->|"method and wall documents"| D5
    PHYS -->|"known result (the target)"| P7["7.0 Import physics as target"]
    D1 -->|"proved structure (the proof obligation)"| P7
    P7 -->|"open import items P1 to P5"| D5
    D1 -->|"published index"| READ
    D3 -->|"published tables"| READ
    KIT -->|"review and approval"| P4
```

## ER diagram

```mermaid
erDiagram
    THEOREM {
        string id PK "T1, T2, RECIP, ONEGEN, ..."
        string statement
        string scope_label "PROVED or CHECKED (supporting only)"
        string section "Clifford, operator algebra, symmetry, Lie campaign, ..."
    }
    PROOF_SCRIPT {
        string path PK "proofs/q4_exact.py, ..."
        string kind "Sage exact, Python, write-up, TeX"
        string engine "SageMath exact, desktop Python, hand algebra"
    }
    TABLE_FILE {
        string path PK "tables/M_eps_diag.csv, ..."
        string validation_status "NC or SC, from the source catalog"
        string manifest_ref FK "tables/MANIFEST.md entry"
    }
    PASCAL_UNIT {
        string path PK "pascal/src/RTheory.Verify.pas, ..."
        string role "operators, theorems, wall, Lie, verifier"
    }
    METHOD_DOC {
        string path PK "docs/METHOD.md, docs/WALL.md, ..."
        string subject "scope labels, wall position, import workflow"
    }
    PHYSICS_IMPORT {
        string id PK "P1 to P5"
        string target "known physics result (the RHS)"
        string status "OPEN; operator side may be PROVED separately"
    }
    THEOREM ||--o{ PROOF_SCRIPT : "verified by"
    THEOREM }o--o{ TABLE_FILE : "supported by"
    THEOREM }o--o{ PASCAL_UNIT : "cross-checked by"
    THEOREM ||--o{ METHOD_DOC : "bounded by"
    PHYSICS_IMPORT }o--o{ THEOREM : "targets proof from"
```

## Grounding notes

- OBSERVED (README.md): the repo holds only PROVED results plus CHECKED items supporting a PROVED theorem; E8, Axiom Zero, SU(8), the Standard Model, dynamic coupling, and theological interpretations are explicitly excluded as not proven.
- OBSERVED (theorems/INDEX.md): 27 PROVED pre-Lie theorems, each with statement, status, and proof script; T4* is carried as strong evidence (2 solutions over ℚ̄, both non-real; 0 real) with the real-algebraic emptiness proof explicitly pending — the gap is stated, not hidden.
- OBSERVED (docs/METHOD.md): scope labels are PROVED and CHECKED only; proof standards are exact-over-numeric, no hidden axioms, failures disclosed, every PROVED theorem reproducible from a script in `proofs/`.
- OBSERVED (docs/WALL.md): the wall is located at the complex-structure gap on ℝ¹⁶ (no real J²=−I in the intertwiner space); SU(8), fermion quantum numbers, a physical projection (Γ), and coupling dynamics are all on the far side.
- OBSERVED (proofs/WORKFLOW.md): the CHECKED→PROVED conversion runs exact SageMath 9.5 on the laptop (WSL `lampy` distro, user root); the execution log records Q4/Q5/Q6a as PROVED 2026-10-07/08 and D2/D3 structure results.
- OBSERVED (git tree, HEAD `7f6eace2`): 191 files — proofs/ 145, tables/ 30 (29 CSV + MANIFEST.md), pascal/ 9, docs/ 3, theorems/ 1 (INDEX.md), plus README.md and LICENSE. **Superseded in part by the 2026-11-07 bullet below (proofs/ count and HEAD).**
- OBSERVED (pascal/README.md, pascal/src/): the Pascal layer re-verifies the corpus in Free Pascal (units for types, operators, theorems, wall, and Lie structure, plus the VerifyProofs program) as an independent second engine.
- OBSERVED (docs/PHYSICS_IMPORT_WORKFLOW.md): physics enters as the known result (the target); the repo's job is the proof. Import items P1–P5 are OPEN; where an operator side is PROVED, the physical identification remains OPEN.
- INFERRED: the Pascal layer is compiled and run on both FPC installs (this VM and laptop WSL) before a result is called verified — the two-engine cross-check pattern; the repo itself states the layer's role but not the dual-run procedure.
- OBSERVED (commits 93c0e26d…83c10989, 2026-11-02 and 2026-11-05; directory listings at HEAD `83c10989`, 2026-11-07): the 2026-11-02 staged sync (13 commits, Kit go-ahead) brought INDEX status fixes, the Lie campaign section, Q7/run150b scripts, and the current Pascal layer sources; the 2026-11-05 sync (83c10989) added 96 new proof scripts/notes — E6 campaign, so(11)/so(16) ladder, blending, Clifford gluing, Kostant-Dirac, pure spinor, fermion miniature — plus 3 corrections. `proofs/` now lists 243 entries (e.g. `blending_odd_solve.py`, `clifford_gluing_toy_check_sage.py` are present at HEAD); `tables/` still lists 30 entries (29 CSV + MANIFEST.md).
- OBSERVED (theorems/INDEX.md at HEAD `83c10989`): the index structure is unchanged — 27 PROVED pre-Lie theorems plus the separately counted Lie campaign (2026-10-12) with per-entry scope labels including CHECKED parts and OPEN items; E104 is carried at its narrowed scope ("PROVED for the five tested patterns only; the universal claim is FALSE — falsified 2026-10-26 by verified counterexample"). The 2026-11 campaign files in `proofs/` are working areas; they do not yet appear as indexed theorems, so no new THEOREM rows are drawn for them.
