# Pascal layer — the proofs as objects

An Object Pascal (Free Pascal) expansion of this codebase: every proof
in `theorems/INDEX.md` becomes a Pascal object, and the core operators
become a Pascal class you can evaluate.

## Layout

| File | Contents |
|------|----------|
| `src/RTheory.Types.pas` | Scope labels as an enum: `PROVED / CHECKED / STRONG-EVIDENCE / OPEN / ASSERTED` |
| `src/RTheory.Operators.pas` | `TOperatorPoint` — srx, sxp, cxp, crx, urx, uxp, sawup, sawdown, FlatWave at a point x; `TCoreOperators` class functions |
| `src/RTheory.Theorems.pas` | `TProofObject` base class + `TProofRegistry`; `RegisterAllProofs` mirrors `theorems/INDEX.md` row by row (31 PROVED + T4\*) |
| `src/RTheory.Verify.pas` | Verifiable proof objects: RECIP, FLIP1/2, INVPAIR, AB, COSRAT, SINPROD, HYPER, CIRCLE, MOBIUS, E77, FW-SIGN, Q7, DURX, HARM — each re-checks its identity numerically at the 8 octant midpoints |
| `programs/VerifyProofs.lpr` | Console program: prints the registry, snapshots the operators, runs all verifiable objects |

## Scope discipline

The Pascal objects mirror the index — they never upgrade a label.

- T4\* is registered as `STRONG-EVIDENCE`, not PROVED.
- A numeric PASS in `VerifyProofs` is a CHECKED-strength confirmation
  of a result whose proof of record remains the exact Sage/algebraic
  derivation named in each object's `ProofSource`.
- SINPROD carries the **corrected** global form
  `urx·uxp = 4·|csc(2x)|` (run 150b). The index line still shows the
  older unsigned form, which holds only where `sin(2x) > 0`; the Pascal
  entry flags this in its statement text.

## Build & run

Requires Free Pascal (`fpc`):

```sh
cd pascal/programs
fpc -Fu../src VerifyProofs.lpr
./VerifyProofs
```

Exit code 0 iff every verifiable object passes.

## Status

v1 (2026-10-13): operators + every `theorems/INDEX.md` row as an
object (31 PROVED + T4\* STRONG-EVIDENCE = 32 objects), with 15
numerically verifiable objects — all 15 PASS under Free Pascal 3.2.2
(worst algebraic residual ~1e-15, derivative residuals < 1e-9).
Note: the index header's "27 PROVED" count is stale against its own
31 rows; the registry follows the rows and flags the mismatch rather
than silently correcting the index.
The Lie campaign entries from the local index (cycle theorem, block
ladder, evening batch) are the natural v2, once that index merge is
published here.
