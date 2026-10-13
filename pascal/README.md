# Pascal layer — the proofs as objects

An Object Pascal (Free Pascal) expansion of this codebase: every proof
in `theorems/INDEX.md` becomes a Pascal object, and the core operators
become a Pascal class you can evaluate.

## Layout

| File | Contents |
|------|----------|
| `src/RTheory.Types.pas` | Scope labels as an enum: `PROVED / CHECKED / STRONG-EVIDENCE / OPEN / ASSERTED` |
| `src/RTheory.Operators.pas` | `TOperatorPoint` — srx, sxp, cxp, crx, urx, uxp, sawup, sawdown, FlatWave at a point x; `TCoreOperators` class functions |
| `src/RTheory.Theorems.pas` | `TProofObject` base class (with verbatim `ScopeLabel` support for mixed-scope entries) + `TProofRegistry`; `RegisterAllProofs` mirrors `theorems/INDEX.md` row by row (31 PROVED + T4\*) |
| `src/RTheory.Verify.pas` | Verifiable proof objects: RECIP, FLIP1/2, INVPAIR, AB, COSRAT, SINPROD, HYPER, CIRCLE, MOBIUS, E77, FW-SIGN, Q7, DURX, HARM — each re-checks its identity numerically at the 8 octant midpoints |
| `src/RTheory.Lie.pas` | `RegisterLieCampaign` — the Lie campaign section of the index as objects (36 entries: foundation, block ladder, variants, evening batch, the scoped E108–E110 wall, and the not-theorems list) |
| `src/RTheory.LieVerify.pas` | A numeric Lie-closure engine (commutator brackets + RREF span basis) and 10 verifiable Lie objects: CYCLE sl(8) dim 63, SO9PATH dim 36 antisymmetric, T(x) span 4 → closure 8, E106, E107, CYCPROD, M1M2 scalar, H3H4H6 deltas, XFAM, SO3EXP θ² |
| `src/RTheory.Wall.pas` | `RegisterWall` — the wall of `docs/WALL.md` as objects: 7 OPEN gaps (complex structure, SU(8), fermion quantum numbers, physical projection Γ, coupling dynamics, E8 host, m₂/m₀) + 1 PROVED no-go (no operator content in so(8) triality) |
| `src/RTheory.WallVerify.pas` | Numeric wall probes: WALL-CHAIN4 (minimal chain closure dim 10 of 15), WALL-WICK (Wick-on-generators closure dim 28, not 63), WALL-COMMUTANT (T_a commutant dim 1 — no commuting complex structure) |
| `programs/VerifyProofs.lpr` | Console program: prints the registry (76 objects), snapshots the operators, runs all 28 verifiable objects |

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

v3 (2026-10-13): the wall is in. Registry = **76 objects** (68 + 7
OPEN wall gaps + 1 PROVED wall no-go); **28 verifiable objects, all
PASS** under Free Pascal 3.2.2 on two independent installs (this
workspace's toolchain and the laptop WSL build) — the 25 v2 checks
plus 3 wall probes: the minimal-chain stall reproduces at closure
dim 10 of 15, Wick-on-generators closes at dim 28 (so(8) conjugate,
not su(8)'s 63), and the sampled T_a commutant is 1-dimensional
(scalars only, so no complex structure commutes with the family).
The wall probes confirm at CHECKED strength; the wall entries
themselves stay OPEN until an exact derivation closes them.

v2 (2026-10-13): the Lie campaign is in. Registry = **68 objects**
(32 pre-Lie + 36 Lie campaign, incl. the wall and the not-theorems
list); **25 verifiable objects, all PASS** under Free Pascal 3.2.2 —
the 15 v1 identities (worst algebraic residual ~1e-15, derivatives
< 1e-9) plus 10 Lie checks: closure dims 63 / 36 / 8 / 8 / 3 come out
exact, basis traces and antisymmetry deviations are 0.0 at Double
precision, and the scalar Lie identities (CYCPROD, M1M2, deltas,
XFAM, θ²) hold to ≤ 1e-14.
Scope notes: mixed PROVED/CHECKED entries carry their exact split in
a verbatim scope label (never rounded to one word); E104's universal
claim stays OPEN/withdrawn; the E108–E110 wall is scoped to the
sampled minimal chain only.
v1 note, still standing: the index header's "27 PROVED" count is
stale against its own 31 pre-Lie rows; the registry follows the rows
and flags the mismatch rather than silently correcting the index.
The Lie section of `theorems/INDEX.md` is merged locally; its proof
files (cycle theorem, so(9) path, sl_n_exact) are on GitHub, the
merged index itself is still local-only.
