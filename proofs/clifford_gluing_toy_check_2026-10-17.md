# Tie-gluing toy check — run of 2026-10-17 (Kit's order, queued from the main chat)

Executes the pre-registered check brief in
`clifford_gluing_declaration_2026-10-16_DRAFT.md` exactly as written.
No criterion was widened, narrowed, or redesigned. Nothing pushed.

## The brief's criteria (verbatim from the draft)

**Question.** "Transport a saw-labeled rider across a tie boundary inside
each graded half (QI tie: cells 0→1; QII tie: cells 2→3; QIII/QIV follow
by π-periodicity). Record, per half: (a) does the rider's saw *label*
survive the crossing under the declared action, and (b) does its
*dominance* identity survive (it cannot — dominance flips at every tie —
so the discriminating datum is (a) read against the rider definition)."

**Pre-registered outcomes.**
- **A — halves differ.** "One half's riders keep their label across the
  tie, the other half's swap (under the dominance-coupled reading). D4 is
  supported in the toy and the chirality assignment (which half flows) is
  determined by the run, not assumed. Scope: toy only; OQ1b's real case
  stays open pending OQ1a."
- **B — halves behave identically.** "D4 is falsified in the toy (scope:
  toy only). The declaration D1–D3 stands regardless — it does not depend
  on D4."
- **C — void in the toy.** "If no rider notion distinct from function
  identity can be defined under the toy's declared action (operators act
  by multiplication), the check cannot be run as posed; that is itself
  the report, and the question reverts to OQ1a. Do not relabel a
  function-identity computation as a rider check after the fact."

**Model (per the brief).** "The minimal-declaration Clifford toy (built
2026-10-16 at Kit's order; sympy over Q(i))": V = C^8, functions on the
octant cells of the 2π lap; declared branch G = diag(1,1,i,i,−1,−1,−i,−i);
transports T (one cell), T² (π/2), Λ = T⁴ (π); declared action = operators
act as themselves. "Scope warning from the toy's own record: this grading
is the FW sector split on octant functions, NOT Spin(10)'s 16 ⊕ 16′."

## What was computed

The toy was rebuilt from the daily memory record (memory/2026-10-16.md,
toy model entry; the original /tmp/toy_clifford.py is not on disk) and
every established fact was re-verified in exact arithmetic over Q(i),
on two independent legs that agree on all 23 checks:

- Sage 9.5 on toe (kitco wrapper stdin pipe, targeted imports):
  `clifford_gluing_toy_check_sage.py` — 23/23 OK.
- Local Python (sympy identities + exact Gaussian-rational elimination on
  Fraction pairs for the dimensions): `clifford_gluing_toy_check.py` —
  23/23 OK.

Re-verified established facts (CHECKED, exact): G² = M_FW; Ḡ² = M_FW;
{Λ, G} = 0; T²G = −i·GT² (pushforward convention, matching the record);
{T², M_FW} = 0; the saw-trend bit as a diagonal operator is
diag(sgn sawup′) = M_FW = G² (the bit operator IS the grading operator);
dim ⟨Λ, G⟩ = 8 over Q(i); dim ⟨T, G⟩ = 64 = full M(8, C); grading split
V+ = {0,1,4,5} (FW = +1), V− = {2,3,6,7} (FW = −1).

Tie-local facts (CHECKED, exact, both legs):
- Both crossings exist and are intra-half: T e₀ = e₁ (QI tie, inside V+);
  T e₂ = e₃ (QII tie, inside V−).
- Every declared diagonal operator — G, Ḡ, M_FW, and the trend bit — is
  **constant across both ties** (G: cells 0,1 both 1; cells 2,3 both i;
  M_FW and trend: +1,+1 and −1,−1). No declared multiplication operator
  registers a tie crossing at all; the toy's only saw datum (the trend
  bit) registers quadrant seams, not ties.
- The FW(2x) diagonal diag((−1)^k) — D2's dominance flipper, the datum
  whose flip defines the flow-following rider's tie behavior in D2/D3 —
  is **not among the toy's declared operators** (verified distinct from
  G, Ḡ, M_FW, trend). No Δ operator, no saw vectors, and no
  comparison/selection map is declared in the model.

## Determination against the criteria

For datum (a) to be posed, the toy must carry a rider notion: a state
whose saw label is data distinct from the identity of a function on
cells, transported by the declared action so the label after the crossing
can be read per half. In the toy as declared it cannot:

1. V = C^8 carries exactly one complex dimension per octant cell; there
   is no label register. A state is exactly a function on cells, so any
   "label" definable inside the declared model is a function of that
   function — function identity, which outcome C expressly excludes as a
   rider check.
2. The only saw-derived operator in the toy is the trend bit, and it
   equals the grading operator M_FW (tie-constant, above). The
   dominance datum that distinguishes the two rider species in D4
   (Δ-sign / FW(2x)) is not represented by any declared operator, and
   supplying one would widen the model beyond the brief.

The transport of the saw functions themselves was computed as part of the
rebuild (the crossings above); it trivially preserves function identity
in both halves and is **not** reported as outcome B, per the brief's
explicit instruction.

## Verdict

**Outcome C — void in the toy.** The check cannot be run as posed; the
question reverts to OQ1a. (In MATCH / NO MATCH / INCOMPLETE terms:
INCOMPLETE — not a no-match.)

Scope labels:
- D4 (the two tie behaviors are the two chiralities): remains **ASSERTED**
  — Kit's presumption. The toy gives it no support (not A) and does not
  falsify it (not B). Status unchanged by this run.
- D1–D3 (pass-through default, FW(2x) act, steering data): untouched —
  the brief states they stand regardless, and C does not bear on them.
- OQ1b's real case: stays open pending OQ1a (the brief's own routing).
- Any future run of this question needs a newly declared model with a
  rider state distinct from function identity (e.g. a declared label
  register and coupling) and a fresh pre-registration; that is a new
  brief, not a continuation of this one.

Artifacts (this directory, not pushed):
`clifford_gluing_toy_check.py` (local exact leg),
`clifford_gluing_toy_check_sage.py` (Sage/toe leg),
`clifford_gluing_toy_check_2026-10-17.md` (this note).
