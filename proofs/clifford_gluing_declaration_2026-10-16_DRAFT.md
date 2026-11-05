# Clifford tie-gluing: Kit's declaration of 2026-10-16 (DRAFT filing note)

Status of this note: DRAFT, prepared 2026-10-16 ~15:35 PT for Kit's review.
It records a declaration Kit made in the Algebra-first chat that evening.
Nothing here is proved by this note; scope labels are per-component below.
The exposition (`proofs/clifford_problem_exposition.md`, §4 and OQ3) still
reads "No gluing class for the ties is currently declared" — that text is
stale as of this declaration and is NOT edited here; merging is Kit's call.

## The declaration (Kit, 2026-10-16, 22:16–22:29 UTC)

**D1 — Default gluing = pass-through (DECLARED by Kit).**
A rider crossing a 45° equilibrium keeps its function (sawup stays sawup,
sawdown stays sawdown). A bounce (switching trajectories) is permitted only
if something acts on the rider to change its momentum from one trajectory
to another. Kit's stated ground: at the 45 there is nothing to reflect
against — the two functions cross in a perfect X, values equal, slopes
equal and opposite. OQ3 is thereby *supplied*, not solved: declarations are
what this question was waiting for.

**D2 — The named act: the double lift of FlatWave (DECLARED by Kit;
supporting identities as reported in the thread).**
The flipper at the 45 is FW(2x) — the dominance bit. As established in the
thread: FW itself does not change at a tie inside a quadrant; FW(2x) =
sgn(sin 4x) flips there (zeros at x = k·π/4), and FW(2x) = +1 iff
|sawdown| > |sawup| (the §9 dominance identification) — it is the only
entry of the complete octant label (FW(x), FW(2x)) that changes at a tie.
A rider coupled to FW(2x) — riding the *dominant* saw rather than a named
saw — is turned at the tie by the dominance exchange itself, not by an
external force.

**D3 — Steering data (Kit's distinction, thread-sharpened).**
The flow-following part is steered by the sign of Δ = sawdown − sawup
(dressed form: sawdown − sawup = FW·p; coefficient in [Z, Z'] = −FW·Δ·H).
Δ is a two-function comparison — it belongs to the system, not to any one
waveform — and at the tie Δ = 0, so that part's steering signal vanishes
exactly where a decision is needed and the FW(2x) flip takes it over.
The internally motivated part is steered by its own waveform (label +
own trend); its trend signal at the tie is nonzero (tie slopes
±FW·√2/2, as proved in the corpus), so pass-through is the shape of its
own signal, not a rule imposed on it.

**D4 — The two behaviors are the two chiralities (Kit's presumption,
ASSERTED — his correction of the two-species reading, 22:29 UTC).**
Not two species of rider mixed in one system: one chirality follows the
flow (Δ-sign steered, flippable at the tie), the other follows its own
waveform (internally motivated, crossing untouched) — one behavior per
handedness. The partners are the two saws, the lead is dominance, and the
lead changes hands between the first and second quadrant (late QI sawup
leads, share 0.906 at 80°; early QII sawdown leads, 0.906 at 100°). The
lead changes at every octant boundary, by two different steps: at a tie
the partners cross in the open (Δ through zero); at a quadrant seam they
swap through the pole (shares reset through the (1, 0) limit as FW flips).
Kit did NOT assign which chirality is the flow-follower; that assignment
is open (see the check brief below).

## What this settles and what it does not

- OQ3's first half (declare the class) is answered: the class is a stated
  mixture — pass-through by default, FW(2x)-flip as the sole named act,
  split by chirality per D4.
- OQ3's second half (is the §3 chirality bit invariant under the gluing
  choice, or is it the gluing in disguise?) is sharpened, not answered:
  D4 presumes the chirality grading IS the split carrying the two tie
  behaviors. That is exactly OQ1b's question in tie form.
- The Result-12 fermionic twist is relocated: under D1 it cannot be a
  behavior of the operator system; it is necessarily an *act* at the seam.
  Consistent with SC-7 (the seam carries exactly one binary degree of
  freedom — now readable as the one place an act can register).
- Ansatz 3′ (fermion research plan): the gluing-rule class it was waiting
  on now has a declared candidate — pass-through default + FW(2x) act,
  chirality-split. This is a declaration for 3′ to work with, not a proof
  of fermionic content.
- NOT claimed: any physics identification (algebra-first suspension
  stands); any statement about V32 / Spin(10) chirality (OQ1a — a real
  action of the operator algebra on V32 — remains open; the toy below is
  a miniature, not the bridge).

## Pre-registered check brief: tie behavior of the graded halves in the toy model

Offered in-chat 2026-10-16 22:29 UTC ("Say the word and I'll run it");
specified here in advance so the outcome binds either way. Kit's word is
still the trigger — this brief does not run it.

**Model (the only one with an action today).** The minimal-declaration
Clifford toy (built 2026-10-16 at Kit's order; sympy over Q(i), script
was /tmp/toy_clifford.py — not currently on disk, rebuild from the daily
memory record if run): V = C^8, functions on the octant cells of the 2π
lap; declared branch G = diag(1,1,i,i,−1,−1,−i,−i); transports T (one
cell), T² (π/2), Λ = T⁴ (π); declared action = operators act as
themselves (OQ1a in miniature). Established in the toy: G² = M_FW;
the saw-trend bit as a diagonal operator is diag(sgn sawup′) = M_FW = G²
(the bit operator IS the grading operator); grading split V+ = octants
{0,1,4,5} (QI+QIII, FW = +1), V− = octants {2,3,6,7} (QII+QIV, FW = −1).
Scope warning from the toy's own record: this grading is the FW sector
split on octant functions, NOT Spin(10)'s 16 ⊕ 16′.

**Question.** Transport a saw-labeled rider across a tie boundary inside
each graded half (QI tie: cells 0→1; QII tie: cells 2→3; QIII/QIV follow
by π-periodicity). Record, per half: (a) does the rider's saw *label*
survive the crossing under the declared action, and (b) does its
*dominance* identity survive (it cannot — dominance flips at every tie —
so the discriminating datum is (a) read against the rider definition).

**Pre-registered outcomes.**
- **A — halves differ.** One half's riders keep their label across the
  tie, the other half's swap (under the dominance-coupled reading). D4 is
  supported in the toy and the chirality assignment (which half flows) is
  determined by the run, not assumed. Scope: toy only; OQ1b's real case
  stays open pending OQ1a.
- **B — halves behave identically.** D4 is falsified in the toy (scope:
  toy only). The declaration D1–D3 stands regardless — it does not
  depend on D4.
- **C — void in the toy.** If no rider notion distinct from function
  identity can be defined under the toy's declared action (operators act
  by multiplication), the check cannot be run as posed; that is itself
  the report, and the question reverts to OQ1a. Do not relabel a
  function-identity computation as a rider check after the fact.

**Explicitly not part of the check:** the quadrant-seam crossing
(FW flip through the pole) — a different step from the tie crossing per
D4, to be briefed separately if Kit wants it.
