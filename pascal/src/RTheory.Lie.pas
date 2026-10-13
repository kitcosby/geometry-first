unit RTheory.Lie;

{ The Lie campaign (2026-10-12) as Pascal objects -- v2 of the Pascal
  layer.

  RegisterLieCampaign mirrors the "Lie campaign" section of
  theorems/INDEX.md row by row, with per-entry scope labels preserved
  exactly, including:
    - the mixed PROVED/CHECKED entries (CYCPROD, M1M2, SO3EXP,
      WICKWALL) -- carried verbatim in ScopeLabel, never rounded to a
      single word;
    - E104, PROVED only for its five tested patterns, with the
      universal claim OPEN (withdrawn 2026-10-12);
    - the E108-E110 wall, scoped to the sampled minimal chain only --
      NOT a general sl(n) wall;
    - the explicitly-not-theorems list (Maxwell dictionary ASSERTED;
      O_prop, the universal parameter theorem, and the Z/W moduli
      OPEN).

  Source note: the Lie section of theorems/INDEX.md is merged locally
  (Kit approved 2026-10-13) and its proof files are on GitHub
  (proofs/sl8_cycle_theorem.tex, proofs/so9_path_theorem.tex,
  proofs/sl_n_exact.tex), but the merged INDEX.md itself is still
  local-only at the time of writing. }

{$mode objfpc}{$H+}

interface

uses
  RTheory.Types, RTheory.Theorems;

procedure RegisterLieCampaign(ARegistry: TProofRegistry);

implementation

procedure RegisterLieCampaign(ARegistry: TProofRegistry);

  procedure P(const AId, AStmt, ASrc: string);
  begin
    ARegistry.Add(TProofObject.Create(AId, AStmt, psProved, ASrc));
  end;

  procedure L(const AId, AStmt, ASrc, ALabel: string);
  var
    O: TProofObject;
  begin
    O := TProofObject.Create(AId, AStmt, psProved, ASrc);
    O.ScopeLabel := ALabel;
    ARegistry.Add(O);
  end;

begin
  { Foundation and main generation theorems }
  P('FOUNDATION',
    'Pythagorean Foundation: operator definitions, reciprocal conjugates, UNA forms, urx*uxp = 4*|csc(2x)|, FlatWave = sgn(sin 2x), |saw_r| + |saw_x| = 1 partition of unity, half-angle forms, octant dominance cycle',
    'proofs/pythagorean_foundation.md');
  P('E77',
    'Double-angle identities: sin(2x)/4 = 1/(urx+uxp), cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp)',
    'analytic (t = tan(x/2)); Foundation');
  P('CYCLE',
    'For any n and any nonzero a_i, X_i = a_i*E_{i,i+1} (indices mod n) generate sl(n); at n = 8 with the eight operators (srx, sxp, cxp, crx, urx, sawup, uxp, sawdown): sl(8), dim 63, rank 7, A7, Killing (35,28,0)',
    'proofs/sl8_cycle_theorem.tex (human proof, Kit)');
  P('SO9PATH',
    'The same eight operators as Y_i = f_i*(E_{i,i+1} - E_{i+1,i}) in 9x9 generate so(9), dim 36, rank 4, B4 -- the algebra is fixed by quiver/symmetry type, not by the values (nonzeroness only)',
    'proofs/so9_path_theorem.tex (human proof, Kit)');
  P('MOBSL6',
    'The six anharmonic transforms of z = srx are units in Q(z) and generate sl(6, Q(z)) by the cycle argument; specializing z = sqrt(3) (x = 60 deg) gives sl(6,R) -- an independent route to sl(6) that does not use the block ladder',
    'recorded 2026-10-12; dependency graph v2');

  { Block ladder -- exact cross-verification (Sage, exact over algebraic closure) }
  P('E70-E72',
    'sl(3,R) via T_a: rank 8, all 28 commutators close, Killing (5,3,0), C2 = 4/9, Cartan 2-dim, A2 roots',
    'proofs/e70_exact.py, proofs/e71_exact_casimir.py, proofs/e72_exact_rest.py');
  P('E88',
    'sl(4,R) via blocks: rank 15, 105 commutators close, Killing (9,6,0)',
    'E88 (see proofs/sl_n_exact.tex)');
  P('E90B-C',
    'sl(5,R): rank 24, 276 commutators close, Killing (14,10,0)',
    'E90b, E90c (see proofs/sl_n_exact.tex)');
  P('E91-E94',
    'sl(6,R): rank 35, 595 commutators close (E91), independent rank-35 re-verification (E94); identification by dimension count. The E92c numeric Killing (18,17) is discarded as unstable; the exact signature is (20,15,0)',
    'E91, E94 (see proofs/sl_n_exact.tex)');
  P('E107',
    'sl(2,R) building block: [[0,srx],[cxp,0]], rank 2->3, closure, Killing (2,1)',
    'proofs/e107_sl2.py (GitHub sha b682e95a8bd0)');
  P('E106',
    'Minimal sl(3,R): [[0,srx,0],[cxp,0,crx],[0,sxp,0]], rank 4->8, closure, Killing (5,3,0); 2 independent functions in 4 positions, via two overlapping sl(2)s (E107)',
    'proofs/e106_minimal.py (GitHub sha 191e285ecdc4)');

  { sl(3) variants and parameter results }
  P('E96', 'FlatWave matrix (traceless-corrected) generates sl(3,R): rank 8, 28 commutators',
    'dependency graph E96');
  P('E97', 'Corner variants (0,0), (sin x, cos x), (sin 2x, cos 2x) all give sl(3,R) -- corners free for the tested cases',
    'dependency graph E97');
  P('E98', 'Sawtooth matrix (diagonal sawup, sawdown, FlatWave) generates sl(3,R): rank 8, closure',
    'dependency graph E98');
  P('E93', 'A specific nonzero v-entry still gives sl(4,R): rank 15, Killing (9,6) -- robustness evidence for that case, not a universal proof',
    'dependency graph E93');
  P('E99', 'The w = -v sign flip preserves sl(4,R), Killing (9,6)',
    'dependency graph E99');
  P('E105',
    'Zero-diagonal sl(3,R): six off-diagonal operators, zero diagonal; rank 4->8, closure, Killing (5,3,0). One of E104''s five tested diagonal patterns',
    'dependency graph E105');
  P('E100',
    'Structural no-go: an 8-dimensional real Lie subalgebra of sl(3,R) equals sl(3,R), so Killing (0,8,0) is unreachable by real traceless matrices. The 100-permutation search illustrates this; it is not itself the proof',
    'proofs/e100_e103_e104_proved.md');
  P('E103',
    'T_a(x) lies in the fixed 4-dimensional span {H, E12, E21, E23-E32} for all x; no sampling of the T_a family exceeds rank 4',
    'proofs/e100_e103_e104_proved.md');
  L('E104',
    'Five specific traceless diagonal patterns -- (urx,-urx,0), (1,-1,0), (srx,cxp,-srx-cxp), (0,0,0), (urx,uxp,-urx-uxp) -- each generate rank 8',
    'proofs/e100_e103_e104_proved.md',
    'PROVED for the five tested patterns only; universal claim OPEN (withdrawn 2026-10-12)');
  P('E101',
    'Wick rotation to su(3): 3 antisymmetric + 5 i*symmetric generators, rank 8, closure, Killing (0,8). Imports i, dagger, and the target real form -- not native to the operators',
    'dependency graph E101');

  { Evening batch (2026-10-12) }
  P('H3H4H6',
    'Deltas: quarter-turn swaps families -- srx(x+pi/2) = crx(x), cxp(x+pi/2) = sxp(x) -- so urx and uxp are the negated quarter-turn deltas of srx and cxp; doubling delta exact per FlatWave branch; reciprocal delta',
    'proofs/harmonics_h3_h4_h6_deltas.md');
  P('XFAM',
    'Cross-family diagonal: (srx+crx) - (sxp+cxp) = 4*cot(2x), the E77 numerator',
    'proofs/cross_family_sum_difference_DRAFT.md');
  P('JACPAIR',
    'Jacobi pair transparency: every complete reciprocal pair contributes exactly 1 to bracket/Jacobi coefficients; yields the Chevalley normalization for sl(8)',
    'proofs/jacobi_pairs.md');
  P('JACXFAM',
    'Cross-family couplings: connector brackets carry srx*cxp as a Cartan-scale coupling (no cancellation), while the connector cross-products telescope to 1 around the full cycle',
    'proofs/jacobi_cross_family.md');
  L('CYCPROD',
    'Cycle product C^8 = I: prod a_i = 1 (four reciprocal pairs), so C is an exact order-8 Coxeter element, characteristic polynomial lambda^8 - 1, spectrum independent of x',
    'proofs/expmap_operators.md',
    'PROVED (items 1-3, 6) / CHECKED (items 4-5)');
  L('M1M2',
    'M1-M2 scalar invariant: tr(M1*M2) = 3 identically, equivalently urx*uxp - srx*crx - sxp*cxp = 2. (Recorded near-miss: det M2 at 30 deg is 1.6188, NOT the golden ratio)',
    'proofs/cross_m1_m2.md',
    'PROVED (trace identity) / CHECKED (entry forms)');
  P('QFLIP',
    'Quadrant flip = perfect-square swap: the pi-shift is w -> i*w, i.e. t -> -1/t; alpha = sgn(Im zeta) selects which Euler square sits upstairs; at block level M1_QIII = diag(1,-1)*M2_QI*diag(1,-1). The mechanism Kit named "self propagation"',
    'proofs/quadrant_flip_polar.md');
  L('SO3EXP',
    'SO(3) exponential: Rodrigues closes exp(Y) for omega = (sxp, cxp, uxp), and theta^2 = 2*(sxp^2 - sxp*cxp + cxp^2), the A2 Eisenstein norm of the cores (since uxp = cxp - sxp)',
    'proofs/expmap_uxp_so3.md',
    'PROVED (structure, theta^2 identity) / CHECKED (numeric matrix at 30 deg)');
  L('WICKWALL',
    'Wick-on-generators wall: Wick-rotating the 8 cycle generators gives a unitary conjugate of so(8), dim 28, NOT su(8). Wick does not commute with generation; the correct su(8) route Wicks the full sl(8) basis (28 keep, 35 x i)',
    'proofs/wick_sl8_wall.md',
    'PROVED (conjugacy identity) / CHECKED (dim-28 closure, compactness)');
  P('MREF4',
    'M_refined generates sl(4,R): raw M_refined closes at dim 16 (gl(4,R)); its traceless correction closes at dim 15, Killing (9,6,0); det factors exactly in t = tan(x/2), vanishing at t^2 = 1/5 in Quadrant I',
    'proofs/m_refined_sl4.md');

  { Walls -- scoped negative results. NOT walls for sl(n) in general:
    walls for the one sampled minimal chain construction only. }
  P('E108-E110',
    'Sampled-minimal-chain wall: the chain reusing the same 2 functions in overlapping sl(2) positions stalls at rank 10 of 15 for sl(4) -- with 4 samples (E108), 8 samples (E109), and Re/Im complex variants (E110) -- because the sampled matrices live in a 4-dim function space (E103). The cycle construction bypasses this wall entirely. NOT a general sl(n) wall',
    'dependency graph E108-E110');

  { Explicitly not theorems }
  ARegistry.Add(TProofObject.Create('MAXWELL-DICT',
    'Maxwell exercise: operator identities PROVED; the E/B identification is a declared dictionary; Maxwell''s equations as spacetime PDEs are NOT derived. A representation witness only',
    psAsserted, 'proofs/maxwell_exercise.md'));
  ARegistry.Add(TProofObject.Create('OPROP',
    'O_prop: no 12x12 candidate passed a leakage test; the proved transport laws model its shape as selection criteria only. The R(theta) = 0 test is NOT executable -- R has no definition in the record (run 233, BLOCKED); pending a definition of R from Kit',
    psOpen, 'proofs/oprop_next_step_correction_DRAFT.md'));
  ARegistry.Add(TProofObject.Create('PARAM-THM',
    'Universal parameter theorem for general traceless 3x3 operator families: OPEN (research order approved 2026-10-12: parameter theorem first, then Z/W)',
    psOpen, 'theorems/INDEX.md, Lie campaign'));
  ARegistry.Add(TProofObject.Create('ZW-MODULI',
    'Z/W complex-structure moduli study: OPEN',
    psOpen, 'theorems/INDEX.md, Lie campaign'));
end;

end.
