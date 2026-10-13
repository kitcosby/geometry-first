unit RTheory.Wall;

{ The wall as Pascal objects -- docs/WALL.md ("Where derivation
  stops") as registry entries.

  RegisterWall records the missing inputs: objects the proved base
  does NOT give. The gaps are registered psOpen -- targets with no
  claim made -- and the unit also carries proved no-go results that
  fix the wall's location. Numeric probes for three of the wall
  statements live in RTheory.WallVerify.

  No duplication: T4* (RTheory.Theorems) already carries the
  X^2 = -I evidence; WICKWALL and E108-E110 (RTheory.Lie) already
  carry two proved walls; OPROP, PARAM-THM and ZW-MODULI are open
  in RTheory.Lie. Entries here are the WALL.md gaps not otherwise
  registered -- the complex structure itself, SU(8), fermion
  quantum numbers, the physical projection, coupling dynamics,
  the E8 host, and m2/m0 -- plus one proved no-go from 2026-10-13
  (no operator content in so(8) triality).

  Statuses are the record's, never upgraded: a psOpen entry stays
  open until an exact derivation in the proof corpus closes it,
  and nothing in this layer upgrades any label. }

{$mode objfpc}{$H+}

interface

uses
  RTheory.Types, RTheory.Theorems;

procedure RegisterWall(ARegistry: TProofRegistry);

implementation

procedure RegisterWall(ARegistry: TProofRegistry);

  procedure O(const AId, AStmt, ASrc: string);
  begin
    ARegistry.Add(TProofObject.Create(AId, AStmt, psOpen, ASrc));
  end;

begin
  { The gaps -- docs/WALL.md: "the base does not give" }
  O('WALL-CPLX',
    'A complex structure J (J^2 = -I) on R^16 compatible with the intertwiner structure V: not in the proved base. T4* is strong evidence that no real X in V squares to -I (2 complex solutions, 0 real); without J there is no C^8. Any complex structure must enter as a labeled hypothesis at the wall',
    'docs/WALL.md; theorems/INDEX.md T4*');
  O('WALL-SU8',
    'SU(8) built from the proved base: requires C^8, hence blocked at WALL-CPLX; the Wick route to it is independently walled (WICKWALL: Wick-rotating the 8 cycle generators gives a unitary conjugate of so(8), dim 28, not su(8))',
    'docs/WALL.md; proofs/wick_sl8_wall.md');
  O('WALL-FERMION',
    'Fermion quantum numbers (10 + 5bar + 1, or (SU(3), SU(2), Y) assignments): nowhere in the proved base; no derivation known',
    'docs/WALL.md');
  O('WALL-GAMMA',
    'A physical projection Gamma (units, scale, spacetime embedding): not derived. The phase-to-radius map is a coordinate choice, not a theorem',
    'docs/WALL.md');
  O('WALL-COUPLING',
    'Dynamics of the coupling g(x,i): not in the proved base',
    'docs/WALL.md');
  O('WALL-E8HOST',
    'E8(-24) as a representation host for the proved base: not derived. Real-form mismatch (E8 contains compact SO(16); Cl(4,4) gives split SO(4,4)), and every operator construction to date yields classical A/B/D algebras only',
    'docs/WALL.md; proofs/so8_triality.md');
  O('WALL-M2M0',
    'The spectral ratio m2/m0: cannot be fixed from scalar trig alone; evaluable only after a field action or a boundary metric contract is declared',
    'certified dependency spine, Volume IV audit');

  { A proved no-go that locates the wall (2026-10-13) }
  ARegistry.Add(TProofObject.Create('WALL-SO8',
    'No operator content in so(8) triality: the Casimir coincidence C_vector = C_spinor = -7*I_8, unique to n = 8, holds for any nonzero path weights (integers included) -- the trig operators contribute nonzeroness only. The claim that the operators MAKE so(8) special is falsified; so(8) specialness itself stands as pure algebra',
    psProved, 'proofs/so8_triality.md'));
end;

end.
