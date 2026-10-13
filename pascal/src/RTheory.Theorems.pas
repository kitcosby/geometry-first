unit RTheory.Theorems;

{ Every theorem in theorems/INDEX.md as a Pascal object.

  TProofObject is the base class: Id, Statement, Status (scope label)
  and ProofSource (the Sage script or algebraic derivation of record).
  TProofRegistry owns the objects and answers the questions the index
  answers: how many are PROVED, where is T1, and so on.

  RegisterAllProofs mirrors theorems/INDEX.md on GitHub row by row
  (31 PROVED rows + T4* strong evidence = 32 objects; the index
  header's "27 PROVED" count is stale against its own rows -- flagged
  here, not silently corrected). One deliberate deviation, flagged in the entry
  itself: SINPROD carries the corrected global form
  urx*uxp = 4*|csc(2x)| (run 150b) -- the GitHub index line still shows
  the older unsigned form, which holds only where sin(2x) > 0. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Classes, RTheory.Types;

type
  TProofObject = class
  private
    FId: string;
    FStatement: string;
    FStatus: TProofStatus;
    FProofSource: string;
  public
    constructor Create(const AId, AStatement: string;
      AStatus: TProofStatus; const AProofSource: string); virtual;
    property Id: string read FId;
    property Statement: string read FStatement;
    property Status: TProofStatus read FStatus;
    property ProofSource: string read FProofSource;
    function StatusName: string;
    function Describe: string; virtual;
    { For objects whose claim can be re-checked numerically in Pascal
      (see RTheory.Verify): returns True and sets AChecked. The base
      implementation makes no numeric claim -- the proof of record is
      the exact script in ProofSource. }
    function VerifyNumeric(out ADetail: string): Boolean; virtual;
  end;

  TProofRegistry = class
  private
    FItems: TStringList; { owns the objects }
    function GetCount: Integer;
    function GetItem(AIndex: Integer): TProofObject;
  public
    constructor Create;
    destructor Destroy; override;
    procedure Add(AProof: TProofObject);
    function FindById(const AId: string): TProofObject;
    function CountByStatus(AStatus: TProofStatus): Integer;
    property Count: Integer read GetCount;
    property Items[AIndex: Integer]: TProofObject read GetItem; default;
  end;

procedure RegisterAllProofs(ARegistry: TProofRegistry);

implementation

{ TProofObject }

constructor TProofObject.Create(const AId, AStatement: string;
  AStatus: TProofStatus; const AProofSource: string);
begin
  inherited Create;
  FId := AId;
  FStatement := AStatement;
  FStatus := AStatus;
  FProofSource := AProofSource;
end;

function TProofObject.StatusName: string;
begin
  Result := ProofStatusName(FStatus);
end;

function TProofObject.Describe: string;
begin
  Result := Format('%s [%s] %s  (proof: %s)',
    [FId, StatusName, FStatement, FProofSource]);
end;

function TProofObject.VerifyNumeric(out ADetail: string): Boolean;
begin
  ADetail := 'not numerically re-checked in Pascal; proof of record: '
    + FProofSource;
  Result := False;
end;

{ TProofRegistry }

constructor TProofRegistry.Create;
begin
  inherited Create;
  FItems := TStringList.Create;
  FItems.OwnsObjects := True;
end;

destructor TProofRegistry.Destroy;
begin
  FItems.Free;
  inherited Destroy;
end;

function TProofRegistry.GetCount: Integer;
begin
  Result := FItems.Count;
end;

function TProofRegistry.GetItem(AIndex: Integer): TProofObject;
begin
  Result := TProofObject(FItems.Objects[AIndex]);
end;

procedure TProofRegistry.Add(AProof: TProofObject);
begin
  FItems.AddObject(AProof.Id, AProof);
end;

function TProofRegistry.FindById(const AId: string): TProofObject;
var
  I: Integer;
begin
  for I := 0 to FItems.Count - 1 do
    if GetItem(I).Id = AId then
      Exit(GetItem(I));
  Result := nil;
end;

function TProofRegistry.CountByStatus(AStatus: TProofStatus): Integer;
var
  I: Integer;
begin
  Result := 0;
  for I := 0 to FItems.Count - 1 do
    if GetItem(I).Status = AStatus then
      Inc(Result);
end;

{ The index, as objects. Statements follow theorems/INDEX.md. }

procedure RegisterAllProofs(ARegistry: TProofRegistry);

  procedure P(const AId, AStmt, ASrc: string);
  begin
    ARegistry.Add(TProofObject.Create(AId, AStmt, psProved, ASrc));
  end;

begin
  { Clifford structure }
  P('T1', 'Cl(4,4) 16-spinor is irreducible (commutant dim 1, exact over Q)',
    'proofs/q4_exact.py');
  P('T2', 'Sector intertwiner space dim = 16 (exact over Q)',
    'proofs/q5_exact.py');
  P('T3', 'Exists X in intertwiner space with X^2 = +I (explicit, over algebraic closure, verified)',
    'proofs/q6a_groebner.py');
  ARegistry.Add(TProofObject.Create('T4*',
    'X^2 = -I: 2 solutions over algebraic closure, both non-real; 0 real found -- strong evidence, NOT a complete real-algebraic proof of emptiness',
    psStrongEvidence, 'proofs/q6b_groebner.py, proofs/q6b_solutions.py'));

  { Operator algebra }
  P('RECIP', 'sawup + sawdown = FlatWave (exact)', 'algebra');
  P('ONEGEN', 'Every operator is a function of z = srx alone', 'algebra');
  P('FLIP1', 'crx = |sec x| - tan x (exact)', 'sec^2 - tan^2 = 1');
  P('FLIP2', 'sxp = |csc x| - cot x (exact)', 'csc^2 - cot^2 = 1');
  P('INVPAIR', 'c1*c4 = 1, c2*c3 = 1 (exact)', 'definitions');

  { Symmetry }
  P('V4ACT', 'V4 acts on (srx,sxp,cxp,crx) via iota=(0 1)(2 3), tau=(0 3)(1 2), iota_d=(0 2)(1 3)',
    'trig + reciprocal');
  P('V4GRP', '{id, iota, tau, iota_d} is Klein V4 (closure verified)',
    'permutation computation');
  P('CHIF', 'FlatWave character chi_F = (+1,-1,-1,+1)', 'FW definition');

  { Mobius and double-angle }
  P('MOBIUS', 'cxp = (srx + eps) / (eps*srx - 1), eps = FlatWave',
    'equivalent to RECIP');
  P('KRED', 'u = (urx,uxp) = K*P, K = [[1,0,0,-1],[0,-1,1,0]]', 'definitions');
  P('AB', 'A = srx - sxp = 2 cot x, B = cxp - crx = 2 tan x', 'FLIP1, FLIP2');
  P('COSRAT', 'K2/D2 = cos(2x), K2 = A - B, D2 = A + B', 'AB + trig');
  P('SINPROD',
    'urx*uxp = 4*|csc(2x)| globally (corrected form, run 150b; +4 csc(2x) on Q1/Q3, -4 csc(2x) on Q2/Q4; the unsigned form in the GitHub index holds only where sin(2x) > 0)',
    'proofs/run150b_doubleangle_product.py');
  P('HYPER', 'D2^2 - K2^2 = 16', 'COSRAT');
  P('CIRCLE', '(sin 2x, cos 2x) = (4/D2, K2/D2)', 'HYPER');

  { Conic geometry }
  P('QCONIC', 'Q_conic = diag(16,4,-1); x^T Q x = 0 iff 4*H_R^2 + V_R^2 = 1/4',
    'algebra');
  P('QHYP', 'Q_hyp = diag(1,-1,-16); y^T Q y = 0 iff Sigma^2 - Delta^2 = 16',
    'algebra');
  P('HYP2ELL', '(Sigma,Delta) -> (1/Sigma, Delta/(2*Sigma)) = (H_R, V_R) maps hyperbola to ellipse; seam Sigma -> infinity maps to (0, 1/2)',
    'algebra');
  P('VERON', 'Veronese nu_2(cos x, sin x) satisfies X0*X2 - X1^2 = 0',
    'algebra');
  P('SPIN1', '<J_z> = cos 2x, <J_x> = sin 2x; <J_z>/2 = V_R, <J_x>/4 = H_R',
    'spin-1 representation');

  { Dynamics }
  P('Q7', 'dS/dx = (srx*|csc x| - crx*|sec x|) * S^2, S = sawup (Bernoulli, exact per-octant)',
    'proofs/q7_dSdx.py');
  P('DURX', 'd(urx)/dx = crx*|sec x| - srx*|csc x| (exact)', 'Q7 + S = 1/urx');
  P('RUN144', '8 (M,P,G) types, uniquely identified; G flips for c1, c4 at octant boundaries',
    'exact lemma + census');
  P('HARM', 'F = (sin 2x / 4, cos 2x / 2), dF/dx = [[0,1],[-4,0]] * F', 'trig');

  { Intertwiner geometry (D2/D3) }
  P('D2CLOSE', 'V not closed under 16x16 multiplication (B0*B0 not in V; I16 not in V)',
    'proofs/d2_structure.py');
  P('D2GRAM', 'Gram G[i,j] = tr(B_i^T B_j) = 4*I16 (orthogonal, uniform norm 2)',
    'proofs/d2_structure.py');
  P('D3FORM', 'Exactly 2-dim Spin(4,4)-invariant symmetric forms, both rank 8; 0 skew',
    'proofs/d3b_forms.py');
  P('D3VOL', 'vol = E0...E7: vol^2 = +I, not scalar, anticommutes with E_mu',
    'proofs/d3b_forms.py');
end;

end.
