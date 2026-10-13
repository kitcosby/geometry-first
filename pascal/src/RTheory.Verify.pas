unit RTheory.Verify;

{ Verifiable proof objects.

  Each class below is a TProofObject whose claim can be re-checked
  numerically in Pascal: VerifyNumeric evaluates the identity at the
  eight octant midpoints x = (2k+1)*pi/8 (k = 0..7), where the |.|
  branches are constant and nothing is singular, and reports the worst
  relative residual. A PASS here is a CHECKED-strength numeric
  confirmation of a result whose proof of record stays the exact
  derivation in ProofSource -- it never upgrades the scope label. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Math, RTheory.Types, RTheory.Operators, RTheory.Theorems;

type
  TVerifiableProof = class(TProofObject)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; virtual; abstract;
    function Tolerance: Double; virtual;
  public
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TRecipProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TFlip1Proof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TFlip2Proof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TInvPairProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TAbProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TCosRatProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TSinProdProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  THyperProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TCircleProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TMobiusProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  end;

  TE77Proof = class(TVerifiableProof) { the double-angle pair, E77 }
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TFlatWaveSignProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TQ7Proof = class(TVerifiableProof) { dS/dx by central differences }
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
    function Tolerance: Double; override;
  end;

  TDurxProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
    function Tolerance: Double; override;
  end;

  THarmProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
    function Tolerance: Double; override;
  end;

{ Convenience constructors mirroring the registry entries. }
function NewRecipProof: TRecipProof;
function NewFlip1Proof: TFlip1Proof;
function NewFlip2Proof: TFlip2Proof;
function NewInvPairProof: TInvPairProof;
function NewAbProof: TAbProof;
function NewCosRatProof: TCosRatProof;
function NewSinProdProof: TSinProdProof;
function NewHyperProof: THyperProof;
function NewCircleProof: TCircleProof;
function NewMobiusProof: TMobiusProof;

implementation

const
  REL_TOL = 1e-9;
  DIFF_TOL = 1e-4; { central differences, h = 1e-5 }

function RelResid(AComputed, AExpected: Double): Double;
var
  Scale: Double;
begin
  Scale := Max(Abs(AExpected), 1.0);
  Result := Abs(AComputed - AExpected) / Scale;
end;

{ TVerifiableProof }

function TVerifiableProof.Tolerance: Double;
begin
  Result := REL_TOL;
end;

function TVerifiableProof.VerifyNumeric(out ADetail: string): Boolean;
var
  K: Integer;
  X, R, Worst: Double;
begin
  Worst := 0.0;
  for K := 0 to 7 do
  begin
    X := (2 * K + 1) * Pi / 8.0; { octant midpoints }
    if not ResidualAt(X, R) then
    begin
      ADetail := Format('evaluation failed at octant %d', [K]);
      Exit(False);
    end;
    if R > Worst then
      Worst := R;
  end;
  ADetail := Format('worst relative residual %.3e over 8 octant midpoints (tol %.0e)',
    [Worst, Tolerance]);
  Result := Worst <= Tolerance;
end;

{ RECIP: sawup + sawdown = FlatWave is the definition of FlatWave here;
  the substantive check is FlatWave = sgn(sin 2x) -- see
  TFlatWaveSignProof. ResidualAt re-evaluates both saws from the raw
  definitions srx/(cxp*srx - 1) and cxp/(cxp*srx - 1). }

function TRecipProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
  SawDownRaw, SawUpRaw: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    SawDownRaw := P.Srx / (P.Cxp * P.Srx - 1.0);
    SawUpRaw := P.Cxp / (P.Cxp * P.Srx - 1.0);
    AResidual := Max(RelResid(SawUpRaw, P.SawUp),
                     RelResid(SawDownRaw, P.SawDown));
    Result := True;
  finally
    P.Free;
  end;
end;

function TFlip1Proof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := RelResid(P.Crx, Abs(1.0 / Cos(AX)) - Sin(AX) / Cos(AX));
    Result := True;
  finally
    P.Free;
  end;
end;

function TFlip2Proof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := RelResid(P.Sxp, Abs(1.0 / Sin(AX)) - Cos(AX) / Sin(AX));
    Result := True;
  finally
    P.Free;
  end;
end;

function TInvPairProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    { c1*c4 = srx*crx is NOT 1 in general; the index pair products are
      srx*sxp = 1 and crx*cxp = 1 (the reciprocal pairs). }
    AResidual := Max(RelResid(P.Srx * P.Sxp, 1.0),
                     RelResid(P.Crx * P.Cxp, 1.0));
    Result := True;
  finally
    P.Free;
  end;
end;

function TAbProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := Max(RelResid(P.A, 2.0 * Cos(AX) / Sin(AX)),
                     RelResid(P.B, 2.0 * Sin(AX) / Cos(AX)));
    Result := True;
  finally
    P.Free;
  end;
end;

function TCosRatProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := RelResid(P.K2 / P.D2, Cos(2.0 * AX));
    Result := True;
  finally
    P.Free;
  end;
end;

function TSinProdProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    { corrected global form (run 150b): urx*uxp = 4*|csc(2x)| }
    AResidual := RelResid(P.Urx * P.Uxp, 4.0 * Abs(1.0 / Sin(2.0 * AX)));
    Result := True;
  finally
    P.Free;
  end;
end;

function THyperProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := RelResid(P.D2 * P.D2 - P.K2 * P.K2, 16.0);
    Result := True;
  finally
    P.Free;
  end;
end;

function TCircleProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := Max(RelResid(4.0 / P.D2, Sin(2.0 * AX)),
                     RelResid(P.K2 / P.D2, Cos(2.0 * AX)));
    Result := True;
  finally
    P.Free;
  end;
end;

function TMobiusProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
  Eps: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    Eps := P.FlatWave; { = +-1 at octant midpoints }
    AResidual := RelResid(P.Cxp, (P.Srx + Eps) / (Eps * P.Srx - 1.0));
    Result := True;
  finally
    P.Free;
  end;
end;

{ E77 }

constructor TE77Proof.Create;
begin
  inherited Create('E77',
    'sin(2x)/4 = 1/(urx+uxp), cos(2x) = (srx+crx-sxp-cxp)/(urx+uxp)',
    psProved, 'analytic (t = tan(x/2)); Foundation');
end;

function TE77Proof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
  SinPart, CosPart: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    SinPart := RelResid(Sin(2.0 * AX) / 4.0, 1.0 / (P.Urx + P.Uxp));
    CosPart := RelResid((P.Srx + P.Crx - P.Sxp - P.Cxp) / (P.Urx + P.Uxp),
      Cos(2.0 * AX));
    AResidual := Max(SinPart, CosPart);
    Result := True;
  finally
    P.Free;
  end;
end;

{ FlatWave = sgn(sin 2x) }

constructor TFlatWaveSignProof.Create;
begin
  inherited Create('FW-SIGN',
    'FlatWave(x) = sgn(sin 2x) (from FOUNDATION: FlatWave = sign of sin 2x)',
    psProved, 'proofs/pythagorean_foundation.md');
end;

function TFlatWaveSignProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := Abs(P.FlatWave - TCoreOperators.FlatWaveSign(AX));
    Result := True;
  finally
    P.Free;
  end;
end;

{ Q7 / DURX / HARM: derivative identities by central differences }

function TQ7Proof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
const
  H = 1e-5;
var
  P: TOperatorPoint;
  DerivNumeric, DerivFormula: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    DerivNumeric := (TCoreOperators.SawUp(AX + H)
      - TCoreOperators.SawUp(AX - H)) / (2.0 * H);
    DerivFormula := (P.Srx * Abs(1.0 / Sin(AX))
      - P.Crx * Abs(1.0 / Cos(AX))) * P.SawUp * P.SawUp;
    AResidual := RelResid(DerivNumeric, DerivFormula);
    Result := True;
  finally
    P.Free;
  end;
end;

function TQ7Proof.Tolerance: Double;
begin
  Result := DIFF_TOL;
end;

function TDurxProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
const
  H = 1e-5;
var
  P: TOperatorPoint;
  DerivNumeric, DerivFormula: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    DerivNumeric := (TCoreOperators.Urx(AX + H)
      - TCoreOperators.Urx(AX - H)) / (2.0 * H);
    DerivFormula := P.Crx * Abs(1.0 / Cos(AX))
      - P.Srx * Abs(1.0 / Sin(AX));
    AResidual := RelResid(DerivNumeric, DerivFormula);
    Result := True;
  finally
    P.Free;
  end;
end;

function TDurxProof.Tolerance: Double;
begin
  Result := DIFF_TOL;
end;

function THarmProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
const
  H = 1e-5;
var
  F1, F2, DF1, DF2: Double;
begin
  F1 := Sin(2.0 * AX) / 4.0;
  F2 := Cos(2.0 * AX) / 2.0;
  DF1 := (Sin(2.0 * (AX + H)) / 4.0 - Sin(2.0 * (AX - H)) / 4.0) / (2.0 * H);
  DF2 := (Cos(2.0 * (AX + H)) / 2.0 - Cos(2.0 * (AX - H)) / 2.0) / (2.0 * H);
  { dF/dx = [[0,1],[-4,0]] * F }
  AResidual := Max(RelResid(DF1, F2), RelResid(DF2, -4.0 * F1));
  Result := True;
end;

function THarmProof.Tolerance: Double;
begin
  Result := DIFF_TOL;
end;

{ Convenience constructors }

function NewRecipProof: TRecipProof;
begin
  Result := TRecipProof.Create('RECIP', 'sawup + sawdown = FlatWave (exact)',
    psProved, 'algebra');
end;

function NewFlip1Proof: TFlip1Proof;
begin
  Result := TFlip1Proof.Create('FLIP1', 'crx = |sec x| - tan x (exact)',
    psProved, 'sec^2 - tan^2 = 1');
end;

function NewFlip2Proof: TFlip2Proof;
begin
  Result := TFlip2Proof.Create('FLIP2', 'sxp = |csc x| - cot x (exact)',
    psProved, 'csc^2 - cot^2 = 1');
end;

function NewInvPairProof: TInvPairProof;
begin
  Result := TInvPairProof.Create('INVPAIR',
    'srx*sxp = 1, crx*cxp = 1 (the reciprocal pairs; index: c1*c4 = 1, c2*c3 = 1)',
    psProved, 'definitions');
end;

function NewAbProof: TAbProof;
begin
  Result := TAbProof.Create('AB', 'A = srx - sxp = 2 cot x, B = cxp - crx = 2 tan x',
    psProved, 'FLIP1, FLIP2');
end;

function NewCosRatProof: TCosRatProof;
begin
  Result := TCosRatProof.Create('COSRAT', 'K2/D2 = cos(2x), K2 = A - B, D2 = A + B',
    psProved, 'AB + trig');
end;

function NewSinProdProof: TSinProdProof;
begin
  Result := TSinProdProof.Create('SINPROD',
    'urx*uxp = 4*|csc(2x)| globally (corrected form, run 150b)',
    psProved, 'proofs/run150b_doubleangle_product.py');
end;

function NewHyperProof: THyperProof;
begin
  Result := THyperProof.Create('HYPER', 'D2^2 - K2^2 = 16',
    psProved, 'COSRAT');
end;

function NewCircleProof: TCircleProof;
begin
  Result := TCircleProof.Create('CIRCLE', '(sin 2x, cos 2x) = (4/D2, K2/D2)',
    psProved, 'HYPER');
end;

function NewMobiusProof: TMobiusProof;
begin
  Result := TMobiusProof.Create('MOBIUS',
    'cxp = (srx + eps) / (eps*srx - 1), eps = FlatWave',
    psProved, 'equivalent to RECIP');
end;

end.
