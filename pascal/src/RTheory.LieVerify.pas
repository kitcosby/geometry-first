unit RTheory.LieVerify;

{ Verifiable Lie-campaign objects.

  A small numeric Lie-closure engine (dense Double matrices, brackets
  by commutator, linear independence by an incremental reduced
  row-echelon basis) plus proof objects that re-check the campaign's
  generation claims numerically at operator values from TOperatorPoint:

    TCycleSl8Proof    cycle generators at x = pi/8 close at dim 63
    TSo9PathProof     path generators close at dim 36, antisymmetric
    TTaFamilyProof    T(x) samples span rank 4 (E103), close at 8 (E70)
    TE106Proof        minimal family spans 4 and closes at 8
    TE107Proof        sl(2) block family spans 2 and closes at 3
    TCycProdProof     prod a_i = 1 and W^8 = I at the octant midpoints
    TM1M2ScalarProof  urx*uxp - srx*crx - sxp*cxp = 2
    TDeltasProof      quarter-turn family swaps (H3/H4)
    TXFamProof        (srx+crx) - (sxp+cxp) = 4 cot(2x)
    TSo3ThetaProof    theta^2 = 2*(sxp^2 - sxp*cxp + cxp^2)

  These are CHECKED-strength numeric confirmations. The proofs of
  record stay the exact/human proofs named in each object's
  ProofSource; a PASS here never upgrades a scope label. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Math, RTheory.Types, RTheory.Operators, RTheory.Theorems,
  RTheory.Verify;

type
  TDynMatrix = array of array of Double;
  TDynVector = array of Double;

  { Incremental RREF basis over flattened matrices. Invariant: every
    stored row has 1 at its own pivot and 0 at every other pivot, so a
    candidate is reduced by one pass over the rows. }
  TSpanBasis = class
  private
    FRows: array of TDynVector;
    FPivots: array of Integer;
    FCount: Integer;
    FLen: Integer;
  public
    constructor Create(ALen: Integer);
    function TryAdd(const AVec: TDynVector): Boolean;
    function RowMatrix(AIndex, AN: Integer): TDynMatrix;
    property Count: Integer read FCount;
  end;

function ZeroMatrix(AN: Integer): TDynMatrix;
function MatMul(const A, B: TDynMatrix): TDynMatrix;
function Commutator(const A, B: TDynMatrix): TDynMatrix;
function FlattenMat(const A: TDynMatrix): TDynVector;
function LieClosureDim(const AGens: array of TDynMatrix; AN: Integer;
  out ABasis: TSpanBasis): Integer;
function SpanDim(const AMats: array of TDynMatrix; AN: Integer): Integer;

type
  TCycleSl8Proof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TSo9PathProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TTaFamilyProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TE106Proof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TE107Proof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TCycProdProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TM1M2ScalarProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TDeltasProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TXFamProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

  TSo3ThetaProof = class(TVerifiableProof)
  protected
    function ResidualAt(AX: Double; out AResidual: Double): Boolean; override;
  public
    constructor Create; reintroduce;
  end;

implementation

{ --- matrix helpers --- }

function ZeroMatrix(AN: Integer): TDynMatrix;
var
  I, J: Integer;
begin
  SetLength(Result, AN, AN);
  for I := 0 to AN - 1 do
    for J := 0 to AN - 1 do
      Result[I][J] := 0.0;
end;

function MatMul(const A, B: TDynMatrix): TDynMatrix;
var
  N, I, J, K: Integer;
  S: Double;
begin
  N := Length(A);
  Result := ZeroMatrix(N);
  for I := 0 to N - 1 do
    for J := 0 to N - 1 do
    begin
      S := 0.0;
      for K := 0 to N - 1 do
        S := S + A[I][K] * B[K][J];
      Result[I][J] := S;
    end;
end;

function Commutator(const A, B: TDynMatrix): TDynMatrix;
var
  N, I, J: Integer;
  AB, BA: TDynMatrix;
begin
  N := Length(A);
  AB := MatMul(A, B);
  BA := MatMul(B, A);
  Result := ZeroMatrix(N);
  for I := 0 to N - 1 do
    for J := 0 to N - 1 do
      Result[I][J] := AB[I][J] - BA[I][J];
end;

function FlattenMat(const A: TDynMatrix): TDynVector;
var
  N, I, J: Integer;
begin
  N := Length(A);
  SetLength(Result, N * N);
  for I := 0 to N - 1 do
    for J := 0 to N - 1 do
      Result[I * N + J] := A[I][J];
end;

function TraceOf(const A: TDynMatrix): Double;
var
  I: Integer;
begin
  Result := 0.0;
  for I := 0 to Length(A) - 1 do
    Result := Result + A[I][I];
end;

{ --- TSpanBasis --- }

constructor TSpanBasis.Create(ALen: Integer);
begin
  inherited Create;
  FLen := ALen;
  FCount := 0;
  SetLength(FRows, 0);
  SetLength(FPivots, 0);
end;

function TSpanBasis.TryAdd(const AVec: TDynVector): Boolean;
var
  V: TDynVector;
  I, J, P: Integer;
  Scale, Best, F: Double;
begin
  SetLength(V, FLen);
  Scale := 0.0;
  for I := 0 to FLen - 1 do
  begin
    V[I] := AVec[I];
    if Abs(AVec[I]) > Scale then
      Scale := Abs(AVec[I]);
  end;
  if Scale = 0.0 then
    Exit(False);
  for I := 0 to FCount - 1 do
  begin
    F := V[FPivots[I]];
    if F <> 0.0 then
      for J := 0 to FLen - 1 do
        V[J] := V[J] - F * FRows[I][J];
  end;
  P := -1;
  Best := Scale * 1e-8;
  for I := 0 to FLen - 1 do
    if Abs(V[I]) > Best then
    begin
      Best := Abs(V[I]);
      P := I;
    end;
  if P < 0 then
    Exit(False);
  F := V[P];
  for J := 0 to FLen - 1 do
    V[J] := V[J] / F;
  for I := 0 to FCount - 1 do
  begin
    F := FRows[I][P];
    if F <> 0.0 then
      for J := 0 to FLen - 1 do
        FRows[I][J] := FRows[I][J] - F * V[J];
  end;
  SetLength(FRows, FCount + 1);
  SetLength(FPivots, FCount + 1);
  FRows[FCount] := V;
  FPivots[FCount] := P;
  Inc(FCount);
  Result := True;
end;

function TSpanBasis.RowMatrix(AIndex, AN: Integer): TDynMatrix;
var
  I, J: Integer;
begin
  Result := ZeroMatrix(AN);
  for I := 0 to AN - 1 do
    for J := 0 to AN - 1 do
      Result[I][J] := FRows[AIndex][I * AN + J];
end;

{ --- span / closure --- }

function SpanDim(const AMats: array of TDynMatrix; AN: Integer): Integer;
var
  B: TSpanBasis;
  I: Integer;
begin
  B := TSpanBasis.Create(AN * AN);
  try
    for I := 0 to High(AMats) do
      B.TryAdd(FlattenMat(AMats[I]));
    Result := B.Count;
  finally
    B.Free;
  end;
end;

function LieClosureDim(const AGens: array of TDynMatrix; AN: Integer;
  out ABasis: TSpanBasis): Integer;
var
  I, J, PassStart: Integer;
  Added: Boolean;
  A, Bm, C: TDynMatrix;
begin
  ABasis := TSpanBasis.Create(AN * AN);
  for I := 0 to High(AGens) do
    ABasis.TryAdd(FlattenMat(AGens[I]));
  repeat
    Added := False;
    PassStart := ABasis.Count;
    for I := 0 to PassStart - 1 do
      for J := I + 1 to PassStart - 1 do
      begin
        A := ABasis.RowMatrix(I, AN);
        Bm := ABasis.RowMatrix(J, AN);
        C := Commutator(A, Bm);
        if ABasis.TryAdd(FlattenMat(C)) then
          Added := True;
      end;
  until not Added;
  Result := ABasis.Count;
end;

{ --- operator values --- }

procedure EightValues(AX: Double; var Vals: array of Double);
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    Vals[0] := P.Srx;  Vals[1] := P.Sxp;
    Vals[2] := P.Cxp;  Vals[3] := P.Crx;
    Vals[4] := P.Urx;  Vals[5] := P.SawUp;
    Vals[6] := P.Uxp;  Vals[7] := P.SawDown;
  finally
    P.Free;
  end;
end;

{ --- TCycleSl8Proof --- }

constructor TCycleSl8Proof.Create;
begin
  inherited Create('CYCLE',
    'Cycle generators X_i = a_i*E_{i,i+1} (mod 8) at x = pi/8 close at dim 63 = sl(8)',
    psProved, 'proofs/sl8_cycle_theorem.tex (human proof, Kit)');
end;

function TCycleSl8Proof.VerifyNumeric(out ADetail: string): Boolean;
var
  Vals: array[0..7] of Double;
  Gens: array of TDynMatrix;
  Basis: TSpanBasis;
  D, I: Integer;
  MaxTrace, Tr: Double;
begin
  EightValues(Pi / 8.0, Vals);
  SetLength(Gens, 8);
  for I := 0 to 7 do
  begin
    Gens[I] := ZeroMatrix(8);
    Gens[I][I][(I + 1) mod 8] := Vals[I];
  end;
  D := LieClosureDim(Gens, 8, Basis);
  try
    MaxTrace := 0.0;
    for I := 0 to Basis.Count - 1 do
    begin
      Tr := Abs(TraceOf(Basis.RowMatrix(I, 8)));
      if Tr > MaxTrace then
        MaxTrace := Tr;
    end;
    ADetail := Format('closure dim %d (expect 63), max |trace| of basis %.3e',
      [D, MaxTrace]);
    Result := (D = 63) and (MaxTrace < 1e-6);
  finally
    Basis.Free;
  end;
end;

{ --- TSo9PathProof --- }

constructor TSo9PathProof.Create;
begin
  inherited Create('SO9PATH',
    'Path generators Y_i = f_i*(E_{i,i+1} - E_{i+1,i}) in 9x9 at x = pi/8 close at dim 36, all basis elements antisymmetric',
    psProved, 'proofs/so9_path_theorem.tex (human proof, Kit)');
end;

function TSo9PathProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Vals: array[0..7] of Double;
  Gens: array of TDynMatrix;
  Basis: TSpanBasis;
  D, I, J, K: Integer;
  M: TDynMatrix;
  MaxDev, Dev: Double;
begin
  EightValues(Pi / 8.0, Vals);
  SetLength(Gens, 8);
  for I := 0 to 7 do
  begin
    Gens[I] := ZeroMatrix(9);
    Gens[I][I][I + 1] := Vals[I];
    Gens[I][I + 1][I] := -Vals[I];
  end;
  D := LieClosureDim(Gens, 9, Basis);
  try
    MaxDev := 0.0;
    for I := 0 to Basis.Count - 1 do
    begin
      M := Basis.RowMatrix(I, 9);
      for J := 0 to 8 do
        for K := 0 to 8 do
        begin
          Dev := Abs(M[J][K] + M[K][J]);
          if Dev > MaxDev then
            MaxDev := Dev;
        end;
    end;
    ADetail := Format('closure dim %d (expect 36), max |M + M^T| %.3e',
      [D, MaxDev]);
    Result := (D = 36) and (MaxDev < 1e-7);
  finally
    Basis.Free;
  end;
end;

{ --- T(x) family: E103 span + E70 closure --- }

function TMatrixAt(AX: Double): TDynMatrix;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    Result := ZeroMatrix(3);
    Result[0][0] := P.Urx;   Result[0][1] := P.Srx;
    Result[1][0] := P.Cxp;   Result[1][1] := -P.Urx; Result[1][2] := P.Uxp;
    Result[2][1] := -P.Uxp;
  finally
    P.Free;
  end;
end;

constructor TTaFamilyProof.Create;
begin
  inherited Create('E103+E70',
    'T(x) samples span rank 4 (E103) and their brackets close at dim 8 = sl(3,R) (E70-E72 generation)',
    psProved, 'proofs/e100_e103_e104_proved.md; proofs/e70_exact.py');
end;

function TTaFamilyProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Samples: array of TDynMatrix;
  Basis: TSpanBasis;
  SD, CD, K: Integer;
begin
  SetLength(Samples, 4);
  for K := 0 to 3 do
    Samples[K] := TMatrixAt((2 * K + 1) * Pi / 8.0);
  SD := SpanDim(Samples, 3);
  CD := LieClosureDim(Samples, 3, Basis);
  Basis.Free;
  ADetail := Format('sample span %d (expect 4), closure dim %d (expect 8)',
    [SD, CD]);
  Result := (SD = 4) and (CD = 8);
end;

{ --- E106 --- }

function M106At(AX: Double): TDynMatrix;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    Result := ZeroMatrix(3);
    Result[0][1] := P.Srx;
    Result[1][0] := P.Cxp;  Result[1][2] := P.Crx;
    Result[2][1] := P.Sxp;
  finally
    P.Free;
  end;
end;

constructor TE106Proof.Create;
begin
  inherited Create('E106',
    'Minimal family [[0,srx,0],[cxp,0,crx],[0,sxp,0]] spans rank 4 and closes at dim 8 = sl(3,R)',
    psProved, 'proofs/e106_minimal.py (GitHub sha 191e285ecdc4)');
end;

function TE106Proof.VerifyNumeric(out ADetail: string): Boolean;
var
  Samples: array of TDynMatrix;
  Basis: TSpanBasis;
  SD, CD, K: Integer;
begin
  SetLength(Samples, 4);
  for K := 0 to 3 do
    Samples[K] := M106At((2 * K + 1) * Pi / 8.0);
  SD := SpanDim(Samples, 3);
  CD := LieClosureDim(Samples, 3, Basis);
  Basis.Free;
  ADetail := Format('sample span %d (expect 4), closure dim %d (expect 8)',
    [SD, CD]);
  Result := (SD = 4) and (CD = 8);
end;

{ --- E107 --- }

function M107At(AX: Double): TDynMatrix;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    Result := ZeroMatrix(2);
    Result[0][1] := P.Srx;
    Result[1][0] := P.Cxp;
  finally
    P.Free;
  end;
end;

constructor TE107Proof.Create;
begin
  inherited Create('E107',
    'sl(2,R) block family [[0,srx],[cxp,0]] spans rank 2 and closes at dim 3 = sl(2,R)',
    psProved, 'proofs/e107_sl2.py (GitHub sha b682e95a8bd0)');
end;

function TE107Proof.VerifyNumeric(out ADetail: string): Boolean;
var
  Samples: array of TDynMatrix;
  Basis: TSpanBasis;
  SD, CD, K: Integer;
begin
  SetLength(Samples, 4);
  for K := 0 to 3 do
    Samples[K] := M107At((2 * K + 1) * Pi / 8.0);
  SD := SpanDim(Samples, 2);
  CD := LieClosureDim(Samples, 2, Basis);
  Basis.Free;
  ADetail := Format('sample span %d (expect 2), closure dim %d (expect 3)',
    [SD, CD]);
  Result := (SD = 2) and (CD = 3);
end;

{ --- CYCPROD --- }

constructor TCycProdProof.Create;
begin
  inherited Create('CYCPROD',
    'Cycle product: prod a_i = 1 and W^8 = I for the weighted cyclic shift W (order-8 Coxeter element)',
    psProved, 'proofs/expmap_operators.md');
  ScopeLabel := 'PROVED (items 1-3, 6) / CHECKED (items 4-5)';
end;

function TCycProdProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  Vals: array[0..7] of Double;
  W, W2, W4, W8: TDynMatrix;
  I, J: Integer;
  Prod, Dev, R1, R2: Double;
begin
  EightValues(AX, Vals);
  Prod := 1.0;
  for I := 0 to 7 do
    Prod := Prod * Vals[I];
  R1 := Abs(Prod - 1.0);
  W := ZeroMatrix(8);
  for I := 0 to 7 do
    W[I][(I + 1) mod 8] := Vals[I];
  W2 := MatMul(W, W);
  W4 := MatMul(W2, W2);
  W8 := MatMul(W4, W4);
  R2 := 0.0;
  for I := 0 to 7 do
    for J := 0 to 7 do
    begin
      if I = J then
        Dev := Abs(W8[I][J] - 1.0)
      else
        Dev := Abs(W8[I][J]);
      if Dev > R2 then
        R2 := Dev;
    end;
  AResidual := Max(R1, R2);
  Result := True;
end;

{ --- M1M2 scalar --- }

constructor TM1M2ScalarProof.Create;
begin
  inherited Create('M1M2',
    'Scalar invariant: urx*uxp - srx*crx - sxp*cxp = 2 (equivalent to tr(M1*M2) = 3)',
    psProved, 'proofs/cross_m1_m2.md');
  ScopeLabel := 'PROVED (trace identity) / CHECKED (entry forms)';
end;

function TM1M2ScalarProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  try
    AResidual := Abs(P.Urx * P.Uxp - P.Srx * P.Crx - P.Sxp * P.Cxp - 2.0) / 2.0;
    Result := True;
  finally
    P.Free;
  end;
end;

{ --- H3/H4 deltas --- }

constructor TDeltasProof.Create;
begin
  inherited Create('H3H4H6',
    'Quarter-turn swaps: srx(x + pi/2) = crx(x), cxp(x + pi/2) = sxp(x)',
    psProved, 'proofs/harmonics_h3_h4_h6_deltas.md');
end;

function TDeltasProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P, Q: TOperatorPoint;
begin
  P := TOperatorPoint.Create(AX);
  Q := TOperatorPoint.Create(AX + Pi / 2.0);
  try
    AResidual := Max(Abs(Q.Srx - P.Crx), Abs(Q.Cxp - P.Sxp));
    Result := True;
  finally
    P.Free;
    Q.Free;
  end;
end;

{ --- XFAM --- }

constructor TXFamProof.Create;
begin
  inherited Create('XFAM',
    'Cross-family diagonal: (srx+crx) - (sxp+cxp) = 4*cot(2x)',
    psProved, 'proofs/cross_family_sum_difference_DRAFT.md');
end;

function TXFamProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
  Expected: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    Expected := 4.0 * Cos(2.0 * AX) / Sin(2.0 * AX);
    AResidual := Abs((P.Srx + P.Crx) - (P.Sxp + P.Cxp) - Expected)
      / Max(Abs(Expected), 1.0);
    Result := True;
  finally
    P.Free;
  end;
end;

{ --- SO3EXP theta^2 --- }

constructor TSo3ThetaProof.Create;
begin
  inherited Create('SO3EXP',
    'theta^2 = sxp^2 + cxp^2 + uxp^2 = 2*(sxp^2 - sxp*cxp + cxp^2) (A2 Eisenstein norm)',
    psProved, 'proofs/expmap_uxp_so3.md');
  ScopeLabel := 'PROVED (structure, theta^2 identity) / CHECKED (numeric matrix at 30 deg)';
end;

function TSo3ThetaProof.ResidualAt(AX: Double; out AResidual: Double): Boolean;
var
  P: TOperatorPoint;
  LHS, RHS: Double;
begin
  P := TOperatorPoint.Create(AX);
  try
    LHS := P.Sxp * P.Sxp + P.Cxp * P.Cxp + P.Uxp * P.Uxp;
    RHS := 2.0 * (P.Sxp * P.Sxp - P.Sxp * P.Cxp + P.Cxp * P.Cxp);
    AResidual := Abs(LHS - RHS) / Max(Abs(RHS), 1.0);
    Result := True;
  finally
    P.Free;
  end;
end;

end.
