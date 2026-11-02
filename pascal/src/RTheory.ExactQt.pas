unit RTheory.ExactQt;

{ EXACT Pascal ports, batch 2 family B -- the Q(t) family.

  TMRefClosureProof  proofs/m_refined_closure.py (Sage/QQbar samples).
                     All sample points lie in Q(sqrt(2), sqrt(3)), so
                     the QQbar computations port EXACTLY onto TField.
  TMRefDetProof      proofs/m_refined_det.py: det M(t) over Q(t) via
                     the RTheory.RatFunc engine, plus the script's
                     4-sample closure re-run.
  TMRefKillingProof  proofs/m_refined_killing.py: structure constants
                     over TField, Killing form, signature by EXACT
                     Sylvester inertia (FSymmetricInertia) -- the
                     Pascal replacement for the script's QQbar
                     eigenvalue-sign count. Inertia is invariant
                     under congruence; eigenvalue signs are not
                     computed and are not needed.
  TQFlipProof        proofs/quadrant_flip_polar_a.py + _b.py, each
                     script's checks replicated AS WRITTEN. Finding
                     (also visible by running the originals): polar_a's
                     three swapped-role assertions are FALSE as
                     written; polar_b's same-name form is the true
                     statement. This object reports both, unharmonized,
                     and therefore FAILS its all-checks conjunction.
  TMaxwellProof      proofs/maxwell_exercise.py, identities only --
                     the E/B dictionary is ASSERTED and stays out.
  TEntropyProof      proofs/entropy_selector_selfconsistency.py,
                     symbolic identities only -- the sign conclusion
                     needs E16 (proved elsewhere) and the physical
                     reading is Kit's ASSERTED postulate; both stay
                     out of this object.

  EXACT arithmetic only; no floating point in any result.

  RESOLVED BLOCKAGE (2026-10-17): the three MREF objects initially
  raised ERationalOverflow inside RTheory.ExactField's bigint
  MagDivMod (multiply-subtract borrow underflow on multi-limb
  operands). Root-caused with a minimal repro and FIXED the same day
  in RTheory.ExactField (disclosed batch-1 infrastructure repair,
  see CONVERSION.md batch 2): the borrow step now loops the
  +BIGBASE correction instead of applying it once. All three MREF
  objects pass from the Pascal engine; batch-1 VerifyProofs still
  passes 31/31 with the fixed unit (regression-verified). The
  exception handlers below remain as honest failure reporting if a
  future infrastructure exception ever fires. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, RTheory.Types, RTheory.Theorems, RTheory.ExactField,
  RTheory.RatFunc;

type
  TMRefClosureProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TMRefDetProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TMRefKillingProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TQFlipProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TMaxwellProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TEntropyProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

implementation

type
  TFMatList = array of TFMatrix;
  TPolyMat = array of array of TPoly;

  TSpan = record
    Pivots: TIntArray;
    Rows: TFieldVectors; { reduced row-echelon basis of the span }
  end;

  TCoreRF = record
    Srx, Sxp, Cxp, Crx, Urx, Uxp, Sawup, Sawdown: TRatFunc;
    Sx, Cx, S2, C2: TRatFunc;
  end;

{ ================= field matrix helpers ================= }

function FComm(const A, B: TFMatrix): TFMatrix;
begin
  Result := FMatSub(FMatMul(A, B), FMatMul(B, A));
end;

function FScaleMat(const S: TField; const M: TFMatrix): TFMatrix;
var
  I, J: Integer;
begin
  Result := FMatZero(Length(M), Length(M[0]));
  for I := 0 to Length(M) - 1 do
    for J := 0 to Length(M[0]) - 1 do
      Result[I][J] := FMul(S, M[I][J]);
end;

function MatAllZero(const M: TFMatrix): Boolean;
var
  I, J: Integer;
begin
  Result := True;
  for I := 0 to Length(M) - 1 do
    for J := 0 to Length(M[0]) - 1 do
      if not FIsZero(M[I][J]) then
        Exit(False);
end;

function MatFromFlat(const V: TFieldArray): TFMatrix;
var
  I, J: Integer;
begin
  Result := FMatZero(4, 4);
  for I := 0 to 3 do
    for J := 0 to 3 do
      Result[I][J] := V[I * 4 + J];
end;

{ ================= exact span over TField ================= }
{ The RREF row space is unique, so this basis is Sage's W.basis(). }

function VecClone(const V: TFieldArray): TFieldArray;
var
  I: Integer;
begin
  SetLength(Result, Length(V));
  for I := 0 to High(V) do
    Result[I] := V[I];
end;

function VecIsZero(const V: TFieldArray): Boolean;
var
  I: Integer;
begin
  Result := True;
  for I := 0 to High(V) do
    if not FIsZero(V[I]) then
      Exit(False);
end;

function SpanReduce(const B: TSpan; const V: TFieldArray): TFieldArray;
var
  I, K: Integer;
  Factor: TField;
begin
  Result := VecClone(V);
  for I := 0 to High(B.Rows) do
    if not FIsZero(Result[B.Pivots[I]]) then
    begin
      Factor := Result[B.Pivots[I]];
      for K := 0 to High(Result) do
        Result[K] := FSub(Result[K], FMul(Factor, B.Rows[I][K]));
    end;
end;

function SpanContains(const B: TSpan; const V: TFieldArray): Boolean;
begin
  Result := VecIsZero(SpanReduce(B, V));
end;

function SpanAdd(var B: TSpan; const V: TFieldArray): Boolean;
var
  W: TFieldArray;
  P, I, K, Pos: Integer;
  Factor, Inv: TField;
  NewRows: TFieldVectors;
  NewPiv: TIntArray;
begin
  W := SpanReduce(B, V);
  P := -1;
  for I := 0 to High(W) do
    if not FIsZero(W[I]) then
    begin
      P := I;
      Break;
    end;
  if P < 0 then
    Exit(False);
  Inv := FInv(W[P]);
  for K := 0 to High(W) do
    W[K] := FMul(W[K], Inv);
  { keep the basis reduced: clear column P from the existing rows }
  for I := 0 to High(B.Rows) do
    if not FIsZero(B.Rows[I][P]) then
    begin
      Factor := B.Rows[I][P];
      for K := 0 to High(W) do
        B.Rows[I][K] := FSub(B.Rows[I][K], FMul(Factor, W[K]));
    end;
  Pos := 0;
  while (Pos < Length(B.Pivots)) and (B.Pivots[Pos] < P) do
    Inc(Pos);
  SetLength(NewRows, Length(B.Rows) + 1);
  SetLength(NewPiv, Length(B.Pivots) + 1);
  for I := 0 to Pos - 1 do
  begin
    NewRows[I] := B.Rows[I];
    NewPiv[I] := B.Pivots[I];
  end;
  NewRows[Pos] := W;
  NewPiv[Pos] := P;
  for I := Pos to High(B.Rows) do
  begin
    NewRows[I + 1] := B.Rows[I];
    NewPiv[I + 1] := B.Pivots[I];
  end;
  B.Rows := NewRows;
  B.Pivots := NewPiv;
  Result := True;
end;

{ ================= m_refined samples over TField ================= }

function SampleTs: TFieldArray;
begin
  SetLength(Result, 5);
  Result[0] := FSub(FldI(2, 1), FSqrt3);              { 2 - sqrt(3) }
  Result[1] := FSub(FSqrt2, FldI(1, 1));             { sqrt(2) - 1 }
  Result[2] := Fld(Rat(0, 1), Rat(0, 1), Rat(1, 3), Rat(0, 1)); { 1/sqrt(3) }
  Result[3] := FldI(1, 2);
  Result[4] := FldI(3, 4);
end;

function SampleMatrix(const T: TField; Corrected: Boolean): TFMatrix;
var
  One, T2, Den, Den2: TField;
  Sx, Cx, S2, C2, Srx, Sxp, Cxp, Crx, Urx, Uxp, Sawup, Sawdown: TField;
  Tr, Corr: TField;
  I: Integer;
begin
  One := FldI(1, 1);
  T2 := FMul(T, T);
  Den := FAdd(One, T2);
  Den2 := FMul(Den, Den);
  Sx := FDiv(FAdd(T, T), Den);
  Cx := FDiv(FSub(One, T2), Den);
  S2 := FDiv(FMul(FldI(4, 1), FMul(T, FSub(One, T2))), Den2);
  C2 := FDiv(FAdd(FSub(One, FMul(FldI(6, 1), T2)), FMul(T2, T2)), Den2);
  Srx := FDiv(One, T);
  Sxp := T;
  Cxp := FDiv(FAdd(One, T), FSub(One, T));
  Crx := FDiv(FSub(One, T), FAdd(One, T));
  Urx := FSub(Srx, Crx);
  Uxp := FSub(Cxp, Sxp);
  Sawup := FInv(Urx);
  Sawdown := FInv(Uxp);
  Result := FMatZero(4, 4);
  Result[0][0] := Sx;   Result[0][1] := Cxp;
  Result[0][2] := FNeg(Sxp); Result[0][3] := Uxp;
  Result[1][0] := Srx;  Result[1][1] := Cx;
  Result[1][2] := C2;   Result[1][3] := Sawup;
  Result[2][0] := FNeg(Crx); Result[2][1] := S2;
  Result[3][0] := Urx;  Result[3][1] := Sawdown;
  if Corrected then
  begin
    Tr := FMatTrace(Result);
    Corr := FDiv(Tr, FldI(4, 1));
    for I := 0 to 3 do
      Result[I][I] := FSub(Result[I][I], Corr);
  end;
end;

function BuildSamples(Count: Integer; Corrected: Boolean): TFMatList;
var
  Ts: TFieldArray;
  I: Integer;
begin
  Ts := SampleTs;
  SetLength(Result, Count);
  for I := 0 to Count - 1 do
    Result[I] := SampleMatrix(Ts[I], Corrected);
end;

{ Commutator closure, script-shaped: G grows only by commutators of
  the current list that lie outside the span; iterate to a fixed
  point, then verify every pair of the final list. Membership and
  dimension are exact rank facts over TField (the span is kept as its
  unique RREF basis, Sage's W). }
procedure ComputeClosure(const Gens: TFMatList; out Dim: Integer;
  out Closed: Boolean; out TracelessAll: Boolean; out Basis: TSpan;
  out GenCount: Integer);
var
  AllM: TFMatList;
  I, J, CurCount: Integer;
  Changed: Boolean;
  Comm: TFMatrix;
begin
  Basis.Pivots := nil;
  Basis.Rows := nil;
  SetLength(AllM, Length(Gens));
  for I := 0 to High(Gens) do
  begin
    AllM[I] := Gens[I];
    SpanAdd(Basis, FFlatten(Gens[I]));
  end;
  repeat
    Changed := False;
    CurCount := Length(AllM);
    for I := 0 to CurCount - 2 do
      for J := I + 1 to CurCount - 1 do
      begin
        Comm := FComm(AllM[I], AllM[J]);
        if not SpanContains(Basis, FFlatten(Comm)) then
        begin
          SetLength(AllM, Length(AllM) + 1);
          AllM[High(AllM)] := Comm;
          SpanAdd(Basis, FFlatten(Comm));
          Changed := True;
        end;
      end;
  until not Changed;
  Closed := True;
  for I := 0 to High(AllM) do
    for J := I + 1 to High(AllM) do
      if not SpanContains(Basis, FFlatten(FComm(AllM[I], AllM[J]))) then
        Closed := False;
  TracelessAll := True;
  for I := 0 to High(AllM) do
    if not FIsZero(FMatTrace(AllM[I])) then
      TracelessAll := False;
  Dim := Length(Basis.Rows);
  GenCount := Length(AllM);
end;

{ ================= rational-function core operators ================= }

function BuildCoreRF: TCoreRF;
var
  T, One, T2, Den, Den2: TRatFunc;
begin
  T := RFVar;
  One := RFOne;
  T2 := RFMul(T, T);
  Den := RFAdd(One, T2);
  Den2 := RFMul(Den, Den);
  with Result do
  begin
    Sx := RFDiv(RFMul(RFInt(2, 1), T), Den);
    Cx := RFDiv(RFSub(One, T2), Den);
    S2 := RFDiv(RFMul(RFInt(4, 1), RFMul(T, RFSub(One, T2))), Den2);
    C2 := RFDiv(RFAdd(RFSub(One, RFMul(RFInt(6, 1), T2)), RFMul(T2, T2)), Den2);
    { numerator 1 - 6t^2 + t^4 }
    Srx := RFDiv(One, T);
    Sxp := T;
    Cxp := RFDiv(RFAdd(One, T), RFSub(One, T));
    Crx := RFDiv(RFSub(One, T), RFAdd(One, T));
    Urx := RFSub(Srx, Crx);
    Uxp := RFSub(Cxp, Sxp);
    Sawup := RFInv(Urx);
    Sawdown := RFInv(Uxp);
  end;
end;

{ D(f) = ((1+t^2)/2) * f'(t), the x-derivative in the t-form. }
function DOp(const F: TRatFunc): TRatFunc;
begin
  Result := RFMul(
    RFMake(PolyMake([Rat(1, 1), Rat(0, 1), Rat(1, 1)]),
      PolyMake([Rat(2, 1)])),
    RFDerivative(F));
end;

function RFHalf(const F: TRatFunc): TRatFunc;
begin
  Result := RFMul(F, RFInt(1, 2));
end;

{ ================= TMRefClosureProof ================= }

constructor TMRefClosureProof.Create;
begin
  inherited Create('MREF-CLOSURE',
    'M_refined closure: 5 exact samples in Q(sqrt2,sqrt3) close at dim 16 (gl(4)); traceless-corrected samples close at dim 15 (sl(4,R))',
    psProved, 'proofs/m_refined_closure.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TMRefClosureProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Dim1, Dim2, GC1, GC2: Integer;
  Closed1, Closed2, TrAll1, TrAll2: Boolean;
  Basis: TSpan;
begin
  try
    ComputeClosure(BuildSamples(5, False), Dim1, Closed1, TrAll1, Basis, GC1);
    ComputeClosure(BuildSamples(5, True), Dim2, Closed2, TrAll2, Basis, GC2);
    ADetail := Format(
      'EXACT (Q(sqrt2,sqrt3), zero-residual): original 5 samples: closure dim %d, closed %s, generating list %d mats, traceless-all-gens computed = %s (the script prints this value; it is False in general since tr(sample) = (-t^2+2t+1)/(t^2+1) <> 0 -- reported as computed, not asserted); traceless-corrected: closure dim %d, closed %s',
      [Dim1, BoolToStr(Closed1, 'yes', 'NO'), GC1,
       BoolToStr(TrAll1, 'True', 'False'), Dim2,
       BoolToStr(Closed2, 'yes', 'NO')]);
    Result := (Dim1 = 16) and Closed1 and (Dim2 = 15) and Closed2;
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT EXCEPTION (infrastructure, not a math mismatch): '
        + E.ClassName + ': ' + E.Message
        + ' -- raised inside RTheory.ExactField bigint rational arithmetic';
      Result := False;
    end;
  end;
end;

{ ================= TMRefDetProof ================= }

constructor TMRefDetProof.Create;
begin
  inherited Create('MREF-DET',
    'det M_refined(t) over Q(t) is exactly -(3t+5)(5t^2-1)/((t+1)(t^2+1)^2), factored as A*B as scripted; zero at t^2 = 1/5; 4-sample closure dims 16 / 15',
    psProved, 'proofs/m_refined_det.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TMRefDetProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Ops: TCoreRF;
  M: TRFMatrix;
  Det, AScr, BScr, ACf, BCf: TRatFunc;
  ExpNum, ExpDen: TPoly;
  NumOk, DenOk, AOk, BOk, ABOk: Boolean;
  Dim1, Dim2, GC1, GC2: Integer;
  Closed1, Closed2, TrA, TrB: Boolean;
  Basis: TSpan;
  ClosNote: string;
begin
  try
    Ops := BuildCoreRF;
    SetLength(M, 4, 4);
    M[0][0] := Ops.Sx;  M[0][1] := Ops.Cxp;
    M[0][2] := RFNeg(Ops.Sxp); M[0][3] := Ops.Uxp;
    M[1][0] := Ops.Srx; M[1][1] := Ops.Cx;
    M[1][2] := Ops.C2;  M[1][3] := Ops.Sawup;
    M[2][0] := RFNeg(Ops.Crx); M[2][1] := Ops.S2;
    M[2][2] := RFZero;  M[2][3] := RFZero;
    M[3][0] := Ops.Urx; M[3][1] := Ops.Sawdown;
    M[3][2] := RFZero;  M[3][3] := RFZero;
    Det := RFMatDet(M);
    ExpNum := PolyMake([Rat(5, 1), Rat(3, 1), Rat(-25, 1), Rat(-15, 1)]);
    ExpDen := PolyMake([Rat(1, 1), Rat(1, 1), Rat(2, 1), Rat(2, 1),
      Rat(1, 1), Rat(1, 1)]);
    NumOk := PolyEq(Det.Num, ExpNum);
    DenOk := PolyEq(Det.Den, ExpDen);
    AScr := RFSub(RFMul(RFNeg(Ops.Sxp), Ops.Sawup), RFMul(Ops.Uxp, Ops.C2));
    BScr := RFSub(RFMul(RFNeg(Ops.Crx), Ops.Sawdown), RFMul(Ops.S2, Ops.Urx));
    ACf := RFMake(PolyMake([Rat(1, 1), Rat(0, 1), Rat(-5, 1)]),
      PolyMake([Rat(-1, 1), Rat(1, 1), Rat(-1, 1), Rat(1, 1)]));
    BCf := RFMake(PolyMake([Rat(-5, 1), Rat(2, 1), Rat(3, 1)]),
      PolyMake([Rat(1, 1), Rat(1, 1), Rat(1, 1), Rat(1, 1)]));
    AOk := RFEq(AScr, ACf);
    BOk := RFEq(BScr, BCf);
    ABOk := RFIsZero(RFSub(RFMul(AScr, BScr), Det));
    { The script's second half (4-sample closure over the field) is
      attempted separately so its infrastructure blockage cannot
      mask the completed Q(t) determinant result above. }
    Dim1 := -1;
    Dim2 := -1;
    Closed1 := False;
    Closed2 := False;
    ClosNote := 'not run';
    try
      ComputeClosure(BuildSamples(4, False), Dim1, Closed1, TrA, Basis, GC1);
      ComputeClosure(BuildSamples(4, True), Dim2, Closed2, TrB, Basis, GC2);
      ClosNote := Format('4-sample closure dim %d closed %s; traceless-corrected dim %d closed %s',
        [Dim1, BoolToStr(Closed1, 'yes', 'NO'),
         Dim2, BoolToStr(Closed2, 'yes', 'NO')]);
    except
      on E: Exception do
        ClosNote := '4-sample closure EXCEPTION (infrastructure, not a math mismatch): '
          + E.ClassName + ': ' + E.Message
          + ' -- raised inside RTheory.ExactField bigint rational arithmetic';
    end;
    ADetail := Format(
      'EXACT (Q(t), zero-residual): det numerator computed = %s (expect -15t^3-25t^2+3t+5: %s); det denominator computed = %s (expect (t+1)(t^2+1)^2 expanded: %s); script A == (1-5t^2)/((t-1)(t^2+1)) %s; script B == (t-1)(3t+5)/((t+1)(t^2+1)) %s; A*B - det == 0 %s; %s',
      [PolyToString(Det.Num), BoolToStr(NumOk, 'yes', 'NO'),
       PolyToString(Det.Den), BoolToStr(DenOk, 'yes', 'NO'),
       BoolToStr(AOk, 'yes', 'NO'), BoolToStr(BOk, 'yes', 'NO'),
       BoolToStr(ABOk, 'yes', 'NO'), ClosNote]);
    Result := NumOk and DenOk and AOk and BOk and ABOk and
      (Dim1 = 16) and Closed1 and (Dim2 = 15) and Closed2;
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT EXCEPTION (infrastructure, not a math mismatch): '
        + E.ClassName + ': ' + E.Message
        + ' -- raised inside RTheory.ExactField bigint rational arithmetic';
      Result := False;
    end;
  end;
end;

{ ================= TMRefKillingProof ================= }

constructor TMRefKillingProof.Create;
begin
  inherited Create('MREF-KILLING',
    'Killing form of the 15-dim traceless M_refined algebra (sl(4,R)): signature (9,6,0) by exact Sylvester inertia over Q(sqrt2,sqrt3)',
    psProved, 'proofs/m_refined_killing.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TMRefKillingProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Dim, GC, I, J, K, L: Integer;
  Closed, TrAll: Boolean;
  Basis: TSpan;
  Bm: TFMatList;
  RowsVec: TFieldVectors;
  Piv: TIntArray;
  SMat, SInv, Km: TFMatrix;
  Stru: array of array of array of TField;
  Comm, Recon: TFMatrix;
  Vec, Vp, Y: TFieldArray;
  S: TField;
  ReconOk: Integer;
  Pos, Neg, Zero: Integer;
begin
  try
    ComputeClosure(BuildSamples(5, True), Dim, Closed, TrAll, Basis, GC);
    if Dim <> 15 then
    begin
      ADetail := Format(
        'EXACT PORT: corrected closure dim %d (expect 15); Killing step not run',
        [Dim]);
      Exit(False);
    end;
    SetLength(Bm, Dim);
    for I := 0 to Dim - 1 do
      Bm[I] := MatFromFlat(Basis.Rows[I]);
    SetLength(RowsVec, Dim);
    for I := 0 to Dim - 1 do
      RowsVec[I] := Basis.Rows[I];
    if FRowsRank(RowsVec, Piv) <> Dim then
    begin
      ADetail := 'EXACT PORT FAILED: basis rank recheck failed';
      Exit(False);
    end;
    { Solve Bmat * y = vec(comm) via the pivot-row restriction:
      S[i][j] = Bmat[Piv[i]][j] = Basis[j][Piv[i]] is invertible
      because projection of the row space onto its RREF pivot
      coordinates is an isomorphism. One inverse serves all solves;
      every solution is then verified by full reconstruction. }
    SMat := FMatZero(Dim, Dim);
    for I := 0 to Dim - 1 do
      for J := 0 to Dim - 1 do
        SMat[I][J] := Basis.Rows[J][Piv[I]];
    SInv := FMatInvert(SMat);
    SetLength(Stru, Dim, Dim, Dim);
    ReconOk := 0;
    for I := 0 to Dim - 1 do
      for J := 0 to Dim - 1 do
      begin
        Comm := FComm(Bm[I], Bm[J]);
        Vec := FFlatten(Comm);
        SetLength(Vp, Dim);
        for K := 0 to Dim - 1 do
          Vp[K] := Vec[Piv[K]];
        SetLength(Y, Dim);
        for K := 0 to Dim - 1 do
        begin
          S := FldI(0, 1);
          for L := 0 to Dim - 1 do
            S := FAdd(S, FMul(SInv[K][L], Vp[L]));
          Y[K] := S;
        end;
        for K := 0 to Dim - 1 do
          Stru[I][J][K] := Y[K];
        Recon := FMatZero(4, 4);
        for K := 0 to Dim - 1 do
          Recon := FMatAdd(Recon, FScaleMat(Y[K], Bm[K]));
        if MatAllZero(FMatSub(Recon, Comm)) then
          Inc(ReconOk);
      end;
    Km := FMatZero(Dim, Dim);
    for I := 0 to Dim - 1 do
      for J := 0 to Dim - 1 do
      begin
        S := FldI(0, 1);
        for K := 0 to Dim - 1 do
          for L := 0 to Dim - 1 do
            S := FAdd(S, FMul(Stru[I][K][L], Stru[J][L][K]));
        Km[I][J] := S;
      end;
    FSymmetricInertia(Km, Pos, Neg, Zero);
    ADetail := Format(
      'EXACT (Q(sqrt2,sqrt3), zero-residual): dim %d, closed %s; structure-constant reconstructions %d/%d exact; Killing inertia computed = (%d, %d, %d) (expected sl(4,R) signature (9,6,0): %s)',
      [Dim, BoolToStr(Closed, 'yes', 'NO'), ReconOk, Dim * Dim,
       Pos, Neg, Zero,
       BoolToStr((Pos = 9) and (Neg = 6) and (Zero = 0), 'MATCH',
         'DISCREPANCY vs (9,6,0)')]);
    Result := Closed and (ReconOk = Dim * Dim) and
      (Pos = 9) and (Neg = 6) and (Zero = 0);
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT EXCEPTION (infrastructure, not a math mismatch): '
        + E.ClassName + ': ' + E.Message
        + ' -- raised inside RTheory.ExactField bigint rational arithmetic';
      Result := False;
    end;
  end;
end;

{ ================= TQFlipProof ================= }

type
  TQOps = record
    V: array of TRatFunc; { srx sxp cxp crx urx uxp sawup sawdown }
  end;

function BuildQ(const Alpha, Beta: TRational): TQOps;
var
  T, One, T2, Srx, Cxp: TRatFunc;
begin
  T := RFVar;
  One := RFOne;
  T2 := RFMul(T, T);
  Srx := RFDiv(
    RFAdd(RFMul(RFConst(Alpha), RFAdd(One, T2)), RFSub(One, T2)),
    RFMul(RFInt(2, 1), T));
  Cxp := RFDiv(
    RFAdd(RFMul(RFConst(Beta), RFAdd(One, T2)), RFMul(RFInt(2, 1), T)),
    RFSub(One, T2));
  SetLength(Result.V, 8);
  Result.V[0] := Srx;
  Result.V[1] := RFInv(Srx);
  Result.V[2] := Cxp;
  Result.V[3] := RFInv(Cxp);
  Result.V[4] := RFSub(Srx, Result.V[3]);
  Result.V[5] := RFSub(Cxp, Result.V[1]);
  Result.V[6] := RFInv(Result.V[4]);
  Result.V[7] := RFInv(Result.V[5]);
end;

function SubNegInvT(const F: TRatFunc): TRatFunc;
begin
  Result := RFMobius(F, Rat(0, 1), Rat(-1, 1), Rat(1, 1), Rat(0, 1));
end;

function RFMat2Mul(const A, B: TRFMatrix): TRFMatrix;
var
  I, J, K: Integer;
  S: TRatFunc;
begin
  SetLength(Result, 2, 2);
  for I := 0 to 1 do
    for J := 0 to 1 do
    begin
      S := RFZero;
      for K := 0 to 1 do
        S := RFAdd(S, RFMul(A[I][K], B[K][J]));
      Result[I][J] := S;
    end;
end;

constructor TQFlipProof.Create;
begin
  inherited Create('QFLIP-EXACT',
    'Quadrant flip (pi-shift t -> -1/t): polar_b same-name QI->QIII form + block flip + Q(i)-Laurent w-identities are exact; polar_a''s swapped-role assertions are checked as written and FAIL (finding -- see VerifyNumeric detail)',
    psProved, 'proofs/quadrant_flip_polar_a.py, proofs/quadrant_flip_polar_b.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TQFlipProof.VerifyNumeric(out ADetail: string): Boolean;
const
  Names: array [0 .. 7] of string =
    ('srx', 'sxp', 'cxp', 'crx', 'urx', 'uxp', 'sawup', 'sawdown');
var
  QIa, QIIIa, QIb, QIIIb: TQOps;
  ChkA: array [0 .. 4] of Boolean;
  ChkAName: array [0 .. 4] of string;
  ChkB: array [0 .. 7] of Boolean;
  M1, M2, Dg, Prod: TRFMatrix;
  BlockOk: Boolean;
  WOk: array [0 .. 3] of Boolean;
  AOkN, BOkN, WOkN, I: Integer;
  Fails: string;
  OneCR, ICR: TComplexRat;
  Wp, Wm, Zp, Zm, TwoL, APlus, AMinus: TLaurent;
begin
  try
    { --- polar_a, as written --- }
    QIa := BuildQ(Rat(1, 1), Rat(1, 1));
    QIIIa := BuildQ(Rat(-1, 1), Rat(-1, 1));
    ChkAName[0] := 'polar_a srx_QI->srx_QIII';
    ChkA[0] := RFEq(SubNegInvT(QIa.V[0]), QIIIa.V[0]);
    ChkAName[1] := 'polar_a urx_QI->uxp_QIII';
    ChkA[1] := RFEq(SubNegInvT(QIa.V[4]), QIIIa.V[5]);
    ChkAName[2] := 'polar_a uxp_QI->urx_QIII';
    ChkA[2] := RFEq(SubNegInvT(QIa.V[5]), QIIIa.V[4]);
    ChkAName[3] := 'polar_a sxp_QI->crx_QIII';
    ChkA[3] := RFEq(SubNegInvT(QIa.V[1]), QIIIa.V[3]);
    ChkAName[4] := 'polar_a cxp_QI->cxp_QIII';
    ChkA[4] := RFEq(SubNegInvT(QIa.V[2]), QIIIa.V[2]);
    { --- polar_b, as written: same-name for all 8 keys --- }
    QIb := BuildQ(Rat(1, 1), Rat(1, 1));
    QIIIb := BuildQ(Rat(-1, 1), Rat(-1, 1));
    for I := 0 to 7 do
      ChkB[I] := RFEq(SubNegInvT(QIb.V[I]), QIIIb.V[I]);
    { --- polar_b block flip --- }
    SetLength(M1, 2, 2);
    M1[0][0] := QIIIb.V[4]; M1[0][1] := QIIIb.V[0];
    M1[1][0] := RFNeg(QIIIb.V[2]); M1[1][1] := RFOne;
    SetLength(M2, 2, 2);
    M2[0][0] := QIb.V[5]; M2[0][1] := QIb.V[1];
    M2[1][0] := RFNeg(QIb.V[3]); M2[1][1] := RFOne;
    SetLength(Dg, 2, 2);
    Dg[0][0] := RFOne; Dg[0][1] := RFZero;
    Dg[1][0] := RFZero; Dg[1][1] := RFInt(-1, 1);
    Prod := RFMat2Mul(RFMat2Mul(Dg, M2), Dg);
    BlockOk := RFEq(M1[0][0], Prod[0][0]) and RFEq(M1[0][1], Prod[0][1])
      and RFEq(M1[1][0], Prod[1][0]) and RFEq(M1[1][1], Prod[1][1]);
    { --- polar_a w-identities over Q(i)-Laurent --- }
    OneCR := CR(Rat(1, 1), Rat(0, 1));
    ICR := CR(Rat(0, 1), Rat(1, 1));
    Wp := LCTerm(1, OneCR);
    Wm := LCTerm(-1, OneCR);
    Zp := LCTerm(2, OneCR);
    Zm := LCTerm(-2, OneCR);
    TwoL := LCTerm(0, CR(Rat(2, 1), Rat(0, 1)));
    APlus := LCAdd(LCAdd(Zp, Zm), TwoL);   { z + z^-1 + 2 }
    AMinus := LCSub(LCAdd(Zp, Zm), TwoL);  { z + z^-1 - 2 }
    WOk[0] := LCRatioEq(LCScale(APlus, ICR), LCSub(Zp, Zm),
      LCScale(LCAdd(Wp, Wm), ICR), LCSub(Wp, Wm));
    WOk[1] := LCEq(LCMul(LCAdd(Wp, Wm), LCAdd(Wp, Wm)), APlus);
    WOk[2] := LCRatioEq(LCScale(AMinus, ICR), LCSub(Zp, Zm),
      LCScale(LCSub(Wp, Wm), ICR), LCAdd(Wp, Wm));
    WOk[3] := LCEq(LCMul(LCSub(Wp, Wm), LCSub(Wp, Wm)), AMinus);
    { --- tally --- }
    AOkN := 0;
    Fails := '';
    for I := 0 to 4 do
      if ChkA[I] then
        Inc(AOkN)
      else
        Fails := Fails + ChkAName[I] + '; ';
    BOkN := 0;
    for I := 0 to 7 do
      if ChkB[I] then
        Inc(BOkN)
      else
        Fails := Fails + 'polar_b same-name ' + Names[I] + '; ';
    if not BlockOk then
      Fails := Fails + 'polar_b block flip; ';
    WOkN := 0;
    for I := 0 to 3 do
      if WOk[I] then
        Inc(WOkN)
      else
        Fails := Fails + 'w-identity ' + IntToStr(I) + '; ';
    if Fails = '' then
      Fails := 'none';
    ADetail := Format(
      'EXACT (Q(t) + Q(i)-Laurent, zero-residual): polar_b same-name checks %d/8; block flip %s; w-identities %d/4 (ratio + perfect-square, both alpha); polar_a t-checks as written %d/5. FAILING: %s -- FINDING: polar_a''s swapped-role assertions are false as written (the original sympy script itself prints False for exactly urx->uxp, uxp->urx, sxp->crx); the true pi-shift statement is polar_b''s same-name form. Not harmonized.',
      [BOkN, BoolToStr(BlockOk, 'yes', 'NO'), WOkN, AOkN, Fails]);
    Result := (AOkN = 5) and (BOkN = 8) and BlockOk and (WOkN = 4);
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

{ ================= TMaxwellProof ================= }

function PolyMatDet3(const M: TPolyMat): TPoly;
var
  I, Counter, X, Digit: Integer;
  Perm: array [0 .. 2] of Integer;
  Used: array [0 .. 2] of Boolean;
  Distinct: Boolean;
  Inversions, A2, B2: Integer;
  Term: TPoly;
begin
  Result := PolyZero;
  for Counter := 0 to 26 do
  begin
    X := Counter;
    for I := 0 to 2 do
      Used[I] := False;
    Distinct := True;
    for I := 0 to 2 do
    begin
      Digit := X mod 3;
      X := X div 3;
      Perm[I] := Digit;
      if Used[Digit] then
        Distinct := False;
      Used[Digit] := True;
    end;
    if not Distinct then
      Continue;
    Inversions := 0;
    for A2 := 0 to 1 do
      for B2 := A2 + 1 to 2 do
        if Perm[A2] > Perm[B2] then
          Inc(Inversions);
    Term := PolyOne;
    for I := 0 to 2 do
      Term := PolyMul(Term, M[I][Perm[I]]);
    if Odd(Inversions) then
      Result := PolySub(Result, Term)
    else
      Result := PolyAdd(Result, Term);
  end;
end;

constructor TMaxwellProof.Create;
begin
  inherited Create('MAXWELL-EXACT',
    'Maxwell exercise identities only (QI, exact over Q(t)): 4 Riccati identities, 2 UNA log-derivative identities, curl-matrix char poly lambda*(lambda^2+a^2+b^2+c^2) by exact grid evaluation. The E/B dictionary is ASSERTED and stays out',
    psProved, 'proofs/maxwell_exercise.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TMaxwellProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Ops: TCoreRF;
  One: TRatFunc;
  C1, C2, C3, C4, C5, C6: Boolean;
  Ia, Ib, Ic, I, J, GridOk, S2Val: Integer;
  PM: TPolyMat;
  Yv: array [0 .. 2, 0 .. 2] of Integer;
  DetP, ExpP: TPoly;
begin
  try
    Ops := BuildCoreRF;
    One := RFOne;
    C1 := RFIsZero(RFAdd(DOp(Ops.Srx),
      RFHalf(RFAdd(One, RFMul(Ops.Srx, Ops.Srx)))));
    C2 := RFIsZero(RFAdd(DOp(Ops.Crx),
      RFHalf(RFAdd(One, RFMul(Ops.Crx, Ops.Crx)))));
    C3 := RFIsZero(RFSub(DOp(Ops.Sxp),
      RFHalf(RFAdd(One, RFMul(Ops.Sxp, Ops.Sxp)))));
    C4 := RFIsZero(RFSub(DOp(Ops.Cxp),
      RFHalf(RFAdd(One, RFMul(Ops.Cxp, Ops.Cxp)))));
    C5 := RFIsZero(RFAdd(RFDiv(DOp(Ops.Urx), Ops.Urx),
      RFHalf(RFAdd(Ops.Srx, Ops.Crx))));
    C6 := RFIsZero(RFSub(RFDiv(DOp(Ops.Uxp), Ops.Uxp),
      RFHalf(RFAdd(Ops.Cxp, Ops.Sxp))));
    { char poly of the curl matrix on the exact 3x3x3 grid (values 0,1,2):
      det(lambda I - Y) with polynomial-in-lambda entries; the
      coefficients are polynomials of degree <= 2 in each of a,b,c,
      so 3 points per variable determine them (interpolation). }
    GridOk := 0;
    for Ia := 0 to 2 do
      for Ib := 0 to 2 do
        for Ic := 0 to 2 do
        begin
          { Y = [[0,-c,b],[c,0,-a],[-b,a,0]] }
          Yv[0][0] := 0;     Yv[0][1] := -Ic;  Yv[0][2] := Ib;
          Yv[1][0] := Ic;    Yv[1][1] := 0;    Yv[1][2] := -Ia;
          Yv[2][0] := -Ib;   Yv[2][1] := Ia;   Yv[2][2] := 0;
          SetLength(PM, 3, 3);
          for I := 0 to 2 do
            for J := 0 to 2 do
              if I = J then
                PM[I][J] := PolySub(PolyVar, PolyConst(Rat(Yv[I][J], 1)))
              else
                PM[I][J] := PolyConst(Rat(-Yv[I][J], 1));
          DetP := PolyMatDet3(PM);
          S2Val := Ia * Ia + Ib * Ib + Ic * Ic;
          ExpP := PolyMake([Rat(0, 1), Rat(S2Val, 1), Rat(0, 1), Rat(1, 1)]);
          if PolyEq(DetP, ExpP) then
            Inc(GridOk);
        end;
    ADetail := Format(
      'EXACT (Q(t), zero-residual): Riccati srx %s, crx %s, sxp %s, cxp %s; UNA log-derivative urx %s, uxp %s; curl char poly det(lambda I - Y) == lambda^3 + (a^2+b^2+c^2) lambda on the exact grid {0,1,2}^3: %d/27 points (grid suffices: coefficients have degree <= 2 in each variable). Excluded per scope: the script''s x-space oscillator (sin(x/2), cos(x/2)) and log-sin items are definition-level trig facts, not part of this exact object; the E/B dictionary is ASSERTED and stays out.',
      [BoolToStr(C1, 'yes', 'NO'), BoolToStr(C2, 'yes', 'NO'),
       BoolToStr(C3, 'yes', 'NO'), BoolToStr(C4, 'yes', 'NO'),
       BoolToStr(C5, 'yes', 'NO'), BoolToStr(C6, 'yes', 'NO'), GridOk]);
    Result := C1 and C2 and C3 and C4 and C5 and C6 and (GridOk = 27);
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

{ ================= TEntropyProof ================= }

constructor TEntropyProof.Create;
begin
  inherited Create('ENTROPY-EXACT',
    'Entropy-selector self-consistency identities only (exact over Q(t)): reciprocal pairs sxp*srx = crx*cxp = 1; mirror uxp(t) = urx((1-t)/(1+t)); derivative identity as scripted in r3. Identities only -- the sign conclusion needs E16 and the physical reading is ASSERTED; both stay out',
    psProved, 'proofs/entropy_selector_selfconsistency.py (+ Pascal exact port, RTheory.ExactQt)');
end;

function TEntropyProof.VerifyNumeric(out ADetail: string): Boolean;
var
  Ops: TCoreRF;
  One: TRatFunc;
  C1, C2, C3, C4: Boolean;
  MirrorArgs: TRatFunc;
begin
  try
    Ops := BuildCoreRF;
    One := RFOne;
    C1 := RFIsZero(RFSub(RFMul(Ops.Sxp, Ops.Srx), One));
    C2 := RFIsZero(RFSub(RFMul(Ops.Crx, Ops.Cxp), One));
    { mirror: uxp(t) == urx(s) at s = (1-t)/(1+t) }
    MirrorArgs := RFMobius(Ops.Urx, Rat(-1, 1), Rat(1, 1), Rat(1, 1),
      Rat(1, 1));
    C3 := RFEq(MirrorArgs, Ops.Uxp);
    { r3 as scripted: D_t(1/uxp)(t) + [D_s(1/urx)](s) evaluated at
      s = (1-t)/(1+t) == 0. The second term is the x-derivative
      function of S = 1/urx evaluated AT the mirror point (D_s =
      ((1+s^2)/2) d/ds); NO dy/dx chain factor is applied -- that
      absence is what makes the identity the content. }
    C4 := RFIsZero(RFAdd(DOp(RFInv(Ops.Uxp)),
      RFMobius(DOp(RFInv(Ops.Urx)), Rat(-1, 1), Rat(1, 1), Rat(1, 1),
        Rat(1, 1))));
    ADetail := Format(
      'EXACT (Q(t), zero-residual): sxp*srx = 1 %s; crx*cxp = 1 %s; mirror uxp(t) = urx((1-t)/(1+t)) %s; r3 derivative identity (no chain factor, as scripted) %s. Scope: identities only; the sign conclusion dS_t/dx <= 0 needs E16 (proved elsewhere); the physical reading is Kit''s ASSERTED postulate and stays out; the original''s numeric corroboration is CHECKED-strength and is not part of this object.',
      [BoolToStr(C1, 'yes', 'NO'), BoolToStr(C2, 'yes', 'NO'),
       BoolToStr(C3, 'yes', 'NO'), BoolToStr(C4, 'yes', 'NO')]);
    Result := C1 and C2 and C3 and C4;
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

end.
