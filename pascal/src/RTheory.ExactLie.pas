unit RTheory.ExactLie;

{ EXACT Pascal ports of the sl(3,R) core scripts e70/e71/e72.

  These are not numeric rechecks: every operation runs in the exact
  biquadratic field Q(sqrt(2), sqrt(3)) (RTheory.ExactField), rationals
  are Int64-normalized, and every equality asserted is componentwise
  EXACT (zero residual). The Sage scripts in ProofSource remain the
  proofs of record; these objects re-derive the same results in this
  repository's own language.

    TE70ExactProof  the 8 generators (4 T_a samples at 15/30/45/60 deg
                    + 3 commutators + 1 double commutator, exactly as
                    e70): all traceless exactly, rank exactly 8, all
                    28 commutators exactly in the span.
    TE71ExactProof  structure constants by the pivot solve, with EVERY
                    commutator reconstructed on all 9 entries exactly
                    (e71 spot-checked one); Killing form; signature
                    (5,3,0) by exact Sylvester inertia -- replaces the
                    QQbar eigenvalue-sign count; C2 from the inverse
                    Killing form equals (4/9)*I exactly.
    TE72ExactProof  centralizer of G0 inside the span: exact nullity 2,
                    its two basis elements commute exactly; the
                    antisymmetric parts have exact rank 3 and a chosen
                    triple closes exactly (so(3)); ad_{H0} has exact
                    rank 6 and its Faddeeva-LeVerrier characteristic
                    polynomial is exactly lambda^2 * (even sextic) --
                    eigenvalues in +/- pairs with exactly two zeros
                    (the A2 root pattern), proved by the exact
                    coefficient form where e72 factored over QQbar. }

{$mode objfpc}{$H+}

interface

uses
  SysUtils, RTheory.Types, RTheory.Theorems, RTheory.ExactField;

type
  TE70ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TE71ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

  TE72ExactProof = class(TProofObject)
  public
    constructor Create; reintroduce;
    function VerifyNumeric(out ADetail: string): Boolean; override;
  end;

implementation

type
  TGens = array of TFMatrix; { 8 generators, each 3x3 }
  TStru = array of array of array of TField; { f[i][j][k] }

{ --- small helpers --- }

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

function FTranspose(const M: TFMatrix): TFMatrix;
var
  I, J: Integer;
begin
  Result := FMatZero(Length(M[0]), Length(M));
  for I := 0 to Length(M) - 1 do
    for J := 0 to Length(M[0]) - 1 do
      Result[J][I] := M[I][J];
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

function FMatVecMul(const M: TFMatrix; const V: TFieldArray): TFieldArray;
var
  I, J: Integer;
  S: TField;
begin
  SetLength(Result, Length(M));
  for I := 0 to Length(M) - 1 do
  begin
    S := FldI(0, 1);
    for J := 0 to Length(V) - 1 do
      S := FAdd(S, FMul(M[I][J], V[J]));
    Result[I] := S;
  end;
end;

{ The 8 generators, built exactly as in e70: sin/cos at 15/30/45/60
  degrees are elements of Q(sqrt(2), sqrt(3)); QI formulas
  srx = (1+c)/s, sxp = 1/srx, cxp = (1+s)/c, crx = 1/cxp,
  urx = srx - crx, uxp = cxp - sxp; T_a = [[urx, srx, 0],
  [cxp, -urx, uxp], [0, -uxp, 0]]. }

procedure AddTa(const S, C: TField; var G: TGens; Idx: Integer);
var
  Srx, Sxp, Cxp, Crx, Urx, Uxp: TField;
begin
  Srx := FDiv(FAdd(FldI(1, 1), C), S);
  Sxp := FInv(Srx);
  Cxp := FDiv(FAdd(FldI(1, 1), S), C);
  Crx := FInv(Cxp);
  Urx := FSub(Srx, Crx);
  Uxp := FSub(Cxp, Sxp);
  G[Idx] := FMatZero(3, 3);
  G[Idx][0][0] := Urx;
  G[Idx][0][1] := Srx;
  G[Idx][1][0] := Cxp;
  G[Idx][1][1] := FNeg(Urx);
  G[Idx][1][2] := Uxp;
  G[Idx][2][1] := FNeg(Uxp);
end;

procedure BuildGenerators(out G: TGens);
var
  H: TRational;
begin
  H := Rat(1, 4);
  SetLength(G, 8);
  { 15 deg: sin = (sqrt6 - sqrt2)/4, cos = (sqrt6 + sqrt2)/4 }
  AddTa(Fld(Rat(0, 1), RatNeg(H), Rat(0, 1), H),
        Fld(Rat(0, 1), H, Rat(0, 1), H), G, 0);
  { 30 deg: sin = 1/2, cos = sqrt3/2 }
  AddTa(FldI(1, 2), Fld(Rat(0, 1), Rat(0, 1), Rat(1, 2), Rat(0, 1)), G, 1);
  { 45 deg: sin = cos = sqrt2/2 }
  AddTa(Fld(Rat(0, 1), Rat(1, 2), Rat(0, 1), Rat(0, 1)),
        Fld(Rat(0, 1), Rat(1, 2), Rat(0, 1), Rat(0, 1)), G, 2);
  { 60 deg: sin = sqrt3/2, cos = 1/2 }
  AddTa(Fld(Rat(0, 1), Rat(0, 1), Rat(1, 2), Rat(0, 1)), FldI(1, 2), G, 3);
  G[4] := FComm(G[0], G[1]);
  G[5] := FComm(G[0], G[2]);
  G[6] := FComm(G[0], G[3]);
  G[7] := FComm(G[4], G[2]);
end;

function TFieldVectorsFromMatrix(const M: TFMatrix): TFieldVectors;
var
  I, J: Integer;
begin
  SetLength(Result, Length(M));
  for I := 0 to Length(M) - 1 do
  begin
    SetLength(Result[I], Length(M[0]));
    for J := 0 to Length(M[0]) - 1 do
      Result[I][J] := M[I][J];
  end;
end;

function GFlatOf(const G: TGens): TFieldVectors;
var
  K: Integer;
begin
  SetLength(Result, 8);
  for K := 0 to 7 do
    Result[K] := FFlatten(G[K]);
end;

{ Structure constants by e71's pivot method: pivots of the 8x9
  generator matrix, N[a][k] = Gflat[k][piv[a]], solve N*f = d on the
  pivot coordinates. Returns False if the pivot matrix is singular. }

function BuildStructure(const G: TGens; out Stru: TStru): Boolean;
var
  Rows: TFieldVectors;
  Piv: TIntArray;
  N, NInv, Comm: TFMatrix;
  D: TFieldArray;
  I, J, K, A: Integer;
begin
  Rows := GFlatOf(G);
  if FRowsRank(Rows, Piv) <> 8 then
    Exit(False);
  N := FMatZero(8, 8);
  for A := 0 to 7 do
    for K := 0 to 7 do
      N[A][K] := Rows[K][Piv[A]];
  try
    NInv := FMatInvert(N);
  except
    Exit(False);
  end;
  SetLength(Stru, 8, 8, 8);
  for I := 0 to 7 do
    for J := 0 to 7 do
    begin
      Comm := FComm(G[I], G[J]);
      SetLength(D, 8);
      for A := 0 to 7 do
        D[A] := FFlatten(Comm)[Piv[A]];
      D := FMatVecMul(NInv, D);
      for K := 0 to 7 do
        Stru[I][J][K] := D[K];
    end;
  Result := True;
end;

{ --- TE70ExactProof --- }

constructor TE70ExactProof.Create;
begin
  inherited Create('E70-EXACT',
    'E70 exact port: 8 generators from T_a(15,30,45,60 deg) are traceless, rank exactly 8, and all 28 commutators lie exactly in the span (sl(3,R))',
    psProved, 'proofs/e70_exact.py (+ Pascal exact port, RTheory.ExactLie)');
end;

function TE70ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  G: TGens;
  Rows, Aug: TFieldVectors;
  Rank, I, J, K, InSpan: Integer;
  TrOk: Boolean;
begin
  try
    BuildGenerators(G);
    TrOk := True;
    for I := 0 to 7 do
      if not FIsZero(FMatTrace(G[I])) then
        TrOk := False;
    Rows := GFlatOf(G);
    Rank := FRowsRankSimple(Rows);
    InSpan := 0;
    for I := 0 to 7 do
      for J := I + 1 to 7 do
      begin
        SetLength(Aug, 9);
        for K := 0 to 7 do
          Aug[K] := Rows[K];
        Aug[8] := FFlatten(FComm(G[I], G[J]));
        if FRowsRankSimple(Aug) = 8 then
          Inc(InSpan);
      end;
    ADetail := Format(
      'EXACT (Q(sqrt2,sqrt3), zero-residual): traceless %s; rank %d/8; commutators in span %d/28',
      [BoolToStr(TrOk, 'yes', 'NO'), Rank, InSpan]);
    Result := TrOk and (Rank = 8) and (InSpan = 28);
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

{ --- TE71ExactProof --- }

constructor TE71ExactProof.Create;
begin
  inherited Create('E71-EXACT',
    'E71 exact port: structure constants reconstruct every commutator exactly; Killing signature (5,3,0) by exact inertia; quadratic Casimir C2 = (4/9)*I exactly',
    psProved, 'proofs/e71_exact_casimir.py (+ Pascal exact port, RTheory.ExactLie)');
end;

function TE71ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  G: TGens;
  Stru: TStru;
  Km, Kinv, C2, Term, Recon, Comm: TFMatrix;
  I, J, K, L, A, B, ReconOk: Integer;
  S, Target: TField;
  Pos, Neg, Zero: Integer;
  C2Ok: Boolean;
begin
  try
    BuildGenerators(G);
    if not BuildStructure(G, Stru) then
    begin
      ADetail := 'EXACT PORT FAILED: structure-constant pivot solve failed';
      Exit(False);
    end;
    { Full reconstruction: every commutator, all 9 entries, exact. }
    ReconOk := 0;
    for I := 0 to 7 do
      for J := 0 to 7 do
      begin
        Recon := FMatZero(3, 3);
        for K := 0 to 7 do
          Recon := FMatAdd(Recon, FScaleMat(Stru[I][J][K], G[K]));
        Comm := FComm(G[I], G[J]);
        if MatAllZero(FMatSub(Recon, Comm)) then
          Inc(ReconOk);
      end;
    { Killing form: K[i][j] = sum_{k,l} f[i][k][l] * f[j][l][k]. }
    Km := FMatZero(8, 8);
    for I := 0 to 7 do
      for J := 0 to 7 do
      begin
        S := FldI(0, 1);
        for K := 0 to 7 do
          for L := 0 to 7 do
            S := FAdd(S, FMul(Stru[I][K][L], Stru[J][L][K]));
        Km[I][J] := S;
      end;
    FSymmetricInertia(Km, Pos, Neg, Zero);
    Kinv := FMatInvert(Km);
    C2 := FMatZero(3, 3);
    for I := 0 to 7 do
      for J := 0 to 7 do
      begin
        Term := FScaleMat(Kinv[I][J], FMatMul(G[I], G[J]));
        C2 := FMatAdd(C2, Term);
      end;
    Target := FldI(4, 9);
    C2Ok := True;
    for A := 0 to 2 do
      for B := 0 to 2 do
        if A = B then
        begin
          if not FEq(C2[A][B], Target) then
            C2Ok := False;
        end
        else if not FIsZero(C2[A][B]) then
          C2Ok := False;
    ADetail := Format(
      'EXACT (Q(sqrt2,sqrt3), zero-residual): commutator reconstructions %d/64; Killing signature (%d,%d,%d) by exact inertia; C2 = (4/9)I %s',
      [ReconOk, Pos, Neg, Zero, BoolToStr(C2Ok, 'exact', 'FAILED')]);
    Result := (ReconOk = 64) and (Pos = 5) and (Neg = 3) and (Zero = 0)
      and C2Ok;
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

{ --- TE72ExactProof --- }

constructor TE72ExactProof.Create;
begin
  inherited Create('E72-EXACT',
    'E72 exact port: centralizer of G0 in the span has exact dim 2 (abelian); antisymmetric parts rank exactly 3 and close (so(3)); ad_{H0} rank 6 with char poly exactly lambda^2 * (even sextic) -- A2 root pairing',
    psProved, 'proofs/e72_exact_rest.py (+ Pascal exact port, RTheory.ExactLie)');
end;

function TE72ExactProof.VerifyNumeric(out ADetail: string): Boolean;
var
  G: TGens;
  Stru: TStru;
  CM: TFMatrix;
  NB: TFieldVectors;
  H: array of TFMatrix;
  AntiM: TGens;
  Anti, ChosenRows, Aug: TFieldVectors;
  ChosenIdx: TIntArray;
  AdM: TFMatrix;
  CP: TFieldArray;
  I, J, K, A, B: Integer;
  AntiRank, AdRank: Integer;
  HComm: TFMatrix;
  HCommOk, CpOk, CloseAll: Boolean;
  S: TField;
begin
  try
    BuildGenerators(G);
    { 1. Centralizer of G0 inside span(G): CM[a][k] = [G0, Gk] flat[a],
      a 9x8 matrix in the 8 unknown coefficients; its nullity (via the
      exact nullspace basis) is the centralizer dimension. }
    CM := FMatZero(9, 8);
    for K := 0 to 7 do
    begin
      HComm := FComm(G[0], G[K]);
      for A := 0 to 8 do
        CM[A][K] := FFlatten(HComm)[A];
    end;
    NB := FMatNullBasis(CM);
    if Length(NB) <> 2 then
    begin
      ADetail := Format(
        'EXACT PORT FAILED: centralizer dim %d (expect 2)', [Length(NB)]);
      Exit(False);
    end;
    SetLength(H, 2);
    for I := 0 to 1 do
    begin
      H[I] := FMatZero(3, 3);
      for K := 0 to 7 do
        H[I] := FMatAdd(H[I], FScaleMat(NB[I][K], G[K]));
    end;
    HComm := FComm(H[0], H[1]);
    HCommOk := True;
    for A := 0 to 2 do
      for B := 0 to 2 do
        if not FIsZero(HComm[A][B]) then
          HCommOk := False;
    { 2. Antisymmetric parts: rank 3; a greedily chosen independent
      triple closes exactly (so(3)). }
    SetLength(AntiM, 8);
    SetLength(Anti, 8);
    for K := 0 to 7 do
    begin
      AntiM[K] := FScaleMat(FldI(1, 2), FMatSub(G[K], FTranspose(G[K])));
      Anti[K] := FFlatten(AntiM[K]);
    end;
    AntiRank := FRowsRankSimple(Anti);
    SetLength(ChosenIdx, 0);
    SetLength(ChosenRows, 0);
    for K := 0 to 7 do
    begin
      SetLength(Aug, Length(ChosenRows) + 1);
      for I := 0 to High(ChosenRows) do
        Aug[I] := ChosenRows[I];
      Aug[High(Aug)] := Anti[K];
      if FRowsRankSimple(Aug) = Length(ChosenRows) + 1 then
      begin
        SetLength(ChosenRows, Length(ChosenRows) + 1);
        ChosenRows[High(ChosenRows)] := Anti[K];
        SetLength(ChosenIdx, Length(ChosenIdx) + 1);
        ChosenIdx[High(ChosenIdx)] := K;
      end;
      if Length(ChosenIdx) = 3 then
        Break;
    end;
    CloseAll := Length(ChosenIdx) = 3;
    if CloseAll then
      for I := 0 to 2 do
        for J := I + 1 to 2 do
        begin
          SetLength(Aug, 4);
          for K := 0 to 2 do
            Aug[K] := ChosenRows[K];
          Aug[3] := FFlatten(FComm(AntiM[ChosenIdx[I]], AntiM[ChosenIdx[J]]));
          if FRowsRankSimple(Aug) <> 3 then
            CloseAll := False;
        end;
    { 3. ad_{H0}: needs the structure constants. (ad_{G_k})[a][b] =
      f[k][b][a], and H0 has coefficients NB[0] in the G basis. }
    if not BuildStructure(G, Stru) then
    begin
      ADetail := 'EXACT PORT FAILED: structure-constant pivot solve failed';
      Exit(False);
    end;
    AdM := FMatZero(8, 8);
    for A := 0 to 7 do
      for B := 0 to 7 do
      begin
        S := FldI(0, 1);
        for K := 0 to 7 do
          S := FAdd(S, FMul(NB[0][K], Stru[K][B][A]));
        AdM[A][B] := S;
      end;
    AdRank := FRowsRankSimple(TFieldVectorsFromMatrix(AdM));
    CP := FMatCharPoly(AdM);
    CpOk := FIsZero(CP[1]) and FIsZero(CP[3]) and FIsZero(CP[5])
      and FIsZero(CP[7]) and FIsZero(CP[8]);
    ADetail := Format(
      'EXACT (Q(sqrt2,sqrt3), zero-residual): centralizer dim %d, basis commutes %s; antisym rank %d, so(3) closure %s; rank(ad_H0) %d; char poly = lambda^2*(even sextic) %s',
      [Length(NB), BoolToStr(HCommOk, 'yes', 'NO'), AntiRank,
       BoolToStr(CloseAll, 'yes', 'NO'), AdRank,
       BoolToStr(CpOk, 'exact', 'FAILED')]);
    Result := (Length(NB) = 2) and HCommOk and (AntiRank = 3)
      and CloseAll and (AdRank = 6) and CpOk;
  except
    on E: Exception do
    begin
      ADetail := 'EXACT PORT FAILED: ' + E.ClassName + ': ' + E.Message;
      Result := False;
    end;
  end;
end;

end.
